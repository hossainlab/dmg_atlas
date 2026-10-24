#!/usr/bin/env python
"""Phase 6 / step 17 — microenvironment subclustering (myeloid, T/NK, vascular, oligodendroglia) and marker annotation.

Input : data/processed/atlas_core.h5ad (obs, obsm['X_scvi'], layers['counts']; read row-wise via h5py, never fully loaded),
        dmg_atlas.states.TME_SUBTYPE_MARKERS (reference marker sets per compartment).
Outputs (never overwritten if present, unless --force):
        results/tables/tme_subclusters.csv              one row per subcluster: size, sample/study breadth, assay mix,
                                                        QC indicators, marker-set label (+ margin), top DE genes
        results/tables/tme_subcluster_markers.csv       Wilcoxon one-vs-rest DE, top 25 genes per subcluster
        results/tables/tme_subcluster_marker_scores.csv mean score_genes score per subcluster x marker set
        data/processed/tme_subcluster_labels.csv.gz     per cell: compartment, subcluster, label
Figures (results/figures/): phase6_tme_umap.png, phase6_tme_marker_scores.png.
Run: nice -n 10 python scripts/17_tme_subclusters.py [--force]

Method
  Compartment membership (qc_pass cells, consensus_call != 'malignant'): lineage_label == c, or lineage_label == 'Unknown'
  and scanvi_label == c with scanvi_maxprob >= 0.9. Oligodendroglia = the 'oligodendrocyte' lineage (normal OPCs are not
  separable from OPC-like tumour cells by label, so they enter only if the consensus call kept them normal).
  Per compartment: kNN graph on the atlas scVI latent (X_scvi, batch = sample_id; n_neighbors 15), Leiden (igraph,
  resolution 0.5, seed 0), UMAP (seed 0). Expression = normalize_total(1e4) + log1p of raw counts. DE: Wilcoxon one-vs-rest.
  Annotation: scanpy score_genes per marker set (genes present only); subcluster label = set with the highest mean score;
  'low_confidence' if that mean is <= 0.05 or beats the runner-up by < 0.05. Flags: 'sample_specific' if one sample holds
  > 50 % of the cells or < 3 samples contribute; 'small' if < 50 cells; 'ambiguous_rich' if > 50 % of cells have
  consensus_call == 'ambiguous' (possible tumour contamination / doublets).
"""
import sys, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py, scipy.sparse as sp  # noqa: E402
import anndata as ad, scanpy as sc  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas.states import TME_SUBTYPE_MARKERS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true")
ap.add_argument("--figures-only", action="store_true", help="redraw phase6_tme_marker_scores.png from the saved score table"); args = ap.parse_args()
CORE = PROCESSED / "atlas_core.h5ad"; TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
log = U.StepLogger(LOGS / "phase6_tme.log", context="17_tme_subclusters")
COMPARTMENTS = ["myeloid", "T_NK", "vascular", "oligodendrocyte"]
RES, N_NEIGH, SEED, MAXPROB, MIN_CELLS, MIN_SCORE, MIN_MARGIN, TOP_DE = 0.5, 15, 0, 0.9, 50, 0.05, 0.05, 25
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
NUCLEI = {"snRNA", "multiome_RNA"}
OUTS = [TABLES / "tme_subclusters.csv", TABLES / "tme_subcluster_markers.csv", TABLES / "tme_subcluster_marker_scores.csv",
        PROCESSED / "tme_subcluster_labels.csv.gz"]


def marker_score_figure(S: pd.DataFrame, path):
    """Mean score_genes score per subcluster x marker set, one panel per compartment (audit MC-09 / MC-20: constrained layout,
    one shared labelled colour bar, no in-panel titles; compartments are identified in the caption, left to right in COMPARTMENTS
    order, and by the subcluster labels on the y axis)."""
    import matplotlib.pyplot as plt
    apply_figure_style()
    mats = {c: S[S.compartment == c].drop(columns="compartment").dropna(axis=1, how="all") for c in COMPARTMENTS}
    lim = float(np.nanmax(np.abs(np.concatenate([m.values.ravel() for m in mats.values()]))))
    fig, axes = plt.subplots(1, len(COMPARTMENTS), figsize=(7.2, 3.4), layout="constrained",
                             gridspec_kw={"width_ratios": [m.shape[1] for m in mats.values()]})
    for ax, (comp, m) in zip(axes, mats.items()):
        im = ax.imshow(m.values, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto", interpolation="nearest")
        ax.set_xticks(range(m.shape[1])); ax.set_xticklabels([c.replace("_", " ") for c in m.columns], rotation=60, ha="right", fontsize=4.5)
        ax.set_yticks(range(m.shape[0])); ax.set_yticklabels([i.replace("_", " ") for i in m.index], fontsize=4.5)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
    cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.01, shrink=0.8); cb.set_label("mean score_genes score (subcluster)")
    fig.savefig(path); plt.close(fig)

if args.figures_only:
    marker_score_figure(pd.read_csv(OUTS[2], index_col=0), FIGS / "phase6_tme_marker_scores.png")
    log.write("phase6_tme_marker_scores.png redrawn from the saved table (--figures-only)"); sys.exit(0)
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)


def read_rows(g, rows: np.ndarray, n_genes: int, block: int = 20000) -> sp.csr_matrix:
    """Sorted row subset of an on-disk CSR group (layers/counts), read in contiguous row blocks via h5py."""
    indptr = g["indptr"][:]; parts = []
    for b0 in range(0, len(indptr) - 1, block):
        r = rows[(rows >= b0) & (rows < b0 + block)]
        if not len(r):
            continue
        b1 = min(b0 + block, len(indptr) - 1); s, e = int(indptr[b0]), int(indptr[b1])
        m = sp.csr_matrix((g["data"][s:e], g["indices"][s:e], indptr[b0:b1 + 1] - s), shape=(b1 - b0, n_genes))
        parts.append(m[r - b0])
    return sp.vstack(parts).tocsr()

COLS = ["sample_id", "study", "assay", "patient_id", "lineage_label", "scanvi_label", "scanvi_maxprob", "consensus_call", "qc_pass",
        "cnv_score", "pct_mito", "doublet_score"]
with log.step("read obs/var/X_scvi"):
    f = h5py.File(CORE, "r")
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in COLS}, index=read_elem(f["obs"]["_index"]))
    var = read_elem(f["var"]); Z = f["obsm"]["X_scvi"][:]
    obs["publication"] = obs.study.astype(str).map(lambda s: PUBLICATION.get(s, s))
    lin, scv = obs.lineage_label.astype(str), obs.scanvi_label.astype(str)
    base = obs.qc_pass.astype(bool) & (obs.consensus_call.astype(str) != "malignant")
    member = {c: np.flatnonzero((base & ((lin == c) | ((lin == "Unknown") & (scv == c) & (obs.scanvi_maxprob >= MAXPROB)))).values)
              for c in COMPARTMENTS}
    log.write(f"atlas {obs.shape[0]} cells x {var.shape[0]} genes; compartment sizes " +
              ", ".join(f"{c} {len(i)} ({int((lin.values[i] == c).sum())} lineage_label + {int((lin.values[i] == 'Unknown').sum())} scANVI)"
                        for c, i in member.items()))

clusters, markers, scores, labels, umaps = [], [], [], [], {}
for comp in COMPARTMENTS:
    rows = member[comp]
    with log.step(f"{comp}: subcluster {len(rows)} cells"):
        X = read_rows(f["layers/counts"], rows, var.shape[0]).astype(np.float32)
        a = ad.AnnData(X=X, obs=obs.iloc[rows].copy(), var=var[[]].copy())
        a.obsm["X_scvi"] = Z[rows]
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        sc.pp.neighbors(a, use_rep="X_scvi", n_neighbors=N_NEIGH, random_state=SEED)
        sc.tl.leiden(a, resolution=RES, flavor="igraph", n_iterations=2, directed=False, random_state=SEED, key_added="sub")
        sc.tl.umap(a, random_state=SEED)
        sc.tl.rank_genes_groups(a, "sub", method="wilcoxon", pts=True)
        sets = {k: [g for g in v if g in a.var_names] for k, v in TME_SUBTYPE_MARKERS[comp].items()}
        for k, g in sets.items():
            sc.tl.score_genes(a, g, score_name=f"ms_{k}", random_state=SEED)
        sub = a.obs["sub"].astype(int)
        M = a.obs[[f"ms_{k}" for k in sets]].groupby(sub.values).mean(); M.columns = list(sets)
        for cl in sorted(sub.unique()):
            idx = sub.values == cl; o = a.obs.loc[idx]
            de = sc.get.rank_genes_groups_df(a, group=str(cl)).head(TOP_DE)
            de.insert(0, "subcluster", f"{comp}_{cl}"); markers.append(de)
            s = M.loc[cl].sort_values(ascending=False)
            low = s.iloc[0] <= MIN_SCORE or (s.iloc[0] - s.iloc[1]) < MIN_MARGIN
            top_sample = o.sample_id.astype(str).value_counts(normalize=True)
            flags = [f for f, c in (("sample_specific", top_sample.iloc[0] > 0.5 or o.sample_id.nunique() < 3), ("small", idx.sum() < MIN_CELLS),
                                    ("ambiguous_rich", (o.consensus_call.astype(str) == "ambiguous").mean() > 0.5)) if c]
            clusters.append(dict(compartment=comp, subcluster=f"{comp}_{cl}", n_cells=int(idx.sum()), frac_compartment=idx.mean(),
                                 label=s.index[0] + ("|low_confidence" if low else ""), best_set=s.index[0], best_score=s.iloc[0],
                                 runner_up=s.index[1], margin=s.iloc[0] - s.iloc[1], n_samples=o.sample_id.nunique(),
                                 n_publications=o.publication.nunique(), frac_top_sample=top_sample.iloc[0],
                                 top_sample=top_sample.index[0], frac_nuclei=o.assay.astype(str).isin(NUCLEI).mean(),
                                 frac_lineage_label=(o.lineage_label.astype(str) == comp).mean(),
                                 frac_consensus_ambiguous=(o.consensus_call.astype(str) == "ambiguous").mean(),
                                 mean_cnv_score=o.cnv_score.astype(float).mean(), mean_pct_mito=o.pct_mito.astype(float).mean(),
                                 mean_doublet_score=o.doublet_score.astype(float).mean(), flags=";".join(flags),
                                 top_genes=";".join(de.names.head(10))))
        M.index = [f"{comp}_{i}" for i in M.index]; M.insert(0, "compartment", comp); scores.append(M)
        lab_map = {c["subcluster"]: c["label"] for c in clusters if c["compartment"] == comp}
        labels.append(pd.DataFrame({"sample_id": a.obs.sample_id.astype(str).values, "compartment": comp,
                                    "subcluster": [f"{comp}_{c}" for c in sub.values]}, index=a.obs_names))
        labels[-1]["label"] = labels[-1].subcluster.map(lab_map)
        umaps[comp] = (a.obsm["X_umap"], sub.values, lab_map)
        log.write(f"{comp}: {sub.nunique()} subclusters\n" + pd.DataFrame([c for c in clusters if c["compartment"] == comp])[
            ["subcluster", "n_cells", "label", "best_score", "margin", "n_samples", "n_publications", "frac_top_sample", "frac_nuclei",
             "frac_consensus_ambiguous", "flags", "top_genes"]].round(3).to_string(index=False))
        del a, X

with log.step("write"):
    pd.DataFrame(clusters).to_csv(OUTS[0], index=False)
    pd.concat(markers, ignore_index=True).to_csv(OUTS[1], index=False)
    pd.concat(scores).rename_axis("subcluster").to_csv(OUTS[2])
    pd.concat(labels).rename_axis("cell_id").to_csv(OUTS[3])

with log.step("figures"):
    apply_figure_style()
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.2))
    for ax, comp in zip(axes, COMPARTMENTS):
        U2, sub, lab_map = umaps[comp]
        cmap = plt.get_cmap("tab20")
        for cl in np.unique(sub):
            m = sub == cl
            ax.scatter(U2[m, 0], U2[m, 1], s=0.3, color=cmap(cl % 20), linewidths=0, rasterized=True)
            cx, cy = np.median(U2[m], axis=0)
            ax.text(cx, cy, str(cl), fontsize=5, ha="center", va="center", weight="bold")
        ax.set_title(f"{comp} (n={len(sub):,})", fontsize=6); ax.set_xticks([]); ax.set_yticks([])
    fig.savefig(FIGS / "phase6_tme_umap.png"); plt.close(fig)
    marker_score_figure(pd.concat(scores), FIGS / "phase6_tme_marker_scores.png")
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
