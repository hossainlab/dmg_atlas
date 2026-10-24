#!/usr/bin/env python
"""Phase 6 / step 22 — cleaned microenvironment subclustering: drop contamination / stress subclusters of step 17, then one scVI
model per compartment on assay-masked features (GPU; run on Modal), re-cluster and re-annotate.

Inputs : data/processed/atlas_core.h5ad (obs, layers['counts'] via h5py), data/processed/tme_subcluster_labels.csv.gz +
         results/tables/tme_subclusters.csv (17_tme_subclusters.py), results/tables/assay_biased_genes.csv (20_assay_masked_metaprograms.py),
         dmg_atlas.states.NORMAL_MARKERS / TME_SUBTYPE_MARKERS, results/tables/genesets/phase6_gene_sets.csv (Stress_HSP).
Outputs (never overwritten if present, unless --force):
  results/tables/tme_clean_exclusions.csv          step-17 subclusters with lineage / stress scores and the exclusion call
  results/tables/tme_subclusters_clean.csv         as tme_subclusters.csv, for the cleaned clustering
  results/tables/tme_subcluster_markers_clean.csv  Wilcoxon one-vs-rest DE, top 25 genes per cleaned subcluster
  results/tables/tme_assay_purity.csv              per compartment: share of cells in assay-pure clusters, step 17 vs cleaned
  data/processed/tme_subcluster_labels_clean.csv.gz, models/scvi_tme_<compartment>/ (saved scVI models)
Figure (results/figures/): phase6_tme_umap_clean.png.
Run: nice -n 10 python scripts/22_tme_scvi_clean.py [--force]     (needs a GPU for reasonable wall time)

Method (fixed before running)
  Exclusion of step-17 subclusters (per compartment): score_genes of every NORMAL_MARKERS lineage and of Stress_HSP on the
  compartment's cells; a subcluster is excluded if (a) another lineage's mean marker score exceeds its own lineage's mean score
  (other-lineage contamination / doublets) or (b) it was flagged 'small' (< 50 cells). A mean Stress_HSP score >= 0.5 is only
  flagged ('stress_hsp'), not excluded: in a first attempt the heat-shock criterion removed mainly the scRNA (dissociated-cell)
  clusters, which would have inflated assay purity. T/NK keeps all subclusters (too few cells to drop any).
  Features: 2,000 seurat_v3 HVGs of the kept cells after removing assay-biased + MT-/ribosomal genes (assay_biased_genes.csv).
  scVI 1.4: batch_key = sample_id, categorical covariate assay_class (cells / nuclei), encode_covariates = True, n_latent 20,
  2 layers, NB likelihood, max_epochs 400 (early stopping on 10 % validation, patience 20), seed 0.
  Leiden (igraph, resolution 0.5, seed 0) on the scVI latent, UMAP, Wilcoxon DE and marker-set annotation exactly as step 17.
  Assay purity: share of cells in clusters with frac_nuclei < 0.1 or > 0.9.
"""
import re, sys, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py, scipy.sparse as sp  # noqa: E402
import anndata as ad, scanpy as sc  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS, ROOT  # noqa: E402
from dmg_atlas.states import NORMAL_MARKERS, TME_SUBTYPE_MARKERS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
CORE = PROCESSED / "atlas_core.h5ad"; TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
MODELS = ROOT / "models"
log = U.StepLogger(LOGS / "phase6_tme.log", context="22_tme_scvi_clean")
COMPARTMENTS = ["myeloid", "T_NK", "vascular", "oligodendrocyte"]
HSP_MAX, RES, SEED, N_HVG, MIN_SCORE, MIN_MARGIN, TOP_DE = 0.5, 0.5, 0, 2000, 0.05, 0.05, 25
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}
OUTS = [TABLES / f for f in ("tme_clean_exclusions.csv", "tme_subclusters_clean.csv", "tme_subcluster_markers_clean.csv", "tme_assay_purity.csv")] + \
       [PROCESSED / "tme_subcluster_labels_clean.csv.gz"]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)
import scvi  # noqa: E402
import torch  # noqa: E402
scvi.settings.seed = SEED


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


def purity(frac_nuclei: pd.Series, sizes: pd.Series) -> float:
    pure = (frac_nuclei < 0.1) | (frac_nuclei > 0.9)
    return float(sizes[pure].sum() / sizes.sum())


COLS = ["sample_id", "study", "assay", "consensus_call", "cnv_score", "pct_mito", "doublet_score"]
with log.step("read inputs"):
    f = h5py.File(CORE, "r")
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in COLS}, index=read_elem(f["obs"]["_index"]))
    for c in ("sample_id", "study", "assay", "consensus_call"):
        obs[c] = obs[c].astype(str)
    obs["publication"] = obs.study.map(lambda s: PUBLICATION.get(s, s)); obs["assay_class"] = obs.assay.map(ASSAY_CLASS)
    var = read_elem(f["var"]); genes = var.index.to_numpy()
    lab17 = pd.read_csv(PROCESSED / "tme_subcluster_labels.csv.gz", index_col=0)
    sub17 = pd.read_csv(TABLES / "tme_subclusters.csv").set_index("subcluster")
    ab = pd.read_csv(TABLES / "assay_biased_genes.csv", index_col=0)
    masked = set(ab.index[ab.call != "none"]) | {g for g in genes if re.match(r"^(MT-|RPS\d|RPL\d|RPS[XY]|RPLP)", g)}
    gs = pd.read_csv(TABLES / "genesets" / "phase6_gene_sets.csv"); hsp = gs.gene[gs.gene_set == "Stress_HSP"].tolist()
    log.write(f"{len(lab17)} step-17 cells; {len(masked)} masked genes; GPU {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none'}")

excl_rows, clusters, markers, labels, umaps, pur = [], [], [], [], {}, []
for comp in COMPARTMENTS:
    cells = lab17.index[lab17.compartment == comp]
    rows = np.sort(obs.index.get_indexer(cells))
    with log.step(f"{comp}: exclusion scores ({len(rows)} cells)"):
        X = read_rows(f["layers/counts"], rows, len(genes)).astype(np.float32)
        a = ad.AnnData(X=X, obs=obs.iloc[rows].copy(), var=pd.DataFrame(index=genes)); del X
        a.obs["sub17"] = lab17.subcluster.reindex(a.obs_names).values
        a.layers["counts"] = a.X.copy()
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        for lin, g in NORMAL_MARKERS.items():
            sc.tl.score_genes(a, [x for x in g if x in a.var_names], score_name=f"lin_{lin}", random_state=SEED)
        sc.tl.score_genes(a, [x for x in hsp if x in a.var_names], score_name="hsp", random_state=SEED)
        M = a.obs.groupby("sub17")[[f"lin_{l}" for l in NORMAL_MARKERS] + ["hsp"]].mean()
        drop = []
        for s17, r in M.iterrows():
            own = r[f"lin_{comp}"]; other = r[[f"lin_{l}" for l in NORMAL_MARKERS if l != comp]]
            reasons = [x for x, c in (("other_lineage:" + other.idxmax()[4:], other.max() > own),
                                      ("small", "small" in str(sub17["flags"].get(s17, "")))) if c]
            excluded = bool(reasons) and comp != "T_NK"
            if r.hsp >= HSP_MAX:
                reasons.append("stress_hsp(flag only)")
            if excluded:
                drop.append(s17)
            excl_rows.append(dict(compartment=comp, subcluster=s17, n_cells=int((a.obs.sub17 == s17).sum()), label_17=sub17.label.get(s17),
                                  own_lineage_score=own, top_other_lineage=other.idxmax()[4:], top_other_score=other.max(), hsp_score=r.hsp,
                                  reasons=";".join(reasons), excluded=excluded))
        pur.append(dict(compartment=comp, n_cells_17=len(rows), purity_17=purity(sub17.loc[sub17.compartment == comp, "frac_nuclei"],
                                                                                 sub17.loc[sub17.compartment == comp, "n_cells"])))
        a = a[~a.obs.sub17.isin(drop)].copy()
        log.write(f"{comp}: excluded {drop} -> {a.n_obs} cells kept")
    with log.step(f"{comp}: scVI + clustering ({a.n_obs} cells)"):
        feat = a[:, [g for g in a.var_names if g not in masked]].copy()
        sc.pp.filter_genes(feat, min_cells=10)
        try:
            sc.pp.highly_variable_genes(feat, flavor="seurat_v3", n_top_genes=N_HVG, layer="counts", subset=True)
        except ImportError:  # scikit-misc missing in the image -> dispersion-based HVGs on the log-normalised data
            log.write("seurat_v3 unavailable (scikit-misc); using flavor='seurat' on log-normalised data")
            sc.pp.highly_variable_genes(feat, flavor="seurat", n_top_genes=N_HVG, subset=True)
        scvi.model.SCVI.setup_anndata(feat, layer="counts", batch_key="sample_id", categorical_covariate_keys=["assay_class"])
        m = scvi.model.SCVI(feat, n_latent=20, n_layers=2, gene_likelihood="nb", encode_covariates=True)
        m.train(max_epochs=400, early_stopping=True, early_stopping_patience=20, train_size=0.9, batch_size=256,
                accelerator="gpu" if torch.cuda.is_available() else "cpu", enable_progress_bar=False)
        a.obsm["X_scvi_tme"] = m.get_latent_representation()
        m.save(MODELS / f"scvi_tme_{comp}", overwrite=True, save_anndata=False)
        sc.pp.neighbors(a, use_rep="X_scvi_tme", n_neighbors=15, random_state=SEED)
        sc.tl.leiden(a, resolution=RES, flavor="igraph", n_iterations=2, directed=False, random_state=SEED, key_added="sub")
        sc.tl.umap(a, random_state=SEED)
        sc.tl.rank_genes_groups(a, "sub", method="wilcoxon", pts=True)
        msets = {k: [g for g in v if g in a.var_names] for k, v in TME_SUBTYPE_MARKERS[comp].items()}
        for k, g in msets.items():
            sc.tl.score_genes(a, g, score_name=f"ms_{k}", random_state=SEED)
        sub = a.obs["sub"].astype(int)
        MS = a.obs[[f"ms_{k}" for k in msets]].groupby(sub.values).mean(); MS.columns = list(msets)
        ep = len(m.history["elbo_train"])
        for cl in sorted(sub.unique()):
            idx = (sub.values == cl); o = a.obs.loc[idx]
            de = sc.get.rank_genes_groups_df(a, group=str(cl)).head(TOP_DE); de.insert(0, "subcluster", f"{comp}_c{cl}"); markers.append(de)
            s = MS.loc[cl].sort_values(ascending=False); low = s.iloc[0] <= MIN_SCORE or (s.iloc[0] - s.iloc[1]) < MIN_MARGIN
            ts = o.sample_id.value_counts(normalize=True)
            clusters.append(dict(compartment=comp, subcluster=f"{comp}_c{cl}", n_cells=int(idx.sum()), frac_compartment=idx.mean(),
                                 label=s.index[0] + ("|low_confidence" if low else ""), best_set=s.index[0], best_score=s.iloc[0],
                                 runner_up=s.index[1], margin=s.iloc[0] - s.iloc[1], n_samples=o.sample_id.nunique(),
                                 n_publications=o.publication.nunique(), frac_top_sample=ts.iloc[0], top_sample=ts.index[0],
                                 frac_nuclei=(o.assay_class == "nuclei").mean(), frac_consensus_ambiguous=(o.consensus_call == "ambiguous").mean(),
                                 mean_cnv_score=o.cnv_score.astype(float).mean(), hsp_score=o.hsp.mean(),
                                 flags=";".join(x for x, c in (("sample_specific", ts.iloc[0] > 0.5 or o.sample_id.nunique() < 3), ("small", idx.sum() < 50)) if c),
                                 top_genes=";".join(de.names.head(10)), scvi_epochs=ep))
        cdf = pd.DataFrame([c for c in clusters if c["compartment"] == comp])
        pur[-1].update(n_cells_clean=int(a.n_obs), purity_clean=purity(cdf.frac_nuclei, cdf.n_cells), n_hvg=int(feat.n_vars), scvi_epochs=ep)
        labels.append(pd.DataFrame({"sample_id": a.obs.sample_id.values, "compartment": comp, "subcluster": [f"{comp}_c{c}" for c in sub.values],
                                    "subcluster_17": a.obs.sub17.values}, index=a.obs_names))
        labels[-1]["label"] = labels[-1].subcluster.map(cdf.set_index("subcluster").label)
        umaps[comp] = (a.obsm["X_umap"], sub.values, (a.obs.assay_class == "nuclei").values)
        log.write(f"{comp}: {sub.nunique()} clean subclusters, scVI {ep} epochs, assay purity {pur[-1]['purity_17']:.2f} -> {pur[-1]['purity_clean']:.2f}\n"
                  + cdf[["subcluster", "n_cells", "label", "best_score", "margin", "n_samples", "n_publications", "frac_nuclei", "flags", "top_genes"]].round(2).to_string(index=False))
        del a, feat, m

with log.step("write + figure"):
    pd.DataFrame(excl_rows).to_csv(OUTS[0], index=False); pd.DataFrame(clusters).to_csv(OUTS[1], index=False)
    pd.concat(markers, ignore_index=True).to_csv(OUTS[2], index=False); pd.DataFrame(pur).to_csv(OUTS[3], index=False)
    pd.concat(labels).rename_axis("cell_id").to_csv(OUTS[4])
    log.write("assay purity:\n" + pd.DataFrame(pur).round(3).to_string(index=False))
    apply_figure_style()
    fig, axes = plt.subplots(2, 4, figsize=(7.2, 3.8))
    for j, comp in enumerate(COMPARTMENTS):
        U2, sub, nuc = umaps[comp]; cmap = plt.get_cmap("tab20")
        for cl in np.unique(sub):
            m_ = sub == cl
            axes[0, j].scatter(U2[m_, 0], U2[m_, 1], s=0.3, color=cmap(cl % 20), linewidths=0, rasterized=True)
            cx, cy = np.median(U2[m_], axis=0); axes[0, j].text(cx, cy, str(cl), fontsize=5, ha="center", va="center", weight="bold")
        axes[1, j].scatter(U2[:, 0], U2[:, 1], s=0.3, c=np.where(nuc, "#2166ac", "#b2182b"), linewidths=0, rasterized=True)
        axes[0, j].set_title(f"{comp} (n={len(sub):,})", fontsize=6)
        for ax in axes[:, j]:
            ax.set_xticks([]); ax.set_yticks([])
    axes[1, 0].set_ylabel("assay: nuclei (blue) / cells (red)", fontsize=5)
    fig.savefig(FIGS / "phase6_tme_umap_clean.png"); plt.close(fig)
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
