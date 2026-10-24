#!/usr/bin/env python
"""Phase 6 / step 19 — per-cell metaprogram scores, metaprogram families, and comparison with the Jessa 2022 author cNMF programmes.

Inputs : data/processed/atlas_core.h5ad (obs, layers['counts'] via h5py), data/processed/malignant_selection_audit.parquet
         (13a_malignant_selection.py), results/tables/metaprograms.csv + metaprogram_genes.csv (16_metaprograms.py),
         data/processed/cnmf_usages_<jessa sample>.parquet (13_cnmf_per_sample.py).
Outputs (never overwritten if present, unless --force):
  data/processed/metaprogram_scores.parquet   per selected malignant cell: UCell score per MP (MPx), within-sample z-score per MP (z_MPx),
                                              family score (family_<f> = mean of member-MP z-scores), 'metaprogram' = assigned family
                                              (the column 13_reference_signatures.py reads with --use-metaprograms), 'max_family_z',
                                              'top_mp' = highest-scoring lineage MP (within-sample z)
  results/tables/metaprogram_families.csv     MP -> family, class, contamination check (median UCell in malignant vs normal lineages,
                                              myeloid / vascular-to-malignant ratios for EVERY MP, flags at 1.5x and 2x)
  results/tables/metaprogram_family_composition.csv   per sample: family proportions (+ study, assay, n); within-sample z assignment
  results/tables/metaprogram_family_composition_assayz.csv   same for the within-assay-class z assignment (used by 18_composition.py)
  results/tables/metaprogram_family_by_assay.csv      family assignment crosstab by assay class / assay (counts and row fractions)
  results/tables/metaprogram_family_reference_classes.csv  per assigned cell: cell_id, family, max_family_z, sample_id, assay_class,
                                              confident (max z >= 1) — the spatial track builds its cell2location reference from this
  results/tables/jessa_cnmf_vs_metaprograms.csv       per Jessa sample x author programme: Spearman rho with every MP score
  results/tables/jessa_cnmf_recovery.csv              per Jessa sample x author programme: best-matching MP and best-matching
                                                      programme of our per-sample cNMF (usage correlation), with rho
Figure (results/figures/): phase6_jessa_cnmf_vs_metaprograms.png.
Run: nice -n 10 python scripts/19_metaprogram_scores.py [--force]

Method (audit corrections 2026-10-19 marked *)
  Cells*: primary selection of 13a_malignant_selection.py (corrected consensus_call == 'malignant', qc_pass, okada__1734_pons excluded);
          no scANVI rescue (MC-04).
  Scores: UCell (dmg_atlas.metaprograms.ucell_scores, max_rank 1500) of each MP gene set (metaprogram_genes.csv) on raw counts, streamed
          from disk in row blocks.
  Contamination check*: UCell also on up to 3,000 cells per normal lineage (lineage_label == scanvi_label, consensus 'normal'); the
          myeloid- and vascular-to-malignant median ratios are reported for every MP (MC-08); an MP is excluded from the families when
          its class is technical/contamination (16_metaprograms.py) or when either ratio >= CONTAM_FOLD (2.0, pre-registered); the
          1.5x flag is reported as sensitivity (at 1.5x the replicated cycling MP would be excluded because its median UCell is ~0.02
          in every compartment, so the ratio is noise-dominated — a reason not to lower the threshold on medians). Neural lineages and T/NK are reported only (lineage-like tumour programmes legitimately score
          higher in their normal counterparts; the T/NK reference is tiny).
  Families*: family from 16_metaprograms.py (published / curated label -> cycling, OPC_like, AC_like, OC_like, MES_like, hypoxia;
          other curated sets keep their name; unlabelled lineage MPs 'novel_<MP>'). Technical programmes (ribosomal, mitochondrial,
          IEG/dissociation stress) and immune contamination never enter a family (MC-07, P7-04).
  Assignment*: every MP's UCell score is z-scored WITHIN SAMPLE (MC-01: pooled z-scoring made the argmax an assay classifier);
          family score = MEAN of the member MPs' z-scores (not max; max favoured families with more members); assignment = argmax,
          'unassigned' if the maximum family score <= 0. Caveat: within-sample standardisation removes between-sample differences in
          family abundance by construction (every sample's cells are spread over the families by relative rank), so the per-sample
          family composition of this column is NOT a composition measurement. A second assignment with z-scores computed WITHIN ASSAY
          CLASS (cells / nuclei; columns metaprogram_assayz, family_assayz_*, metaprogram_family_composition_assayz.csv) keeps
          between-sample variation inside each assay class; 18_composition.py uses that one for the mp_family composition set.
  Jessa comparison: author columns author_cNMF_program_malignant_<i> are per-sample programme usages (numbering is sample-specific),
          so all correlations are computed within sample over cells with a non-missing author value (>= 50 cells); recovery is
          summarised by assay (MC-11); the heatmap key is the full sample_id (MC-10).
"""
import sys, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py, scipy.sparse as sp  # noqa: E402
import scipy.stats as ss  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
CORE = PROCESSED / "atlas_core.h5ad"; TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
log = U.StepLogger(LOGS / "fix_metaprograms.log", context="19_metaprogram_scores")
NORMAL = ["myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular"]
N_NORMAL, MIN_AUTHOR_CELLS, SEED = 3000, 50, 0
CONTAM_LINEAGES, CONTAM_FOLD, CONTAM_FOLD_SENS, CONFIDENT_Z = ("myeloid", "vascular"), 2.0, 1.5, 1.0
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}
OUTS = [PROCESSED / "metaprogram_scores.parquet", TABLES / "metaprogram_families.csv", TABLES / "metaprogram_family_composition.csv",
        TABLES / "jessa_cnmf_vs_metaprograms.csv", TABLES / "jessa_cnmf_recovery.csv", TABLES / "metaprogram_family_by_assay.csv",
        TABLES / "metaprogram_family_reference_classes.csv"]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)


def ucell_stream(g, rows: np.ndarray, genes: np.ndarray, sets: dict, block: int = 20000) -> pd.DataFrame:
    """UCell scores for a sorted row subset of an on-disk CSR group, read block-wise (bounded RSS)."""
    indptr = g["indptr"][:]; parts = []
    for b0 in range(0, len(indptr) - 1, block):
        r = rows[(rows >= b0) & (rows < b0 + block)]
        if not len(r):
            continue
        b1 = min(b0 + block, len(indptr) - 1); s, e = int(indptr[b0]), int(indptr[b1])
        m = sp.csr_matrix((g["data"][s:e], g["indices"][s:e], indptr[b0:b1 + 1] - s), shape=(b1 - b0, len(genes)))[r - b0]
        parts.append(MP.ucell_scores(m, genes, sets))
    return pd.concat(parts, ignore_index=True)


COLS = ["sample_id", "study", "assay", "lineage_label", "scanvi_label", "consensus_call", "qc_pass"]
with log.step("read obs + select cells"):
    f = h5py.File(CORE, "r")
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in COLS}, index=read_elem(f["obs"]["_index"]))
    for c in COLS:
        if c != "qc_pass":
            obs[c] = obs[c].astype(str)
    obs["assay_class"] = obs.assay.map(ASSAY_CLASS)
    var = read_elem(f["var"]); genes = var.index.to_numpy()
    selp = pd.read_parquet(PROCESSED / "malignant_selection_audit.parquet", columns=["primary", "cnv_confidence"]).reindex(obs.index)
    mal = np.flatnonzero(selp.primary.fillna(False).astype(bool).values)
    rng = np.random.default_rng(SEED); norm = {}
    for lin in NORMAL:
        pool = np.flatnonzero(((obs.lineage_label == lin) & (obs.scanvi_label == lin) & (obs.consensus_call == "normal") & obs.qc_pass.astype(bool)).values)
        norm[lin] = np.sort(rng.choice(pool, size=min(N_NORMAL, len(pool)), replace=False))
    log.write(f"{len(mal)} malignant cells (primary selection; by assay class {obs.assay_class.iloc[mal].value_counts().to_dict()}); "
              f"normal reference cells " + str({k: len(v) for k, v in norm.items()}))

mpt = pd.read_csv(TABLES / "metaprograms.csv")
mpg = pd.read_csv(TABLES / "metaprogram_genes.csv")
sets = {m: list(g.gene) for m, g in mpg.groupby("metaprogram", sort=False)}
sets = {m: sets[m] for m in mpt.metaprogram if m in sets}

with log.step("UCell on malignant + normal reference cells (streamed)"):
    rows = np.sort(np.concatenate([mal] + list(norm.values())))
    S = ucell_stream(f["layers/counts"], rows, genes, sets); S.index = obs.index[rows]
    Smal = S.loc[obs.index[mal]]
    log.write(f"UCell done: {S.shape}; peak RSS {U.peak_rss_gb():.2f} GB")

with log.step("families + contamination check"):
    fam = []
    for r in mpt.itertuples():
        med_m = float(Smal[r.metaprogram].median())
        meds = {lin: float(S.loc[obs.index[idx], r.metaprogram].median()) for lin, idx in norm.items()}
        ratios = {lin: meds[lin] / med_m if med_m > 0 else np.inf for lin in CONTAM_LINEAGES}
        top_lin = max(CONTAM_LINEAGES, key=ratios.get)
        contam15, contam20 = ratios[top_lin] >= CONTAM_FOLD_SENS, ratios[top_lin] >= CONTAM_FOLD
        technical = isinstance(r.technical_class, str)
        family = None if (technical or contam20) else (r.family if isinstance(r.family, str) else f"novel_{r.metaprogram}")
        fam.append(dict(metaprogram=r.metaprogram, label=r.label, label_provenance=r.label_provenance, technical_class=r.technical_class,
                        replicated=r.replicated, replicated_strict=r.replicated_strict, n_samples_cells=r.n_samples_cells,
                        n_samples_nuclei=r.n_samples_nuclei, family=family,
                        excluded_reason=("technical/contamination class" if technical else (f"{top_lin} ratio >= {CONTAM_FOLD}" if contam20 else None)),
                        median_malignant=med_m, ratio_myeloid=ratios["myeloid"], ratio_vascular=ratios["vascular"],
                        top_normal_lineage=top_lin, median_top_normal=meds[top_lin],
                        contamination_flag_1p5x=contam15, contamination_flag_2x=contam20,
                        **{f"median_{k}": v for k, v in meds.items()}))
    fam = pd.DataFrame(fam)
    fam.to_csv(OUTS[1], index=False)
    log.write("\n" + fam[["metaprogram", "label", "family", "excluded_reason", "median_malignant", "ratio_myeloid", "ratio_vascular",
                         "contamination_flag_1p5x", "contamination_flag_2x"]].round(3).to_string(index=False))

with log.step("within-sample z-scores, family assignment"):
    keep = fam.dropna(subset=["family"])
    sid = obs.sample_id.values[mal]
    Zm = Smal.groupby(sid).transform(lambda x: (x - x.mean()) / (x.std(ddof=0) if x.std(ddof=0) > 0 else 1.0))
    F = pd.DataFrame({fm: Zm[g.metaprogram.tolist()].mean(axis=1) for fm, g in keep.groupby("family")})
    assign = F.idxmax(axis=1).where(F.max(axis=1) > 0, "unassigned")
    acl = obs.assay_class.values[mal]
    Za = Smal.groupby(acl).transform(lambda x: (x - x.mean()) / (x.std(ddof=0) if x.std(ddof=0) > 0 else 1.0))
    Fa = pd.DataFrame({fm: Za[g.metaprogram.tolist()].mean(axis=1) for fm, g in keep.groupby("family")})
    assign_a = Fa.idxmax(axis=1).where(Fa.max(axis=1) > 0, "unassigned")
    out = Smal.copy()
    out = out.join(Zm.add_prefix("z_")).join(F.add_prefix("family_")).join(Fa.add_prefix("family_assayz_"))
    out["max_family_z"] = F.max(axis=1); out["metaprogram"] = assign; out["top_mp"] = Zm[keep.metaprogram.tolist()].idxmax(axis=1)
    out["max_family_z_assayz"] = Fa.max(axis=1); out["metaprogram_assayz"] = assign_a
    out.insert(0, "sample_id", sid); out.insert(1, "assay_class", acl)
    out.to_parquet(OUTS[0])
    meta = obs.iloc[mal].groupby("sample_id").agg(study=("study", "first"), assay=("assay", "first"), assay_class=("assay_class", "first"), n_cells=("assay", "size"))
    pd.crosstab(out.sample_id, out.metaprogram, normalize="index").join(meta).to_csv(OUTS[2])
    pd.crosstab(out.sample_id, out.metaprogram_assayz, normalize="index").join(meta).to_csv(TABLES / "metaprogram_family_composition_assayz.csv")
    bas = []
    for col, zs in (("metaprogram", "within_sample"), ("metaprogram_assayz", "within_assay_class")):
        ct = pd.crosstab(out.assay_class, out[col]); ct_a = pd.crosstab(obs.assay.values[mal], out[col].values)
        ba = pd.concat([ct.assign(level="assay_class"), ct_a.assign(level="assay")]); ba["z_scoring"] = zs
        frac = ba.drop(columns=["level", "z_scoring"]).div(ba.drop(columns=["level", "z_scoring"]).sum(1), axis=0).add_suffix("_frac")
        bas.append(pd.concat([ba, frac], axis=1))
    pd.concat(bas).rename_axis("group").to_csv(OUTS[5])
    ref = out.loc[out.metaprogram != "unassigned", ["sample_id", "assay_class", "metaprogram", "max_family_z", "metaprogram_assayz", "max_family_z_assayz"]] \
        .rename(columns={"metaprogram": "family", "metaprogram_assayz": "family_assayz"})
    ref["confident"] = ref.max_family_z >= CONFIDENT_Z
    ref.rename_axis("cell_id").to_csv(OUTS[6])
    ct = pd.crosstab(out.assay_class, out.metaprogram); cta = pd.crosstab(out.assay_class, out.metaprogram_assayz)
    log.write(f"family assignment (within-sample z): {assign.value_counts().to_dict()}\nby assay class (row fractions):\n" + ct.div(ct.sum(1), axis=0).round(3).to_string()
              + f"\nfamily assignment (within-assay-class z): {assign_a.value_counts().to_dict()}\nby assay class (row fractions):\n" + cta.div(cta.sum(1), axis=0).round(3).to_string()
              + f"\nreference classes: {len(ref)} assigned cells, {int(ref.confident.sum())} with max z >= {CONFIDENT_Z}")

with log.step("Jessa author cNMF comparison"):
    auth_cols = [f"author_cNMF_program_malignant_{i}" for i in range(1, 10)]
    A = pd.DataFrame({c: read_elem(f["obs"][c]) for c in auth_cols}, index=obs.index).iloc[mal]
    A.columns = [c.rsplit("_", 1)[1] for c in auth_cols]
    long, rec = [], []
    smeta = obs.iloc[mal].groupby("sample_id")[["study", "assay", "assay_class"]].first()
    for s_, idx in out.groupby("sample_id").groups.items():
        a = A.loc[idx]
        if a.notna().sum().max() < MIN_AUTHOR_CELLS:
            continue
        upath = PROCESSED / f"cnmf_usages_{s_}.parquet"
        Uu = pd.read_parquet(upath).reindex(idx) if upath.exists() else None
        for p in a.columns[a.notna().sum() >= MIN_AUTHOR_CELLS]:
            m = a[p].notna()
            rho = {mp: ss.spearmanr(a.loc[m, p], out.loc[idx[m.values], mp]).statistic for mp in sets}
            for mp, r in rho.items():
                long.append(dict(sample_id=s_, author_program=int(p), n_cells=int(m.sum()), metaprogram=mp, rho=r))
            best = max(rho, key=lambda k: -np.inf if np.isnan(rho[k]) else rho[k])
            row = dict(sample_id=s_, study=smeta.study[s_], assay=smeta.assay[s_], assay_class=smeta.assay_class[s_],
                       author_program=int(p), n_cells=int(m.sum()), best_mp=best, best_mp_rho=rho[best],
                       best_mp_label=fam.set_index("metaprogram").label.get(best), best_mp_family=fam.set_index("metaprogram").family.get(best))
            if Uu is not None:
                uu = Uu.loc[idx[m.values]].dropna(how="all")
                if len(uu) >= MIN_AUTHOR_CELLS:
                    rr = {q: ss.spearmanr(a.loc[uu.index, p], uu[q]).statistic for q in uu.columns}
                    bq = max(rr, key=lambda k: -np.inf if np.isnan(rr[k]) else rr[k])
                    row.update(best_own_program=bq, best_own_program_rho=rr[bq])
            rec.append(row)
    long = pd.DataFrame(long); rec = pd.DataFrame(rec)
    long.to_csv(OUTS[3], index=False); rec.to_csv(OUTS[4], index=False)
    by_assay = rec.groupby("assay").agg(n=("best_own_program_rho", "size"), own_median_rho=("best_own_program_rho", "median"),
                                        own_frac_ge_0p5=("best_own_program_rho", lambda x: (x >= 0.5).mean()),
                                        mp_median_rho=("best_mp_rho", "median"), mp_frac_ge_0p5=("best_mp_rho", lambda x: (x >= 0.5).mean()))
    log.write(f"{rec.sample_id.nunique()} Jessa samples, {len(rec)} author programmes; recovery by assay:\n" + by_assay.round(3).to_string()
              + "\nbest-MP family:\n" + rec.groupby("best_mp_family", dropna=False).agg(n=("best_mp_rho", "size"), median_rho=("best_mp_rho", "median")).round(2).to_string())

with log.step("figure"):
    apply_figure_style()
    M = long.assign(ap=long.sample_id + "|p" + long.author_program.astype(str)).pivot_table(index="ap", columns="metaprogram", values="rho")
    M = M[[m for m in sets if m in M.columns]]
    order = M.idxmax(axis=1).map({m: i for i, m in enumerate(M.columns)}).sort_values().index
    fig, ax = plt.subplots(figsize=(4.6, max(3.0, 0.06 * len(M) + 1.0)))
    im = ax.imshow(M.loc[order].values, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, rotation=90, fontsize=5)
    ax.set_yticks([]); ax.set_ylabel(f"{len(M)} Jessa author cNMF programmes (sample | programme)")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03); cb.set_label("within-sample Spearman ρ")
    fig.savefig(FIGS / "phase6_jessa_cnmf_vs_metaprograms.png"); plt.close(fig)
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
