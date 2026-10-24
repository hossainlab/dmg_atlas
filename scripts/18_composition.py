#!/usr/bin/env python
"""Phase 6 / step 18 — per-sample composition models (cell types, malignant states, TME subclusters) vs clinical covariates.

Inputs : data/processed/atlas_core.h5ad (obs; layers['counts'] rows of selected malignant cells, read via h5py),
         data/processed/tme_subcluster_labels.csv.gz + results/tables/tme_subclusters.csv (17_tme_subclusters.py),
         dmg_atlas.states.STATE_GENE_SETS.
Outputs (results/tables/, never overwritten if present, unless --force):
         composition_per_sample.csv    long: sample x composition set x category -> k, n, proportion, covariates
         composition_results.csv       one row per composition set x category x covariate contrast (primary + secondary model)
         composition_estimability.csv  per covariate: samples / publications per level, publications containing both levels
         malignant_state_assignment.csv.gz  per selected malignant cell: state scores + assigned state
Figure (results/figures/): phase6_composition_effects.png.
Run: nice -n 10 python scripts/18_composition.py [--force]

Audit corrections (2026-10-19; MC-04, MC-06, MC-23, MC-19):
  * malignant cells = primary selection of 13a_malignant_selection.py (corrected consensus 'malignant', qc_pass, okada__1734_pons
    excluded; no scANVI rescue). Scoring is done per sample (score_genes within sample, streamed from disk) and the state is the
    argmax of the WITHIN-SAMPLE z-scored scores ('state'); the pre-audit raw-score argmax is kept as 'state_raw' for comparison.
  * new composition set mp_family: metaprogram-family assignment of 19_metaprogram_scores.py with WITHIN-ASSAY-CLASS z-scoring
    (column metaprogram_assayz; the within-sample column equalises families within every sample and carries no composition signal;
    lineage families only; 'unassigned' kept in the denominator), MIN_N 200.
  * every contrast reports n_patients per arm; the H3.1 vs H3.3 contrast is additionally fitted with site as a covariate
    (coef_site_adj, p_site_adj, q_site_adj) and restricted to pons samples (coef_pons, p_pons, n_pons_ref/test), because all H3.1
    samples are pons nuclei (MC-06). The primary model of every contrast is unchanged otherwise.
  * figure: colour scale covers the full coefficient range (MC-19), human-readable contrast names, caption in the notes.
  * composition_results_vs_pre_audit.csv: coefficient / q comparison with the pre-audit results for celltype and malignant_state.

Composition sets (samples = atlas samples; a sample enters a set if its denominator n >= MIN_N)
  celltype       : qc_pass cells by scanvi_label (malignant + 6 normal lineages); MIN_N 200.
  malignant_state: selected malignant cells (consensus_call == 'malignant' OR scanvi malignant with maxprob >= 0.9) assigned
                   to OPC/AC/OC/MES-like/cycling by scanpy score_genes argmax on normalize_total(1e4)+log1p counts, scored
                   on all selected cells at once (cells with max score <= 0 -> 'unassigned', kept in the denominator);
                   same rule as the cell2location reference fallback (13_reference_signatures.py). MIN_N 200.
  tme_<compartment>: subcluster labels (17_tme_subclusters.py) within each compartment; subclusters flagged 'small' are
                   pooled into 'other'. MIN_N 30.
Contrasts (reference level first): site thalamus vs pons; age_group adult vs pediatric; h3_variant H3.1 vs H3.3;
  treatment ONC201 + standard care vs standard care without imipridone; tissue_timing autopsy vs resection;
  assay_class nuclei vs cells. Samples with other / unknown levels are dropped for that contrast.
Models, per category, response = empirical logit log((k + 0.5) / (n - k + 0.5)):
  primary  : OLS  y ~ contrast + publication (fixed effects) [+ assay_class when it varies within a publication], cluster-robust
             SE by patient_id. The contrast is estimated only from within-publication variation, so it is reported as
             'not_estimable' when no publication contains both levels (e.g. autopsy vs resection: autopsy = Okada only).
  secondary: MixedLM y ~ contrast [+ assay_class], random intercept per publication (REML); reported with its convergence flag.
  Publications: jessa_scrna and jessa_multiome = Jessa 2022 (shared patients). BH within composition set x contrast.
  scCODA was not used (not installed in the project image); the per-category logit models are the plan's stated alternative.
"""
import sys, argparse, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py, scipy.sparse as sp  # noqa: E402
import scipy.stats as ss  # noqa: E402
import anndata as ad, scanpy as sc  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402
from dmg_atlas.states import STATE_GENE_SETS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
CORE = PROCESSED / "atlas_core.h5ad"; TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
log = U.StepLogger(LOGS / "fix_metaprograms.log", context="18_composition")
SEED = 0
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}
CONTRASTS = {  # name: (column, reference level, test level)
    "site_thalamus_vs_pons": ("site", "pons", "thalamus"),
    "age_adult_vs_pediatric": ("age_group", "pediatric", "adult"),
    "h3_H3.1_vs_H3.3": ("h3_variant", "H3.3", "H3.1"),
    "treatment_ONC201_vs_standard": ("treatment", "standard care, no imipridone", "imipridone (ONC201) + standard care"),
    "timing_autopsy_vs_resection": ("tissue_timing", "resection", "autopsy"),
    "assay_nuclei_vs_cells": ("assay_class", "cells", "nuclei"),
}
OUTS = [TABLES / f for f in ("composition_per_sample.csv", "composition_results.csv", "composition_estimability.csv",
                             "malignant_state_assignment.csv.gz", "composition_results_vs_pre_audit.csv")]
from dmg_atlas.figstyle import CONTRAST_NAMES  # noqa: E402
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

COLS = ["sample_id", "study", "assay", "patient_id", "site", "age_group", "h3_variant", "treatment", "tissue_timing",
        "scanvi_label", "scanvi_maxprob", "consensus_call", "qc_pass"]  # scanvi columns kept for the record; not used for selection
with log.step("read obs"):
    f = h5py.File(CORE, "r")
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in COLS}, index=read_elem(f["obs"]["_index"]))
    for c in COLS:
        if c not in ("scanvi_maxprob", "qc_pass"):
            obs[c] = obs[c].astype(str)
    var = read_elem(f["var"])
    obs["publication"] = obs.study.map(lambda s: PUBLICATION.get(s, s)); obs["assay_class"] = obs.assay.map(ASSAY_CLASS)
    smeta = obs.groupby("sample_id")[["study", "publication", "assay", "assay_class", "patient_id"] + [v[0] for v in CONTRASTS.values() if v[0] != "assay_class"]].first()
    log.write(f"{len(obs)} cells, {len(smeta)} samples, {smeta.publication.nunique()} publications")

with log.step("malignant state assignment (per sample, within-sample z)"):
    selp = pd.read_parquet(PROCESSED / "malignant_selection_audit.parquet", columns=["primary"]).reindex(obs.index)
    sel = selp.primary.fillna(False).astype(bool).values
    rows = np.flatnonzero(sel)
    parts = []
    for s_, idx in pd.Series(rows).groupby(obs.sample_id.values[rows]):
        r_ = np.sort(idx.values)
        X = read_rows(f["layers/counts"], r_, var.shape[0]).astype(np.float32)
        a = ad.AnnData(X=X, obs=obs.iloc[r_][["sample_id"]].copy(), var=var[[]].copy()); del X
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        for st, genes in STATE_GENE_SETS.items():
            sc.tl.score_genes(a, [g for g in genes if g in a.var_names], score_name=f"score_{st}", random_state=SEED)
        parts.append(a.obs[[f"score_{s}" for s in STATE_GENE_SETS]].copy()); del a
    S = pd.concat(parts); S = S.loc[obs.index[rows]]
    raw = S.idxmax(axis=1).str.replace("score_", "", regex=False); raw[S.max(axis=1) <= 0] = "unassigned"
    Z = S.groupby(obs.sample_id.values[rows]).transform(lambda x: (x - x.mean()) / (x.std(ddof=0) if x.std(ddof=0) > 0 else 1.0))
    st = Z.idxmax(axis=1).str.replace("score_", "", regex=False); st[Z.max(axis=1) <= 0] = "unassigned"
    S = S.join(Z.add_prefix("z_")); S["state_raw"] = raw.values; S["state"] = st.values
    S.insert(0, "sample_id", obs.sample_id.values[rows]); S.insert(1, "assay_class", obs.assay_class.values[rows])
    S.to_csv(OUTS[3])
    log.write(f"{len(rows)} selected malignant cells; states (within-sample z) {S.state.value_counts().to_dict()}; raw argmax {S.state_raw.value_counts().to_dict()}\n"
              + "state by assay class (row fractions):\n" + pd.crosstab(S.assay_class, S.state, normalize="index").round(3).to_string())

with log.step("per-sample compositions"):
    long = []
    def add(set_name, lab: pd.Series, min_n):
        tab = pd.crosstab(obs.sample_id.reindex(lab.index).values, lab.values)  # lab is indexed by cell id
        n = tab.sum(1)
        tab = tab[n >= min_n]
        for s, r in tab.iterrows():
            for c, k in r.items():
                long.append(dict(set=set_name, sample_id=s, category=c, k=int(k), n=int(r.sum())))
    qc = obs[obs.qc_pass.astype(bool)]
    add("celltype", pd.Series(qc.scanvi_label.values, index=qc.index), 200)
    add("malignant_state", pd.Series(S.state.values, index=S.index), 200)
    mps = PROCESSED / "metaprogram_scores.parquet"
    if mps.exists():
        fm = pd.read_parquet(mps, columns=["metaprogram_assayz"]).metaprogram_assayz  # within-assay-class z (see 19 docstring)
        add("mp_family", pd.Series(fm.values, index=fm.index), 200)
    else:
        log.write("metaprogram_scores.parquet missing: mp_family set skipped")
    tl = pd.read_csv(PROCESSED / "tme_subcluster_labels.csv.gz", index_col=0)
    tsub = pd.read_csv(TABLES / "tme_subclusters.csv").set_index("subcluster")
    small = set(tsub.index[tsub["flags"].fillna("").str.contains("small")])
    for comp, d in tl.groupby("compartment"):
        cat = np.where(d.subcluster.isin(small), "other", d.subcluster + ":" + d.label.str.replace("|low_confidence", "?", regex=False))
        add(f"tme_{comp}", pd.Series(cat, index=d.index), 30)
    long = pd.DataFrame(long)
    long["proportion"] = long.k / long.n
    long = long.join(smeta, on="sample_id")
    long.to_csv(OUTS[0], index=False)
    log.write("samples per set: " + long.groupby("set").sample_id.nunique().to_dict().__repr__())

with log.step("estimability"):
    est = []
    for name, (col, ref, test) in CONTRASTS.items():
        sm = smeta[smeta[col].isin([ref, test])]
        both = [p for p, d in sm.groupby("publication") if d[col].nunique() == 2]
        est.append(dict(contrast=name, column=col, reference=ref, test=test, n_ref=int((sm[col] == ref).sum()), n_test=int((sm[col] == test).sum()),
                        n_patients_ref=int(sm.patient_id[sm[col] == ref].nunique()), n_patients_test=int(sm.patient_id[sm[col] == test].nunique()),
                        sites_ref=";".join(f"{k}:{v}" for k, v in sm.site[sm[col] == ref].value_counts().items()),
                        sites_test=";".join(f"{k}:{v}" for k, v in sm.site[sm[col] == test].value_counts().items()),
                        assay_ref=";".join(f"{k}:{v}" for k, v in sm.assay_class[sm[col] == ref].value_counts().items()),
                        assay_test=";".join(f"{k}:{v}" for k, v in sm.assay_class[sm[col] == test].value_counts().items()),
                        publications_ref=";".join(sorted(sm.publication[sm[col] == ref].unique())),
                        publications_test=";".join(sorted(sm.publication[sm[col] == test].unique())),
                        publications_with_both=";".join(both), estimable=bool(both)))
    est = pd.DataFrame(est); est.to_csv(OUTS[2], index=False)
    log.write("\n" + est.to_string(index=False))

with log.step("models"):
    res = []
    for (set_name, cat), d in long.groupby(["set", "category"]):
        d = d.copy(); d["y"] = np.log((d.k + 0.5) / (d.n - d.k + 0.5))
        for name, (col, ref, test) in CONTRASTS.items():
            e = est.set_index("contrast").loc[name]
            dd = d[d[col].isin([ref, test])].copy()
            dd["x"] = (dd[col] == test).astype(int)
            row = dict(set=set_name, category=cat, contrast=name, n_samples=len(dd), n_ref=int((dd.x == 0).sum()), n_test=int((dd.x == 1).sum()),
                       median_prop_ref=dd.proportion[dd.x == 0].median(), median_prop_test=dd.proportion[dd.x == 1].median(),
                       n_patients_ref=int(dd.patient_id[dd.x == 0].nunique()), n_patients_test=int(dd.patient_id[dd.x == 1].nunique()),
                       estimable=bool(e.estimable))
            within = [p for p, g in dd.groupby("publication") if g.x.nunique() == 2]
            if not within or dd.x.sum() < 2 or (dd.x == 0).sum() < 2:
                row["estimable"] = False; res.append(row); continue
            assay_term = col != "assay_class" and any(g.assay_class.nunique() == 2 for _, g in dd.groupby("publication"))
            rhs = "x" + (" + C(assay_class)" if assay_term else "")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m = smf.ols(f"y ~ {rhs} + C(publication)", dd).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(dd.patient_id)[0]})
                row.update(coef=m.params["x"], se=m.bse["x"], p=m.pvalues["x"], ci_lo=m.conf_int().loc["x", 0], ci_hi=m.conf_int().loc["x", 1],
                           assay_adjusted=assay_term, publications_within=";".join(within), n_patients=dd.patient_id.nunique())
                if col == "h3_variant":  # MC-06: site covariate and pons-only restriction
                    if dd.site.nunique() >= 2:
                        ms = smf.ols(f"y ~ {rhs} + C(site) + C(publication)", dd).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(dd.patient_id)[0]})
                        row.update(coef_site_adj=ms.params["x"], se_site_adj=ms.bse["x"], p_site_adj=ms.pvalues["x"])
                    dp = dd[dd.site == "pons"]
                    row.update(n_pons_ref=int((dp.x == 0).sum()), n_pons_test=int((dp.x == 1).sum()))
                    if dp.x.sum() >= 2 and (dp.x == 0).sum() >= 2 and any(g.x.nunique() == 2 for _, g in dp.groupby("publication")):
                        rhs_p = "x" + (" + C(assay_class)" if any(g.assay_class.nunique() == 2 for _, g in dp.groupby("publication")) else "")
                        mp_ = smf.ols(f"y ~ {rhs_p} + C(publication)", dp).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(dp.patient_id)[0]})
                        row.update(coef_pons=mp_.params["x"], se_pons=mp_.bse["x"], p_pons=mp_.pvalues["x"])
                try:
                    mm = smf.mixedlm(f"y ~ {rhs}", dd, groups=dd.publication).fit(reml=True)
                    row.update(coef_mixed=mm.params["x"], p_mixed=mm.pvalues["x"], mixed_converged=bool(mm.converged))
                except Exception as ex:  # singular random-effect fits with very few publications
                    row.update(coef_mixed=np.nan, p_mixed=np.nan, mixed_converged=False, mixed_error=type(ex).__name__)
            res.append(row)
    res = pd.DataFrame(res)
    res["q_bh"] = np.nan
    for _, g in res[res.p.notna()].groupby(["set", "contrast"]):
        res.loc[g.index, "q_bh"] = ss.false_discovery_control(g.p.values, method="bh")
    for c in ("p_site_adj", "p_pons"):
        if c in res:
            res[c.replace("p_", "q_")] = np.nan
            for _, g in res[res[c].notna()].groupby(["set", "contrast"]):
                res.loc[g.index, c.replace("p_", "q_")] = ss.false_discovery_control(g[c].values, method="bh")
    res["note"] = np.where((res.set == "mp_family") & (res.contrast == "assay_nuclei_vs_cells"),
                           "not interpretable: family scores standardised within assay class (19_metaprogram_scores.py)", "")
    res.to_csv(OUTS[1], index=False)
    pre = TABLES / "composition_results.pre_audit.csv"
    if pre.exists():
        old = pd.read_csv(pre); old = old[old.set.isin(["celltype", "malignant_state"])]; key = ["set", "category", "contrast"]
        cmp_ = res[res.set.isin(["celltype", "malignant_state"])][key + ["n_ref", "n_test", "coef", "q_bh"]].merge(
            old[key + ["n_ref", "n_test", "coef", "q_bh"]], on=key, how="outer", suffixes=("_new", "_pre_audit"))
        cmp_["sig_pre_audit"] = cmp_.q_bh_pre_audit < 0.05; cmp_["sig_new"] = cmp_.q_bh_new < 0.05
        cmp_.to_csv(OUTS[4], index=False)
        log.write(f"vs pre-audit: q<0.05 pre {int(cmp_.sig_pre_audit.sum())}, new {int(cmp_.sig_new.sum())}, lost "
                  f"{int((cmp_.sig_pre_audit & ~cmp_.sig_new).sum())}, gained {int((~cmp_.sig_pre_audit & cmp_.sig_new).sum())}")
    h3 = res[(res.contrast == "h3_H3.1_vs_H3.3") & res.p.notna() & res.set.isin(["celltype", "malignant_state", "mp_family"])]
    log.write("H3.1 vs H3.3 (primary / site-adjusted / pons-only):\n" + h3[["set", "category", "n_ref", "n_test", "n_patients_ref", "n_patients_test",
              "coef", "q_bh", "coef_site_adj", "q_site_adj", "coef_pons", "q_pons", "n_pons_ref", "n_pons_test"]].round(3).to_string(index=False))
    top = res[res.q_bh < 0.1].sort_values("q_bh")
    log.write(f"{len(res)} tests, {int(res.p.notna().sum())} estimable, {len(top)} with q < 0.1:\n" +
              top[["set", "category", "contrast", "n_ref", "n_test", "median_prop_ref", "median_prop_test", "coef", "ci_lo", "ci_hi", "q_bh",
                   "coef_mixed", "p_mixed"]].round(3).to_string(index=False))

with log.step("figure"):
    apply_figure_style()
    sets = ["celltype", "malignant_state", "mp_family"]
    d = res[res.set.isin(sets) & res.p.notna()].copy()
    d["row"] = d.set.map({"celltype": "cell type", "malignant_state": "state", "mp_family": "MP family"}) + ": " + d.category.str.replace("_", "-")
    M = d.pivot_table(index="row", columns="contrast", values="coef"); Q = d.pivot_table(index="row", columns="contrast", values="q_bh")
    M = M[[c for c in CONTRAST_NAMES if c in M.columns]]; Q = Q[M.columns]
    if "assay_nuclei_vs_cells" in M.columns:  # by construction uninterpretable (within-assay-class standardisation)
        M.loc[M.index.str.startswith("MP family"), "assay_nuclei_vs_cells"] = np.nan; Q.loc[Q.index.str.startswith("MP family"), "assay_nuclei_vs_cells"] = np.nan
    lim = float(np.ceil(np.nanmax(np.abs(M.values)) * 2) / 2)
    fig, ax = plt.subplots(figsize=(3.8, 0.17 * len(M) + 1.3))
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            q = Q.iat[i, j]
            if pd.notna(q) and q < 0.1:
                ax.text(j, i, "**" if q < 0.01 else ("*" if q < 0.05 else "·"), ha="center", va="center", fontsize=6)
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels([CONTRAST_NAMES.get(c, c) for c in M.columns], rotation=45, ha="right", fontsize=5)
    ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=5)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03); cb.set_label("Δ logit proportion (within-publication)")
    fig.savefig(FIGS / "phase6_composition_effects.png"); plt.close(fig)
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
