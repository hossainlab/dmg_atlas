#!/usr/bin/env python
"""Phase 9 quality assessment / step 34 — sensitivity of the cell-type composition results to the labels of ambiguous cells.

18_composition.py builds the 'celltype' composition set from scanvi_label on every qc_pass cell, so the 29 % of cells whose
corrected consensus call is ambiguous (lineage_label 'Unknown') are placed in a lineage by the classifier alone. (The
metaprogram / malignant-state track already uses consensus calls only, 13a_malignant_selection.py.) Here the celltype set is
rebuilt under four labelings and refitted with the 18_composition.py primary model, unchanged:
  scanvi            : scanvi_label, all qc_pass cells (= 18_composition.py; reproduction asserted against composition_per_sample.csv)
  scanvi_confident  : scanvi_label, qc_pass cells with scanvi_maxprob >= 0.9 only
  consensus_only    : corrected lineage_label (reference_labels_corrected.csv.gz), ambiguous ('Unknown') cells dropped
  consensus_ambig   : corrected lineage_label with 'ambiguous' kept as its own category in the denominator
Model (18_composition.py primary): per category, y = log((k + 0.5) / (n - k + 0.5)); OLS y ~ x [+ C(assay_class) when it varies
within a publication] + C(publication), cluster-robust SE by patient_id; samples with n >= 200; BH within labeling x contrast.
Inputs : data/processed/atlas_core.h5ad (obs only), results/tables/reflabels/reference_labels_corrected.csv.gz,
         results/tables/composition_per_sample.csv, results/tables/composition_results.csv
Outputs (results/tables/quality/, skipped if present unless --force):
         composition_label_sensitivity_per_sample.csv  labeling x sample x category: k, n, proportion
         composition_label_sensitivity_results.csv     labeling x category x contrast: coef, se, p, q_bh, n_samples
         composition_label_sensitivity_compare.csv     per category x contrast: coef / q under every labeling, sign agreement,
                                                       q < 0.05 under scanvi and under consensus_only
         composition_label_sensitivity_inputs.json     md5 of every input
Run: nice -n 10 python scripts/34_composition_label_sensitivity.py [--force]
Log: logs/phase9_quality.log
"""
import sys, json, hashlib, argparse, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py  # noqa: E402
import scipy.stats as ss  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
log = U.StepLogger(LOGS / "phase9_quality.log", context="34_composition_label_sensitivity")
CORE = PROCESSED / "atlas_core.h5ad"; TABLES = RESULTS / "tables"; OUT = TABLES / "quality"; OUT.mkdir(parents=True, exist_ok=True)
LABELS = TABLES / "reflabels/reference_labels_corrected.csv.gz"
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}
CONTRASTS = {  # = 18_composition.py
    "site_thalamus_vs_pons": ("site", "pons", "thalamus"),
    "age_adult_vs_pediatric": ("age_group", "pediatric", "adult"),
    "h3_H3.1_vs_H3.3": ("h3_variant", "H3.3", "H3.1"),
    "treatment_ONC201_vs_standard": ("treatment", "standard care, no imipridone", "imipridone (ONC201) + standard care"),
    "timing_autopsy_vs_resection": ("tissue_timing", "resection", "autopsy"),
    "assay_nuclei_vs_cells": ("assay_class", "cells", "nuclei"),
}
MIN_N, MAXP = 200, 0.9
OUTS = [OUT / f"composition_label_sensitivity_{n}" for n in ("per_sample.csv", "results.csv", "compare.csv", "inputs.json")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fit(d, col, ref, test):
    """18_composition.py primary model for one category x contrast; None when not estimable."""
    dd = d[d[col].isin([ref, test])].copy(); dd["x"] = (dd[col] == test).astype(int)
    if not any(g.x.nunique() == 2 for _, g in dd.groupby("publication")) or dd.x.sum() < 2 or (dd.x == 0).sum() < 2:
        return None
    rhs = "x" + (" + C(assay_class)" if col != "assay_class" and any(g.assay_class.nunique() == 2 for _, g in dd.groupby("publication")) else "")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = smf.ols(f"y ~ {rhs} + C(publication)", dd).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(dd.patient_id)[0]})
    return dict(n_samples=len(dd), n_ref=int((dd.x == 0).sum()), n_test=int((dd.x == 1).sum()), coef=m.params["x"], se=m.bse["x"], p=m.pvalues["x"])


with log.step("read obs + corrected labels"):
    cols = ["sample_id", "study", "assay", "patient_id", "site", "age_group", "h3_variant", "treatment", "tissue_timing",
            "scanvi_label", "scanvi_maxprob", "qc_pass"]
    with h5py.File(CORE, "r") as f:
        obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in cols}, index=read_elem(f["obs"]["_index"]))
    for c in cols:
        if c not in ("scanvi_maxprob", "qc_pass"):
            obs[c] = obs[c].astype(str)
    lab = pd.read_csv(LABELS, index_col="cell_id", usecols=["cell_id", "lineage_label"]).lineage_label.astype(str)
    assert obs.index.isin(lab.index).all(), "atlas cells without a corrected label"
    obs["lineage"] = lab.reindex(obs.index).values
    obs["publication"] = obs.study.map(lambda s: PUBLICATION.get(s, s)); obs["assay_class"] = obs.assay.map(ASSAY_CLASS)
    smeta = obs.groupby("sample_id")[["study", "publication", "assay", "assay_class", "patient_id"] + [v[0] for v in CONTRASTS.values() if v[0] != "assay_class"]].first()
    qc = obs[obs.qc_pass.astype(bool)]
    log.write(f"{len(obs)} cells, {len(qc)} qc_pass; corrected lineage among qc_pass {qc.lineage.value_counts().to_dict()}")

with log.step("per-sample compositions"):
    labelings = {"scanvi": qc.scanvi_label, "scanvi_confident": qc.scanvi_label[qc.scanvi_maxprob >= MAXP],
                 "consensus_only": qc.lineage[qc.lineage != "Unknown"], "consensus_ambig": qc.lineage.replace("Unknown", "ambiguous")}
    long = []
    for name, l in labelings.items():
        tab = pd.crosstab(qc.sample_id.reindex(l.index).values, l.values)
        tab = tab[tab.sum(1) >= MIN_N]
        t = tab.stack().rename("k").reset_index(); t.columns = ["sample_id", "category", "k"]
        t["n"] = t.sample_id.map(tab.sum(1)); t.insert(0, "labeling", name); long.append(t)
    long = pd.concat(long, ignore_index=True); long["proportion"] = long.k / long.n
    ref = pd.read_csv(TABLES / "composition_per_sample.csv"); ref = ref[ref.set == "celltype"]
    chk = long[long.labeling == "scanvi"].merge(ref[["sample_id", "category", "k", "n"]], on=["sample_id", "category"], how="outer", suffixes=("", "_18"))
    assert chk.k.notna().all() and chk.k_18.notna().all() and (chk.k == chk.k_18).all() and (chk.n == chk.n_18).all(), \
        "scanvi labeling does not reproduce composition_per_sample.csv (celltype)"
    long = long.join(smeta, on="sample_id")
    long[["labeling", "sample_id", "category", "k", "n", "proportion", "publication", "assay_class", "patient_id"]].to_csv(OUTS[0], index=False)
    mf = long[long.category == "malignant"].pivot(index="sample_id", columns="labeling", values="proportion")
    log.write("samples per labeling " + long.groupby("labeling").sample_id.nunique().to_dict().__repr__()
              + "\nper-sample malignant fraction, Spearman vs scanvi: " + mf.corr("spearman")["scanvi"].round(3).to_dict().__repr__()
              + "\nmedian malignant fraction: " + mf.median().round(3).to_dict().__repr__())

with log.step("models"):
    res = []
    for (name, cat), d in long.groupby(["labeling", "category"]):
        d = d.copy(); d["y"] = np.log((d.k + 0.5) / (d.n - d.k + 0.5))
        for cn, (col, r, t) in CONTRASTS.items():
            row = fit(d, col, r, t)
            res.append(dict(labeling=name, category=cat, contrast=cn, estimable=row is not None, **(row or {})))
    res = pd.DataFrame(res); res["q_bh"] = np.nan
    for _, g in res[res.p.notna()].groupby(["labeling", "contrast"]):
        res.loc[g.index, "q_bh"] = ss.false_discovery_control(g.p.values, method="bh")
    r18 = pd.read_csv(TABLES / "composition_results.csv"); r18 = r18[(r18.set == "celltype") & r18.coef.notna()]
    chk = res[(res.labeling == "scanvi") & res.coef.notna()].merge(r18[["category", "contrast", "coef", "q_bh"]], on=["category", "contrast"], suffixes=("", "_18"))
    assert len(chk) == len(r18) and np.allclose(chk.coef, chk.coef_18) and np.allclose(chk.q_bh, chk.q_bh_18), "scanvi model does not reproduce composition_results.csv"
    res.to_csv(OUTS[1], index=False)

with log.step("compare + write"):
    est = res[res.estimable]
    w = est.pivot_table(index=["category", "contrast"], columns="labeling", values=["coef", "q_bh"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]; w = w.reset_index()
    alts = [k for k in labelings if k != "scanvi"]
    for a in alts:
        w[f"same_sign_{a}"] = np.sign(w[f"coef_{a}"]) == np.sign(w["coef_scanvi"])
    w["sig_scanvi"] = w.q_bh_scanvi < 0.05
    w["sig_consensus_only"] = w.q_bh_consensus_only < 0.05
    w.to_csv(OUTS[2], index=False)
    json.dump(dict(script="scripts/34_composition_label_sensitivity.py", params=dict(min_n=MIN_N, maxprob=MAXP),
                   inputs={str(p.relative_to(RESULTS.parent)): md5(p) for p in [LABELS, TABLES / "composition_per_sample.csv", TABLES / "composition_results.csv"]}
                   | {"data/processed/atlas_core.h5ad": "obs only; md5 not recomputed (5.3 GB); see results/tables/release_manifest.csv"}),
              open(OUTS[3], "w"), indent=1)
    sig = w[w.sig_scanvi]
    log.write(f"{len(w)} estimable category x contrast cells; {len(sig)} with q < 0.05 under scanvi labels; of these, same sign and q < 0.05 "
              f"under consensus_only: {int((sig.same_sign_consensus_only & sig.sig_consensus_only).sum())}; same sign under consensus_only: "
              f"{int(sig.same_sign_consensus_only.sum())}; new q < 0.05 under consensus_only only: {int((w.sig_consensus_only & ~w.sig_scanvi).sum())}\n"
              + w[w.sig_scanvi | w.sig_consensus_only].round(4).to_string(index=False))
