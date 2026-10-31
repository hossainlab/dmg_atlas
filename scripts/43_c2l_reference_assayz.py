#!/usr/bin/env python
"""Phase 9 science fix / step 43 — cell2location reference with malignant classes from the within-assay-class assignment.

40_family_assignment_check.py showed that metaprogram families assigned by scores standardised WITHIN SAMPLE (the classes of
the published reference, 13b_c2l_reference_mpfam.py) do not track whether a sample's own factorisation found the programme,
whereas families standardised WITHIN ASSAY CLASS do for 5 of 7 families. This script rebuilds the reference exactly as 13b
(same normal classes, exclusions, cap of 60,000 cells and <= 8,000 per class stratified by sample, seed 0, same gene set)
but takes malignant classes from metaprogram_family_reference_classes.csv[family_assayz] with max_family_z_assayz >= 1.
The gene set does not depend on the classes, so the existing Visium inputs (modal/inputs/visium_<cohort>.h5ad) remain valid.
Outputs: modal/inputs/c2l_reference_assayz.h5ad, results/tables/quality/c2l_reference_classes_assayz.csv
Run: nice -n 10 python scripts/43_c2l_reference_assayz.py
Next: modal/phase6_c2l.py signatures / deconvolve on Modal, then scripts/44_spatial_assayz.py.
"""
from __future__ import annotations
import sys, json, importlib.util
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402
from dmg_atlas.states import all_marker_genes  # noqa: E402
import h5py, numpy as np, pandas as pd, anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402

spec = importlib.util.spec_from_file_location("ref13b", ROOT / "scripts/13b_c2l_reference_mpfam.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)  # reuse stratified(), read_rows() and constants
LOG = StepLogger(ROOT / "logs/phase9_quality.log", context="43_c2l_reference_assayz")
OUT_H5AD = ROOT / "modal/inputs/c2l_reference_assayz.h5ad"; OUT_TAB = ROOT / "results/tables/quality/c2l_reference_classes_assayz.csv"


def main():
    rng = np.random.default_rng(B.SEED)
    f = h5py.File(B.ATLAS, "r"); obs = read_elem(f["obs"]); var = read_elem(f["var"])
    usable = (~obs.sample_id.astype(str).isin(B.EXCLUDE_SAMPLES) & ~obs.study.astype(str).str.lower().str.startswith(B.EXCLUDE_STUDY_PREFIX)).values
    cls = pd.Series(pd.NA, index=obs.index, dtype="object")
    agree = obs.lineage_label.astype(str) == obs.scanvi_label.astype(str)
    for lin in B.NORMAL:
        m = (obs.scanvi_label.astype(str) == lin) & (agree if lin != "T_NK" else (obs.lineage_label.astype(str) == "T_NK"))
        cls[m.values & usable] = lin
    mp = pd.read_csv(B.MPFAM)
    mpc = mp[(mp.max_family_z_assayz >= 1) & mp.family_assayz.notna() & (mp.family_assayz != "unassigned")].set_index("cell_id")
    fam = mpc["family_assayz"].reindex(obs.index)
    mal = fam.notna().values & usable
    cls[mal] = "malignant_" + fam[mal].astype(str)
    LOG.write("class sizes (all): " + json.dumps(cls.value_counts().to_dict()))
    full = cls.dropna(); sizes = full.value_counts(); cap = B.CLASS_CAP
    while np.minimum(sizes, cap).sum() > B.TOTAL_CAP:
        cap -= 100
    pos = pd.Series(np.arange(obs.shape[0]), index=obs.index); keep = []
    for c in sizes.index:
        idx = np.sort(pos[full.index[full == c]].values)
        keep.append(B.stratified(idx, obs.sample_id.astype(str).values[idx], min(cap, len(idx)), rng, floor=20))
    keep = np.sort(np.concatenate(keep)); labels = full.reindex(obs.index[keep]).values
    T = ROOT / "results/tables"
    hvg = pd.read_csv(T / "atlas_input_hvg.csv"); hvg = set(hvg.symbol[hvg.highly_variable])
    mhvg = pd.read_csv(T / "atlas_malignant_hvg.csv"); mhvg = set(mhvg.symbol[mhvg.highly_variable])
    genes = sorted((hvg | mhvg | set(all_marker_genes())) & set(var.index)); gmask = var.index.isin(genes)
    published = ad.read_h5ad(ROOT / "modal/inputs/c2l_reference.h5ad", backed="r")
    assert list(published.var_names) == list(var.index[gmask]), "gene set differs from the published reference"
    with LOG.step(f"read {len(keep)} reference rows"):
        X = B.read_rows(f, keep)[:, np.where(gmask)[0]].tocsr(); X.sort_indices()
    o = obs.iloc[keep][["sample_id", "patient_id", "study", "assay", "pct_mito", "lineage_label", "scanvi_label"]].copy()
    for c in o.columns:
        if c != "pct_mito":
            o[c] = o[c].astype(str)
    o["pct_mito"] = o["pct_mito"].astype(np.float32); o["c2l_class"] = labels.astype(str)
    ref = ad.AnnData(X=X.astype(np.int32), obs=o, var=pd.DataFrame(index=var.index[gmask]))
    ref.uns["phase9"] = {"state_source": "metaprogram_family_reference_classes.csv[family_assayz], max_family_z_assayz >= 1",
                         "class_cap": int(cap), "seed": B.SEED, "n_cells": int(ref.n_obs)}
    ref.write_h5ad(OUT_H5AD, compression="gzip", compression_opts=4)
    tab = o.groupby("c2l_class").agg(n_cells=("sample_id", "size"), n_samples=("sample_id", "nunique"),
                                     n_patients=("patient_id", "nunique"), n_studies=("study", "nunique")).reset_index()
    tab.to_csv(OUT_TAB, index=False)
    LOG.write(f"wrote {OUT_H5AD} {ref.shape}; class cap {cap}\n" + tab.to_string(index=False) + f"\npeak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
