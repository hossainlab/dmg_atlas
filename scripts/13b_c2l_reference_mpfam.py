#!/usr/bin/env python
"""Phase 6 / step 13b (audit remediation 2026-10-20, findings SP-06, SP-10) — cell2location reference built from
the CORRECTED metaprogram families instead of the score_genes-argmax fallback.

Malignant classes : results/tables/metaprogram_family_reference_classes.csv, confident cells only, column `family`
                    (within-sample z assignment; the within-assay-class variant `family_assayz` is written out as a
                    sensitivity class table but is not used for the reference).
Normal classes    : refreshed scanvi_label (atlas_core.h5ad, wave-1/2 retrain) in agreement with the consensus
                    lineage_label; T_NK from lineage_label alone (scANVI T_NK recall 0.467, Phase 4 notes §5).
Exclusions        : sample okada__1734_pons (low_confidence_sample, Phase 4 wave-2 decision) and any ScPCA study.
Subsample         : <= 60,000 cells, <= 8,000 per class, stratified by sample_id within class (seed 0).
Genes             : atlas HVG union malignant HVG union curated markers (dmg_atlas.states).
Outputs           : modal/inputs/c2l_reference.h5ad  (old -> .pre_audit.h5ad)
                    results/tables/c2l_reference_classes.csv          (old -> .pre_audit.csv)
                    results/tables/c2l_reference_class_composition.csv (class x study / class x assay cell counts)
"""
from __future__ import annotations

import sys, json, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402
from dmg_atlas.states import all_marker_genes  # noqa: E402

import h5py, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_spatial.log", context="13b_c2l_reference_mpfam")
ATLAS = ROOT / "data/processed/atlas_core.h5ad"
MPFAM = ROOT / "results/tables/metaprogram_family_reference_classes.csv"
MODAL_IN = ROOT / "modal/inputs"; TABLES = ROOT / "results/tables"
NORMAL = ["myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular"]
EXCLUDE_SAMPLES = {"okada__1734_pons"}
EXCLUDE_STUDY_PREFIX = ("scpca",)
SEED, TOTAL_CAP, CLASS_CAP = 0, 60_000, 8_000


def stratified(idx, strata, n, rng, floor=20):
    if len(idx) <= n:
        return idx
    s = pd.Series(strata, index=np.arange(len(idx)))
    sizes = s.value_counts()
    alloc = (sizes / sizes.sum() * n).round().astype(int)
    alloc = np.minimum(alloc.clip(lower=np.minimum(floor, sizes)), sizes)
    while alloc.sum() > n:
        alloc[alloc.idxmax()] -= 1
    out = []
    for k, m in alloc.items():
        pool = s.index[s == k].values
        out.append(rng.choice(pool, size=int(m), replace=False) if m < len(pool) else pool)
    return np.sort(idx[np.concatenate(out)])


def read_rows(f, rows, block=20_000):
    g = f["layers/counts"]; indptr = g["indptr"][:]; n_genes = int(g.attrs["shape"][1]); parts = []
    for b0 in range(0, len(indptr) - 1, block):
        r = rows[(rows >= b0) & (rows < b0 + block)]
        if not len(r):
            continue
        b1 = min(b0 + block, len(indptr) - 1); s, e = int(indptr[b0]), int(indptr[b1])
        m = sp.csr_matrix((g["data"][s:e], g["indices"][s:e], indptr[b0:b1 + 1] - s), shape=(b1 - b0, n_genes))
        parts.append(m[r - b0])
    return sp.vstack(parts).tocsr()


def main():
    rng = np.random.default_rng(SEED)
    f = h5py.File(ATLAS, "r")
    obs = read_elem(f["obs"]); var = read_elem(f["var"])
    LOG.write(f"atlas_core: {obs.shape[0]} cells x {var.shape[0]} genes")

    keep_sample = ~obs.sample_id.astype(str).isin(EXCLUDE_SAMPLES)
    keep_study = ~obs.study.astype(str).str.lower().str.startswith(EXCLUDE_STUDY_PREFIX)
    LOG.write(f"exclusions: okada__1734_pons {int((~keep_sample).sum())} cells; ScPCA {int((~keep_study).sum())} cells "
              f"(no ScPCA study present in atlas_core)")
    usable = (keep_sample & keep_study).values

    # normal classes from the refreshed scanvi_label
    cls = pd.Series(pd.NA, index=obs.index, dtype="object")
    agree = obs.lineage_label.astype(str) == obs.scanvi_label.astype(str)
    for lin in NORMAL:
        m = (obs.scanvi_label.astype(str) == lin) & (agree if lin != "T_NK" else (obs.lineage_label.astype(str) == "T_NK"))
        cls[m.values & usable] = lin
    LOG.write("normal class sizes (post-exclusion): " + json.dumps(cls.value_counts().to_dict()))

    # malignant classes from the corrected metaprogram families (confident cells only)
    mp = pd.read_csv(MPFAM)
    LOG.write(f"metaprogram_family_reference_classes.csv: {mp.shape}, confident {int(mp.confident.sum())}")
    mpc = mp[mp.confident.astype(bool)].set_index("cell_id")
    fam = mpc["family"].reindex(obs.index)
    mal = fam.notna().values & usable
    cls[mal] = "malignant_" + fam[mal].astype(str)
    LOG.write("class sizes (all): " + json.dumps(cls.value_counts().to_dict()))

    full = cls.dropna()
    sizes = full.value_counts()
    cap = CLASS_CAP
    while np.minimum(sizes, cap).sum() > TOTAL_CAP:
        cap -= 100
    LOG.write(f"per-class cap {cap} (requested <= {CLASS_CAP}); total {int(np.minimum(sizes, cap).sum())}")
    pos = pd.Series(np.arange(obs.shape[0]), index=obs.index)
    keep = []
    for c in sizes.index:
        idx = np.sort(pos[full.index[full == c]].values)
        keep.append(stratified(idx, obs.sample_id.astype(str).values[idx], min(cap, len(idx)), rng, floor=20))
    keep = np.sort(np.concatenate(keep))
    labels = full.reindex(obs.index[keep]).values

    hvg = pd.read_csv(TABLES / "atlas_input_hvg.csv"); hvg = set(hvg.symbol[hvg.highly_variable])
    mhvg = pd.read_csv(TABLES / "atlas_malignant_hvg.csv"); mhvg = set(mhvg.symbol[mhvg.highly_variable])
    markers = set(all_marker_genes())
    genes = sorted((hvg | mhvg | markers) & set(var.index))
    gmask = var.index.isin(genes)
    LOG.write(f"gene set: HVG {len(hvg)} u malignant HVG {len(mhvg)} u markers {len(markers)} -> {len(genes)} genes")

    with LOG.step(f"read {len(keep)} reference rows"):
        X = read_rows(f, keep)[:, np.where(gmask)[0]].tocsr(); X.sort_indices()
    o = obs.iloc[keep][["sample_id", "patient_id", "study", "assay", "pct_mito", "lineage_label", "scanvi_label"]].copy()
    for c in o.columns:
        if c != "pct_mito":
            o[c] = o[c].astype(str)
    o["pct_mito"] = o["pct_mito"].astype(np.float32)
    o["c2l_class"] = labels.astype(str)
    ref = ad.AnnData(X=X.astype(np.int32), obs=o, var=pd.DataFrame(index=var.index[gmask]))
    ref.var["in_reference_hvg"] = ref.var_names.isin(hvg)
    ref.var["in_malignant_hvg"] = ref.var_names.isin(mhvg)
    ref.var["is_marker"] = ref.var_names.isin(markers)
    ref.uns["phase6"] = {"state_source": "metaprogram_family_reference_classes.csv[family], confident==True",
                         "normal_source": "atlas_core scanvi_label (wave-2 retrain) & lineage_label agreement; T_NK from lineage_label",
                         "excluded_samples": sorted(EXCLUDE_SAMPLES), "seed": SEED, "class_cap": int(cap),
                         "n_cells": int(ref.n_obs), "audit": "2026-10-20 SP-06/SP-10"}
    out = MODAL_IN / "c2l_reference.h5ad"
    if out.exists() and not (MODAL_IN / "c2l_reference.pre_audit.h5ad").exists():
        shutil.move(str(out), str(MODAL_IN / "c2l_reference.pre_audit.h5ad"))
    ref.write_h5ad(out, compression="gzip", compression_opts=4)
    LOG.write(f"wrote {out} ({out.stat().st_size/1e6:.1f} MB, {ref.shape}, nnz {X.nnz})")

    tab = o.groupby("c2l_class").agg(n_cells=("sample_id", "size"), n_samples=("sample_id", "nunique"),
                                     n_patients=("patient_id", "nunique"), n_studies=("study", "nunique")).reset_index()
    tab["state_source"] = "metaprogram families (confident)"
    p = TABLES / "c2l_reference_classes.csv"
    if p.exists() and not (TABLES / "c2l_reference_classes.pre_audit.csv").exists():
        shutil.copy(p, TABLES / "c2l_reference_classes.pre_audit.csv")
    tab.to_csv(p, index=False)
    comp = (pd.crosstab(o.c2l_class, o.study).join(pd.crosstab(o.c2l_class, o.assay), rsuffix="_assay")
            .rename_axis("c2l_class").reset_index())
    comp.to_csv(TABLES / "c2l_reference_class_composition.csv", index=False)
    LOG.write("reference classes:\n" + tab.to_string(index=False))
    LOG.write("composition:\n" + comp.to_string(index=False))
    LOG.write(f"done | peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
