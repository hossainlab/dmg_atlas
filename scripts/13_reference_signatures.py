#!/usr/bin/env python
"""Phase 6 / step 13 — cell2location reference signatures from atlas_core.h5ad.

Local part (`prepare`): build the per-cell class label, draw a stratified reference subsample (<= 60,000 cells,
<= CAP per class, stratified by sample_id within class), restrict to HVG + marker genes and write a compressed h5ad
(modal/inputs/c2l_reference.h5ad) that is shipped to Modal as a job input.
Remote part: modal/phase6_c2l.py signatures  (RegressionModel, 250 epochs, A10G) -> out/signatures.csv.
Local part (`collect`): copy the harvested signature matrix to results/tables/c2l_signatures.csv and write the class table.
Options: --use-metaprograms takes the malignant classes from data/processed/metaprogram_scores.parquet (19_metaprogram_scores.py,
column 'metaprogram' = family; 'unassigned' cells and cells whose winning z-scored family score is < 1 are dropped) instead of the fallback; --tag <t> suffixes every output
(modal/inputs/c2l_reference_<t>.h5ad, results/tables/c2l_reference_classes_<t>.csv, results/tables/c2l_signatures_<t>.csv) so the
original fallback outputs are left untouched. The first (fallback) run used neither option.

Class definition
  normal lineages  : cells whose consensus `lineage_label` equals the scANVI `scanvi_label`
                     (myeloid, oligodendrocyte, astrocyte, neuron, vascular); T_NK uses `lineage_label == T_NK`
                     alone because the scANVI classifier's T_NK recall is 0.467 (Phase 4 notes §5).
  malignant states : consensus-malignant AND scANVI-malignant cells, split by state.
                     With --use-metaprograms the metaprogram-family assignment is used; otherwise (fallback, recorded in the output) scanpy.score_genes with the gene sets in
                     dmg_atlas.states (Filbin 2018 OPC-/AC-/OC-like + cycling, Liu 2022 / Neftel 2019 MES-like), argmax.
"""
from __future__ import annotations

import sys
import json
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402  (sets thread caps)
from dmg_atlas.states import STATE_GENE_SETS, NORMAL_MARKERS, all_marker_genes  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import anndata as ad  # noqa: E402
import scanpy as sc  # noqa: E402
from anndata.io import read_elem  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase6_spatial.log", context="13_reference_signatures")
ATLAS = ROOT / "data/processed/atlas_core.h5ad"
METAPROG = ROOT / "data/processed/metaprogram_scores.parquet"
MODAL_IN = ROOT / "modal/inputs"
TABLES = ROOT / "results/tables"
NORMAL = ["myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular"]
STATES = list(STATE_GENE_SETS)
SEED = 0
TOTAL_CAP = 60_000
CLASS_CAP = 8_000
N_MALIGNANT_CANDIDATES = 90_000
MP_MIN_Z = 1.0   # --use-metaprograms only; set after seeing the z-argmax assignment of 19_metaprogram_scores.py spread cells thinly


def stratified(idx: np.ndarray, strata: np.ndarray, n: int, rng, floor: int = 20) -> np.ndarray:
    """Draw n of idx, proportional to stratum size with a floor per stratum (all cells if stratum <= floor)."""
    if len(idx) <= n:
        return idx
    s = pd.Series(strata, index=np.arange(len(idx)))
    sizes = s.value_counts()
    alloc = (sizes / sizes.sum() * n).round().astype(int)
    alloc = alloc.clip(lower=np.minimum(floor, sizes))
    alloc = np.minimum(alloc, sizes)
    # trim/boost proportionally to hit n
    while alloc.sum() > n:
        big = alloc.idxmax(); alloc[big] -= 1
    out = []
    for k, m in alloc.items():
        pool = s.index[s == k].values
        out.append(rng.choice(pool, size=int(m), replace=False) if m < len(pool) else pool)
    return np.sort(idx[np.concatenate(out)])


def read_rows(f: h5py.File, rows: np.ndarray, block: int = 20_000) -> sp.csr_matrix:
    """Row-subset of layers['counts'] read from disk in contiguous row blocks via h5py (rows must be sorted). Replaces the
    anndata sparse_dataset row indexing, which fails with the scipy/anndata pair in the Modal image; same result."""
    g = f["layers/counts"]; indptr = g["indptr"][:]; n_genes = int(g.attrs["shape"][1]); parts = []
    for b0 in range(0, len(indptr) - 1, block):
        r = rows[(rows >= b0) & (rows < b0 + block)]
        if not len(r):
            continue
        b1 = min(b0 + block, len(indptr) - 1); s, e = int(indptr[b0]), int(indptr[b1])
        m = sp.csr_matrix((g["data"][s:e], g["indices"][s:e], indptr[b0:b1 + 1] - s), shape=(b1 - b0, n_genes))
        parts.append(m[r - b0])
    return sp.vstack(parts).tocsr()


def prepare(use_mp: bool = False, tag: str = ""):
    rng = np.random.default_rng(SEED)
    with LOG.step("read obs/var"):
        f = h5py.File(ATLAS, "r")
        obs = read_elem(f["obs"]); var = read_elem(f["var"])
        n_cells = obs.shape[0]
        LOG.write(f"atlas_core: {n_cells} cells x {var.shape[0]} genes")

    # ---------------- normal classes
    cls = pd.Series(pd.NA, index=obs.index, dtype="object")
    agree = obs.lineage_label.astype(str) == obs.scanvi_label.astype(str)
    for lin in NORMAL:
        m = (obs.lineage_label.astype(str) == lin) & (agree if lin != "T_NK" else True)
        cls[m] = lin
    malig = (obs.lineage_label.astype(str) == "malignant") & (obs.scanvi_label.astype(str) == "malignant")
    LOG.write("normal class sizes: " + json.dumps(cls.value_counts().to_dict()) + f"; malignant candidates {int(malig.sum())}")

    # ---------------- malignant states
    state_source = None
    mal_idx_all = np.where(malig.values)[0]
    if use_mp:
        mp = pd.read_parquet(METAPROG)
        LOG.write(f"metaprogram_scores.parquet found: {mp.shape}, columns {list(mp.columns)[:12]}")
        # expected: index = cell id, a column naming the assigned metaprogram
        col = next((c for c in ("state", "metaprogram", "assignment", "mp_assignment", "program") if c in mp.columns), None)
        if col is None:
            raise SystemExit("metaprogram_scores.parquet has no assignment column; cannot use it")
        assign = mp[col].reindex(obs.index[mal_idx_all]).astype("object")
        assign[assign == "unassigned"] = pd.NA
        zc = [c for c in mp.columns if c.startswith("z_")]
        if zc:  # reference purity: only cells whose winning z-scored family score is >= MP_MIN_Z (typical members of the family)
            zmax = mp[zc].max(axis=1).reindex(obs.index[mal_idx_all])
            assign[(zmax < MP_MIN_Z).values] = pd.NA
            LOG.write(f"metaprogram reference: kept cells with max family z >= {MP_MIN_Z}: {int(assign.notna().sum())} of {len(assign)}")
        state_source = f"metaprogram_scores.parquet[{col}]"
        cand_idx = mal_idx_all
    else:
        state_source = "fallback: scanpy.score_genes with dmg_atlas.states gene sets (Filbin 2018 / Liu 2022), argmax"
        LOG.write("metaprogram_scores.parquet NOT found -> " + state_source)
        cand_idx = stratified(mal_idx_all, obs.sample_id.astype(str).values[mal_idx_all], N_MALIGNANT_CANDIDATES, rng, floor=50)
        with LOG.step(f"read {len(cand_idx)} malignant candidate rows + score states"):
            X = read_rows(f, cand_idx)
            a = ad.AnnData(X=X.astype(np.float32), obs=obs.iloc[cand_idx][["sample_id", "study", "assay"]].copy(), var=var[[]].copy())
            sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
            scores = {}
            for st, genes in STATE_GENE_SETS.items():
                g = [x for x in genes if x in a.var_names]
                sc.tl.score_genes(a, g, score_name=f"score_{st}", random_state=SEED)
                scores[st] = a.obs[f"score_{st}"].values
                LOG.write(f"  {st}: {len(g)}/{len(set(genes))} genes present")
            S = pd.DataFrame(scores, index=a.obs_names)
            assign = S.idxmax(axis=1).astype("object")
            assign[S.max(axis=1) <= 0] = pd.NA
            S["state"] = assign.values
            S.to_csv(TABLES / "c2l_reference_state_scores.csv.gz")
            del a, X
    counts = assign.value_counts(dropna=False)
    LOG.write("malignant state assignment: " + json.dumps({str(k): int(v) for k, v in counts.items()}))

    # ---------------- subsample per class
    cls_mal = pd.Series(assign.values, index=obs.index[cand_idx])
    cls_mal = cls_mal.dropna()
    full = cls.dropna()
    full = pd.concat([full, cls_mal.map(lambda s: f"malignant_{s}")])
    sizes = full.value_counts()
    # cap per class so the total stays <= TOTAL_CAP: shrink the cap until it fits
    cap = CLASS_CAP
    while np.minimum(sizes, cap).sum() > TOTAL_CAP:
        cap -= 100
    LOG.write(f"per-class cap {cap} (requested <= {CLASS_CAP}), total {int(np.minimum(sizes, cap).sum())}")
    keep = []
    pos = pd.Series(np.arange(n_cells), index=obs.index)
    for c, n in sizes.items():
        ids = full.index[full == c]
        idx = np.sort(pos[ids].values)
        idx = stratified(idx, obs.sample_id.astype(str).values[idx], min(cap, len(idx)), rng, floor=20)
        keep.append(idx)
    keep = np.sort(np.concatenate(keep))
    labels = full.reindex(obs.index[keep]).values

    # ---------------- gene set
    hvg = pd.read_csv(TABLES / "atlas_input_hvg.csv"); hvg = set(hvg.symbol[hvg.highly_variable])
    mhvg = pd.read_csv(TABLES / "atlas_malignant_hvg.csv"); mhvg = set(mhvg.symbol[mhvg.highly_variable])
    markers = set(all_marker_genes())
    genes = sorted((hvg | mhvg | markers) & set(var.index))
    gmask = var.index.isin(genes)
    LOG.write(f"gene set: HVG {len(hvg)} ∪ malignant HVG {len(mhvg)} ∪ markers {len(markers)} -> {len(genes)} genes "
              f"({len(markers - set(var.index))} marker genes absent from the atlas: {sorted(markers - set(var.index))})")

    with LOG.step(f"read {len(keep)} reference rows and write h5ad"):
        X = read_rows(f, keep)[:, np.where(gmask)[0]].tocsr()
        X.sort_indices()
        o = obs.iloc[keep][["sample_id", "patient_id", "study", "assay", "pct_mito", "lineage_label", "scanvi_label"]].copy()
        for c in o.columns:
            if c != "pct_mito":
                o[c] = o[c].astype(str)
        o["c2l_class"] = labels.astype(str)
        o["pct_mito"] = o["pct_mito"].astype(np.float32)
        ref = ad.AnnData(X=X.astype(np.int32), obs=o, var=pd.DataFrame(index=var.index[gmask]))
        ref.var["in_reference_hvg"] = ref.var_names.isin(hvg)
        ref.var["in_malignant_hvg"] = ref.var_names.isin(mhvg)
        ref.var["is_marker"] = ref.var_names.isin(markers)
        ref.uns["phase6"] = {"state_source": state_source, "seed": SEED, "class_cap": int(cap), "n_cells": int(ref.n_obs)}
        MODAL_IN.mkdir(parents=True, exist_ok=True)
        out = MODAL_IN / f"c2l_reference{tag}.h5ad"
        ref.write_h5ad(out, compression="gzip", compression_opts=4)
        mb = out.stat().st_size / 1e6
        LOG.write(f"wrote {out} ({mb:.1f} MB, {ref.shape}, nnz {X.nnz}, max count {int(X.max())})")
    tab = o.groupby("c2l_class").agg(n_cells=("sample_id", "size"), n_samples=("sample_id", "nunique"),
                                     n_studies=("study", "nunique")).reset_index()
    tab["state_source"] = state_source
    tab.to_csv(TABLES / f"c2l_reference_classes{tag}.csv", index=False)
    LOG.write("reference classes:\n" + tab.to_string(index=False))
    LOG.write(f"done | peak RSS {peak_rss_gb():.2f} GB")


def collect(harvest_dir: str, tag: str = ""):
    h = Path(harvest_dir)
    sig = pd.read_csv(h / "signatures.csv", index_col=0)
    sig.to_csv(TABLES / f"c2l_signatures{tag}.csv")
    run = json.load(open(h / "signatures_run.json")) if (h / "signatures_run.json").exists() else {}
    LOG.write(f"collected signatures {sig.shape} -> results/tables/c2l_signatures{tag}.csv; run {json.dumps(run)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["prepare", "collect"])
    ap.add_argument("--harvest", default=None)
    ap.add_argument("--use-metaprograms", action="store_true")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    tag = f"_{a.tag}" if a.tag else ""
    if a.use_metaprograms and not METAPROG.exists():
        raise SystemExit(f"--use-metaprograms but {METAPROG} is missing (run 19_metaprogram_scores.py)")
    prepare(a.use_metaprograms, tag) if a.mode == "prepare" else collect(a.harvest, tag)
