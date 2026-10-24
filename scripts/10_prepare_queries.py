#!/usr/bin/env python
"""Phase 5 / step 10 — build scANVI query files restricted to the 3,000 reference HVGs.

Queries
  filbin   Smart-seq2, linear TPM            -> pseudo-counts (see below)
  labelle  Smart-seq2 + 10x (P9) log2(TPM/10+1) -> pseudo-counts
  liu      Smart-seq2, RSEM expected counts  -> round()
  braun    Braun 2023 first-trimester hindbrain/pons/medulla, 10x UMI counts (subsampled, stratified)
  siletti  Siletti 2023 adult pons (6 cellxgene dissections), 10x snRNA UMI counts (subsampled, stratified)

Pseudo-counts for TPM cohorts (documented deviation, these cohorts are projection-only):
  linear TPM per cell is rescaled so that the cell's total over all harmonised genes equals
  n_detected_genes(cell) x m, where m is the median (counts / detected genes) of the matching
  counts-based cohort in this project: Liu RSEM Smart-seq2 for Smart-seq2 tables, Jessa 10x scRNA
  for the LaBelle P9 10x CPM table. Values are then rounded to integers. The HVG subset is taken after
  rescaling. No MT genes exist in the Filbin/LaBelle tables, so pct_mito = 0 for those cells.

Outputs: data/processed/queries/<query>.hvg.h5ad (layers['counts'] int32 CSR, X = same),
         results/tables/query_preparation.csv, logs/phase5_mapping.log
"""
from __future__ import annotations

import sys
import json
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402  (sets thread caps)

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase5_mapping.log", context="10_prepare_queries")
OUT = ROOT / "data/processed/queries"
OUT.mkdir(parents=True, exist_ok=True)
HVG = pd.read_csv(ROOT / "results/tables/atlas_input_hvg.csv")
HVG = HVG.loc[HVG.highly_variable, "symbol"].tolist()
assert len(HVG) == 3000
# reference var order (the model validates var_names order)
with h5py.File(ROOT / "data/processed/atlas_input_hvg.train.h5ad") as f:
    REF_VAR = read_elem(f["var"]).index.tolist()
assert set(REF_VAR) == set(HVG)
UNIVERSE = pd.read_csv(ROOT / "data/interim/gene_universe.csv")
SYM2ENS = dict(zip(UNIVERSE.canonical_symbol, UNIVERSE.ensembl_id))
MAX_NORMAL = 100_000
RNG = np.random.default_rng(0)
records = []


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def to_ref_order(X: sp.csr_matrix, var_names: list[str]) -> tuple[sp.csr_matrix, int]:
    """Column-align X (cells x var_names) to REF_VAR, zero-filling missing genes. Returns (matrix, n_missing)."""
    pos = {g: i for i, g in enumerate(var_names)}
    cols, present = [], []
    for g in REF_VAR:
        if g in pos:
            cols.append(pos[g]); present.append(True)
        else:
            present.append(False)
    present = np.array(present)
    Xs = X[:, cols].tocsr()
    if present.all():
        return Xs.astype(np.int32), 0
    Xc = Xs.tocoo()
    newcol = np.flatnonzero(present)[Xc.col]
    full = sp.csr_matrix((Xc.data.astype(np.int32), (Xc.row, newcol)), shape=(X.shape[0], len(REF_VAR)))
    return full, int((~present).sum())


def finish(name: str, X: sp.csr_matrix, obs: pd.DataFrame, n_missing: int, notes: dict, src: Path | list[Path]):
    X = sp.csr_matrix(X, dtype=np.int32)
    X.eliminate_zeros()
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(REF_VAR, name="symbol")))
    adata.layers["counts"] = adata.X.copy()
    adata.obs["lineage_label"] = "Unknown"
    for c in ["sample_id", "study", "assay", "lineage_label"]:
        adata.obs[c] = adata.obs[c].astype(str).astype("category")
    adata.obs["pct_mito"] = adata.obs["pct_mito"].astype(np.float32)
    adata.obs["query_n_counts_hvg"] = np.asarray(X.sum(1)).ravel()
    adata.obs["query_n_genes_hvg"] = np.asarray((X > 0).sum(1)).ravel()
    srcs = src if isinstance(src, list) else [src]
    adata.uns["phase5_query"] = {"query": name, "n_missing_hvg": n_missing, "sources": [str(s) for s in srcs], **notes}
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].astype(str).replace("nan", np.nan).astype("category")
    out = OUT / f"{name}.hvg.h5ad"
    adata.write_h5ad(out, compression="gzip")
    rec = {"query": name, "cells": adata.n_obs, "samples": adata.obs.sample_id.nunique(), "n_missing_hvg": n_missing,
           "median_hvg_counts_per_cell": float(np.median(adata.obs.query_n_counts_hvg)),
           "median_hvg_genes_per_cell": float(np.median(adata.obs.query_n_genes_hvg)),
           "assays": ";".join(sorted(adata.obs.assay.unique())), "file": str(out.relative_to(ROOT)),
           "file_md5": md5(out), **{k: v for k, v in notes.items() if not isinstance(v, (dict, list))}}
    records.append(rec)
    LOG.write(f"{name}: {adata.n_obs} cells, {adata.obs.sample_id.nunique()} samples, missing HVGs {n_missing}, "
              f"median HVG counts/cell {rec['median_hvg_counts_per_cell']:.0f}")
    return adata


def median_counts_per_gene(study: str, assay: str) -> float:
    with h5py.File(ROOT / f"data/interim/{study}.qc.h5ad") as f:
        obs = read_elem(f["obs"])
    q = obs[(obs.qc_pass.astype(bool)) & (obs.assay.astype(str) == assay)]
    return float((q.n_counts / q.n_genes).median())


# ----------------------------------------------------------------------------------------------- TPM cohorts
def tpm_cohort(name: str):
    path = ROOT / f"data/interim/{name}.qc.h5ad"
    a = ad.read_h5ad(path)
    a = a[a.obs.qc_pass.astype(bool)].copy()
    layer = "tpm" if "tpm" in a.layers else "log2_tpm10"  # audit DQ-01: LaBelle layer renamed to its true scale
    T = sp.csr_matrix(a.layers[layer], dtype=np.float64)
    rowmax = np.asarray(T.max(1).todense()).ravel()
    is_log = rowmax < 30  # Phase 3 rule: log2(TPM/10+1) tables have per-cell maxima 9-12.4
    notes = {"n_cells_log2_scale": int(is_log.sum()), "n_cells_linear_scale": int((~is_log).sum())}
    if is_log.any():
        Tl = T[is_log].copy()
        Tl.data = 10.0 * (np.power(2.0, Tl.data) - 1.0)
        T = sp.vstack([Tl, T[~is_log]]).tocsr() if (~is_log).any() else Tl
        order = np.concatenate([np.flatnonzero(is_log), np.flatnonzero(~is_log)])
        a = a[order].copy()
    n_det = np.asarray((T > 0).sum(1)).ravel().astype(float)
    rowsum = np.asarray(T.sum(1)).ravel()
    m_ss2 = median_counts_per_gene("liu", "smartseq2")
    m_10x = median_counts_per_gene("jessa_scrna", "scRNA")
    assay = a.obs.assay.astype(str).values
    m = np.where(assay == "smartseq2", m_ss2, m_10x)
    target = n_det * m
    scale = np.where(rowsum > 0, target / np.maximum(rowsum, 1e-12), 0.0)
    C = sp.diags(scale) @ T
    C.data = np.round(C.data)
    C.eliminate_zeros()
    C = sp.csr_matrix(C, dtype=np.int32)
    notes.update({"pseudocount_rule": "round(linear_TPM * n_detected_genes * m / rowsum)",
                  "m_smartseq2_counts_per_gene_liu": round(m_ss2, 3), "m_10x_counts_per_gene_jessa_scrna": round(m_10x, 3),
                  "median_pseudocount_total_all_genes": float(np.median(np.asarray(C.sum(1)).ravel())),
                  "pct_mito_rule": "0 (no MT- genes in deposited TPM tables)"})
    Xh, n_missing = to_ref_order(C, a.var_names.tolist())
    not_all = int((~a.var.loc[[g for g in REF_VAR if g in a.var_names], "measured_in_all_samples"].astype(bool)).sum())
    notes["n_hvg_not_measured_in_all_samples"] = not_all
    obs = a.obs.copy()
    obs["pct_mito"] = 0.0
    obs["query"] = name
    finish(name, Xh, obs, n_missing, notes, path)


# ----------------------------------------------------------------------------------------------- Liu
def liu():
    path = ROOT / "data/interim/liu.qc.h5ad"
    a = ad.read_h5ad(path)
    a = a[a.obs.qc_pass.astype(bool)].copy()
    C = sp.csr_matrix(a.layers["counts"], dtype=np.float64)
    frac_nonint = float(np.mean(C.data != np.round(C.data)))
    C.data = np.round(C.data)
    C.eliminate_zeros()
    C = sp.csr_matrix(C, dtype=np.int32)
    Xh, n_missing = to_ref_order(C, a.var_names.tolist())
    not_all = int((~a.var.loc[[g for g in REF_VAR if g in a.var_names], "measured_in_all_samples"].astype(bool)).sum())
    obs = a.obs.copy()
    obs["query"] = "liu"
    notes = {"counts_rule": "round(RSEM expected counts)", "frac_nonint_values_before_rounding": round(frac_nonint, 4),
             "n_hvg_not_measured_in_all_samples": not_all, "pct_mito_rule": "Phase 2 obs pct_mito (RSEM counts)"}
    finish("liu", Xh, obs, n_missing, notes, path)


# ----------------------------------------------------------------------------------------------- normal references
def stratified_subsample(obs: pd.DataFrame, strata_cols: list[str], n_max: int, floor: int = 100) -> np.ndarray:
    """Proportional allocation per stratum with a per-stratum floor (min(size, floor)); returns sorted row indices."""
    n = len(obs)
    if n <= n_max:
        return np.arange(n)
    key = obs[strata_cols].astype(str).agg("|".join, axis=1).values
    groups = pd.Series(np.arange(n)).groupby(key)
    sizes = groups.size()
    alloc = np.floor(sizes * n_max / n).astype(int)
    alloc = np.maximum(alloc, np.minimum(sizes, floor))
    # if floors pushed the total over n_max, shrink the strata above the floor proportionally
    excess = int(alloc.sum() - n_max)
    if excess > 0:
        room = (alloc - floor).clip(lower=0)
        cut = np.floor(room * excess / max(room.sum(), 1)).astype(int)
        alloc = alloc - cut
        while alloc.sum() > n_max:
            alloc[alloc.idxmax()] -= 1
    picks = []
    for k, idx in groups:
        idx = idx.values
        take = int(alloc[k])
        picks.append(idx if take >= len(idx) else RNG.choice(idx, take, replace=False))
    return np.sort(np.concatenate(picks))


def read_rows_hvg(path: Path, rows: np.ndarray, chunk: int = 20_000):
    """Read selected rows of a CSR h5ad X in row chunks; return (HVG-aligned csr, total counts, MT counts, n_missing)."""
    with h5py.File(path) as f:
        var = read_elem(f["var"])
        ens = pd.Index(var.index.astype(str).str.replace(r"\.\d+$", "", regex=True))
        sym = var["Gene"].astype(str).values if "Gene" in var else var["feature_name"].astype(str).values
        # HVG -> column: Ensembl ID first, symbol second
        col_of = {}
        sympos = {}
        for i, s in enumerate(sym):
            sympos.setdefault(s, i)
        for g in REF_VAR:
            e = SYM2ENS.get(g)
            j = ens.get_loc(e) if (e is not None and e in ens) else None
            if j is None and g in sympos:
                j = sympos[g]
            if j is not None:
                col_of[g] = j
        n_missing = len(REF_VAR) - len(col_of)
        hv_cols = np.array([col_of[g] for g in REF_VAR if g in col_of])
        hv_target = np.array([i for i, g in enumerate(REF_VAR) if g in col_of])
        mt_cols = np.flatnonzero(pd.Series(sym).str.startswith("MT-").values)
        shape = tuple(f["X"].attrs["shape"])
        indptr = f["X/indptr"][:]
        data_ds, idx_ds = f["X/data"], f["X/indices"]
        blocks, tot, mt = [], [], []
        rows = np.asarray(rows)
        for start in range(0, shape[0], chunk):
            stop = min(start + chunk, shape[0])
            sel = rows[(rows >= start) & (rows < stop)]
            if len(sel) == 0:
                continue
            a, b = indptr[start], indptr[stop]
            ip = indptr[start:stop + 1] - a
            M = sp.csr_matrix((data_ds[a:b], idx_ds[a:b], ip), shape=(stop - start, shape[1]))
            M = M[sel - start]
            tot.append(np.asarray(M.sum(1)).ravel())
            mt.append(np.asarray(M[:, mt_cols].sum(1)).ravel() if len(mt_cols) else np.zeros(M.shape[0]))
            H = M[:, hv_cols].tocoo()
            blocks.append(sp.csr_matrix((H.data, (H.row, hv_target[H.col])), shape=(M.shape[0], len(REF_VAR))))
            del M
        X = sp.vstack(blocks).tocsr()
    return X, np.concatenate(tot), np.concatenate(mt), n_missing, len(mt_cols)


def braun():
    path = ROOT / "data/external/braun2023/braun2023_hindbrain_pons_medulla.h5ad"
    with h5py.File(path) as f:
        obs = read_elem(f["obs"])
    if "CellID" not in obs.columns:
        obs["CellID"] = obs.index.astype(str)
    rows = stratified_subsample(obs, ["CellClass", "donor_id"], MAX_NORMAL)
    LOG.write(f"braun: {len(obs)} cells -> {len(rows)} after stratified subsampling (CellClass x donor_id, floor 100)")
    X, tot, mt, n_missing, n_mt = read_rows_hvg(path, rows)
    o = obs.iloc[rows].copy()
    o.index = ("braun__" + o["donor_id"].astype(str) + ":" + o["CellID"].astype(str)).values
    o["original_cell_id"] = obs.iloc[rows]["CellID"].astype(str).values
    o["sample_id"] = ("braun__" + o["donor_id"].astype(str)).values
    o["study"] = "braun2023"
    o["assay"] = "scRNA"
    o["query"] = "braun"
    o["n_counts"] = tot
    o["n_genes"] = np.nan
    o["pct_mito"] = np.where(tot > 0, 100.0 * mt / np.maximum(tot, 1), 0.0)
    o["author_CellClass"] = o["CellClass"].astype(str)
    o["author_cell_type"] = o["cell_type"].astype(str)
    o["author_Region"] = o["Region"].astype(str)
    o["author_Age_pcw"] = pd.to_numeric(o["Age"], errors="coerce")
    o["author_Subregion"] = o["Subregion"].astype(str)
    keep = ["original_cell_id", "sample_id", "study", "assay", "query", "n_counts", "pct_mito", "donor_id", "author_CellClass",
            "author_cell_type", "author_Region", "author_Subregion", "author_Age_pcw", "development_stage", "sex"]
    o = o[keep]
    notes = {"subsample": f"{len(obs)} -> {len(rows)} stratified by CellClass x donor_id (proportional, floor 100), seed 0",
             "counts_rule": "raw UMI counts (X is integer-valued)", "pct_mito_rule": f"100*MT-/total over all {n_mt} MT- genes",
             "gene_matching": "Ensembl ID (2020-A universe) first, then symbol"}
    finish("braun", X, o, n_missing, notes, path)


def siletti():
    paths = sorted((ROOT / "data/external/siletti2023").glob("siletti2023_pons_*.h5ad"))
    obs_all = []
    for p in paths:
        with h5py.File(p) as f:
            o = read_elem(f["obs"])
        if "CellID" not in o.columns:
            o["CellID"] = o.index.astype(str)
        o["_file"] = p.name
        o["_row"] = np.arange(len(o))
        obs_all.append(o)
    obs = pd.concat(obs_all)
    obs.index = np.arange(len(obs))
    rows = stratified_subsample(obs, ["supercluster_term", "donor_id", "_file"], MAX_NORMAL, floor=50)
    LOG.write(f"siletti: {len(obs)} nuclei -> {len(rows)} after stratified subsampling (supercluster_term x donor_id x dissection, floor 50)")
    sel = obs.iloc[rows]
    Xs, tots, mts, nmiss, nmt = [], [], [], None, None
    for p in paths:
        r = sel.loc[sel._file == p.name, "_row"].values
        X, tot, mt, n_missing, n_mt = read_rows_hvg(p, np.sort(r))
        Xs.append(X); tots.append(tot); mts.append(mt); nmiss, nmt = n_missing, n_mt
    # re-order sel to match the per-file sorted row order
    sel = pd.concat([sel[sel._file == p.name].sort_values("_row") for p in paths])
    X = sp.vstack(Xs).tocsr(); tot = np.concatenate(tots); mt = np.concatenate(mts)
    o = sel.copy()
    o.index = ("siletti__" + o["donor_id"].astype(str) + ":" + o["CellID"].astype(str)).values
    o["original_cell_id"] = sel["CellID"].astype(str).values
    o["sample_id"] = ("siletti__" + o["donor_id"].astype(str) + "__" + o["roi"].astype(str).str.replace("Human ", "")).values
    o["study"] = "siletti2023"
    o["assay"] = "snRNA"
    o["query"] = "siletti"
    o["n_counts"] = tot
    o["pct_mito"] = np.where(tot > 0, 100.0 * mt / np.maximum(tot, 1), 0.0)
    o["author_supercluster_term"] = o["supercluster_term"].astype(str)
    o["author_cell_type"] = o["cell_type"].astype(str)
    o["author_roi"] = o["roi"].astype(str)
    o["author_tissue"] = o["tissue"].astype(str)
    keep = ["original_cell_id", "sample_id", "study", "assay", "query", "n_counts", "pct_mito", "donor_id",
            "author_supercluster_term", "author_cell_type", "author_roi", "author_tissue", "development_stage", "sex"]
    o = o[keep]
    notes = {"subsample": f"{len(obs)} -> {len(rows)} stratified by supercluster_term x donor_id x dissection (proportional, floor 50), seed 0",
             "counts_rule": "raw UMI counts (X is integer-valued)", "pct_mito_rule": f"100*MT-/total over all {nmt} MT- genes",
             "gene_matching": "Ensembl ID (2020-A universe) first, then symbol",
             "sample_id_rule": "donor x dissection (roi); donors are the biological replicates"}
    finish("siletti", X, o, nmiss, notes, paths)


if __name__ == "__main__":
    which = sys.argv[1:] or ["filbin", "labelle", "liu", "braun", "siletti"]
    for q in which:
        with LOG.step(f"prepare {q}"):
            {"filbin": lambda: tpm_cohort("filbin"), "labelle": lambda: tpm_cohort("labelle"), "liu": liu,
             "braun": braun, "siletti": siletti}[q]()
    tab = ROOT / "results/tables/query_preparation.csv"
    df = pd.DataFrame(records)
    if tab.exists():
        old = pd.read_csv(tab)
        df = pd.concat([old[~old["query"].isin(df["query"])], df])
    df.to_csv(tab, index=False)
    LOG.write(f"wrote {tab} | peak RSS {peak_rss_gb():.2f} GB")
