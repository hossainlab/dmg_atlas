#!/usr/bin/env python
"""Phase 4, step 1 — assemble the integrated UMI reference cohort.

Cells: Phase-3 working set (qc_pass; Jessa also author_qc_pass; the four
cytoplasm-contaminated nuclei samples use the relaxed 20 % mito rule), minus
predicted_doublet_fixed; scpca excluded entirely (441 unassigned low-quality
nuclei). Genes: the 18,882 symbols measured in every sample of all UMI
studies. Writes
  data/processed/atlas_input.h5ad      (all shared genes; X = layers['counts'] = int32 CSR)
  data/processed/atlas_input_hvg.h5ad  (3,000 seurat_v3 HVGs, batch_key = sample_id)
  results/tables/atlas_input_cells_per_sample.csv, atlas_input_hvg.csv
Counts are read from the Phase-2 files row-chunk by row-chunk through h5py,
so peak RSS is one copy of the final matrix plus one chunk.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path("/Users/jubayer/Projects/dmg_atlas")
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger  # noqa: E402  (sets thread caps first)

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402
import scanpy as sc  # noqa: E402
from dmg_atlas.cnv import working_mask  # noqa: E402

STUDIES = ["jessa_scrna", "jessa_multiome", "okada", "collot", "zang", "desisto"]
EXCLUDED = {"scpca": "441 unassigned low-quality nuclei (SCPCS000637, median 362 genes); SCPCS000636 has 0 qc_pass cells"}
CHUNK = 20_000
N_HVG = 3000
LOG = StepLogger(ROOT / "logs/phase4_integration.log", context="07_assemble")
OUT = ROOT / "data/processed"
TAB = ROOT / "results/tables"
TAB.mkdir(parents=True, exist_ok=True)

CALL_COLS = [
    "cnv_score", "cnv_call", "cnv_confidence", "marker_call", "consensus_call", "consensus_call_eventtrust", "consensus_celltype_coarse",
    "author_call", "marker_type", "cnv_cluster", "high_mito_nuclei", "low_confidence_sample",
    "score_myeloid", "score_T_NK", "score_oligodendrocyte", "score_vascular", "score_neuron",
    "score_astrocyte", "score_OPC_like", "score_AC_like", "score_MES_like", "score_cycling",
    # audit correction CNV-03 (2026-10-19): chr4q12_gain is no longer a hallmark event column; arm_4q12 and gene_dev_PDGFRA/KIT are exploratory
    "chr1q_gain", "chr2_gain", "chr7_gain", "chr13q_loss", "chr14q_loss", "chr10q_loss", "chr17p_loss",
    "arm_1q", "arm_2", "arm_7", "arm_4q12", "arm_13q", "arm_14q", "arm_10q", "arm_17p", "gene_dev_PDGFRA", "gene_dev_KIT",
]
LINEAGES = ["myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular", "other"]


def lineage_label(obs: pd.DataFrame) -> pd.Series:
    lab = pd.Series("Unknown", index=obs.index, dtype=object)
    lab[obs["consensus_call"].eq("malignant")] = "malignant"
    normal = obs["consensus_call"].eq("normal")
    ct = obs["consensus_celltype_coarse"].astype(str)
    lab[normal & ct.isin(LINEAGES)] = ct[normal & ct.isin(LINEAGES)]
    lab[normal & ~ct.isin(LINEAGES)] = "other"
    return lab


def categorise(obs: pd.DataFrame) -> pd.DataFrame:
    for c in obs.columns:
        if obs[c].dtype == object:
            obs[c] = obs[c].astype(str).where(obs[c].notna(), None).astype("category")
    return obs


def main():
    shared = pd.read_csv(ROOT / "results/qc/genes_shared_all_umi_studies.csv")["symbol"].tolist()
    LOG.write(f"shared UMI gene list: {len(shared)} symbols")

    with LOG.step("load malignant_calls"):
        calls = pd.read_csv(OUT / "malignant_calls.csv", index_col="cell_id", usecols=["cell_id"] + CALL_COLS, low_memory=False)

    # ---- pass 0: obs, masks, var alignment -------------------------------------------------
    plans = []
    obs_parts = []
    gene_set = set(shared)
    for s in STUDIES:
        path = ROOT / f"data/interim/{s}.qc.h5ad"
        with h5py.File(path) as f:
            obs = read_elem(f["obs"])
            var = read_elem(f["var"])
        mask, hm = working_mask(obs, s)
        n_work = int(mask.sum())
        dbl = obs["predicted_doublet_fixed"].fillna(False).astype(bool)
        mask = mask & ~dbl
        in_calls = obs.index.isin(calls.index)
        miss = int((mask & ~in_calls).sum())
        if miss:
            LOG.write(f"WARNING {s}: {miss} selected cells have no malignant_calls row -> dropped")
            mask = mask & in_calls
        rows = np.flatnonzero(mask.values)
        sub = obs.iloc[rows].copy()
        sub["high_mito_nuclei"] = hm.iloc[rows].values
        LOG.write(f"{s}: barcodes {len(obs)} | phase3 working set {n_work} | doublets removed {int((dbl & working_mask(obs, s)[0]).sum())} | kept {len(rows)}")
        gene_set &= set(var.index)
        plans.append((s, path, rows, var.index.to_numpy()))
        obs_parts.append(sub)
    genes = [g for g in shared if g in gene_set]
    LOG.write(f"genes present in all six studies: {len(genes)}")
    gpos = {g: i for i, g in enumerate(genes)}

    # ---- pass 1: upper bound on nnz from indptr only --------------------------------------
    nnz_ub = 0
    for s, path, rows, _ in plans:
        with h5py.File(path) as f:
            indptr = f["layers/counts/indptr"][:]
        nnz_ub += int((indptr[rows + 1] - indptr[rows]).sum())
    n_cells = sum(len(p[2]) for p in plans)
    LOG.write(f"cells {n_cells}, nnz upper bound {nnz_ub} ({nnz_ub * 8 / 1e9:.2f} GB)")

    # ---- pass 2: fill ---------------------------------------------------------------------
    data = np.empty(nnz_ub, dtype=np.int32)
    indices = np.empty(nnz_ub, dtype=np.int32)
    indptr_out = np.zeros(n_cells + 1, dtype=np.int64)
    pos = 0
    row0 = 0
    with LOG.step("read counts (row chunks)"):
        for s, path, rows, varnames in plans:
            colmap = np.full(len(varnames), -1, dtype=np.int64)
            for j, g in enumerate(varnames):
                k = gpos.get(g)
                if k is not None:
                    colmap[j] = k
            with h5py.File(path) as f:
                g = f["layers/counts"]
                indptr = g["indptr"][:]
                for a in range(0, len(rows), CHUNK):
                    r = rows[a:a + CHUNK]
                    lo, hi = int(indptr[r[0]]), int(indptr[r[-1] + 1])
                    d = g["data"][lo:hi]
                    ix = g["indices"][lo:hi]
                    ip = indptr[r[0]:r[-1] + 2] - lo
                    block = sp.csr_matrix((d, ix, ip), shape=(len(ip) - 1, len(varnames)))
                    block = block[r - r[0]]                      # drop unselected rows in the span
                    block = block[:, colmap >= 0]                # keep measured shared genes (study order)
                    block.indices = colmap[colmap >= 0][block.indices].astype(np.int32)
                    block.sort_indices()
                    n = block.nnz
                    data[pos:pos + n] = block.data
                    indices[pos:pos + n] = block.indices
                    indptr_out[row0 + 1:row0 + block.shape[0] + 1] = pos + block.indptr[1:]
                    pos += n
                    row0 += block.shape[0]
            LOG.write(f"{s}: filled rows -> {row0}, nnz so far {pos}")
    X = sp.csr_matrix((data[:pos], indices[:pos], indptr_out), shape=(n_cells, len(genes)))
    del data, indices
    LOG.write(f"matrix {X.shape}, nnz {X.nnz}, density {X.nnz / X.shape[0] / X.shape[1]:.4f}")

    # ---- obs ---------------------------------------------------------------------------
    obs = pd.concat(obs_parts, axis=0)
    assert obs.index.is_unique
    obs = obs.drop(columns=[c for c in CALL_COLS if c in obs.columns and c != "high_mito_nuclei"], errors="ignore")
    joined = calls.loc[obs.index, [c for c in CALL_COLS if c != "high_mito_nuclei"]]
    obs = pd.concat([obs, joined], axis=1)
    obs["lineage_label"] = lineage_label(obs)
    obs["pct_mito"] = obs["pct_mito"].astype(np.float32)
    obs = categorise(obs)
    var = pd.DataFrame(index=pd.Index(genes, name="symbol"))

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = adata.X
    adata.uns["phase4"] = {
        "studies": STUDIES, "excluded": EXCLUDED, "n_genes_shared": len(genes),
        "hvg_flavor": "seurat_v3_paper", "cell_rule": "qc_pass (Jessa: & author_qc_pass; 4 high-mito nuclei samples: mito<20%) & ~predicted_doublet_fixed",
    }

    # ---- HVGs ----------------------------------------------------------------------------
    with LOG.step("HVG seurat_v3 batch_key=sample_id"):
        span_used = None
        for span in (0.3, 0.5):
            try:
                # scanpy's flavor='seurat_v3' ranks genes by median per-batch rank FIRST, so with 57 batches
                # (some of 150-300 cells) ~25 % of the 3,000 genes were HVG in only 1-2 samples (e.g. FOXA1, OPN4).
                # 'seurat_v3_paper' uses the same per-batch statistics but ranks by the number of batches first
                # (the Seurat/scvi-tools behaviour); both rankings are kept in var for the record.
                sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=N_HVG, batch_key="sample_id", layer="counts", span=span)
                adata.var["hvg_seurat_v3_medianrank_first"] = adata.var["highly_variable"].values
                sc.pp.highly_variable_genes(adata, flavor="seurat_v3_paper", n_top_genes=N_HVG, batch_key="sample_id", layer="counts", span=span)
                ov = int((adata.var["hvg_seurat_v3_medianrank_first"] & adata.var["highly_variable"]).sum())
                LOG.write(f"HVG flavour seurat_v3_paper (nbatches first) selected; overlap with scanpy seurat_v3 (median-rank first) = {ov}/{N_HVG}; "
                          f"nbatches min/median of selected = {int(adata.var.loc[adata.var.highly_variable, 'highly_variable_nbatches'].min())}/"
                          f"{float(adata.var.loc[adata.var.highly_variable, 'highly_variable_nbatches'].median())}")
                span_used = span
                break
            except Exception as e:  # loess failures on small batches
                LOG.write(f"HVG span={span} failed: {type(e).__name__}: {str(e)[:200]}")
        if span_used is None:
            raise RuntimeError("HVG selection failed")
        adata.uns["phase4"]["hvg_span"] = span_used
        LOG.write(f"HVG span used {span_used}; n HVG {int(adata.var.highly_variable.sum())}")
    adata.var[["highly_variable", "highly_variable_rank", "means", "variances", "variances_norm", "highly_variable_nbatches", "hvg_seurat_v3_medianrank_first"]].to_csv(TAB / "atlas_input_hvg.csv")

    # ---- reports ---------------------------------------------------------------------------
    per_sample = obs.groupby(["study", "sample_id", "assay"], observed=True).size().rename("n_cells").reset_index()
    per_sample.to_csv(TAB / "atlas_input_cells_per_sample.csv", index=False)
    LOG.write("cells per study: " + obs.groupby("study", observed=True).size().to_dict().__repr__())
    LOG.write("cells per assay: " + obs.groupby("assay", observed=True).size().to_dict().__repr__())
    LOG.write("lineage_label counts: " + obs["lineage_label"].value_counts().to_dict().__repr__())
    LOG.write(f"total cells {adata.n_obs}, samples {obs['sample_id'].nunique()}")

    # ---- write ----------------------------------------------------------------------------
    with LOG.step("write atlas_input_hvg.h5ad"):
        hv = adata[:, adata.var.highly_variable].copy()
        hv.layers["counts"] = hv.X
        hv.write_h5ad(OUT / "atlas_input_hvg.h5ad")
    with LOG.step("write atlas_input.h5ad"):
        adata.write_h5ad(OUT / "atlas_input.h5ad")
    LOG.write("done")


if __name__ == "__main__":
    main()
