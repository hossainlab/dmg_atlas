#!/usr/bin/env python
"""Phase 4 — slim, compressed training files for the Modal Volume.

  data/processed/atlas_input_hvg.train.h5ad            all reference cells x 3,000 HVGs (counts), slim obs
  data/processed/atlas_input_malignant_hvg.h5ad        consensus malignant cells x 3,000 HVGs recomputed on that subset
Slim obs = sample_id, patient_id, study, assay, pct_mito, lineage_label, consensus_call, consensus_celltype_coarse.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path("/Users/jubayer/Projects/dmg_atlas")
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402
import scanpy as sc  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase4_integration.log", context="07b_training_files")
OUT = ROOT / "data/processed"
SLIM = ["sample_id", "patient_id", "study", "assay", "pct_mito", "lineage_label", "consensus_call", "consensus_celltype_coarse"]
N_HVG = 3000
MIN_BATCH = 100


def slim_obs(obs):
    o = obs[SLIM].copy()
    for c in o.columns:
        if c != "pct_mito":
            o[c] = o[c].astype(str).astype("category")
    o["pct_mito"] = o["pct_mito"].astype(np.float32)
    return o


def read_rows_csr(path, rows, n_vars, chunk=20000):
    parts = []
    with h5py.File(path) as f:
        g = f["layers/counts"]
        indptr = g["indptr"][:]
        for a in range(0, len(rows), chunk):
            r = rows[a:a + chunk]
            lo, hi = int(indptr[r[0]]), int(indptr[r[-1] + 1])
            block = sp.csr_matrix((g["data"][lo:hi], g["indices"][lo:hi], indptr[r[0]:r[-1] + 2] - lo), shape=(len(r) if r[-1] - r[0] + 1 == len(r) else r[-1] - r[0] + 1, n_vars))
            parts.append(block[r - r[0]])
    return sp.vstack(parts, format="csr")


def main():
    with LOG.step("slim reference training file"):
        a = ad.read_h5ad(OUT / "atlas_input_hvg.h5ad")
        t = ad.AnnData(X=a.layers["counts"], obs=slim_obs(a.obs), var=a.var[["highly_variable_rank", "highly_variable_nbatches"]].copy())
        t.layers["counts"] = t.X
        t.write_h5ad(OUT / "atlas_input_hvg.train.h5ad", compression="gzip", compression_opts=4)
        LOG.write(f"atlas_input_hvg.train.h5ad: {t.shape}, {(OUT / 'atlas_input_hvg.train.h5ad').stat().st_size / 1e6:.0f} MB")
        del a, t

    with LOG.step("malignant subset + HVG"):
        with h5py.File(OUT / "atlas_input.h5ad") as f:
            obs = read_elem(f["obs"])
            var = read_elem(f["var"])
        rows = np.flatnonzero(obs["consensus_call"].astype(str).eq("malignant").values)
        X = read_rows_csr(OUT / "atlas_input.h5ad", rows, len(var))
        m = ad.AnnData(X=X, obs=slim_obs(obs.iloc[rows]), var=pd.DataFrame(index=var.index))
        m.layers["counts"] = m.X
        n_per = m.obs["sample_id"].value_counts()
        small = n_per[n_per < MIN_BATCH].index
        hb = m.obs["sample_id"].astype(str).where(~m.obs["sample_id"].isin(small), "pooled_small_samples")
        m.obs["hvg_batch"] = hb.astype("category")
        LOG.write(f"malignant cells {m.n_obs} in {n_per.size} samples; {len(small)} samples with < {MIN_BATCH} malignant cells pooled for HVG "
                  f"({int(n_per[small].sum())} cells): {list(small)}")
        span_used = None
        for span in (0.3, 0.5):
            try:
                sc.pp.highly_variable_genes(m, flavor="seurat_v3_paper", n_top_genes=N_HVG, batch_key="hvg_batch", layer="counts", span=span)
                span_used = span
                break
            except Exception as e:
                LOG.write(f"malignant HVG span={span} failed: {type(e).__name__}: {str(e)[:200]}")
        LOG.write(f"malignant HVG span {span_used}; nbatches min/median {int(m.var.loc[m.var.highly_variable, 'highly_variable_nbatches'].min())}/"
                  f"{float(m.var.loc[m.var.highly_variable, 'highly_variable_nbatches'].median())}")
        m.var[["highly_variable", "highly_variable_rank", "highly_variable_nbatches", "means", "variances_norm"]].to_csv(ROOT / "results/tables/atlas_malignant_hvg.csv")
        mh = m[:, m.var.highly_variable].copy()
        mh.layers["counts"] = mh.X
        mh.uns["phase4"] = {"hvg_flavor": "seurat_v3_paper", "hvg_batch": f"sample_id (samples < {MIN_BATCH} malignant cells pooled)", "span": span_used}
        mh.write_h5ad(OUT / "atlas_input_malignant_hvg.h5ad", compression="gzip", compression_opts=4)
        LOG.write(f"atlas_input_malignant_hvg.h5ad: {mh.shape}, {(OUT / 'atlas_input_malignant_hvg.h5ad').stat().st_size / 1e6:.0f} MB")


if __name__ == "__main__":
    main()
