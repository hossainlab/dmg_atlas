#!/usr/bin/env python
"""Audit remediation (reference-labels track, 19 Oct 2026) — rebuild data/processed/atlas_input_malignant_hvg.h5ad
from the CORRECTED consensus_call (malignant_calls.csv regenerated in the CNV track), excluding the
low-confidence sample okada__1734_pons (sample_qc_pass False, DQ-03; 1,163 genes / 22 CNV windows, 64.5 % score-only
CNV calls in its reference cells, CNV-09). Same recipe as scripts/07b_make_training_files.py (malignant subset part):
3,000 seurat_v3_paper HVGs with batch_key = sample_id (samples < 100 malignant cells pooled), layer counts, span 0.5
(07b tried 0.3 first and fell back to 0.5 on failure; here 0.5 is used directly and recorded).
Previous file kept as atlas_input_malignant_hvg.pre_audit.h5ad; previous HVG table as atlas_malignant_hvg.pre_audit.csv.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402
import scanpy as sc  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_reflabels.log", context="07d_rebuild_malignant_input")
OUT = ROOT / "data/processed"
LABELS = ROOT / "results/tables/reflabels/reference_labels_corrected.csv.gz"
SLIM = ["sample_id", "patient_id", "study", "assay", "pct_mito", "lineage_label", "consensus_call", "consensus_celltype_coarse"]
N_HVG = 3000
MIN_BATCH = 100
SPAN = 0.5
EXCLUDE_SAMPLES = ["okada__1734_pons"]


def slim_obs(obs):
    o = obs[SLIM].copy()
    for c in o.columns:
        if c != "pct_mito":
            o[c] = o[c].astype(str).astype("category")
    o["pct_mito"] = o["pct_mito"].astype(np.float32)
    return o


def read_rows_csr(path, rows, n_vars, chunk=20000):
    """Row-span CSR reader (verbatim from 07b)."""
    parts = []
    with h5py.File(path) as f:
        g = f["layers/counts"]
        indptr = g["indptr"][:]
        for a in range(0, len(rows), chunk):
            r = rows[a:a + chunk]
            lo, hi = int(indptr[r[0]]), int(indptr[r[-1] + 1])
            block = sp.csr_matrix((g["data"][lo:hi], g["indices"][lo:hi], indptr[r[0]:r[-1] + 2] - lo), shape=(r[-1] - r[0] + 1, n_vars))
            parts.append(block[r - r[0]])
    return sp.vstack(parts, format="csr")


def main():
    for src, dst in [(OUT / "atlas_input_malignant_hvg.h5ad", OUT / "atlas_input_malignant_hvg.pre_audit.h5ad"),
                     (ROOT / "results/tables/atlas_malignant_hvg.csv", ROOT / "results/tables/atlas_malignant_hvg.pre_audit.csv")]:
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)
            LOG.write(f"pre-audit copy kept: {dst.name}")

    with LOG.step("select corrected consensus-malignant reference cells"):
        with h5py.File(OUT / "atlas_input.h5ad") as f:
            obs = read_elem(f["obs"])
            var = read_elem(f["var"])
        lab = pd.read_csv(LABELS, index_col="cell_id", usecols=["cell_id", "consensus_call", "consensus_call_pre_audit", "lineage_label", "consensus_celltype_coarse", "sample_id"])
        lab = lab.loc[obs.index]
        obs["consensus_call_pre_audit"] = obs["consensus_call"].astype(str).values
        obs["consensus_call"] = lab["consensus_call"].astype(str).values
        obs["lineage_label"] = lab["lineage_label"].astype(str).values
        obs["consensus_celltype_coarse"] = lab["consensus_celltype_coarse"].astype(str).values
        mal = obs["consensus_call"].eq("malignant").values
        excl = obs["sample_id"].astype(str).isin(EXCLUDE_SAMPLES).values
        LOG.write(f"corrected consensus malignant {int(mal.sum()):,} (pre-audit {int(obs['consensus_call_pre_audit'].eq('malignant').sum()):,}); "
                  f"excluded low-confidence sample cells among them {int((mal & excl).sum()):,} ({EXCLUDE_SAMPLES})")
        rows = np.flatnonzero(mal & ~excl)
        X = read_rows_csr(OUT / "atlas_input.h5ad", rows, len(var))
        m = ad.AnnData(X=X, obs=slim_obs(obs.iloc[rows]), var=pd.DataFrame(index=var.index))
        m.layers["counts"] = m.X
        LOG.write(f"malignant-only input: {m.shape}; peak RSS {peak_rss_gb():.2f} GB")

    with LOG.step(f"HVG seurat_v3_paper batch=sample_id (pooled < {MIN_BATCH}), span {SPAN}"):
        n_per = m.obs["sample_id"].value_counts()
        small = n_per[n_per < MIN_BATCH].index
        hb = m.obs["sample_id"].astype(str).where(~m.obs["sample_id"].isin(small), "pooled_small_samples")
        m.obs["hvg_batch"] = hb.astype("category")
        LOG.write(f"{n_per.size} samples; {len(small)} with < {MIN_BATCH} malignant cells pooled ({int(n_per[small].sum())} cells): {list(small)}")
        sc.pp.highly_variable_genes(m, flavor="seurat_v3_paper", n_top_genes=N_HVG, batch_key="hvg_batch", layer="counts", span=SPAN)
        old = pd.read_csv(ROOT / "results/tables/atlas_malignant_hvg.pre_audit.csv", index_col=0)
        ov = int((old["highly_variable"].reindex(m.var_names).fillna(False).astype(bool) & m.var["highly_variable"]).sum())
        LOG.write(f"HVG nbatches min/median {int(m.var.loc[m.var.highly_variable, 'highly_variable_nbatches'].min())}/"
                  f"{float(m.var.loc[m.var.highly_variable, 'highly_variable_nbatches'].median())}; overlap with pre-audit HVG set {ov}/{N_HVG}")
        m.var[["highly_variable", "highly_variable_rank", "highly_variable_nbatches", "means", "variances_norm"]].to_csv(ROOT / "results/tables/atlas_malignant_hvg.csv")
        mh = m[:, m.var.highly_variable].copy()
        mh.layers["counts"] = mh.X
        mh.uns["phase4"] = {"hvg_flavor": "seurat_v3_paper", "hvg_batch": f"sample_id (samples < {MIN_BATCH} malignant cells pooled)", "span": SPAN,
                            "audit_fix_2026_10_19": f"corrected consensus_call; excluded {EXCLUDE_SAMPLES}; pre-audit HVG overlap {ov}/{N_HVG}"}
        mh.write_h5ad(OUT / "atlas_input_malignant_hvg.h5ad", compression="gzip", compression_opts=4)
        LOG.write(f"atlas_input_malignant_hvg.h5ad: {mh.shape}, {(OUT / 'atlas_input_malignant_hvg.h5ad').stat().st_size / 1e6:.0f} MB; peak RSS {peak_rss_gb():.2f} GB")
    LOG.write("done")


if __name__ == "__main__":
    main()
