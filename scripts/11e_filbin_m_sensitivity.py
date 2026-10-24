#!/usr/bin/env python
"""Phase 5 audit fix IM-12 (2026-10-19) — sensitivity of the Filbin (TPM) projection to the pseudo-count scale m.

The TPM cohorts are converted to pseudo-counts as round(linear TPM x n_detected x m / rowsum) with m = 33.82 (Liu RSEM median counts per
detected gene). Because the reference encoder is applied zero-shot to log1p(counts) (IM-02), m directly sets where TPM cells land.
This script rebuilds the Filbin pseudo-counts at m x {0.25, 0.5, 1, 2, 4}, encodes each with the reference scANVI model through
SCANVI.load_query_data WITHOUT training (zero-shot projection), and reports for every scale: median 15-NN distance to the reference
(pynndescent index on all 311,618 reference cells), median kNN-uncertainty percentile, kNN-lineage / reference-classifier label agreement
with the m x 1 deliverable, and the predicted-malignant fraction under both labels.
Writes results/tables/phase5_filbin_m_sensitivity.csv. Log: logs/fix_phase5.log
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import anndata as ad  # noqa: E402
import h5py  # noqa: E402
from anndata.io import read_elem  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_phase5.log", context="11e_filbin_m_sensitivity")
REF_DIR = ROOT / "models/scanvi"
SCALES = [0.25, 0.5, 1.0, 2.0, 4.0]
K = 15


def main():
    import torch
    import scvi
    import pynndescent
    from scvi.model import SCANVI
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "8")))
    scvi.settings.seed = 1
    scvi.settings.dl_num_workers = 0
    base = ad.read_h5ad(ROOT / "data/processed/queries/filbin.hvg.h5ad")  # m x 1 deliverable (var in reference order)
    q = ad.read_h5ad(ROOT / "data/interim/filbin.qc.h5ad")
    q = q[base.obs_names].copy()
    T = sp.csr_matrix(q.layers["tpm"], dtype=np.float64)
    assert float(T.max()) >= 30, "Filbin TPM table expected on linear scale"
    n_det = np.asarray((T > 0).sum(1)).ravel().astype(float)
    rowsum = np.asarray(T.sum(1)).ravel()
    m1 = 33.817  # results/tables/query_preparation.csv (Liu RSEM median counts per detected gene)
    col = pd.Series(np.arange(q.n_vars), index=q.var_names)
    keep = [g for g in base.var_names if g in col.index]
    src_idx = col[keep].values
    dst_idx = pd.Series(np.arange(base.n_vars), index=base.var_names)[keep].values

    with LOG.step("reference kNN index (pynndescent, 311,618 cells)"):
        Zref = np.load(REF_DIR / "latent_f16.npy").astype(np.float32)
        ref_names = np.array((REF_DIR / "obs_names.txt").read_text().split("\n"))
        index = pynndescent.NNDescent(Zref, n_neighbors=30, metric="euclidean", random_state=0, n_jobs=8)
        index.prepare()
        ref_lab = pd.read_csv(ROOT / "results/tables/atlas_core_obs_labels.csv", index_col=0).loc[ref_names, "scanvi_label"].astype(str).values
        ref_knn_sorted = np.sort(pd.read_csv(ROOT / "results/tables/phase5_reference_knn_distances.csv").knn_mean_dist.values)
        lab_list = sorted(set(ref_lab))

    with h5py.File(ROOT / "data/processed/query_mapping.h5ad") as f:
        obs = read_elem(f["obs"])
    deliv = obs.loc[base.obs_names]
    rows, lat = [], {}
    for s in SCALES:
        with LOG.step(f"scale m x {s}"):
            scale = np.where(rowsum > 0, n_det * m1 * s / np.maximum(rowsum, 1e-12), 0.0)
            C = sp.diags(scale) @ T
            C.data = np.round(C.data)
            C.eliminate_zeros()
            X = sp.lil_matrix((base.n_obs, base.n_vars), dtype=np.float32)
            X[:, dst_idx] = C[:, src_idx]
            a = ad.AnnData(X=sp.csr_matrix(X), obs=base.obs.copy(), var=base.var.copy())
            a.layers["counts"] = a.X.copy()
            if s == 1.0:
                d = np.abs(a.X - base.X)
                LOG.write(f"m x 1 rebuild vs filbin.hvg.h5ad: max |diff| {d.max():.0f}, mean |diff| {d.mean():.4f}")
            SCANVI.prepare_query_anndata(a, str(REF_DIR))
            model = SCANVI.load_query_data(a, str(REF_DIR), accelerator="cpu", device="auto")
            model.is_trained_ = True  # zero-shot: no optimisation step; module holds the reference parameters
            model.module.eval()
            Z = model.get_latent_representation()
            P = model.predict(soft=True)
            lat[s] = Z
            nn_idx, dist = index.query(Z, k=K)
            votes = ref_lab[nn_idx]
            counts = np.stack([(votes == l).sum(1) for l in lab_list], 1)
            knn_lin = np.array(lab_list)[counts.argmax(1)]
            unc = np.searchsorted(ref_knn_sorted, dist.mean(1)) / len(ref_knn_sorted)
            pred = P.idxmax(1).values
            rows.append({"m_scale": s, "m": round(m1 * s, 3), "median_total_pseudocount_hvg": float(np.median(np.asarray(a.X.sum(1)).ravel())),
                         "median_knn_mean_dist": float(np.median(dist.mean(1))), "median_uncertainty_knn": float(np.median(unc)),
                         "frac_knn_malignant": float((knn_lin == "malignant").mean()), "frac_pred_malignant": float((pred == "malignant").mean()),
                         "knn_lineage_agree_with_deliverable": float((knn_lin == deliv.knn_lineage.astype(str).values).mean()),
                         "pred_lineage_agree_with_deliverable": float((pred == deliv.predicted_lineage.astype(str).values).mean()),
                         "consensus_malignant_to_knn_malignant": float((knn_lin[deliv.consensus_call.astype(str).values == "malignant"] == "malignant").mean()),
                         "consensus_malignant_to_pred_malignant": float((pred[deliv.consensus_call.astype(str).values == "malignant"] == "malignant").mean())})
            LOG.write(str({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()}))
    for r in rows:
        r["median_latent_shift_vs_m1"] = float(np.median(np.linalg.norm(lat[r["m_scale"]] - lat[1.0], axis=1)))
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/tables/phase5_filbin_m_sensitivity.csv", index=False)
    LOG.write("\n" + out.round(4).to_string() + f"\ndone | peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
