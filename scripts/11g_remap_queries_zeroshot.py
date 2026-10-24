#!/usr/bin/env python
"""Audit remediation (reference-labels track, 19 Oct 2026) — re-map the five query HVG files with the RETRAINED reference
scANVI (models/scanvi) as a zero-shot projection (IM-02: the reference encoder ignores batch covariates, so no query training
is performed; SCANVI.load_query_data() is called only to align the query to the model and is NOT trained).

  encode      per query: latent (encoder mean), soft class probabilities and argmax label from the reference classifier ->
              models/scanvi_query_<q>/ (same file layout as Phase 5; previous surgery-trained directories moved to
              models/scanvi_query_<q>_pre_audit). Verification: for Filbin every z_encoder / classifier tensor of the query model
              is compared with models/scanvi/model.pt (max |diff| must be 0: zero-shot = reference weights).
  uncertainty scripts/11_map_queries.py uncertainty() against the NEW reference latent (kNN index, 15-NN percentiles, UMAP
              refit); then query_mapping.h5ad gets the previous per-cell results back as *_pre_audit (original Phase 5 frozen
              classifier + kNN on the old latent) and *_fix1 (wave-1 reference classifier on the old latent) for cells present
              in both files, plus obs['liu_prep'] (IM-24). The previous file is kept as query_mapping.pre_audit.h5ad (wave-1 state).
Log: logs/fix_reflabels.log
"""
from __future__ import annotations

import json
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_reflabels.log", context="11g_remap_queries_zeroshot")
REF_DIR = ROOT / "models/scanvi"
QDIR = ROOT / "data/processed/queries"
QUERIES = ["filbin", "labelle", "liu", "braun", "siletti"]
QMAP = ROOT / "data/processed/query_mapping.h5ad"
CARRY = ["predicted_lineage", "scanvi_maxprob", "uncertainty_prob", "classifier_knn_agree", "knn_lineage", "knn_vote_frac", "knn_mean_dist", "uncertainty_knn"]


def model_dir(q):
    return ROOT / f"models/scanvi_query_{q}"


def encode():
    import torch
    import scvi
    from scvi.model import SCANVI
    torch.set_num_threads(8)
    scvi.settings.seed = 1
    scvi.settings.dl_num_workers = 0
    for q in QUERIES:
        with LOG.step(f"zero-shot encode {q}"):
            t0 = time.time()
            adata = ad.read_h5ad(QDIR / f"{q}.hvg.h5ad")
            SCANVI.prepare_query_anndata(adata, str(REF_DIR))
            model = SCANVI.load_query_data(adata, str(REF_DIR), accelerator="cpu", device="auto")
            model.module.eval()
            model.is_trained = True  # no training: the reference weights are used as they are
            latent = model.get_latent_representation(batch_size=4096)
            probs = model.predict(soft=True, batch_size=4096)
            pred = probs.idxmax(axis=1).values.astype(str)
            maxprob = probs.max(axis=1).values.astype(np.float32)
            wall = time.time() - t0
            verify = {}
            if q == "filbin":
                # zero-shot check: every encoder / classifier tensor of the query model equals the saved reference model.pt
                ref_sd = torch.load(REF_DIR / "model.pt", map_location="cpu", weights_only=False)["model_state_dict"]
                qsd = model.module.state_dict()
                diffs = {k: float((qsd[k].float() - ref_sd[k].float()).abs().max()) for k in ref_sd
                         if (k.startswith("z_encoder") or k.startswith("classifier")) and k in qsd and qsd[k].shape == ref_sd[k].shape}
                verify = {"encoder_classifier_tensors_compared": len(diffs), "max_abs_diff_vs_reference_model_pt": max(diffs.values())}
                LOG.write(f"filbin verification: {verify}")
            out = model_dir(q)
            pre = ROOT / f"models/scanvi_query_{q}_pre_audit"
            if out.exists() and not pre.exists():
                shutil.move(str(out), str(pre))
            elif out.exists():
                shutil.rmtree(out)
            out.mkdir(parents=True)
            np.save(out / "latent_f16.npy", latent.astype(np.float16))
            np.save(out / "soft_probs_f16.npy", probs.values.astype(np.float16))
            np.save(out / "maxprob.npy", maxprob)
            (out / "soft_prob_columns.txt").write_text("\n".join(probs.columns))
            (out / "predicted_labels.txt").write_text("\n".join(pred))
            (out / "obs_names.txt").write_text("\n".join(adata.obs_names))
            (out / "history.json").write_text(json.dumps({}))
            counts = pd.Series(pred).value_counts().to_dict()
            run = {"query": q, "n_cells": int(adata.n_obs), "mode": "zero-shot projection with the retrained reference scANVI (models/scanvi); no query training",
                   "max_epochs": 0, "epochs_run": 0, "batch_size": 4096, "s_per_epoch": 0.0, "trainable_params": 0,
                   "total_params": int(sum(p.numel() for p in model.module.parameters())), "wall_s": round(wall, 1), "device": "cpu", "seed": 1,
                   "median_maxprob": float(np.median(maxprob)), "frac_maxprob_ge_0.8": float(np.mean(maxprob >= 0.8)),
                   "predicted_counts": {k: int(v) for k, v in counts.items()}, "reference_model": "models/scanvi (audit retrain 2026-10-19)",
                   "query_file_md5": None, "scvi_tools": scvi.__version__, "torch": torch.__version__, **verify}
            (out / "run.json").write_text(json.dumps(run, indent=1, default=str))
            LOG.write(f"{q}: {adata.n_obs} cells, wall {wall:.0f}s, median maxprob {run['median_maxprob']:.3f}, predicted {counts}; peak RSS {peak_rss_gb():.2f} GB")
            del model, adata


def carry_over():
    """Add *_pre_audit (original Phase 5) and *_fix1 (wave-1) per-cell columns from the previous query_mapping file, plus liu_prep."""
    pre = QMAP.with_name("query_mapping.pre_audit.h5ad")
    with LOG.step("carry over previous per-cell results as *_pre_audit / *_fix1"):
        with h5py.File(pre) as f:
            old = read_elem(f["obs"])
            old_probs = f["obsm/scanvi_probs"][:]
            old_probs_pre = f["obsm/scanvi_probs_pre_audit"][:] if "scanvi_probs_pre_audit" in f["obsm"] else None
            old_cols = list(f["uns/phase5/scanvi_prob_columns"][:].astype(str))
        A = ad.read_h5ad(QMAP)
        common = A.obs_names.intersection(old.index)
        LOG.write(f"cells in new file {A.n_obs:,}; in previous {len(old):,}; common {len(common):,} "
                  f"(new-only {A.n_obs - len(common):,}: {A.obs.loc[~A.obs_names.isin(common), 'query'].value_counts().to_dict()}; "
                  f"dropped {len(old) - len(common):,}: {old.loc[~old.index.isin(common), 'query'].value_counts().to_dict()})")
        for c in CARRY:
            orig = f"{c}_pre_audit" if f"{c}_pre_audit" in old.columns else c  # original Phase 5 value
            A.obs[f"{c}_pre_audit"] = old[orig].reindex(A.obs_names).values
            if f"{c}_pre_audit" in old.columns:  # wave-1 value differs from the original only for the classifier columns
                A.obs[f"{c}_fix1"] = old[c].reindex(A.obs_names).values
        for c in A.obs.columns:
            if A.obs[c].dtype == object:
                A.obs[c] = A.obs[c].astype(str).replace({"nan": np.nan, "None": np.nan}).astype("category")
        pos = old.index.get_indexer(A.obs_names)
        P = np.full((A.n_obs, len(old_cols)), np.nan, dtype=np.float32)
        P[pos >= 0] = old_probs[pos[pos >= 0]]
        A.obsm["scanvi_probs_fix1"] = P
        if old_probs_pre is not None:
            P0 = np.full_like(P, np.nan)
            P0[pos >= 0] = old_probs_pre[pos[pos >= 0]]
            A.obsm["scanvi_probs_pre_audit"] = P0
        prep = A.obs["author_prep"].astype(object).map({"fresh_cells": "fresh cell", "frozen_nuclei": "frozen nucleus"}).fillna("unannotated")
        A.obs["liu_prep"] = pd.Categorical(np.where(A.obs["query"].astype(str) == "liu", prep.astype(str), "n/a"))
        A.uns["phase5"]["audit_fix_2026_10_19_reflabels"] = {
            "reference_model": "models/scanvi retrained with corrected consensus labels (patient-grouped validation, early stopping); zero-shot projection, no query training",
            "columns": "predicted_lineage/scanvi_maxprob/uncertainty_prob/knn_*/uncertainty_knn/classifier_knn_agree recomputed on the new latent; *_pre_audit = original Phase 5 "
                       "(frozen-classifier labels, kNN on the pre-audit latent); *_fix1 = wave-1 reference classifier on the pre-audit latent (IM-01); NaN for cells absent from the previous file",
            "liu_subset": f"Liu restricted to the {int((A.obs['query'] == 'liu').sum()):,} current qc_pass cells (DQ-04)",
            "previous_file": "data/processed/query_mapping.pre_audit.h5ad (wave-1 state)",
        }
        A.obs["classifier_knn_agree"] = A.obs["classifier_knn_agree"].astype(bool)
        A.write_h5ad(QMAP, compression="gzip")
        LOG.write(f"query_mapping.h5ad rewritten: obs {A.obs.shape[1]} cols, obsm {list(A.obsm)}; peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "encode":
        encode()
    elif cmd == "uncertainty":
        pre = QMAP.with_name("query_mapping.pre_audit.h5ad")
        if QMAP.exists() and not pre.exists():
            shutil.copy2(QMAP, pre)
            LOG.write("previous query_mapping.h5ad kept as query_mapping.pre_audit.h5ad (wave-1 state)")
        sys.path.insert(0, str(ROOT / "scripts"))
        import importlib
        m = importlib.import_module("11_map_queries")
        m.uncertainty()
        carry_over()
