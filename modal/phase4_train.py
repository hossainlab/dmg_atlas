#!/usr/bin/env python
"""Phase 4 training jobs for the DMG atlas, run on Modal (A10G).

Usage (inside the sandbox; data on the Volume mounted at /data, outputs under ./out):
  python phase4_train.py scvi      --data /data/atlas_input_hvg.h5ad --seeds 0 1 2
  python phase4_train.py controls  --data /data/atlas_input_hvg.h5ad            # PCA + Harmony
  python phase4_train.py scanvi    --data /data/atlas_input_hvg.h5ad --scvi-dir scvi_seed0 --max-epochs 25
  python phase4_train.py scvi      --data /data/atlas_input_malignant_hvg.h5ad --seeds 0 --tag scvi_malignant
Each run writes out/<tag>/{model.pt, latent_f16.npy, obs_names.txt, run.json} and appends to out/runs.jsonl.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import time
from pathlib import Path

import numpy as np

OUT = Path("out")
OUT.mkdir(exist_ok=True)
BATCH = "sample_id"
CAT_COV = ["study", "assay"]
CONT_COV = ["pct_mito"]


def gpu_info():
    import torch
    if not torch.cuda.is_available():
        return {"device": "cpu"}
    return {"device": torch.cuda.get_device_name(0),
            "max_mem_alloc_gb": round(torch.cuda.max_memory_allocated() / 1e9, 3),
            "max_mem_reserved_gb": round(torch.cuda.max_memory_reserved() / 1e9, 3)}


def peak_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024 / 1e9


def record(tag, d):
    d = {"tag": tag, **d}
    (OUT / tag).mkdir(parents=True, exist_ok=True)
    json.dump(d, open(OUT / tag / "run.json", "w"), indent=1, default=str)
    with open(OUT / "runs.jsonl", "a") as fh:
        fh.write(json.dumps(d, default=str) + "\n")
    print("RUN", json.dumps(d, default=str), flush=True)


def load(path):
    import anndata as ad
    t0 = time.time()
    adata = ad.read_h5ad(path)
    for c in CAT_COV + [BATCH]:
        adata.obs[c] = adata.obs[c].astype(str).astype("category")
    adata.obs[CONT_COV[0]] = adata.obs[CONT_COV[0]].astype(np.float32)
    print(f"loaded {path}: {adata.shape}, nnz {adata.layers['counts'].nnz}, {time.time() - t0:.0f}s", flush=True)
    return adata


def save_outputs(tag, model, adata, extra=None):
    d = OUT / tag
    d.mkdir(parents=True, exist_ok=True)
    model.save(str(d), overwrite=True, save_anndata=False)
    z = model.get_latent_representation(batch_size=4096).astype(np.float16)
    np.save(d / "latent_f16.npy", z)
    with open(d / "obs_names.txt", "w") as fh:
        fh.write("\n".join(adata.obs_names))
    hist = {k: v.iloc[:, 0].tolist() for k, v in model.history.items() if hasattr(v, "iloc")}
    json.dump(hist, open(d / "history.json", "w"))
    if extra:
        for k, v in extra.items():
            np.save(d / f"{k}.npy", v)
    # also persist models on the Volume so a later job can load them without re-upload
    vol_models = Path("/data/models")
    try:
        vol_models.mkdir(parents=True, exist_ok=True)
        if (vol_models / tag).exists():
            shutil.rmtree(vol_models / tag)
        shutil.copytree(d, vol_models / tag)
    except Exception as e:  # read-only mount etc.
        print("volume copy skipped:", e)
    return z


def run_scvi(args):
    import scvi
    import torch
    scvi.settings.num_threads = 8
    adata = load(args.data)
    n = adata.n_obs
    max_epochs = args.max_epochs or scvi.model._utils.get_max_epochs_heuristic(n)
    print(f"n_cells {n}; max_epochs heuristic {max_epochs}", flush=True)
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=BATCH,
                                  categorical_covariate_keys=CAT_COV, continuous_covariate_keys=CONT_COV)
    for seed in args.seeds:
        tag = args.tag or f"scvi_seed{seed}"
        scvi.settings.seed = seed
        torch.cuda.reset_peak_memory_stats()
        model = scvi.model.SCVI(adata, n_latent=args.n_latent, n_layers=2, gene_likelihood="nb")
        t0 = time.time()
        model.train(max_epochs=max_epochs, batch_size=1024, early_stopping=True,
                    early_stopping_monitor="elbo_validation", early_stopping_patience=args.patience,
                    early_stopping_min_delta=0.0, check_val_every_n_epoch=1, train_size=0.9,
                    accelerator="gpu", devices=1)
        wall = time.time() - t0
        epochs_run = len(model.history["elbo_train"])
        save_outputs(tag, model, adata)
        ev = model.history["elbo_validation"].iloc[:, 0]
        record(tag, {"model": "SCVI", "seed": seed, "n_cells": n, "n_genes": adata.n_vars, "max_epochs": int(max_epochs),
                     "epochs_run": int(epochs_run), "early_stopped": bool(epochs_run < max_epochs),
                     "batch_size": 1024, "n_latent": args.n_latent, "n_layers": 2, "gene_likelihood": "nb",
                     "batch_key": BATCH, "categorical_covariates": CAT_COV, "continuous_covariates": CONT_COV,
                     "elbo_validation_final": float(ev.iloc[-1]), "elbo_validation_best": float(ev.min()),
                     "wall_s": round(wall, 1), "s_per_epoch": round(wall / epochs_run, 2), "gpu": gpu_info(),
                     "peak_rss_gb": round(peak_rss_gb(), 2), "data": args.data})
        del model
        torch.cuda.empty_cache()


def run_scanvi(args):
    import scvi
    import torch
    scvi.settings.num_threads = 8
    adata = load(args.data)
    adata.obs["lineage_label"] = adata.obs["lineage_label"].astype(str).astype("category")
    scvi.settings.seed = args.seed
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=BATCH,
                                  categorical_covariate_keys=CAT_COV, continuous_covariate_keys=CONT_COV)
    vae = scvi.model.SCVI.load(args.scvi_dir, adata=adata)
    torch.cuda.reset_peak_memory_stats()
    lvae = scvi.model.SCANVI.from_scvi_model(vae, labels_key="lineage_label", unlabeled_category="Unknown")
    t0 = time.time()
    lvae.train(max_epochs=args.max_epochs, n_samples_per_label=args.n_samples_per_label, batch_size=1024, accelerator="gpu", devices=1,
               check_val_every_n_epoch=1)
    wall = time.time() - t0
    tag = args.scanvi_tag
    soft = lvae.predict(soft=True, batch_size=4096)
    pred = soft.idxmax(axis=1).values
    maxp = soft.max(axis=1).values.astype(np.float32)
    z = save_outputs(tag, lvae, adata, extra={"maxprob": maxp})
    np.save(OUT / tag / "soft_probs_f16.npy", soft.values.astype(np.float16))
    with open(OUT / tag / "soft_prob_columns.txt", "w") as fh:
        fh.write("\n".join(map(str, soft.columns)))
    with open(OUT / tag / "predicted_labels.txt", "w") as fh:
        fh.write("\n".join(map(str, pred)))
    lab = adata.obs["lineage_label"].astype(str).values
    labelled = lab != "Unknown"
    agree = float((pred[labelled] == lab[labelled]).mean())
    unk_conf = float((maxp[~labelled] >= 0.8).mean())
    epochs_run = len(lvae.history["elbo_train"])
    record(tag, {"model": "SCANVI", "seed": args.seed, "from": args.scvi_dir, "n_cells": adata.n_obs, "max_epochs": args.max_epochs,
                 "epochs_run": int(epochs_run), "n_samples_per_label": args.n_samples_per_label, "batch_size": 1024,
                 "agreement_labelled": agree, "n_labelled": int(labelled.sum()), "n_unknown": int((~labelled).sum()),
                 "frac_unknown_maxprob_ge_0.8": unk_conf,
                 "wall_s": round(wall, 1), "s_per_epoch": round(wall / epochs_run, 2), "gpu": gpu_info(),
                 "peak_rss_gb": round(peak_rss_gb(), 2), "data": args.data})


def run_controls(args):
    import scanpy as sc
    import torch
    sc.settings.n_jobs = 8
    adata = load(args.data)
    adata.X = adata.layers["counts"].copy()
    t0 = time.time()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=50, svd_solver="arpack", random_state=0)
    pca_wall = time.time() - t0
    tag = "pca"
    (OUT / tag).mkdir(exist_ok=True)
    np.save(OUT / tag / "latent_f16.npy", adata.obsm["X_pca"].astype(np.float16))
    np.save(OUT / tag / "variance_ratio.npy", adata.uns["pca"]["variance_ratio"])
    with open(OUT / tag / "obs_names.txt", "w") as fh:
        fh.write("\n".join(adata.obs_names))
    record(tag, {"model": "PCA", "n_comps": 50, "input": "log1p(normalize_total 1e4) on HVG counts, unscaled",
                 "wall_s": round(pca_wall, 1), "peak_rss_gb": round(peak_rss_gb(), 2)})
    # Harmony
    try:
        import harmonypy as hm
    except ImportError:
        import subprocess, sys
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "harmonypy"], check=True)
        import harmonypy as hm
    print("harmonypy", getattr(hm, "__version__", "?"), flush=True)
    t0 = time.time()
    torch.cuda.reset_peak_memory_stats()
    ho = hm.run_harmony(adata.obsm["X_pca"].astype(np.float32), adata.obs, BATCH, max_iter_harmony=20, random_state=0)
    Z = ho.Z_corr
    Z = Z.cpu().numpy() if hasattr(Z, "cpu") else np.asarray(Z)
    if Z.shape[0] != adata.n_obs:
        Z = Z.T
    wall = time.time() - t0
    tag = "harmony"
    (OUT / tag).mkdir(exist_ok=True)
    np.save(OUT / tag / "latent_f16.npy", Z.astype(np.float16))
    with open(OUT / tag / "obs_names.txt", "w") as fh:
        fh.write("\n".join(adata.obs_names))
    record(tag, {"model": "Harmony", "harmonypy_version": getattr(hm, "__version__", "?"), "batch": BATCH, "max_iter_harmony": 20,
                 "wall_s": round(wall, 1), "gpu": gpu_info(), "peak_rss_gb": round(peak_rss_gb(), 2)})


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["scvi", "scanvi", "controls"])
    p.add_argument("--data", required=True)
    p.add_argument("--seeds", type=int, nargs="*", default=[0])
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--tag", default=None)
    p.add_argument("--n-latent", type=int, default=30)
    p.add_argument("--max-epochs", type=int, default=None)
    p.add_argument("--patience", type=int, default=5)
    p.add_argument("--scvi-dir", default=None)
    p.add_argument("--n-samples-per-label", type=int, default=None, help="scANVI labelled subsample per class per epoch; None = all labelled cells")
    p.add_argument("--scanvi-tag", default="scanvi")
    args = p.parse_args()
    print("env:", {k: os.environ.get(k) for k in ("CUDA_VISIBLE_DEVICES",)}, flush=True)
    {"scvi": run_scvi, "scanvi": run_scanvi, "controls": run_controls}[args.mode](args)
