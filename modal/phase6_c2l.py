#!/usr/bin/env python
"""Phase 6 cell2location jobs for the DMG atlas, run on Modal (A10G; image im-gl8KES0WsH7PoFjf83hK6v,
cell2location 0.1.5, scvi-tools 1.4.2, torch 2.5.1+cu124).

  python phase6_c2l.py signatures --ref c2l_reference.h5ad --epochs 250
      -> out/signatures.csv (genes x classes, posterior mean of means_per_cluster_mu_fg), out/signatures_run.json,
         out/signatures_history.csv, out/regression_model/ (saved RegressionModel)
  python phase6_c2l.py deconvolve --vis visium_<cohort>.h5ad --sig signatures.csv --cohort <cohort>
         --target-epochs 20000 --budget-s 2400
      -> out/<section>_abundance.csv (q05 + mean cell abundances per spot, one row per spot),
         out/<cohort>_run.json, out/<cohort>_history.csv
Training length for deconvolution is adaptive: a 150-epoch probe on a throw-away model measures s/epoch, then
max_epochs = clip(budget_s / s_per_epoch, 2000, target). Everything written under ./out/ is harvested.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("out"); OUT.mkdir(exist_ok=True)


def gpu_info():
    import torch
    if not torch.cuda.is_available():
        return {"device": "cpu"}
    return {"device": torch.cuda.get_device_name(0),
            "max_mem_alloc_gb": round(torch.cuda.max_memory_allocated() / 1e9, 3)}


def peak_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024 / 1e9


def history_df(mod):
    h = mod.history
    if isinstance(h, dict):
        parts = [v.rename(columns=lambda c: c) for v in h.values()]
        return pd.concat(parts, axis=1)
    return pd.DataFrame(h)


def signatures(args):
    import anndata as ad
    import scvi
    from cell2location.models import RegressionModel
    scvi.settings.seed = 0
    t0 = time.time()
    ref = ad.read_h5ad(args.ref)
    for c in ("sample_id", "study", "assay", "c2l_class"):
        ref.obs[c] = ref.obs[c].astype(str).astype("category")
    print(f"reference {ref.shape}, classes {ref.obs.c2l_class.nunique()}, batches {ref.obs.sample_id.nunique()}", flush=True)
    RegressionModel.setup_anndata(adata=ref, batch_key="sample_id", labels_key="c2l_class",
                                  categorical_covariate_keys=["study", "assay"])
    mod = RegressionModel(ref)
    t1 = time.time()
    mod.train(max_epochs=args.epochs, batch_size=2500, train_size=1, lr=0.002, accelerator="gpu")
    train_s = time.time() - t1
    print(f"trained {args.epochs} epochs in {train_s:.0f}s", flush=True)
    hist = history_df(mod); hist.to_csv(OUT / "signatures_history.csv")
    ref = mod.export_posterior(ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500})
    names = list(ref.uns["mod"]["factor_names"])
    key = "means_per_cluster_mu_fg"
    if key in ref.varm.keys():
        sig = pd.DataFrame(np.asarray(ref.varm[key]), index=ref.var_names,
                           columns=[c for c in ref.varm[key].columns] if hasattr(ref.varm[key], "columns") else names)
        sig.columns = [c.replace(f"{key}_", "") for c in sig.columns]
    else:
        sig = ref.var[[f"{key}_{i}" for i in names]].copy(); sig.columns = names
    sig = sig[names]
    sig.to_csv(OUT / "signatures.csv")
    run = {"n_cells": int(ref.n_obs), "n_genes": int(ref.n_vars), "classes": names, "epochs": args.epochs,
           "train_s": round(train_s, 1), "total_s": round(time.time() - t0, 1),
           "final_elbo": float(hist.iloc[-1, 0]) if hist.shape[1] else None, **gpu_info(),
           "peak_rss_gb": round(peak_rss_gb(), 2), "cell2location": __import__("cell2location").__version__,
           "scvi": scvi.__version__}
    json.dump(run, open(OUT / "signatures_run.json", "w"), indent=1)
    print("RUN", json.dumps(run), flush=True)
    mod.save(str(OUT / "regression_model"), overwrite=True, save_anndata=False)


def deconvolve(args):
    import anndata as ad
    import scvi
    import torch
    from cell2location.models import Cell2location
    scvi.settings.seed = 0
    t0 = time.time()
    vis = ad.read_h5ad(args.vis)
    sig = pd.read_csv(args.sig, index_col=0)
    genes = [g for g in sig.index if g in set(vis.var_names)]
    vis = vis[:, genes].copy(); sig = sig.loc[genes]
    vis.obs["section_id"] = vis.obs["section_id"].astype(str).astype("category")
    print(f"visium {vis.shape}, sections {vis.obs.section_id.nunique()}, signature genes used {len(genes)}/{sig.shape[0]}", flush=True)
    Cell2location.setup_anndata(adata=vis, batch_key="section_id")
    kw = dict(cell_state_df=sig, N_cells_per_location=args.n_cells, detection_alpha=args.detection_alpha)

    # timing probe
    probe = Cell2location(vis, **kw)
    tp = time.time()
    probe.train(max_epochs=args.probe_epochs, batch_size=None, train_size=1, accelerator="gpu")
    s_per_epoch = (time.time() - tp) / args.probe_epochs
    del probe; torch.cuda.empty_cache()
    max_epochs = int(np.clip(args.budget_s / s_per_epoch, 2000, args.target_epochs))
    print(f"probe: {s_per_epoch*1000:.1f} ms/epoch -> max_epochs {max_epochs} (target {args.target_epochs}, budget {args.budget_s}s)", flush=True)

    mod = Cell2location(vis, **kw)
    t1 = time.time()
    mod.train(max_epochs=max_epochs, batch_size=None, train_size=1, accelerator="gpu")
    train_s = time.time() - t1
    hist = history_df(mod); hist.to_csv(OUT / f"{args.cohort}_history.csv")
    print(f"trained {max_epochs} epochs in {train_s:.0f}s; ELBO last {float(hist.iloc[-1, 0]):.1f}", flush=True)
    t2 = time.time()
    vis = mod.export_posterior(vis, sample_kwargs={"num_samples": args.num_samples, "batch_size": min(vis.n_obs, 4000)})
    post_s = time.time() - t2
    classes = list(sig.columns)
    q05 = pd.DataFrame(np.asarray(vis.obsm["q05_cell_abundance_w_sf"]), index=vis.obs_names,
                       columns=[f"q05_{c}" for c in classes])
    mean = pd.DataFrame(np.asarray(vis.obsm["means_cell_abundance_w_sf"]), index=vis.obs_names,
                        columns=[f"mean_{c}" for c in classes])
    keep_obs = [c for c in ("section_id", "cohort", "array_row", "array_col", "pxl_row", "pxl_col", "n_counts", "n_genes") if c in vis.obs]
    df = pd.concat([vis.obs[keep_obs], q05, mean], axis=1)
    for sec, sub in df.groupby("section_id", observed=True):
        sub.to_csv(OUT / f"{sec}_abundance.csv")
    run = {"cohort": args.cohort, "n_spots": int(vis.n_obs), "n_genes": int(vis.n_vars), "n_sections": int(vis.obs.section_id.nunique()),
           "sections": sorted(vis.obs.section_id.unique().tolist()), "classes": classes, "N_cells_per_location": args.n_cells,
           "detection_alpha": args.detection_alpha, "probe_ms_per_epoch": round(s_per_epoch * 1000, 1), "max_epochs": max_epochs,
           "target_epochs": args.target_epochs, "train_s": round(train_s, 1), "posterior_s": round(post_s, 1),
           "num_samples": args.num_samples, "total_s": round(time.time() - t0, 1), "final_elbo": float(hist.iloc[-1, 0]),
           **gpu_info(), "peak_rss_gb": round(peak_rss_gb(), 2)}
    json.dump(run, open(OUT / f"{args.cohort}_run.json", "w"), indent=1)
    print("RUN", json.dumps(run), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    s = sub.add_parser("signatures"); s.add_argument("--ref", required=True); s.add_argument("--epochs", type=int, default=250)
    d = sub.add_parser("deconvolve"); d.add_argument("--vis", required=True); d.add_argument("--sig", required=True)
    d.add_argument("--cohort", required=True); d.add_argument("--target-epochs", type=int, default=20000)
    d.add_argument("--budget-s", type=float, default=2400); d.add_argument("--probe-epochs", type=int, default=150)
    d.add_argument("--n-cells", type=float, default=5); d.add_argument("--detection-alpha", type=float, default=20)
    d.add_argument("--num-samples", type=int, default=1000)
    a = ap.parse_args()
    signatures(a) if a.mode == "signatures" else deconvolve(a)
