#!/usr/bin/env python
"""Phase 5 / step 11 — scArches-style scANVI query mapping.

Sub-commands
  train --query <name> [--max-epochs 100] [--batch-size auto]
        SCANVI.load_query_data(query, models/scanvi) with the scArches defaults (only the new batch /
        covariate category weights are trainable: freeze_expression, freeze_decoder_first_layer,
        freeze_batchnorm_encoder, freeze_classifier all True), train, predict labels + soft
        probabilities, latent. Writes models/scanvi_query_<name>/ (model.pt, latent_f16.npy,
        predicted_labels.txt, soft_probs_f16.npy, maxprob.npy, obs_names.txt, run.json, history.json).
  uncertainty
        kNN uncertainty of every mapped query cell against the reference X_scanvi (k = 15, pynndescent):
        mean distance to the 15 nearest reference cells expressed as the percentile of the reference's own
        15-NN mean-distance distribution (50,000-cell reference subsample, self excluded) — the
        scArches / HLCA convention. Also projects query cells onto a UMAP refitted on a 100,000-cell
        reference subsample of X_scanvi (umap-learn; the stored scanpy UMAP has no transform model).
        Writes data/processed/query_mapping.h5ad and results/tables/mapping_uncertainty.csv.
Log: logs/phase5_mapping.log
"""
from __future__ import annotations

import sys
import json
import time
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402  (thread caps first)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import anndata as ad  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase5_mapping.log", context="11_map_queries")
REF_DIR = ROOT / "models/scanvi"
QDIR = ROOT / "data/processed/queries"
QUERIES = ["filbin", "labelle", "liu", "braun", "siletti"]
K = 15


def model_dir(q: str) -> Path:
    return ROOT / f"models/scanvi_query_{q}"


# ------------------------------------------------------------------------------------------------------ train
def train(q: str, max_epochs: int, batch_size: str, seed: int = 1):
    import torch
    import scvi
    from scvi.model import SCANVI
    torch.set_num_threads(int(__import__("os").environ.get("OMP_NUM_THREADS", "8")))
    scvi.settings.seed = seed
    scvi.settings.dl_num_workers = 0
    adata = ad.read_h5ad(QDIR / f"{q}.hvg.h5ad")
    n = adata.n_obs
    bs = int(batch_size) if batch_size != "auto" else int(min(1024, max(128, 2 ** int(np.log2(max(n // 100, 1))))))
    with LOG.step(f"train {q} (n={n}, batch_size={bs}, max_epochs={max_epochs})"):
        t0 = time.time()
        SCANVI.prepare_query_anndata(adata, str(REF_DIR))  # asserts var alignment (already aligned by step 10)
        model = SCANVI.load_query_data(adata, str(REF_DIR), accelerator="cpu", device="auto")
        trainable = sum(p.numel() for p in model.module.parameters() if p.requires_grad)
        total = sum(p.numel() for p in model.module.parameters())
        LOG.write(f"{q}: trainable params {trainable} / {total}")
        model.train(max_epochs=max_epochs, batch_size=bs, accelerator="cpu", devices=1,
                    check_val_every_n_epoch=10, plan_kwargs={"weight_decay": 0.0}, early_stopping=False,
                    n_samples_per_label=None)
        wall_train = time.time() - t0
        epochs_run = len(model.history["train_loss_epoch"]) if "train_loss_epoch" in model.history else max_epochs
        latent = model.get_latent_representation()
        probs = model.predict(soft=True)
        pred = probs.idxmax(axis=1).values
        maxprob = probs.max(axis=1).values
        out = model_dir(q)
        out.mkdir(parents=True, exist_ok=True)
        model.save(str(out), overwrite=True, save_anndata=False)
        np.save(out / "latent_f16.npy", latent.astype(np.float16))
        np.save(out / "soft_probs_f16.npy", probs.values.astype(np.float16))
        np.save(out / "maxprob.npy", maxprob.astype(np.float32))
        (out / "soft_prob_columns.txt").write_text("\n".join(probs.columns))
        (out / "predicted_labels.txt").write_text("\n".join(pred))
        (out / "obs_names.txt").write_text("\n".join(adata.obs_names))
        hist = {k: v.iloc[:, 0].tolist() for k, v in model.history.items()}
        (out / "history.json").write_text(json.dumps(hist))
        pred_counts = pd.Series(pred).value_counts().to_dict()
        run = {"query": q, "n_cells": n, "max_epochs": max_epochs, "epochs_run": int(epochs_run), "batch_size": bs,
               "seed": seed, "wall_s": round(wall_train, 1), "s_per_epoch": round(wall_train / max(epochs_run, 1), 2),
               "device": "cpu", "trainable_params": int(trainable), "total_params": int(total),
               "freeze": {"freeze_expression": True, "freeze_decoder_first_layer": True, "freeze_batchnorm_encoder": True,
                          "freeze_batchnorm_decoder": False, "freeze_classifier": True, "freeze_dropout": False},
               "plan_kwargs": {"weight_decay": 0.0}, "check_val_every_n_epoch": 10,
               "median_maxprob": float(np.median(maxprob)), "frac_maxprob_ge_0.8": float(np.mean(maxprob >= 0.8)),
               "predicted_counts": {k: int(v) for k, v in pred_counts.items()},
               "train_elbo_first_last": [hist.get("elbo_train", [None])[0], hist.get("elbo_train", [None])[-1]],
               "peak_rss_gb": round(peak_rss_gb(), 2), "scvi_tools": scvi.__version__, "torch": torch.__version__}
        (out / "run.json").write_text(json.dumps(run, indent=1, default=str))
        LOG.write(f"{q}: epochs {epochs_run}, wall {wall_train:.0f}s ({run['s_per_epoch']} s/epoch), median maxprob "
                  f"{run['median_maxprob']:.3f}, predicted {pred_counts}")


# ------------------------------------------------------------------------------------------- uncertainty / UMAP
def load_ref_latent():
    Z = np.load(REF_DIR / "latent_f16.npy").astype(np.float32)
    names = (REF_DIR / "obs_names.txt").read_text().split("\n")
    assert len(names) == Z.shape[0]
    return Z, np.array(names)


def load_query(q: str):
    d = model_dir(q)
    a = ad.read_h5ad(QDIR / f"{q}.hvg.h5ad")
    names = (d / "obs_names.txt").read_text().split("\n")
    assert list(a.obs_names) == names, q
    a.obsm["X_scanvi"] = np.load(d / "latent_f16.npy").astype(np.float32)
    cols = (d / "soft_prob_columns.txt").read_text().split("\n")
    probs = np.load(d / "soft_probs_f16.npy").astype(np.float32)
    a.obsm["scanvi_probs"] = pd.DataFrame(probs, index=a.obs_names, columns=cols)
    a.obs["predicted_lineage"] = pd.Categorical((d / "predicted_labels.txt").read_text().split("\n"))
    a.obs["scanvi_maxprob"] = np.load(d / "maxprob.npy").astype(np.float32)
    a.obs["uncertainty_prob"] = 1.0 - a.obs["scanvi_maxprob"]
    run = json.loads((d / "run.json").read_text())
    a.uns["phase5_mapping"] = {k: run[k] for k in ["epochs_run", "wall_s", "batch_size", "s_per_epoch", "trainable_params"]}
    return a


def uncertainty(seed: int = 0):
    import pynndescent
    import umap
    rng = np.random.default_rng(seed)
    Zref, ref_names = load_ref_latent()
    with LOG.step("build reference kNN index (pynndescent, X_scanvi, 311,618 cells)"):
        index = pynndescent.NNDescent(Zref, n_neighbors=30, metric="euclidean", random_state=seed, n_jobs=8)
        index.prepare()
    with LOG.step("reference 15-NN distance distribution (50,000-cell subsample)"):
        sub = rng.choice(Zref.shape[0], 50_000, replace=False)
        _, d = index.query(Zref[sub], k=K + 1)
        ref_knn = d[:, 1:].mean(1)  # drop self
        ref_knn_sorted = np.sort(ref_knn)
        pd.DataFrame({"cell_id": ref_names[sub], "knn_mean_dist": ref_knn}).to_csv(
            ROOT / "results/tables/phase5_reference_knn_distances.csv", index=False)
    with LOG.step("fit UMAP on 100,000-cell reference subsample (umap-learn, n_neighbors 30, min_dist 0.5, seed 0)"):
        usub = rng.choice(Zref.shape[0], 100_000, replace=False)
        reducer = umap.UMAP(n_neighbors=30, min_dist=0.5, n_components=2, random_state=seed, n_jobs=1)
        ref_umap = reducer.fit_transform(Zref[usub])
        pd.DataFrame({"cell_id": ref_names[usub], "umap_1": ref_umap[:, 0], "umap_2": ref_umap[:, 1]}).to_csv(
            ROOT / "results/tables/phase5_umap_reference_fit.csv", index=False)
    ref_lab = pd.read_csv(ROOT / "results/tables/atlas_core_obs_labels.csv", index_col=0).loc[ref_names, "scanvi_label"].astype(str).values
    calls = pd.read_csv(ROOT / "data/processed/malignant_calls.csv", usecols=["cell_id", "cnv_call", "consensus_call",
                                                                             "consensus_celltype_coarse", "marker_call"],
                        index_col="cell_id")
    parts, unc_rows = [], []
    for q in QUERIES:
        with LOG.step(f"uncertainty + UMAP projection {q}"):
            a = load_query(q)
            nn_idx, d = index.query(a.obsm["X_scanvi"], k=K)
            a.obs["knn_mean_dist"] = d.mean(1).astype(np.float32)
            votes = ref_lab[nn_idx]  # kNN label transfer (HLCA-style): majority scANVI label of the 15 nearest reference cells
            lab_list = sorted(set(ref_lab))
            counts = np.stack([(votes == l).sum(1) for l in lab_list], 1)
            a.obs["knn_lineage"] = pd.Categorical(np.array(lab_list)[counts.argmax(1)])
            a.obs["knn_vote_frac"] = (counts.max(1) / K).astype(np.float32)
            a.obs["classifier_knn_agree"] = (a.obs["knn_lineage"].astype(str).values == a.obs["predicted_lineage"].astype(str).values)
            a.obs["uncertainty_knn"] = (np.searchsorted(ref_knn_sorted, a.obs["knn_mean_dist"].values) / len(ref_knn_sorted)).astype(np.float32)
            a.obsm["X_umap_ref"] = reducer.transform(a.obsm["X_scanvi"]).astype(np.float32)
            for c in ["cnv_call", "consensus_call", "consensus_celltype_coarse", "marker_call"]:
                a.obs[c] = calls[c].reindex(a.obs_names).astype(object).values if q in ("filbin", "labelle", "liu") else np.nan
            a.obs["query"] = q
            parts.append(a)
            for s, g in a.obs.groupby("sample_id", observed=True):
                unc_rows.append({"query": q, "sample_id": s, "n_cells": len(g),
                                 "median_uncertainty_prob": g.uncertainty_prob.median(), "p90_uncertainty_prob": g.uncertainty_prob.quantile(0.9),
                                 "frac_uncertainty_prob_gt_0.5": (g.uncertainty_prob > 0.5).mean(),
                                 "median_uncertainty_knn": g.uncertainty_knn.median(), "p90_uncertainty_knn": g.uncertainty_knn.quantile(0.9),
                                 "frac_uncertainty_knn_gt_0.5": (g.uncertainty_knn > 0.5).mean(),
                                 "median_knn_mean_dist": g.knn_mean_dist.median(), "median_maxprob": g.scanvi_maxprob.median(),
                                 "frac_maxprob_ge_0.8": (g.scanvi_maxprob >= 0.8).mean(),
                                 "frac_classifier_knn_agree": g.classifier_knn_agree.mean()})
    with LOG.step("assemble data/processed/query_mapping.h5ad"):
        for a in parts:  # harmonise obs columns
            for c in ["author_annotation", "author_prep", "author_location", "author_age", "author_broad_annot",
                      "author_detailed_annot", "author_CellClass", "author_cell_type", "author_Region", "author_Subregion",
                      "author_Age_pcw", "author_supercluster_term", "author_roi", "author_tissue", "patient_id", "site",
                      "age_group", "h3_variant", "donor_id", "labelle_platform", "labelle_table", "original_cell_id", "n_counts", "n_genes"]:
                if c not in a.obs:
                    a.obs[c] = np.nan
        keep = ["query", "study", "sample_id", "patient_id", "donor_id", "assay", "site", "age_group", "h3_variant", "original_cell_id",
                "n_counts", "n_genes", "pct_mito", "query_n_counts_hvg", "query_n_genes_hvg", "labelle_platform", "labelle_table",
                "author_annotation", "author_prep", "author_location", "author_age", "author_broad_annot", "author_detailed_annot",
                "author_CellClass", "author_cell_type", "author_Region", "author_Subregion", "author_Age_pcw",
                "author_supercluster_term", "author_roi", "author_tissue",
                "cnv_call", "consensus_call", "consensus_celltype_coarse", "marker_call",
                "predicted_lineage", "scanvi_maxprob", "uncertainty_prob", "knn_mean_dist", "uncertainty_knn",
                "knn_lineage", "knn_vote_frac", "classifier_knn_agree"]
        for a in parts:
            a.obs = a.obs[keep].copy()
            for c in a.obs.columns:
                if a.obs[c].dtype == object or isinstance(a.obs[c].dtype, pd.CategoricalDtype):
                    a.obs[c] = a.obs[c].astype(str).replace({"nan": np.nan, "None": np.nan}).astype("category")
            a.obsm["scanvi_probs"] = a.obsm["scanvi_probs"].values.astype(np.float32)
            del a.layers["counts"]
            a.uns = {}
        cols = (model_dir("filbin") / "soft_prob_columns.txt").read_text().split("\n")
        M = ad.concat(parts, join="outer", merge="same", label=None, index_unique=None)
        M.var = parts[0].var.copy()
        M.layers["counts"] = M.X
        M.uns["phase5"] = {"scanvi_prob_columns": cols, "reference_model": "models/scanvi", "k": K,
                           "knn_lineage": "majority scANVI label (atlas_core scanvi_label) among the 15 nearest reference cells in X_scanvi; "
                                          "robust alternative to the frozen classifier for off-manifold queries (see notes)",
                           "uncertainty_knn": "percentile of mean 15-NN distance in the reference's own 15-NN distribution (50k subsample)",
                           "X_umap_ref": "umap-learn transform of query X_scanvi with a UMAP fitted on a 100,000-cell reference subsample "
                                         "(results/tables/phase5_umap_reference_fit.csv); coordinates differ from atlas_core X_umap",
                           "per_query_training": {q: json.loads((model_dir(q) / "run.json").read_text()) for q in QUERIES},
                           "counts": "layers['counts'] = X: reference-HVG counts (UMI for Braun/Siletti, round(RSEM) for Liu, "
                                     "TPM-derived pseudo-counts for Filbin/LaBelle; see docs/phase5_mapping_notes.md)"}
        out = ROOT / "data/processed/query_mapping.h5ad"
        M.write_h5ad(out, compression="gzip")
        LOG.write(f"wrote {out}: {M.n_obs} cells x {M.n_vars} genes; obs {M.obs.shape[1]} cols; obsm {list(M.obsm)}")
    unc = pd.DataFrame(unc_rows)
    unc.to_csv(ROOT / "results/tables/mapping_uncertainty.csv", index=False)
    LOG.write(f"wrote results/tables/mapping_uncertainty.csv ({len(unc)} sample rows) | peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["train", "uncertainty"])
    ap.add_argument("--query", choices=QUERIES)
    ap.add_argument("--max-epochs", type=int, default=100)
    ap.add_argument("--batch-size", default="auto")
    args = ap.parse_args()
    if args.cmd == "train":
        train(args.query, args.max_epochs, args.batch_size)
    else:
        uncertainty()
