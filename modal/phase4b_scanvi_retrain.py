#!/usr/bin/env python
"""Audit remediation (reference-labels track, 19 Oct 2026): retrain scANVI from models/scvi_seed1 with the
CORRECTED lineage_label (consensus_call regenerated in the CNV track) and a PATIENT-GROUPED validation split.

Run on Modal (A10G), data on the Volume mounted at /data:
  python phase4b_scanvi_retrain.py --data /data/atlas_input_hvg.h5ad --labels reference_labels_corrected.csv.gz \
      --scvi-dir /data/models/scvi_seed1 --pre-audit-scanvi /data/models/scanvi --max-epochs 25 --seed 1 --holdout-frac 0.10

Differences from the Phase 4 run (modal/phase4_train.py run_scanvi), all deliberate (finding IM-03):
  * obs['lineage_label'] is overwritten from --labels (corrected consensus calls); the pre-audit label is kept in
    obs['lineage_label_pre_audit'] for the change report.
  * validation set = ALL cells of ~10 % of PATIENTS (patient_id), chosen with the run seed such that every labelled
    class is present in the validation patients; passed through datasplitter_kwargs={'external_indexing': [...]}.
    Phase 4 used scvi-tools' random 90/10 CELL split (every validation patient was also a training patient).
  * early stopping on validation_classification_loss (patience 5, min_delta 0) and the best-epoch weights are
    restored at the end (SaveCheckpoint(load_best_on_end=True)). Phase 4 ran a fixed 25 epochs without early stopping.
  Unchanged: from_scvi_model(scvi_seed1), labels_key lineage_label, unlabeled 'Unknown', n_samples_per_label None
  (all labelled cells), batch_size 1024, max_epochs 25, seed 1.
Outputs under out/scanvi/: model.pt, latent_f16.npy, soft_probs_f16.npy, soft_prob_columns.txt, predicted_labels.txt,
  maxprob.npy, obs_names.txt, history.json, run.json, split.json, validation_metrics.json, validation_confusion.csv,
  per_class_recall.csv, pre_audit_comparison.json
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
import pandas as pd

OUT = Path("out")
OUT.mkdir(exist_ok=True)
BATCH = "sample_id"
CAT_COV = ["study", "assay"]
CONT_COV = ["pct_mito"]
UNLAB = "Unknown"


def gpu_info():
    import torch
    if not torch.cuda.is_available():
        return {"device": "cpu"}
    return {"device": torch.cuda.get_device_name(0),
            "max_mem_alloc_gb": round(torch.cuda.max_memory_allocated() / 1e9, 3),
            "max_mem_reserved_gb": round(torch.cuda.max_memory_reserved() / 1e9, 3)}


def peak_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024 / 1e9


def load(path, labels_path):
    import anndata as ad
    t0 = time.time()
    adata = ad.read_h5ad(path)
    for c in CAT_COV + [BATCH, "patient_id"]:
        adata.obs[c] = adata.obs[c].astype(str).astype("category")
    adata.obs[CONT_COV[0]] = adata.obs[CONT_COV[0]].astype(np.float32)
    lab = pd.read_csv(labels_path, index_col="cell_id", usecols=["cell_id", "lineage_label", "lineage_label_pre_audit", "okada_1734_pons_flag"])
    assert lab.index.equals(adata.obs_names) or set(lab.index) == set(adata.obs_names), "label file / training file cell mismatch"
    lab = lab.loc[adata.obs_names]
    old = adata.obs["lineage_label"].astype(str).values
    assert (old == lab["lineage_label_pre_audit"].astype(str).values).all(), "training file lineage_label != label file pre-audit column"
    adata.obs["lineage_label_pre_audit"] = pd.Categorical(old)
    adata.obs["lineage_label"] = pd.Categorical(lab["lineage_label"].astype(str).values)
    adata.obs["okada_1734_pons_flag"] = lab["okada_1734_pons_flag"].astype(bool).values
    changed = int((old != adata.obs["lineage_label"].astype(str).values).sum())
    print(f"loaded {path}: {adata.shape}, nnz {adata.layers['counts'].nnz}, {time.time() - t0:.0f}s; lineage_label changed for {changed} cells", flush=True)
    return adata


def patient_grouped_split(obs, frac, seed, max_tries=200):
    """Hold out whole patients until >= frac of LABELLED cells are held out; every labelled class must appear in the hold-out."""
    lab = obs["lineage_label"].astype(str).values
    pat = obs["patient_id"].astype(str).values
    labelled = lab != UNLAB
    classes = sorted(set(lab[labelled]))
    patients = np.array(sorted(set(pat)))
    n_lab_total = labelled.sum()
    for t in range(max_tries):
        rng = np.random.default_rng(seed * 1000 + t)
        order = rng.permutation(patients)
        held, n_held = [], 0
        for p in order:
            m = pat == p
            held.append(p)
            n_held += int((m & labelled).sum())
            if n_held >= frac * n_lab_total:
                break
        vm = np.isin(pat, held)
        present = {c: int(((lab == c) & vm).sum()) for c in classes}
        if all(v > 0 for v in present.values()) and n_held <= 0.16 * n_lab_total:
            return held, vm, present, t
    raise RuntimeError("no patient split found with all classes present")


def per_class(y_true, y_pred, classes):
    rows = []
    for c in classes:
        m = y_true == c
        tp = int((m & (y_pred == c)).sum())
        pp = int((y_pred == c).sum())
        rows.append({"class": c, "n_true": int(m.sum()), "n_pred": pp, "recall": tp / max(m.sum(), 1), "precision": tp / max(pp, 1)})
    return pd.DataFrame(rows)


def main(args):
    import scvi
    import torch
    from scvi.train import SaveCheckpoint
    scvi.settings.num_threads = 8
    scvi.settings.seed = args.seed
    adata = load(args.data, args.labels)
    lab = adata.obs["lineage_label"].astype(str).values
    labelled = lab != UNLAB
    classes = sorted(set(lab[labelled]))

    held, vm, present, tries = patient_grouped_split(adata.obs, args.holdout_frac, args.seed)
    val_idx = np.flatnonzero(vm)
    train_idx = np.flatnonzero(~vm)
    split = {"validation_patients": list(map(str, held)), "n_validation_cells": int(vm.sum()), "n_validation_labelled": int((vm & labelled).sum()),
             "frac_validation_labelled": float((vm & labelled).sum() / labelled.sum()), "n_train_cells": int((~vm).sum()),
             "n_train_labelled": int((~vm & labelled).sum()), "validation_class_counts": present, "seed_tries": tries,
             "validation_samples": sorted(set(adata.obs["sample_id"].astype(str).values[vm])),
             "validation_studies": adata.obs["study"].astype(str)[vm].value_counts().to_dict(),
             "rule": f"whole patients (patient_id) added in a seeded random order until >= {args.holdout_frac:.0%} of labelled cells are held out; "
                     "all 7 labelled classes must be present; unlabelled (Unknown) cells of those patients are also in the validation set"}
    print("SPLIT", json.dumps(split, default=str), flush=True)
    d = OUT / args.tag
    d.mkdir(parents=True, exist_ok=True)
    json.dump(split, open(d / "split.json", "w"), indent=1, default=str)
    pd.DataFrame({"cell_id": adata.obs_names, "set": np.where(vm, "validation", "train")}).to_csv(d / "split_cells.csv.gz", index=False)

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=BATCH, categorical_covariate_keys=CAT_COV, continuous_covariate_keys=CONT_COV)
    vae = scvi.model.SCVI.load(args.scvi_dir, adata=adata)
    torch.cuda.reset_peak_memory_stats()
    lvae = scvi.model.SCANVI.from_scvi_model(vae, labels_key="lineage_label", unlabeled_category=UNLAB)
    ckpt_dir = Path("ckpt_scanvi")
    ckpt = SaveCheckpoint(dirpath=str(ckpt_dir), monitor="validation_classification_loss", mode="min", load_best_on_end=True)
    t0 = time.time()
    lvae.train(max_epochs=args.max_epochs, n_samples_per_label=None, batch_size=1024, accelerator="gpu", devices=1,
               check_val_every_n_epoch=1,
               datasplitter_kwargs={"external_indexing": [train_idx, val_idx, np.array([], dtype=int)]},
               early_stopping=True, early_stopping_monitor="validation_classification_loss", early_stopping_mode="min",
               early_stopping_patience=args.patience, early_stopping_min_delta=0.0, callbacks=[ckpt])
    wall = time.time() - t0
    hist = {k: v.iloc[:, 0].tolist() for k, v in lvae.history.items() if hasattr(v, "iloc")}
    epochs_run = len(hist["elbo_train"])
    vcl = hist["validation_classification_loss"]
    best_epoch = int(np.argmin(vcl)) + 1
    print(f"epochs run {epochs_run}; best validation_classification_loss {min(vcl):.4f} at epoch {best_epoch}; "
          f"val acc by epoch {np.round(hist['validation_accuracy'], 4).tolist()}", flush=True)
    shutil.rmtree(ckpt_dir, ignore_errors=True)

    # ---- outputs (same layout as Phase 4) ----
    lvae.save(str(d), overwrite=True, save_anndata=False)
    z = lvae.get_latent_representation(batch_size=4096).astype(np.float16)
    np.save(d / "latent_f16.npy", z)
    (d / "obs_names.txt").write_text("\n".join(adata.obs_names))
    json.dump(hist, open(d / "history.json", "w"))
    soft = lvae.predict(soft=True, batch_size=4096)
    pred = soft.idxmax(axis=1).values.astype(str)
    maxp = soft.max(axis=1).values.astype(np.float32)
    np.save(d / "maxprob.npy", maxp)
    np.save(d / "soft_probs_f16.npy", soft.values.astype(np.float16))
    (d / "soft_prob_columns.txt").write_text("\n".join(map(str, soft.columns)))
    (d / "predicted_labels.txt").write_text("\n".join(pred))

    # ---- metrics ----
    tr_l = ~vm & labelled
    va_l = vm & labelled
    acc_train = float((pred[tr_l] == lab[tr_l]).mean())
    acc_val = float((pred[va_l] == lab[va_l]).mean())
    pc_val = per_class(lab[va_l], pred[va_l], classes)
    pc_train = per_class(lab[tr_l], pred[tr_l], classes)
    pc = pc_val.merge(pc_train, on="class", suffixes=("_validation", "_train"))
    pc.to_csv(d / "per_class_recall.csv", index=False)
    conf = pd.crosstab(pd.Series(lab[va_l], name="true"), pd.Series(pred[va_l], name="predicted"))
    conf.to_csv(d / "validation_confusion.csv")
    bal_val = float(pc_val["recall"].mean())
    # per held-out patient accuracy
    pat = adata.obs["patient_id"].astype(str).values
    per_pat = pd.DataFrame({"patient_id": pat[va_l], "correct": pred[va_l] == lab[va_l]}).groupby("patient_id")["correct"].agg(["size", "mean"])
    per_pat.to_csv(d / "validation_per_patient_accuracy.csv")
    unk_conf = float((maxp[~labelled] >= 0.8).mean())
    metrics = {"train_accuracy": acc_train, "validation_accuracy": acc_val, "validation_balanced_accuracy": bal_val,
               "validation_recall": dict(zip(pc_val["class"], pc_val["recall"].round(4))),
               "validation_precision": dict(zip(pc_val["class"], pc_val["precision"].round(4))),
               "train_recall": dict(zip(pc_train["class"], pc_train["recall"].round(4))),
               "best_epoch": best_epoch, "epochs_run": int(epochs_run), "validation_classification_loss_by_epoch": [round(v, 5) for v in vcl],
               "validation_accuracy_by_epoch": [round(v, 5) for v in hist["validation_accuracy"]],
               "frac_unknown_maxprob_ge_0.8": unk_conf, "predicted_counts_all": pd.Series(pred).value_counts().to_dict(),
               "predicted_counts_unknown": pd.Series(pred[~labelled]).value_counts().to_dict(),
               "per_patient_validation_accuracy": per_pat["mean"].round(4).to_dict()}
    print("METRICS", json.dumps(metrics, default=str), flush=True)
    json.dump(metrics, open(d / "validation_metrics.json", "w"), indent=1, default=str)

    # ---- pre-audit classifier on the same held-out patients (NOTE: it was trained on those patients) ----
    cmp = {}
    if args.pre_audit_scanvi:
        try:
            adata_old = adata.copy()
            adata_old.obs["lineage_label"] = adata_old.obs["lineage_label_pre_audit"]
            old = scvi.model.SCANVI.load(args.pre_audit_scanvi, adata=adata_old)
            soft_old = old.predict(soft=True, batch_size=4096)
            pred_old = soft_old.idxmax(axis=1).values.astype(str)
            saved = (Path(args.pre_audit_scanvi) / "predicted_labels.txt").read_text().split("\n")
            repro = float((np.array(saved) == pred_old).mean()) if len(saved) == len(pred_old) else None
            lab_old = adata.obs["lineage_label_pre_audit"].astype(str).values
            va_l_old = vm & (lab_old != UNLAB)
            both = va_l & va_l_old
            cmp = {"pre_audit_model_reproduces_saved_labels": repro,
                   "note": "the pre-audit model saw every validation patient during training (random cell split); its numbers are NOT held-out estimates",
                   "pre_audit_acc_vs_corrected_labels_heldout": float((pred_old[va_l] == lab[va_l]).mean()),
                   "pre_audit_acc_vs_pre_audit_labels_heldout": float((pred_old[va_l_old] == lab_old[va_l_old]).mean()),
                   "new_acc_vs_corrected_labels_heldout": acc_val,
                   "new_acc_vs_pre_audit_labels_heldout": float((pred[va_l_old] == lab_old[va_l_old]).mean()),
                   "new_vs_pre_audit_prediction_agreement_all_cells": float((pred == pred_old).mean()),
                   "new_vs_pre_audit_prediction_agreement_heldout": float((pred[vm] == pred_old[vm]).mean()),
                   "pre_audit_recall_vs_corrected_heldout": dict(zip(classes, per_class(lab[va_l], pred_old[va_l], classes)["recall"].round(4))),
                   "n_heldout_labelled_both": int(both.sum()),
                   "crosstab_new_vs_pre_audit_all": pd.crosstab(pd.Series(pred_old, name="pre_audit"), pd.Series(pred, name="new")).to_dict()}
            (d / "pre_audit_predicted_labels.txt").write_text("\n".join(pred_old))
            del old, adata_old
        except Exception as e:
            cmp = {"error": f"{type(e).__name__}: {str(e)[:300]}"}
        print("PRE_AUDIT_CMP", json.dumps(cmp, default=str), flush=True)
        json.dump(cmp, open(d / "pre_audit_comparison.json", "w"), indent=1, default=str)

    run = {"tag": args.tag, "model": "SCANVI", "seed": args.seed, "from": args.scvi_dir, "n_cells": adata.n_obs, "max_epochs": args.max_epochs,
           "epochs_run": int(epochs_run), "best_epoch": best_epoch, "early_stopping": {"monitor": "validation_classification_loss", "patience": args.patience, "min_delta": 0.0, "load_best_on_end": True},
           "split": "patient-grouped (see split.json)", "n_samples_per_label": None, "batch_size": 1024,
           "agreement_labelled_all": float((pred[labelled] == lab[labelled]).mean()), "train_accuracy": acc_train, "validation_accuracy": acc_val,
           "n_labelled": int(labelled.sum()), "n_unknown": int((~labelled).sum()), "frac_unknown_maxprob_ge_0.8": unk_conf,
           "labels": "corrected consensus_call (CNV track, 19 Oct 2026) via reference_labels_corrected.csv.gz",
           "wall_s": round(wall, 1), "s_per_epoch": round(wall / epochs_run, 2), "gpu": gpu_info(), "peak_rss_gb": round(peak_rss_gb(), 2),
           "data": args.data, "scvi_tools": scvi.__version__, "torch": torch.__version__}
    json.dump(run, open(d / "run.json", "w"), indent=1, default=str)
    with open(OUT / "runs.jsonl", "a") as fh:
        fh.write(json.dumps(run, default=str) + "\n")
    print("RUN", json.dumps(run, default=str), flush=True)
    # persist on the Volume too (pre-audit copy kept)
    try:
        vol = Path("/data/models")
        if (vol / "scanvi").exists() and not (vol / "scanvi_pre_audit").exists():
            shutil.copytree(vol / "scanvi", vol / "scanvi_pre_audit")
        if (vol / args.tag).exists():
            shutil.rmtree(vol / args.tag)
        shutil.copytree(d, vol / args.tag)
        print("volume copy done", flush=True)
    except Exception as e:
        print("volume copy skipped:", e, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--data", required=True)
    p.add_argument("--labels", required=True)
    p.add_argument("--scvi-dir", required=True)
    p.add_argument("--pre-audit-scanvi", default=None)
    p.add_argument("--max-epochs", type=int, default=25)
    p.add_argument("--patience", type=int, default=5)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--holdout-frac", type=float, default=0.10)
    p.add_argument("--tag", default="scanvi")
    args = p.parse_args()
    print("env:", {k: os.environ.get(k) for k in ("CUDA_VISIBLE_DEVICES",)}, flush=True)
    main(args)
