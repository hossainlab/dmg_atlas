#!/usr/bin/env python
"""Phase 4/5 audit fix IM-03 (2026-10-19) — held-out accuracy of the scANVI reference classifier and a patient-grouped bound.

1. Train / validation split of the saved reference model (models/scanvi/model.pt attr_dict train_indices_ / validation_indices_, a random
   90/10 split of cells made by scvi-tools) evaluated with the saved predicted_labels.txt against lineage_label (atlas_input_hvg.h5ad):
   accuracy, macro recall, per-class recall / precision -> results/tables/scanvi_validation_split.csv; the training history (per-epoch
   train / validation accuracy and classification loss) -> results/tables/scanvi_training_history.csv.
2. Patient-grouped 5-fold cross-validation (GroupKFold on patient_id, 52 patients) of two light classifiers trained on the UNSUPERVISED
   scVI latent (models/scvi_seed1, X_scvi) of the 220,309 labelled cells: multinomial logistic regression and a 128-128 MLP (the scANVI
   classifier's width/depth), each after per-fold standardisation; the same classifiers under a random cell-level 5-fold split for
   comparison -> results/tables/scanvi_patient_grouped_cv.csv (per fold) and scanvi_patient_grouped_cv_summary.csv (mean, sd).
   This is NOT a cross-validation of scANVI itself (retraining scANVI five times was out of the local CPU budget); it bounds how well a
   classifier of comparable capacity on the same information generalises to unseen patients. The scVI latent was fitted on all cells
   (unsupervised, no labels), so the patient hold-out is exact for the classifier but transductive for the embedding.
Log: logs/fix_phase5.log
"""
from __future__ import annotations

import os
import sys
import json
import time
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import h5py  # noqa: E402
import torch  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.neural_network import MLPClassifier  # noqa: E402
from sklearn.model_selection import GroupKFold, StratifiedKFold  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402
from sklearn.metrics import recall_score, precision_score, f1_score  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_phase5.log", context="11d_scanvi_validation")
T = ROOT / "results/tables"
CLASSES = ["malignant", "myeloid", "oligodendrocyte", "astrocyte", "neuron", "vascular", "T_NK"]


def main():
    ref = torch.load(ROOT / "models/scanvi/model.pt", map_location="cpu", weights_only=False)
    attr = ref["attr_dict"]
    tr_idx, va_idx = np.asarray(attr["train_indices_"]), np.asarray(attr["validation_indices_"])
    pred = np.array((ROOT / "models/scanvi/predicted_labels.txt").read_text().split("\n"))
    names = np.array((ROOT / "models/scanvi/obs_names.txt").read_text().split("\n"))
    with h5py.File(ROOT / "data/processed/atlas_input_hvg.h5ad") as f:
        oi = read_elem(f["obs"]["_index"]).astype(str)
        lab = read_elem(f["obs"]["lineage_label"]).astype(str)
        pid = read_elem(f["obs"]["patient_id"]).astype(str)
    assert (oi == names).all() and len(pred) == len(lab)
    labelled = lab != "Unknown"
    is_tr = np.zeros(len(lab), bool); is_tr[tr_idx] = True
    is_va = np.zeros(len(lab), bool); is_va[va_idx] = True

    # ------------------------------------------------------------------ 1. saved split
    with LOG.step("scANVI train / validation split accuracy"):
        rows = []
        for name, m in [("all_labelled", labelled), ("train_split", labelled & is_tr), ("validation_split", labelled & is_va)]:
            y, p = lab[m], pred[m]
            r = {"split": name, "n_labelled_cells": int(m.sum()), "n_patients": int(len(np.unique(pid[m]))), "accuracy": round(float((y == p).mean()), 4),
                 "macro_recall": round(recall_score(y, p, labels=CLASSES, average="macro"), 4),
                 "macro_f1": round(f1_score(y, p, labels=CLASSES, average="macro"), 4)}
            for c, rc, pc in zip(CLASSES, recall_score(y, p, labels=CLASSES, average=None), precision_score(y, p, labels=CLASSES, average=None, zero_division=0)):
                r[f"recall_{c}"] = round(float(rc), 4); r[f"precision_{c}"] = round(float(pc), 4); r[f"n_{c}"] = int((y == c).sum())
            rows.append(r)
        val = pd.DataFrame(rows)
        val["split_type"] = "random 90/10 cell split by scvi-tools (every validation patient also in training)"
        val.to_csv(T / "scanvi_validation_split.csv", index=False)
        LOG.write("\n" + val[["split", "n_labelled_cells", "accuracy", "macro_recall"] + [f"recall_{c}" for c in CLASSES]].to_string())
        hist = json.loads((ROOT / "models/scanvi/history.json").read_text())
        H = pd.DataFrame({k: hist[k] for k in ["train_accuracy", "validation_accuracy", "train_classification_loss", "validation_classification_loss",
                                                 "elbo_train", "elbo_validation"]})
        H.index.name = "epoch"
        H.to_csv(T / "scanvi_training_history.csv")
        LOG.write(f"history: validation_accuracy max {H.validation_accuracy.max():.4f} at epoch {int(H.validation_accuracy.idxmax())}, final {H.validation_accuracy.iloc[-1]:.4f}; "
                  f"validation_classification_loss min {H.validation_classification_loss.min():.4f} at epoch {int(H.validation_classification_loss.idxmin())}, "
                  f"final {H.validation_classification_loss.iloc[-1]:.4f}; train_accuracy final {H.train_accuracy.iloc[-1]:.4f}")

    # ------------------------------------------------------------------ 2. patient-grouped CV on X_scvi
    with LOG.step("patient-grouped 5-fold CV of light classifiers on X_scvi (labelled cells)"):
        X = np.load(ROOT / "models/scvi_seed1/latent_f16.npy").astype(np.float32)
        XL, yL, gL = X[labelled], lab[labelled], pid[labelled]
        scls = np.array(sorted(CLASSES)); yI = np.searchsorted(scls, yL)
        makers = {"logreg": lambda: LogisticRegression(max_iter=2000, C=1.0),
                  "mlp128x2": lambda: MLPClassifier(hidden_layer_sizes=(128, 128), max_iter=60, early_stopping=True, n_iter_no_change=5, random_state=0, batch_size=1024)}
        out = []
        for cname, mk in makers.items():
            for scheme, splitter, groups in [("patient_grouped_5fold", GroupKFold(n_splits=5), gL),
                                             ("random_cell_5fold", StratifiedKFold(n_splits=5, shuffle=True, random_state=0), None)]:
                t0 = time.time()
                for k, (tr, te) in enumerate(splitter.split(XL, yI, groups)):
                    sc = StandardScaler().fit(XL[tr])
                    p = scls[mk().fit(sc.transform(XL[tr]), yI[tr]).predict(sc.transform(XL[te]))]
                    r = {"classifier": cname, "scheme": scheme, "fold": k, "n_test": len(te), "n_test_patients": len(np.unique(gL[te])),
                         "accuracy": float((yL[te] == p).mean()), "macro_recall": float(recall_score(yL[te], p, labels=CLASSES, average="macro"))}
                    for c, rc in zip(CLASSES, recall_score(yL[te], p, labels=CLASSES, average=None)):
                        r[f"recall_{c}"] = float(rc); r[f"n_test_{c}"] = int((yL[te] == c).sum())
                    out.append(r)
                LOG.write(f"{cname} {scheme}: {time.time() - t0:.0f}s")
        cv = pd.DataFrame(out)
        cv.to_csv(T / "scanvi_patient_grouped_cv.csv", index=False)
        metrics = ["accuracy", "macro_recall"] + [f"recall_{c}" for c in CLASSES]
        summ = cv.groupby(["classifier", "scheme"])[metrics].agg(["mean", "std"]).round(4)
        summ.columns = [f"{a}_{b}" for a, b in summ.columns]
        summ = summ.reset_index()
        summ["note"] = ("bound on patient-level generalisation of a classifier on the unsupervised scVI latent (X_scvi, 30 dims, per-fold standardised); "
                        "not a cross-validation of scANVI")
        summ.to_csv(T / "scanvi_patient_grouped_cv_summary.csv", index=False)
        LOG.write("\n" + summ[["classifier", "scheme", "accuracy_mean", "accuracy_std", "macro_recall_mean"] + [f"recall_{c}_mean" for c in CLASSES]].to_string())
    LOG.write(f"done | peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
