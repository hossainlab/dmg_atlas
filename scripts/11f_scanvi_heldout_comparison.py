#!/usr/bin/env python
"""Audit remediation (reference-labels track) — compare the retrained scANVI (models/scanvi, patient-grouped validation,
early stopping) with the pre-audit classifier (models/scanvi_pre_audit, random cell split, 25 fixed epochs) on the SAME
held-out patients (models/scanvi/split.json). Caveat recorded in the output: the pre-audit model was trained on every
held-out patient, so its numbers are training-set numbers, not held-out estimates.
Writes results/tables/reflabels/scanvi_heldout_comparison.csv, scanvi_heldout_per_class.csv, scanvi_label_change_crosstab.csv
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_reflabels.log", context="11f_scanvi_heldout_comparison")
OUT = ROOT / "results/tables/reflabels"
NEW, OLD = ROOT / "models/scanvi", ROOT / "models/scanvi_pre_audit"


def read_pred(d):
    names = pd.Index((d / "obs_names.txt").read_text().split("\n"))
    pred = pd.Series((d / "predicted_labels.txt").read_text().split("\n"), index=names)
    maxp = pd.Series(np.load(d / "maxprob.npy"), index=names)
    return pred, maxp


def recall_table(y, p, classes):
    return pd.Series({c: float((p[y == c] == c).mean()) if (y == c).any() else np.nan for c in classes})


def main():
    lab = pd.read_csv(OUT / "reference_labels_corrected.csv.gz", index_col="cell_id")
    split = pd.read_csv(NEW / "split_cells.csv.gz", index_col="cell_id")["set"]
    pn, mn = read_pred(NEW)
    po, mo = read_pred(OLD)
    idx = lab.index
    pn, mn, po, mo, split = pn.loc[idx], mn.loc[idx], po.loc[idx], mo.loc[idx], split.loc[idx]
    val = (split == "validation").values
    new_lab = lab["lineage_label"].astype(str).values
    old_lab = lab["lineage_label_pre_audit"].astype(str).values
    classes = sorted(set(new_lab) - {"Unknown"})
    rows = []
    for name, pred in [("new_scanvi", pn.values), ("pre_audit_scanvi", po.values)]:
        for labname, y in [("corrected", new_lab), ("pre_audit", old_lab)]:
            for setname, m in [("heldout_patients", val), ("train_patients", ~val), ("all", np.ones(len(y), bool))]:
                mm = m & (y != "Unknown")
                rows.append({"classifier": name, "labels": labname, "cells": setname, "n_labelled": int(mm.sum()),
                             "accuracy": float((pred[mm] == y[mm]).mean()),
                             "balanced_accuracy": float(recall_table(y[mm], pred[mm], classes).mean())})
    cmp = pd.DataFrame(rows)
    cmp["note"] = np.where(cmp.classifier.eq("pre_audit_scanvi") & cmp.cells.eq("heldout_patients"),
                           "pre-audit model trained on these patients (random cell split): NOT a held-out estimate", "")
    cmp.to_csv(OUT / "scanvi_heldout_comparison.csv", index=False)
    LOG.write("comparison:\n" + cmp.to_string())
    mm = val & (new_lab != "Unknown")
    pc = pd.DataFrame({"new_scanvi_recall_heldout_vs_corrected": recall_table(new_lab[mm], pn.values[mm], classes),
                       "pre_audit_scanvi_recall_heldout_vs_corrected": recall_table(new_lab[mm], po.values[mm], classes),
                       "n_heldout_labelled": pd.Series({c: int((new_lab[mm] == c).sum()) for c in classes})})
    mm2 = val & (old_lab != "Unknown")
    pc["pre_audit_scanvi_recall_heldout_vs_pre_audit_labels"] = recall_table(old_lab[mm2], po.values[mm2], classes)
    pc.index.name = "class"
    pc.to_csv(OUT / "scanvi_heldout_per_class.csv")
    LOG.write("per-class (held-out patients):\n" + pc.round(4).to_string())
    ct = pd.crosstab(pd.Series(po.values, name="scanvi_label_pre_audit"), pd.Series(pn.values, name="scanvi_label_new"), margins=True)
    ct.to_csv(OUT / "scanvi_label_change_crosstab.csv")
    agree = float((po.values == pn.values).mean())
    LOG.write(f"scanvi_label agreement new vs pre-audit over all 311,618 cells: {agree:.4f}; changed {(po.values != pn.values).sum():,}\n" + ct.to_string())
    LOG.write(f"malignant scanvi_label: pre-audit {int((po.values == 'malignant').sum()):,} -> new {int((pn.values == 'malignant').sum()):,}; "
              f"median maxprob pre {float(mo.median()):.4f} new {float(mn.median()):.4f}; frac maxprob>=0.8 pre {float((mo >= 0.8).mean()):.4f} new {float((mn >= 0.8).mean()):.4f}")
    json.dump({"agreement_all": agree, "n_changed": int((po.values != pn.values).sum())}, open(OUT / "scanvi_label_change_summary.json", "w"))


if __name__ == "__main__":
    main()
