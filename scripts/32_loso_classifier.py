#!/usr/bin/env python
"""Phase 9 quality assessment / step 32 — score the leave-one-publication-out (LOPO) reference-classifier test.

Inputs : models/loso_scanvi/<publication>/{predictions.csv.gz, split.json, run.json} from modal/phase9_loso_scanvi.py
         (scVI + scANVI refitted without the publication, the publication mapped by scArches surgery, its labels hidden),
         results/tables/reflabels/scanvi_heldout_comparison.csv (held-out-PATIENT baseline of models/scanvi).
Truth = corrected lineage_label (reference_labels_corrected.csv.gz); 'Unknown' (ambiguous consensus) cells are not scored but
their predicted-label distribution is reported. Balanced accuracy = mean recall over classes with >= MIN_CLASS true cells in
the held-out publication (classes below that are listed, not averaged).
Outputs (results/tables/quality/, skipped if present unless --force):
         loso_classifier_summary.csv      one row per fold + the held-out-patient baseline: accuracy, balanced accuracy, macro F1,
                                          malignant precision / recall, coverage and accuracy at maxprob >= 0.9, per-patient range
         loso_classifier_per_class.csv    fold x class: n_true, n_pred, recall, precision, F1
         loso_classifier_per_patient.csv  fold x patient: n labelled, accuracy
         loso_classifier_confusion.csv    fold x true x predicted counts (labelled cells)
         loso_classifier_unknown.csv      fold x predicted class for ambiguous ('Unknown') cells: n, fraction, median maxprob
         loso_classifier_inputs.json      md5 of every input
Run: nice -n 10 python scripts/32_loso_classifier.py [--force]
Log: logs/phase9_quality.log
"""
import sys, json, hashlib, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd  # noqa: E402
from dmg_atlas import RESULTS, LOGS  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
log = U.StepLogger(LOGS / "phase9_quality.log", context="32_loso_classifier")
ROOT = RESULTS.parent; MODELS = ROOT / "models/loso_scanvi"; TABLES = RESULTS / "tables"; OUT = TABLES / "quality"; OUT.mkdir(parents=True, exist_ok=True)
BASE = TABLES / "reflabels/scanvi_heldout_comparison.csv"
FOLDS = ["okada", "jessa2022", "collot", "zang", "desisto"]
UNLAB, MIN_CLASS, MAXP = "Unknown", 20, 0.9
OUTS = [OUT / f"loso_classifier_{n}" for n in ("summary.csv", "per_class.csv", "per_patient.csv", "confusion.csv", "unknown.csv", "inputs.json")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


with log.step("read folds"):
    missing = [f for f in FOLDS if not (MODELS / f / "predictions.csv.gz").exists()]
    assert not missing, f"folds without predictions (run modal/phase9_loso_scanvi.py): {missing}"
    P = {f: pd.read_csv(MODELS / f / "predictions.csv.gz") for f in FOLDS}
    runs = {f: json.load(open(MODELS / f / "run.json")) for f in FOLDS}
    inputs = {str(p.relative_to(ROOT)): md5(p) for f in FOLDS for p in (MODELS / f / "predictions.csv.gz", MODELS / f / "split.json", MODELS / f / "run.json")}
    inputs[str(BASE.relative_to(ROOT))] = md5(BASE)
    log.write("cells per fold " + {f: len(d) for f, d in P.items()}.__repr__())

with log.step("score"):
    summ, pcls, ppat, conf, unk = [], [], [], [], []
    for f, d in P.items():
        lab = d[d.true != UNLAB]; y, p = lab.true.values, lab.pred.values
        classes = sorted(set(y))
        pc = []
        for c in classes:
            tp = int(((y == c) & (p == c)).sum()); nt = int((y == c).sum()); npred = int((p == c).sum())
            rec, prec = tp / nt, (tp / npred if npred else np.nan)
            pc.append(dict(fold=f, cls=c, n_true=nt, n_pred=npred, recall=rec, precision=prec,
                           f1=(2 * rec * prec / (rec + prec) if npred and rec + prec else 0.0), averaged=nt >= MIN_CLASS))
        pc = pd.DataFrame(pc); pcls.append(pc); avg = pc[pc.averaged]
        pp = lab.assign(ok=y == p).groupby("patient_id").ok.agg(["size", "mean"]).reset_index()
        pp.columns = ["patient_id", "n_labelled", "accuracy"]; pp.insert(0, "fold", f); ppat.append(pp)
        conf.append(pd.crosstab(pd.Series(y, name="true"), pd.Series(p, name="pred")).stack().rename("n").reset_index().assign(fold=f))
        hi = lab.maxprob.values >= MAXP
        mal = pc.set_index("cls").reindex(["malignant"])
        r = runs[f]
        summ.append(dict(fold=f, test="leave-one-publication-out (scVI + scANVI refitted; scArches surgery)", n_cells=len(d), n_labelled=len(lab),
                         n_patients=d.patient_id.nunique(), frac_unknown=round(float((d.true == UNLAB).mean()), 4),
                         accuracy=float((y == p).mean()), balanced_accuracy=float(avg.recall.mean()), macro_f1=float(avg.f1.mean()),
                         classes_averaged=";".join(avg.cls), classes_below_min=";".join(pc.cls[~pc.averaged]),
                         malignant_recall=float(mal.recall.iloc[0]), malignant_precision=float(mal.precision.iloc[0]),
                         worst_class=avg.sort_values("recall").cls.iloc[0], worst_class_recall=float(avg.recall.min()),
                         coverage_maxprob_ge_0_9=float(hi.mean()), accuracy_maxprob_ge_0_9=float((y[hi] == p[hi]).mean()) if hi.any() else np.nan,
                         patient_accuracy_min=float(pp.accuracy.min()), patient_accuracy_median=float(pp.accuracy.median()),
                         scanvi_inner_validation_accuracy=r["scanvi_inner_validation_accuracy"], scvi_epochs=r["scvi_epochs"],
                         scanvi_epochs=r["scanvi_epochs"], scanvi_best_epoch=r["scanvi_best_epoch"]))
        u = d[d.true == UNLAB]
        if len(u):
            ug = u.groupby("pred").maxprob.agg(["size", "median"]).reset_index(); ug.columns = ["pred", "n", "median_maxprob"]
            ug["fraction"] = ug.n / len(u); ug.insert(0, "fold", f); unk.append(ug)
    summ = pd.DataFrame(summ)
    b = pd.read_csv(BASE).query("classifier == 'new_scanvi' and labels == 'corrected' and cells == 'heldout_patients'").iloc[0]
    summ = pd.concat([summ, pd.DataFrame([dict(fold="baseline_heldout_patients", test="models/scanvi: 5 held-out patients, scVI trained on all cells",
                                               n_labelled=int(b.n_labelled), accuracy=b.accuracy, balanced_accuracy=b.balanced_accuracy)])], ignore_index=True)
    pcls = pd.concat(pcls); ppat = pd.concat(ppat); conf = pd.concat(conf)[["fold", "true", "pred", "n"]]; unk = pd.concat(unk)
    lab_all = pd.concat([d[d.true != UNLAB] for d in P.values()])
    pooled = dict(fold="pooled_all_folds", test="all held-out publications pooled (cell-weighted)", n_cells=sum(len(d) for d in P.values()),
                  n_labelled=len(lab_all), accuracy=float((lab_all.true == lab_all.pred).mean()),
                  balanced_accuracy=float(np.mean([((lab_all.true == c) & (lab_all.pred == c)).sum() / (lab_all.true == c).sum() for c in sorted(set(lab_all.true))])))
    summ = pd.concat([summ, pd.DataFrame([pooled])], ignore_index=True)

with log.step("write"):
    for df, p in zip([summ, pcls, ppat, conf, unk], OUTS[:5]):
        df.to_csv(p, index=False)
    json.dump(dict(script="scripts/32_loso_classifier.py", params=dict(min_class=MIN_CLASS, maxprob=MAXP), inputs=inputs), open(OUTS[5], "w"), indent=1)
    log.write("summary:\n" + summ.drop(columns=["test"]).round(4).to_string(index=False)
              + "\nper-class recall (fold x class):\n" + pcls.pivot(index="cls", columns="fold", values="recall").round(3).to_string())
