#!/usr/bin/env python
"""Phase 9 science check / step 40 — do per-cell metaprogram-family assignments measure per-sample composition?

Two assignments exist (19_metaprogram_scores.py): 'metaprogram' (UCell scores standardised WITHIN SAMPLE; used for Fig. 3d and
the cell2location reference) and 'metaprogram_assayz' (standardised WITHIN ASSAY CLASS; used for the mp_family composition
models). An assignment that measures composition should give a higher per-sample fraction of a family in samples whose own
cNMF found a programme of that family. For each lineage family and assignment: AUC of the per-sample family fraction for
discriminating samples with vs without a member programme of that family (per-sample cNMF, 16_metaprograms.py), and the
fraction of samples with >= 5 % of cells assigned to the family.
Inputs : data/processed/metaprogram_scores.parquet, results/tables/metaprogram_members.csv, results/tables/metaprograms.csv
Outputs (results/tables/quality/): family_assignment_check.csv, family_assignment_per_sample.csv, family_assignment_check_inputs.json
Run: nice -n 10 python scripts/40_family_assignment_check.py [--force]
Log: logs/phase9_quality.log
"""
import sys, json, hashlib, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
log = U.StepLogger(LOGS / "phase9_quality.log", context="40_family_assignment_check")
T = RESULTS / "tables"; OUT = T / "quality"
OUTS = [OUT / f for f in ("family_assignment_check.csv", "family_assignment_per_sample.csv", "family_assignment_check_inputs.json")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)
INPUTS = [PROCESSED / "metaprogram_scores.parquet", T / "metaprogram_members.csv", T / "metaprograms.csv"]
MIN_FRAC = 0.05


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


with log.step("check"):
    sc = pd.read_parquet(INPUTS[0], columns=["sample_id", "assay_class", "metaprogram", "metaprogram_assayz"])
    mem = pd.read_csv(INPUTS[1]); mps = pd.read_csv(INPUTS[2]).set_index("metaprogram")
    mem = mem[mem.metaprogram != "unassigned"].copy()
    mem["sample_id"] = mem.program.str.split("|").str[0]
    mem["family"] = mem.metaprogram.map(mps.family)
    found = mem.dropna(subset=["family"]).groupby(["sample_id", "family"]).size().gt(0).unstack(fill_value=False)
    samples = sorted(set(sc.sample_id) & set(mem.sample_id))
    fams = [f for f in mps.family.dropna().unique() if f in set(sc.metaprogram)]
    per, rows = [], []
    for col, name in [("metaprogram", "within_sample"), ("metaprogram_assayz", "within_assay_class")]:
        frac = pd.crosstab(sc.sample_id, sc[col], normalize="index").reindex(samples).fillna(0)
        for f in fams:
            y = found.reindex(samples)[f].fillna(False).astype(bool).values if f in found else np.zeros(len(samples), bool)
            x = frac[f].values if f in frac else np.zeros(len(samples))
            auc = roc_auc_score(y, x) if 0 < y.sum() < len(y) else np.nan
            rows.append(dict(assignment=name, family=f, n_samples=len(samples), n_samples_programme_found=int(y.sum()), auc=auc,
                             frac_samples_ge_5pct=float((x >= MIN_FRAC).mean()), median_fraction=float(np.median(x))))
        per.append(frac.assign(assignment=name).reset_index())
    res = pd.DataFrame(rows); res.to_csv(OUTS[0], index=False); pd.concat(per).to_csv(OUTS[1], index=False)
    log.write(res.round(3).to_string(index=False))
    json.dump(dict(script="scripts/40_family_assignment_check.py", params=dict(min_frac=MIN_FRAC),
                   inputs={str(p.relative_to(RESULTS.parent)): md5(p) for p in INPUTS}), open(OUTS[2], "w"), indent=1)
