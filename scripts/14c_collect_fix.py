#!/usr/bin/env python
"""Phase 6 / step 14c (audit remediation 2026-10-20; SP-03) — harvest the 20,000-epoch cell2location runs.

Moves the pre-audit per-cohort abundance directories to data/processed/spatial_pre_audit/<cohort>/, copies the new
q05/mean abundance CSVs into data/processed/spatial/<cohort>/, copies every *_history.csv into
results/tables/c2l_histories/, writes results/tables/c2l_deconvolution_runs.json, the plateau statistic table
results/tables/c2l_elbo_plateau.csv (SP-03: the claim 'all traces had plateaued' is now backed by the harvested
histories) and, once all four cohorts are present, results/tables/spatial_abundance_long.csv (old -> .pre_audit.csv).

usage: python 14c_collect_fix.py <cohort>=<harvest_dir> [...]
"""
from __future__ import annotations
import sys, json, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger  # noqa: E402
import numpy as np, pandas as pd  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_spatial.log", context="14c_collect_fix")
TABLES = ROOT / "results/tables"; SPATIAL = ROOT / "data/processed/spatial"
PRE = ROOT / "data/processed/spatial_pre_audit"; HIST = TABLES / "c2l_histories"
COHORTS = ["brainstem", "collot_st", "huchede", "ren"]


def plateau(hist: pd.DataFrame):
    e = hist.iloc[:, 0].values.astype(float)
    tot = e[0] - e[-1]
    last = e[int(len(e) * 0.95)] - e[-1]
    return len(e), float(e[0]), float(e[-1]), float(100 * last / tot)


def main(harvest: dict[str, str]):
    PRE.mkdir(parents=True, exist_ok=True); HIST.mkdir(parents=True, exist_ok=True)
    runs = json.load(open(TABLES / "c2l_deconvolution_runs.json")) if (TABLES / "c2l_deconvolution_runs.json").exists() else {}
    if (TABLES / "c2l_deconvolution_runs.json").exists() and not (TABLES / "c2l_deconvolution_runs.pre_audit.json").exists():
        shutil.copy(TABLES / "c2l_deconvolution_runs.json", TABLES / "c2l_deconvolution_runs.pre_audit.json")
        runs = {}
    rows = []
    for cohort, hdir in harvest.items():
        h = Path(hdir)
        src = SPATIAL / cohort
        if src.exists() and not (PRE / cohort).exists():
            shutil.move(str(src), str(PRE / cohort))
            LOG.write(f"{cohort}: pre-audit abundances -> {PRE/cohort}")
        src.mkdir(parents=True, exist_ok=True)
        n = 0
        for f in sorted(h.glob("*_abundance.csv")):
            shutil.copy(f, src / f.name); n += 1
        hf = h / f"{cohort}_history.csv"
        shutil.copy(hf, HIST / hf.name)
        hist = pd.read_csv(hf, index_col=0)
        ne, e0, e1, pct = plateau(hist)
        runs[cohort] = json.load(open(h / f"{cohort}_run.json"))
        rows.append(dict(cohort=cohort, n_sections=n, epochs=ne, elbo_first=e0, elbo_last=e1,
                         pct_of_total_reduction_in_last_5pct=pct, history_file=f"results/tables/c2l_histories/{hf.name}"))
        LOG.write(f"{cohort}: {n} sections, {ne} epochs, ELBO {e0:.4g} -> {e1:.4g}, last 5% = {pct:.4f}% of reduction")
    json.dump(runs, open(TABLES / "c2l_deconvolution_runs.json", "w"), indent=1)
    pl = pd.DataFrame(rows)
    p = TABLES / "c2l_elbo_plateau.csv"
    if p.exists():
        old = pd.read_csv(p); pl = pd.concat([old[~old.cohort.isin(pl.cohort)], pl], ignore_index=True)
    pl.sort_values("cohort").to_csv(p, index=False)
    have = [c for c in COHORTS if (SPATIAL / c).exists() and list((SPATIAL / c).glob("*_abundance.csv"))]
    LOG.write(f"cohorts with new abundances: {have}")
    if set(have) == set(COHORTS):
        sec = pd.read_csv(ROOT / "results/qc/spatial_sections.csv"); sec = sec[sec.species == "human"]
        long = []
        for cohort in COHORTS:
            for f in sorted((SPATIAL / cohort).glob("*_abundance.csv")):
                df = pd.read_csv(f, index_col=0)
                q = [c for c in df.columns if c.startswith("q05_")]; m = [c for c in df.columns if c.startswith("mean_")]
                a = df[q].copy(); a.columns = [c[4:] for c in q]
                b = df[m].copy(); b.columns = [c[5:] for c in m]
                a = a.stack().rename("q05").reset_index(); a.columns = ["spot", "class", "q05"]
                b = b.stack().rename("mean").reset_index(); b.columns = ["spot", "class", "mean"]
                l = a.merge(b, on=["spot", "class"])
                l.insert(0, "section_id", df["section_id"].iloc[0]); l.insert(0, "cohort", cohort)
                l = l.merge(df[["array_row", "array_col", "pxl_row", "pxl_col"]].rename_axis("spot").reset_index(), on="spot")
                long.append(l)
        long = pd.concat(long, ignore_index=True)
        long = long.merge(sec[["section_id", "region_or_diagnosis"]], on="section_id", how="left")
        out = TABLES / "spatial_abundance_long.csv"
        if out.exists() and not (TABLES / "spatial_abundance_long.pre_audit.csv").exists():
            shutil.move(str(out), str(TABLES / "spatial_abundance_long.pre_audit.csv"))
        long.to_csv(out, index=False)
        LOG.write(f"spatial_abundance_long.csv: {long.shape}, sections {long.section_id.nunique()}, "
                  f"spots {long.spot.nunique()}, classes {long['class'].nunique()}, "
                  f"NaN {int(long.q05.isna().sum())}, negative {int((long.q05 < 0).sum())}")
    else:
        LOG.write("spatial_abundance_long.csv NOT written (waiting for: " + ",".join(set(COHORTS) - set(have)) + ")")


if __name__ == "__main__":
    main(dict(x.split("=", 1) for x in sys.argv[1:]))
