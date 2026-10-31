#!/usr/bin/env python
"""Phase 9 science fix / step 44 — collect the cell2location runs with the within-assay-class reference and re-run the
spatial statistics.

Deconvolution: modal/phase6_c2l.py on Modal A10G with signatures from modal/inputs/c2l_reference_assayz.h5ad
(43_c2l_reference_assayz.py; 250 epochs) and the published per-cohort Visium inputs, same settings as the published runs
(5 cells per spot, detection alpha 20, target 20,000 epochs within a 2,400 s GPU budget). Harvested to models/c2l_assayz/<cohort>/.
This script builds results/tables/spatial_abundance_long_assayz.csv exactly as 14c_collect_fix.py builds the published table
(q05 and mean abundance per spot and class, spot coordinates, section annotation), records the runs and ELBO plateau
statistics, then runs 15b_spatial_tests_fix.py --tag _assayz (all spatial tables with the suffix _assayz; the published
tables are not modified) and 39_spatial_robustness.py --tag _assayz (Hartung-Knapp-Sidik-Jonkman re-pooling).
Outputs: results/tables/spatial_abundance_long_assayz.csv, results/tables/quality/c2l_assayz_runs.json,
         results/tables/quality/c2l_assayz_elbo_plateau.csv, and the _assayz tables written by 15b / 39
Run: nice -n 10 python scripts/44_spatial_assayz.py
Log: logs/phase9_quality.log
"""
from __future__ import annotations
import sys, json, subprocess
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger  # noqa: E402
import pandas as pd  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase9_quality.log", context="44_spatial_assayz")
H = ROOT / "models/c2l_assayz"; T = ROOT / "results/tables"; Q = T / "quality"
COHORTS = ["brainstem", "collot_st", "huchede", "ren"]


def plateau(hist: pd.DataFrame):
    e = hist.iloc[:, 0].values; n = len(e); tail = max(1, n // 20)
    return n, float(e[0]), float(e[-1]), float(100 * (e[-tail - 1] - e[-1]) / (e[0] - e[-1])) if e[0] != e[-1] else 0.0


with LOG.step("collect"):
    missing = [c for c in COHORTS if not list((H / c).glob("*_abundance.csv"))]
    assert not missing, f"cohorts without harvested abundances: {missing}"
    sec = pd.read_csv(ROOT / "results/qc/spatial_sections.csv"); sec = sec[sec.species == "human"]
    long, runs, pl = [], {}, []
    for c in COHORTS:
        runs[c] = json.load(open(H / c / f"{c}_run.json"))
        n, e0, e1, pct = plateau(pd.read_csv(H / c / f"{c}_history.csv", index_col=0))
        pl.append(dict(cohort=c, epochs=n, elbo_first=e0, elbo_last=e1, pct_of_total_reduction_in_last_5pct=pct))
        for f in sorted((H / c).glob("*_abundance.csv")):
            df = pd.read_csv(f, index_col=0)
            q = [x for x in df.columns if x.startswith("q05_")]; m = [x for x in df.columns if x.startswith("mean_")]
            a = df[q].copy(); a.columns = [x[4:] for x in q]; b = df[m].copy(); b.columns = [x[5:] for x in m]
            a = a.stack().rename("q05").reset_index(); a.columns = ["spot", "class", "q05"]
            b = b.stack().rename("mean").reset_index(); b.columns = ["spot", "class", "mean"]
            l = a.merge(b, on=["spot", "class"]); l.insert(0, "section_id", df["section_id"].iloc[0]); l.insert(0, "cohort", c)
            l = l.merge(df[["array_row", "array_col", "pxl_row", "pxl_col"]].rename_axis("spot").reset_index(), on="spot")
            long.append(l)
    long = pd.concat(long, ignore_index=True).merge(sec[["section_id", "region_or_diagnosis"]], on="section_id", how="left")
    pub = pd.read_csv(T / "spatial_abundance_long.csv", usecols=["section_id", "spot"]).drop_duplicates()
    assert set(map(tuple, long[["section_id", "spot"]].drop_duplicates().values)) == set(map(tuple, pub.values)), "spots differ from the published run"
    long.to_csv(T / "spatial_abundance_long_assayz.csv", index=False)
    json.dump(runs, open(Q / "c2l_assayz_runs.json", "w"), indent=1); pd.DataFrame(pl).to_csv(Q / "c2l_assayz_elbo_plateau.csv", index=False)
    LOG.write(f"spatial_abundance_long_assayz.csv {long.shape}; sections {long.section_id.nunique()}; epochs "
              + str({p['cohort']: p['epochs'] for p in pl}))

for script in ("15b_spatial_tests_fix.py", "39_spatial_robustness.py"):
    with LOG.step(f"run {script} --tag _assayz"):
        r = subprocess.run([sys.executable, str(ROOT / "scripts" / script), "--tag", "_assayz", "--force"], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
