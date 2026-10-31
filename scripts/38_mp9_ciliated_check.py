#!/usr/bin/env python
"""Phase 9 science check / step 38 — is the ciliated metaprogram (MP9) tumour or contaminating normal ependyma?

Normal ependymal cells express the same motile-cilia genes as MP9 (ADGB, CFAP43, DNAH9) but carry no tumour copy-number
alterations. Event flags are inherited from the copy-number CLUSTER, so they cannot decide this; per-cell arm-level values can.
Per sample with >= MIN_HIGH MP9-high malignant cells (within-sample z of the MP9 UCell score >= Z_HIGH), and >= MIN_REF normal
reference cells:
  tumour axis   = mean arm profile (8 arm values) of the sample's other malignant cells minus that of its normal reference cells
  tumour index  = per-cell projection onto the tumour axis, scaled so that normal reference = 0 and other malignant = 1
  compared      : MP9-high malignant cells vs the sample's other malignant cells vs its normal reference cells
                  (Mann-Whitney MP9-high vs normal reference; median tumour index; per-cell cnv_score)
Also: author cell-type labels of MP9-high cells where the source study provides them (Jessa 2022 consensus, Okada cellType).
Inputs : data/processed/metaprogram_scores.parquet (19_metaprogram_scores.py), data/processed/malignant_calls.csv (06)
Outputs (results/tables/quality/): mp9_tumour_index_per_sample.csv, mp9_author_labels.csv, mp9_ciliated_check_inputs.json
Run: nice -n 10 python scripts/38_mp9_ciliated_check.py [--force]
Log: logs/phase9_quality.log
"""
import sys, json, hashlib, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd  # noqa: E402
import scipy.stats as ss  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
log = U.StepLogger(LOGS / "phase9_quality.log", context="38_mp9_ciliated_check")
OUT = RESULTS / "tables/quality"; OUT.mkdir(parents=True, exist_ok=True)
OUTS = [OUT / f for f in ("mp9_tumour_index_per_sample.csv", "mp9_author_labels.csv", "mp9_ciliated_check_inputs.json")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)
ARMS = ["arm_1q", "arm_2", "arm_7", "arm_13q", "arm_14q", "arm_10q", "arm_17p", "arm_4q12"][:7]
Z_HIGH, MIN_HIGH, MIN_REF = 2.0, 20, 30
AUTHOR = ["author_Cell_type_consensus_Jessa2022", "author_cellType", "author_celltype", "author_annotation"]
SCORES, CALLS = PROCESSED / "metaprogram_scores.parquet", PROCESSED / "malignant_calls.csv"


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


with log.step("read"):
    sc = pd.read_parquet(SCORES, columns=["sample_id", "assay_class", "z_MP9", "metaprogram"])
    mc = pd.read_csv(CALLS, usecols=["cell_id", "sample_id", "study", "cnv_score", "consensus_call", "is_reference_final"] + ARMS + AUTHOR,
                     index_col="cell_id", low_memory=False)
    mc = mc[mc.sample_id.isin(set(sc.sample_id))]
    mal = sc.join(mc.drop(columns="sample_id"), how="left")
    assert mal.cnv_score.notna().all(), "metaprogram cells without copy-number calls"
    ref = mc[mc.is_reference_final.astype(str).str.lower().eq("true") & mc.consensus_call.eq("normal")]
    log.write(f"{len(mal)} selected malignant cells in {mal.sample_id.nunique()} samples; MP9-high (z >= {Z_HIGH}): {int((mal.z_MP9 >= Z_HIGH).sum())}; "
              f"normal reference cells in these samples: {len(ref)}")

with log.step("tumour index per sample"):
    rows, cells = [], []
    for s, d in mal.groupby("sample_id"):
        hi = d[d.z_MP9 >= Z_HIGH]; lo = d[d.z_MP9 < Z_HIGH]; r = ref[ref.sample_id == s]
        if len(hi) < MIN_HIGH or len(r) < MIN_REF:
            continue
        mu_t, mu_n = lo[ARMS].mean().values, r[ARMS].mean().values
        axis = mu_t - mu_n; nrm = float(axis @ axis)
        if nrm == 0:
            continue
        ti = lambda x: ((x[ARMS].values - mu_n) @ axis) / nrm  # noqa: E731
        t_hi, t_lo, t_ref = ti(hi), ti(lo), ti(r)
        u = ss.mannwhitneyu(t_hi, t_ref, alternative="greater")
        rows.append(dict(sample_id=s, study=d.study.iloc[0], assay_class=d.assay_class.iloc[0], n_mp9_high=len(hi), n_other_malignant=len(lo),
                         n_reference=len(r), tumour_index_mp9_high=np.median(t_hi), tumour_index_other_malignant=np.median(t_lo),
                         tumour_index_reference=np.median(t_ref), frac_mp9_high_above_ref_p95=float((t_hi > np.quantile(t_ref, 0.95)).mean()),
                         cnv_score_mp9_high=hi.cnv_score.median(), cnv_score_other_malignant=lo.cnv_score.median(),
                         cnv_score_reference=r.cnv_score.median(), p_mp9_vs_reference=u.pvalue))
    res = pd.DataFrame(rows)
    res["q_mp9_vs_reference"] = ss.false_discovery_control(res.p_mp9_vs_reference.values, method="bh") if len(res) else []
    res.to_csv(OUTS[0], index=False)
    log.write(f"{len(res)} samples evaluable\n" + res.round(3).to_string(index=False))

with log.step("author labels"):
    hi = mal[mal.z_MP9 >= Z_HIGH]
    lab = pd.concat([hi[[c]].dropna().rename(columns={c: "author_label"}).assign(source=c) for c in AUTHOR if c in hi])
    tab = (lab.groupby(["source", "author_label"]).size().rename("n_mp9_high_cells").reset_index()
           .sort_values(["source", "n_mp9_high_cells"], ascending=[True, False]))
    tab["fraction"] = tab.n_mp9_high_cells / tab.groupby("source").n_mp9_high_cells.transform("sum")
    tab.to_csv(OUTS[1], index=False)
    log.write(tab.round(3).to_string(index=False))
    json.dump(dict(script="scripts/38_mp9_ciliated_check.py", params=dict(z_high=Z_HIGH, min_high=MIN_HIGH, min_ref=MIN_REF, arms=ARMS),
                   inputs={str(p.relative_to(RESULTS.parent)): md5(p) for p in (SCORES, CALLS)}), open(OUTS[2], "w"), indent=1)
