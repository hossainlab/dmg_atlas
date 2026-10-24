#!/usr/bin/env python
"""Phase 6 / step 13a (audit correction, 2026-10-19) — scripted malignant-cell selection for the metaprogram track.

Replaces dmg_atlas.metaprograms.select_malignant (consensus 'malignant' OR scANVI malignant with maxprob >= 0.9; audit MC-04:
29 % of the input were scANVI rescues at assay-dependent rates) with the corrected Phase-3 consensus calls
(data/processed/malignant_calls.csv, regenerated in audit wave 1).

Selections written per cell (data/processed/malignant_selection_audit.parquet):
  primary    : consensus_call == 'malignant' (event_backed or score_only CNV confidence) AND qc_pass,
               excluding okada__1734_pons (sample_qc_pass False, 1,163 genes) and the ScPCA samples (not in atlas_core).
  eventtrust : consensus_call_eventtrust == 'malignant' AND qc_pass, same exclusions (sensitivity set, CNV-05).
  old_sel    : the pre-audit rule, for the change tally only.
Also writes results/tables/malignant_selection_by_assay.csv (per-assay breakdown of both selections, with the
pre-audit rule and the cnv_confidence split) and results/tables/cnmf/sample_rerun_plan.csv (per sample: old vs new counts,
fraction of the malignant set changed, analysed_now, rerun flag; rerun when the set changed by > 5 % of the pre-audit set).
Run: python scripts/13a_malignant_selection.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, h5py  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402

log = U.StepLogger(LOGS / "fix_metaprograms.log", context="13a_malignant_selection")
CORE = PROCESSED / "atlas_core.h5ad"; CALLS = PROCESSED / "malignant_calls.csv"
EXCLUDE_SAMPLES = {"okada__1734_pons"}          # sample_qc_pass False (DQ-03); ScPCA samples are not in atlas_core
EXCLUDE_STUDIES = {"scpca"}
MIN_CELLS, RERUN_FRAC = 200, 0.05
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}

with log.step("read obs + corrected calls"):
    with h5py.File(CORE, "r") as f:
        cols = ["sample_id", "study", "assay", "patient_id", "scanvi_label", "scanvi_maxprob", "consensus_call", "qc_pass"]
        obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in cols}, index=read_elem(f["obs"]["_index"]))
    for c in cols:
        if c not in ("scanvi_maxprob", "qc_pass"):
            obs[c] = obs[c].astype(str)
    obs["qc_pass"] = obs.qc_pass.astype(bool)
    mc = pd.read_csv(CALLS, usecols=["cell_id", "consensus_call", "cnv_confidence", "consensus_call_eventtrust", "low_confidence_sample"],
                     index_col="cell_id", low_memory=False)
    missing = ~obs.index.isin(mc.index)
    assert not missing.any(), f"{missing.sum()} atlas cells without a corrected call"
    obs = obs.join(mc.add_prefix("new_"))
    excl = obs.sample_id.isin(EXCLUDE_SAMPLES) | obs.study.isin(EXCLUDE_STUDIES)
    obs["old_sel"] = ((obs.consensus_call == "malignant") | ((obs.scanvi_label == "malignant") & (obs.scanvi_maxprob >= 0.9))) & obs.qc_pass
    obs["primary"] = (obs.new_consensus_call == "malignant") & obs.qc_pass & ~excl
    obs["eventtrust"] = (obs.new_consensus_call_eventtrust == "malignant") & obs.qc_pass & ~excl
    obs["assay_class"] = obs.assay.map(ASSAY_CLASS)
    log.write(f"old rule {int(obs.old_sel.sum())}, primary {int(obs.primary.sum())}, eventtrust {int(obs.eventtrust.sum())} cells; "
              f"excluded samples {sorted(set(obs.sample_id[excl]))}")

with log.step("write selection"):
    out = obs[["sample_id", "study", "assay", "assay_class", "patient_id", "qc_pass", "old_sel", "primary", "eventtrust"]].copy()
    out["cnv_confidence"] = obs.new_cnv_confidence.astype(str); out["consensus_call"] = obs.new_consensus_call.astype(str)
    out["consensus_call_eventtrust"] = obs.new_consensus_call_eventtrust.astype(str)
    out.index.name = "cell_id"
    out.to_parquet(PROCESSED / "malignant_selection_audit.parquet")

with log.step("per-assay breakdown"):
    rows = []
    for ac, d in obs[~excl & obs.qc_pass].groupby("assay_class"):
        rows.append(dict(assay_class=ac, n_qc_pass=len(d), n_samples=d.sample_id.nunique(),
                         n_old_rule=int(d.old_sel.sum()), n_old_scanvi_rescue=int((d.old_sel & (d.consensus_call != "malignant")).sum()),
                         n_primary=int(d.primary.sum()),
                         n_primary_event_backed=int((d.primary & (d.new_cnv_confidence == "event_backed")).sum()),
                         n_primary_score_only=int((d.primary & (d.new_cnv_confidence == "score_only")).sum()),
                         n_eventtrust=int(d.eventtrust.sum()), n_eventtrust_added=int((d.eventtrust & ~d.primary).sum()),
                         n_old_not_primary=int((d.old_sel & ~d.primary).sum()), n_primary_not_old=int((d.primary & ~d.old_sel).sum())))
    ba = pd.DataFrame(rows)
    tot = ba.drop(columns="assay_class").sum().to_dict(); tot["assay_class"] = "all"; ba = pd.concat([ba, pd.DataFrame([tot])], ignore_index=True)
    ba["frac_primary_score_only"] = ba.n_primary_score_only / ba.n_primary
    ba["frac_old_scanvi_rescue"] = ba.n_old_scanvi_rescue / ba.n_old_rule
    ba.to_csv(RESULTS / "tables" / "malignant_selection_by_assay.csv", index=False)
    log.write("\n" + ba.round(3).to_string(index=False))

with log.step("rerun plan"):
    samp_prev = pd.read_csv(RESULTS / "tables" / "cnmf.pre_audit" / "samples_analysed.csv", index_col=0) \
        if (RESULTS / "tables" / "cnmf.pre_audit" / "samples_analysed.csv").exists() else \
        pd.read_csv(RESULTS / "tables" / "cnmf" / "samples_analysed.csv", index_col=0)
    g = obs.groupby("sample_id")
    per = pd.DataFrame({"n_cells": g.size(), "study": g.study.first(), "assay": g.assay.first(), "assay_class": g.assay_class.first(),
                        "patient_id": g.patient_id.first(), "n_old": g.old_sel.sum(), "n_primary": g.primary.sum(), "n_eventtrust": g.eventtrust.sum(),
                        "n_old_not_primary": (obs.old_sel & ~obs.primary).groupby(obs.sample_id).sum(),
                        "n_primary_not_old": (obs.primary & ~obs.old_sel).groupby(obs.sample_id).sum(),
                        "n_primary_event_backed": (obs.primary & (obs.new_cnv_confidence == "event_backed")).groupby(obs.sample_id).sum(),
                        "n_primary_score_only": (obs.primary & (obs.new_cnv_confidence == "score_only")).groupby(obs.sample_id).sum()})
    per["frac_changed"] = (per.n_old_not_primary + per.n_primary_not_old) / per.n_old.replace(0, np.nan)
    per["analysed_pre_audit"] = per.index.isin(samp_prev.index[samp_prev.analysed])
    per["analysed_now"] = (per.n_primary >= MIN_CELLS) & ~per.index.isin(EXCLUDE_SAMPLES) & ~per.study.isin(EXCLUDE_STUDIES)
    per["analysed_eventtrust"] = (per.n_eventtrust >= MIN_CELLS) & ~per.index.isin(EXCLUDE_SAMPLES) & ~per.study.isin(EXCLUDE_STUDIES)
    per["rerun"] = per.analysed_now & ((per.frac_changed.fillna(1.0) > RERUN_FRAC) | ~per.analysed_pre_audit)
    per["reuse_pre_audit"] = per.analysed_now & ~per.rerun
    per["dropped"] = per.analysed_pre_audit & ~per.analysed_now
    (RESULTS / "tables" / "cnmf").mkdir(parents=True, exist_ok=True)
    per.to_csv(RESULTS / "tables" / "cnmf" / "sample_rerun_plan.csv")
    log.write(f"analysed now {int(per.analysed_now.sum())} (cells {int((per.analysed_now & (per.assay_class == 'cells')).sum())}, "
              f"nuclei {int((per.analysed_now & (per.assay_class == 'nuclei')).sum())}); rerun {int(per.rerun.sum())}; reuse "
              f"{sorted(per.index[per.reuse_pre_audit])}; dropped {sorted(per.index[per.dropped])}; "
              f"eventtrust-analysable {int(per.analysed_eventtrust.sum())}")
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
