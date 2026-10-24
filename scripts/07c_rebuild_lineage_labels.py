#!/usr/bin/env python
"""Audit remediation (reference-labels track, 19 Oct 2026) — rebuild obs['lineage_label'] for the
311,618 reference cells from the corrected data/processed/malignant_calls.csv, exactly as
scripts/07_assemble_atlas_input.py lineage_label() did (malignant = consensus_call malignant; normal
lineages from consensus_celltype_coarse; else Unknown).

Outputs
  results/tables/reflabels/reference_labels_corrected.csv.gz   cell_id, sample_id, patient_id, study, assay,
        lineage_label_pre_audit, lineage_label, consensus_call_pre_audit, consensus_call,
        consensus_call_eventtrust, cnv_confidence, low_confidence_sample, okada_1734_pons_flag
  results/tables/reflabels/lineage_label_change_crosstab.csv   pre-audit x corrected lineage_label
  results/tables/reflabels/lineage_label_change_by_sample.csv  per-sample counts of changed cells
  results/tables/reflabels/lineage_label_counts.csv            counts old vs new
The Modal scANVI job consumes reference_labels_corrected.csv.gz (the training h5ad on the Volume keeps
its pre-audit obs; the job overwrites obs['lineage_label'] from this file).
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from anndata.io import read_elem  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_reflabels.log", context="07c_rebuild_lineage_labels")
OUT = ROOT / "results/tables/reflabels"
OUT.mkdir(parents=True, exist_ok=True)
LINEAGES = ["myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular", "other"]
LOW_CONF_SAMPLE = "okada__1734_pons"   # sample_qc_pass False (DQ-03), 90 % score-only CNV calls (CNV-09)


def lineage_label(obs: pd.DataFrame) -> pd.Series:
    """Verbatim copy of scripts/07_assemble_atlas_input.py::lineage_label."""
    lab = pd.Series("Unknown", index=obs.index, dtype=object)
    lab[obs["consensus_call"].eq("malignant")] = "malignant"
    normal = obs["consensus_call"].eq("normal")
    ct = obs["consensus_celltype_coarse"].astype(str)
    lab[normal & ct.isin(LINEAGES)] = ct[normal & ct.isin(LINEAGES)]
    lab[normal & ~ct.isin(LINEAGES)] = "other"
    return lab


def main():
    with LOG.step("read reference obs (atlas_input_hvg.train.h5ad)"):
        with h5py.File(ROOT / "data/processed/atlas_input_hvg.train.h5ad") as f:
            obs = read_elem(f["obs"])
        obs = obs.rename(columns={"lineage_label": "lineage_label_pre_audit", "consensus_call": "consensus_call_pre_audit",
                                  "consensus_celltype_coarse": "consensus_celltype_coarse_pre_audit"})
        for c in obs.columns:
            if c != "pct_mito":
                obs[c] = obs[c].astype(str)
        LOG.write(f"reference cells {len(obs):,}; samples {obs.sample_id.nunique()}; patients {obs.patient_id.nunique()}")

    with LOG.step("read corrected malignant_calls.csv"):
        cols = ["cell_id", "consensus_call", "consensus_call_eventtrust", "consensus_celltype_coarse", "cnv_confidence", "low_confidence_sample", "cnv_call", "marker_call"]
        calls = pd.read_csv(ROOT / "data/processed/malignant_calls.csv", usecols=cols, index_col="cell_id", low_memory=False)
        missing = obs.index.difference(calls.index)
        assert len(missing) == 0, f"{len(missing)} reference cells missing from malignant_calls.csv"
        calls = calls.loc[obs.index]
        LOG.write(f"calls rows matched {len(calls):,}; consensus_call counts {calls.consensus_call.value_counts().to_dict()}")

    with LOG.step("rebuild lineage_label"):
        obs = pd.concat([obs, calls], axis=1)
        obs["lineage_label"] = lineage_label(obs)
        # sanity: re-applying the rule to the pre-audit columns must reproduce the pre-audit label
        chk = lineage_label(obs[["consensus_call_pre_audit", "consensus_celltype_coarse_pre_audit"]]
                            .rename(columns={"consensus_call_pre_audit": "consensus_call", "consensus_celltype_coarse_pre_audit": "consensus_celltype_coarse"}))
        repro = float((chk.values == obs["lineage_label_pre_audit"].values).mean())
        LOG.write(f"rule applied to pre-audit columns reproduces pre-audit lineage_label: {repro:.6f}")
        assert repro == 1.0
        obs["okada_1734_pons_flag"] = obs["sample_id"].eq(LOW_CONF_SAMPLE)
        changed = obs["lineage_label"] != obs["lineage_label_pre_audit"]
        LOG.write(f"lineage_label changed for {int(changed.sum()):,} cells ({changed.mean():.4f})")
        ct = pd.crosstab(obs["lineage_label_pre_audit"], obs["lineage_label"], margins=True)
        ct.to_csv(OUT / "lineage_label_change_crosstab.csv")
        LOG.write("crosstab (pre-audit rows x corrected cols):\n" + ct.to_string())
        counts = pd.DataFrame({"pre_audit": obs["lineage_label_pre_audit"].value_counts(), "corrected": obs["lineage_label"].value_counts()}).fillna(0).astype(int)
        counts["delta"] = counts["corrected"] - counts["pre_audit"]
        counts.to_csv(OUT / "lineage_label_counts.csv")
        LOG.write("counts:\n" + counts.to_string())
        bys = (obs.assign(changed=changed).groupby(["study", "sample_id"], observed=True)
               .agg(n_cells=("changed", "size"), n_changed=("changed", "sum"),
                    n_malignant_pre=("lineage_label_pre_audit", lambda s: int((s == "malignant").sum())),
                    n_malignant_new=("lineage_label", lambda s: int((s == "malignant").sum())),
                    n_unknown_pre=("lineage_label_pre_audit", lambda s: int((s == "Unknown").sum())),
                    n_unknown_new=("lineage_label", lambda s: int((s == "Unknown").sum())),
                    frac_score_only=("cnv_confidence", lambda s: float((s == "score_only").mean())),
                    low_confidence_sample=("low_confidence_sample", "first")).reset_index())
        bys["frac_changed"] = bys["n_changed"] / bys["n_cells"]
        bys.sort_values("n_changed", ascending=False).to_csv(OUT / "lineage_label_change_by_sample.csv", index=False)
        LOG.write("top changed samples:\n" + bys.sort_values("n_changed", ascending=False).head(12).to_string())
        o = obs[obs.okada_1734_pons_flag]
        LOG.write(f"{LOW_CONF_SAMPLE}: {len(o):,} cells; score_only {float((o.cnv_confidence == 'score_only').mean()):.3f}; "
                  f"lineage new {o.lineage_label.value_counts().to_dict()}; pre {o.lineage_label_pre_audit.value_counts().to_dict()}")

    with LOG.step("write reference_labels_corrected.csv.gz"):
        keep = ["sample_id", "patient_id", "study", "assay", "lineage_label_pre_audit", "lineage_label", "consensus_call_pre_audit", "consensus_call",
                "consensus_call_eventtrust", "consensus_celltype_coarse", "cnv_confidence", "cnv_call", "marker_call", "low_confidence_sample", "okada_1734_pons_flag"]
        out = obs[keep].copy()
        out.index.name = "cell_id"
        out.to_csv(OUT / "reference_labels_corrected.csv.gz", compression="gzip")
        LOG.write(f"written {OUT / 'reference_labels_corrected.csv.gz'} ({(OUT / 'reference_labels_corrected.csv.gz').stat().st_size / 1e6:.1f} MB)")
    LOG.write("done")


if __name__ == "__main__":
    main()
