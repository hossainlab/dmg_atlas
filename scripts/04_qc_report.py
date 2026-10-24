#!/usr/bin/env python
"""Phase 2 report: merge per-study QC tables, run sanity checks, write results/qc/qc_report.md + qc_per_sample.csv
and optional violin figures (results/qc/qc_violin_<study>.png). Reads only obs/var from the .qc.h5ad files (no matrices)."""
import sys
import re
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas import utils  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import h5py  # noqa: E402
import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from anndata.io import read_elem  # noqa: E402

from dmg_atlas import INTERIM, RESULTS, LOGS  # noqa: E402

QC = RESULTS / "qc"
STUDIES = ["jessa_scrna", "jessa_multiome", "collot", "desisto", "zang", "scpca", "okada", "filbin", "liu", "labelle"]
UMI_STUDIES = ["jessa_scrna", "jessa_multiome", "collot", "desisto", "zang", "scpca", "okada"]


def read_obs_var(path):
    with h5py.File(path) as f:
        obs = read_elem(f["obs"]); var = read_elem(f["var"])
        layers = list(f["layers"].keys())
    return obs, var, layers


def main():
    log = utils.StepLogger(LOGS / "phase2_qc.log", context="report")
    t0 = time.time()
    per_sample = pd.concat([pd.read_csv(QC / f"qc_per_sample_{s}.csv") for s in STUDIES if (QC / f"qc_per_sample_{s}.csv").exists()], ignore_index=True)
    per_sample.to_csv(QC / "qc_per_sample.csv", index=False)
    harm = pd.read_csv(QC / "gene_harmonisation.csv")

    obs_all, var_sets, layer_info = {}, {}, {}
    for s in STUDIES:
        for suffix in ("", "outgroup."):
            p = INTERIM / f"{s}.{suffix}qc.h5ad"
            if p.exists():
                obs, var, layers = read_obs_var(p)
                obs_all[(s, suffix)] = obs
                layer_info[(s, suffix)] = (len(var), layers, p.stat().st_size / 1e9)
                if suffix == "":
                    var_sets[s] = set(var.index[var.get("measured_in_all_samples", pd.Series(True, index=var.index)).astype(bool)])

    # ---- sanity (c): genes shared by all UMI studies
    umi_present = [s for s in UMI_STUDIES if s in var_sets]
    shared_umi = set.intersection(*[var_sets[s] for s in umi_present])
    shared_all = set.intersection(*[var_sets[s] for s in var_sets])
    log.write(f"harmonised genes shared by all {len(umi_present)} UMI studies: {len(shared_umi)}; by all {len(var_sets)} studies: {len(shared_all)}")
    pd.Series(sorted(shared_umi), name="symbol").to_csv(QC / "genes_shared_all_umi_studies.csv", index=False)

    # ---- sanity (a): Jessa malignant fraction passing QC
    jessa_lines = []
    for s in ("jessa_scrna", "jessa_multiome"):
        for suffix in ("", "outgroup."):
            o = obs_all.get((s, suffix))
            if o is None:
                continue
            mal = o["author_Malignant_normal_consensus_Jessa2022"].astype(str) == "Malignant"
            nor = o["author_Malignant_normal_consensus_Jessa2022"].astype(str) == "Normal"
            aq = o["author_qc_pass"].astype(bool)
            jessa_lines.append({
                "file": f"{s}.{suffix}qc.h5ad", "author_malignant": int(mal.sum()), "malignant_pass_qc": int((mal & o.qc_pass).sum()),
                "frac_malignant_pass": float((mal & o.qc_pass).sum() / max(mal.sum(), 1)),
                "author_normal": int(nor.sum()), "frac_normal_pass": float((nor & o.qc_pass).sum() / max(nor.sum(), 1)),
                "author_qc_cells": int(aq.sum()), "frac_author_qc_pass_ours": float((aq & o.qc_pass).sum() / max(aq.sum(), 1)),
                "cells_not_in_author_metadata": int((~aq).sum()), "frac_non_author_cells_passing_ours": float(((~aq) & o.qc_pass).sum() / max((~aq).sum(), 1)),
                "malignant_flagged_doublet": float(o.loc[mal & o.qc_pass, "predicted_doublet"].astype("boolean").mean()) if (mal & o.qc_pass).any() else np.nan,
            })
    jessa_df = pd.DataFrame(jessa_lines)
    jessa_df.to_csv(QC / "sanity_jessa_malignant_qc.csv", index=False)

    # ---- sample-level flags for every study (audit DQ-03 / DQ-04), derived from obs so that studies not re-ingested get them too:
    #      mt_genes_present = any cell of the sample has pct_mito > 0 (annotation carries MT- genes with signal)
    #      sample_qc_pass   = median n_genes of qc_pass cells >= configs/qc.yaml sample_gate.min_median_genes_pass (flag only)
    import yaml
    cfg = yaml.safe_load(open(ROOT / "configs" / "qc.yaml"))
    gate_thr = float(cfg.get("sample_gate", {}).get("min_median_genes_pass", 500))
    flags = []
    for (s_, suffix), o in obs_all.items():
        g = o.groupby("sample_id", observed=True)
        mt = g.pct_mito.max().fillna(0) > 0
        med = o[o.qc_pass.astype(bool)].groupby("sample_id", observed=True).n_genes.median()
        rule = g.qc_mito_rule.first() if "qc_mito_rule" in o else None
        for sid in g.size().index:
            flags.append({"sample_id": sid, "mt_genes_present": bool(mt.get(sid, False)), "sample_median_genes_pass": float(med.get(sid, np.nan)),
                          "sample_qc_pass": bool(med.get(sid, np.nan) >= gate_thr), "qc_mito_rule": rule.get(sid, "") if rule is not None else ""})
    flags = pd.DataFrame(flags)
    for c in ["mt_genes_present", "sample_qc_pass", "qc_mito_rule", "sample_gate_min_median_genes"]:
        if c in per_sample:
            per_sample = per_sample.drop(columns=c)
    per_sample = per_sample.merge(flags, on="sample_id", how="left")
    per_sample["sample_gate_min_median_genes"] = gate_thr
    per_sample.to_csv(QC / "qc_per_sample.csv", index=False)
    flagged = per_sample[~per_sample.sample_qc_pass.astype(bool)].sort_values("median_genes")
    no_mt = per_sample[~per_sample.mt_genes_present.astype(bool)]
    log.write(f"sample gate (median genes of passing cells >= {gate_thr:.0f}): {len(flagged)} flagged: {flagged.sample_id.tolist()}; samples without MT genes: {len(no_mt)}")

    # ---- sanity (b): assay label vs median pct_mito (atlas + outgroup samples; samples without MT- genes are not informative)
    ps = per_sample.copy()
    contradictions = []
    for _, r in ps[ps.mt_genes_present].iterrows():
        if r.qc_assay_class == "nuclei" and r.median_pct_mito_before > 5:
            contradictions.append(f"{r.sample_id} labelled {r.assay} but median pct_mito before QC = {r.median_pct_mito_before:.2f}% (>5%)")
        if r.qc_assay_class == "cells" and r.median_pct_mito_before < 1:
            contradictions.append(f"{r.sample_id} labelled {r.assay} but median pct_mito before QC = {r.median_pct_mito_before:.2f}% (<1%, nucleus-like)")
    assay_tab = ps[ps.mt_genes_present].groupby("qc_assay_class").median_pct_mito_before.describe()[["count", "min", "25%", "50%", "75%", "max"]]

    # ---- per-study totals
    tot = ps.groupby(["study", "include_in_atlas"]).agg(samples=("sample_id", "count"), cells_before=("cells_before_qc", "sum"),
                                                        cells_after=("cells_after_qc", "sum")).reset_index()
    tot["pct_retained"] = 100 * tot.cells_after / tot.cells_before
    # weighted median genes: use per-sample medians -> report median of medians plus pooled from obs
    pooled = []
    for (s, suffix), o in obs_all.items():
        p = o[o.qc_pass]
        pooled.append({"study": s, "include_in_atlas": suffix == "", "median_genes_pooled": float(np.median(p.n_genes)) if len(p) else np.nan,
                       "median_counts_pooled": float(np.median(p.n_counts)) if len(p) else np.nan,
                       "median_pct_mito_pooled": float(np.nanmedian(p.pct_mito)) if len(p) else np.nan,
                       "doublet_rate_pooled": float(p.predicted_doublet.astype("boolean").mean()) if p.predicted_doublet.notna().any() else np.nan,
                       "doublet_rate_fixed_pooled": float(p.predicted_doublet_fixed.astype("boolean").mean()) if "predicted_doublet_fixed" in p and p.predicted_doublet_fixed.notna().any() else np.nan,
                       "n_genes_var": layer_info[(s, suffix)][0], "layers": ",".join(layer_info[(s, suffix)][1]), "file_gb": round(layer_info[(s, suffix)][2], 2)})
    tot = tot.merge(pd.DataFrame(pooled), on=["study", "include_in_atlas"], how="left")
    tot.to_csv(QC / "qc_per_study.csv", index=False)

    # ---- violin figures (plain, no titles or text in panels; identification in the report caption).
    #      audit DQ-13: TPM / log-scale cohorts (Filbin, LaBelle) get the genes panel only (counts and % mito carry no information
    #      there); the per-sample min/max gene thresholds are drawn; studies with > 24 samples are split over two rows.
    for (s, suffix), o in obs_all.items():
        if suffix:
            continue
        layers = layer_info[(s, suffix)][1]
        is_tpm = any(L in layers for L in ("tpm", "log2_tpm10"))
        cols = [("n_genes", "genes detected")] if is_tpm else [("n_genes", "genes detected"), ("n_counts", "counts"), ("pct_mito", "% mitochondrial")]
        samples = sorted(o.sample_id.unique())
        n_rows = 2 if len(samples) > 24 else 1
        half = len(samples) // 2 + len(samples) % 2
        chunks = [samples[:half], samples[half:]] if n_rows == 2 else [samples]
        fig, axes = plt.subplots(n_rows, len(cols), figsize=(max(6, 0.35 * max(len(c) for c in chunks) + 3) * (1 if is_tpm else len(cols) / 1.5), 3.2 * n_rows), squeeze=False)
        for r, smp in enumerate(chunks):
            for ax, (col, ylab) in zip(axes[r], cols):
                data = [o.loc[o.sample_id == sm, col].dropna().values for sm in smp]
                data = [d if len(d) else np.array([np.nan]) for d in data]
                ax.violinplot([np.log10(d + 1) if col != "pct_mito" else d for d in data], showmedians=True, widths=0.8)
                if col == "n_genes":
                    lo = [o.loc[o.sample_id == sm, "qc_min_genes"].iloc[0] for sm in smp]
                    hi = [o.loc[o.sample_id == sm, "qc_max_genes"].iloc[0] for sm in smp]
                    ax.hlines(np.log10(np.array(lo, dtype=float) + 1), np.arange(1, len(smp) + 1) - 0.4, np.arange(1, len(smp) + 1) + 0.4, colors="0.3", linestyles="--", linewidth=0.8)
                    ax.hlines(np.log10(np.array(hi, dtype=float) + 1), np.arange(1, len(smp) + 1) - 0.4, np.arange(1, len(smp) + 1) + 0.4, colors="0.3", linestyles=":", linewidth=0.8)
                if col == "pct_mito":
                    th = [o.loc[o.sample_id == sm, "qc_pct_mito_max"].iloc[0] for sm in smp]
                    ax.hlines(th, np.arange(1, len(smp) + 1) - 0.4, np.arange(1, len(smp) + 1) + 0.4, colors="0.3", linestyles="--", linewidth=0.8)
                ax.set_ylabel(f"log10({ylab}+1)" if col != "pct_mito" else ylab)
                ax.set_xticks(range(1, len(smp) + 1)); ax.set_xticklabels([sm.split("__", 1)[1] for sm in smp], rotation=90, fontsize=6)
                ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout(); fig.savefig(QC / f"qc_violin_{s}.png", dpi=300); plt.close(fig)

    # ---- log-derived peak RSS and wall time
    peak, wall = {}, {}
    for line in open(LOGS / "phase2_qc.log"):
        m = re.search(r"\[(\w+)\] DONE \w+ \| total wall ([\d.]+) min \| peak RSS ([\d.]+) GB", line)
        if m:
            wall[m.group(1)] = float(m.group(2)); peak[m.group(1)] = float(m.group(3))

    # ---- report
    L = ["# Phase 2 QC report — DMG atlas", "", f"Generated {time.strftime('%Y-%m-%d %H:%M')} by scripts/04_qc_report.py from data/interim/*.qc.h5ad. "
         "Thresholds: configs/qc.yaml (nuclei: pct_mito < 5 %, cells: < 20 %, min genes 300 for 10x / 1000 for Smart-seq2, per-sample max genes = 99.5th percentile; "
         "scrublet expected rate 0.06 on 10x samples, flagged only). Cells are flagged (`obs['qc_pass']`), none removed. Doublets are not removed.", ""]
    L += ["## Per-study totals", "", "| study | set | samples | cells before | cells after | % retained | median genes (pooled) | median counts | median % mito | doublet rate (auto thr) | doublet rate (score>0.25) | genes in var | layers | file GB | wall min | peak RSS GB |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in tot.sort_values(["include_in_atlas", "study"], ascending=[False, True]).iterrows():
        tpm_study = any(L_ in str(r.layers) for L_ in ("tpm", "log2_tpm10"))
        L.append(f"| {r.study} | {'atlas' if r.include_in_atlas else 'out-group'} | {r.samples} | {r.cells_before:,} | {r.cells_after:,} | {r.pct_retained:.1f} | {r.median_genes_pooled:.0f} | {'n/a' if tpm_study else f'{r.median_counts_pooled:.0f}'} | {'n/a' if tpm_study else f'{r.median_pct_mito_pooled:.2f}'} | "
                 f"{'' if pd.isna(r.doublet_rate_pooled) else f'{r.doublet_rate_pooled:.3f}'} | {'' if pd.isna(r.doublet_rate_fixed_pooled) else f'{r.doublet_rate_fixed_pooled:.3f}'} | {r.n_genes_var:,} | {r.layers} | {r.file_gb} | {wall.get(r.study, float('nan')):.1f} | {peak.get(r.study, float('nan')):.2f} |")
    atlas = tot[tot.include_in_atlas]
    L += ["", f"Atlas (included) samples: **{int(atlas.samples.sum())} samples, {int(atlas.cells_before.sum()):,} cells before QC, {int(atlas.cells_after.sum()):,} after QC ({100 * atlas.cells_after.sum() / atlas.cells_before.sum():.1f} %)**. "
          f"Out-group files: {int(tot[~tot.include_in_atlas].cells_before.sum()):,} cells before / {int(tot[~tot.include_in_atlas].cells_after.sum()):,} after.", ""]
    L += ["## Sample-level flags (audit corrections 2026-10-19)", "",
          f"`sample_qc_pass` (flag only, nothing dropped): median detected genes of the sample's qc_pass cells >= {gate_thr:.0f} (configs/qc.yaml `sample_gate`). "
          f"**{len(flagged)} sample(s) flagged:**", ""]
    L += [f"- {r.sample_id} ({r.study}, {'atlas' if r.include_in_atlas else 'out-group'}): {r.cells_after_qc:,} passing cells, median {r.median_genes:.0f} genes / {r.median_counts:.0f} counts" for _, r in flagged.iterrows()] or ["- none"]
    L += ["", f"`mt_genes_present` = the sample's feature annotation carries MT- genes with signal; where False the mitochondrial filter is inoperative. "
          f"**{len(no_mt)} sample(s) without MT- genes** ({', '.join(sorted(no_mt.study.unique()))}): " + ", ".join(no_mt.sample_id.tolist()), ""]
    if "qc_mito_rule" in per_sample and per_sample.qc_mito_rule.astype(str).str.contains("author_prep_aware").any():
        apa = per_sample[per_sample.qc_mito_rule.astype(str).str.contains("author_prep_aware")]
        L += [f"Author-prep-aware mito rule (nuclei threshold {cfg['thresholds']['nuclei']['pct_mito_max']} % applied to Smart-seq2 frozen-nuclei samples with MT- genes): "
              f"{len(apa)} samples: " + ", ".join(apa.sample_id.tolist()), ""]
    L += ["Value scales: UMI studies and Liu = counts; Filbin = linear TPM (`layers['tpm']`, per-cell sum 1e6); LaBelle = log2(TPM/10+1) as deposited "
          "(`layers['log2_tpm10']`, `uns['value_scale']`). 'median counts' is therefore not meaningful for Filbin/LaBelle and is shown as n/a.", ""]
    L += ["## Sanity checks", "", f"**(a) Jessa author-labelled malignant cells passing our QC** (`results/qc/sanity_jessa_malignant_qc.csv`):", "",
          "| file | author malignant | pass QC | fraction | author normal pass fraction | author-QC cells passing ours | barcodes not in author metadata | of which pass ours | malignant flagged doublet |", "|---|---|---|---|---|---|---|---|---|"]
    for _, r in jessa_df.iterrows():
        L.append(f"| {r.file} | {r.author_malignant:,} | {r.malignant_pass_qc:,} | {r.frac_malignant_pass:.3f} | {r.frac_normal_pass:.3f} | {r.frac_author_qc_pass_ours:.3f} | {r.cells_not_in_author_metadata:,} | {r.frac_non_author_cells_passing_ours:.3f} | {r.malignant_flagged_doublet:.3f} |")
    L += ["", "**(b) Assay label vs median pct_mito before QC** (UMI studies; TPM tables of Filbin/LaBelle carry no MT- genes so pct_mito is 0 there and not informative):", "", "```", assay_tab.round(2).to_string(), "```", ""]
    if contradictions:
        L += [f"Samples whose median pct_mito contradicts the manifest assay label ({len(contradictions)}):", ""] + [f"- {c}" for c in contradictions]
    else:
        L += ["No sample contradicts its assay label."]
    L += ["", f"**(c) Gene overlap:** {len(shared_umi):,} harmonised symbols are measured in every sample of all {len(umi_present)} UMI studies ({', '.join(umi_present)}); "
          f"{len(shared_all):,} are shared by all {len(var_sets)} studies including the Smart-seq2/TPM cohorts. Gene universe = GRCh38-2020-A symbols (36,591 unique). Per-study losses: `results/qc/gene_harmonisation.csv`.", ""]
    L += ["## Gene harmonisation per study", "", "| study | input genes | mapped by Ensembl | exact symbol | alias | ambiguous dropped | unmapped dropped | duplicates collapsed | output genes |", "|---|---|---|---|---|---|---|---|---|"]
    hs = harm.groupby("study").first().reset_index()
    for _, r in hs.iterrows():
        L.append(f"| {r.study} | {r.n_input_genes:,} | {r.mapped_ensembl:,} | {r.mapped_exact_symbol:,} | {r.mapped_alias:,} | {r.ambiguous_dropped} | {r.unmapped_dropped:,} | {r.duplicate_symbols_collapsed} | {r.n_output_genes:,} |")
    L += ["", "## Per-sample QC (`results/qc/qc_per_sample.csv`)", "", "| sample | study | assay | set | before | after | median genes | median counts | median % mito | doublet rate (auto thr) | doublet rate (score>0.25) | fail low genes | fail high genes | fail mito | sample gate | MT genes |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in ps.sort_values(["include_in_atlas", "study", "sample_id"], ascending=[False, True, True]).iterrows():
        dr = "" if pd.isna(r.doublet_rate_flagged) else f"{r.doublet_rate_flagged:.3f}"
        df_ = "" if pd.isna(r.get("doublet_rate_fixed025", np.nan)) else f"{r.doublet_rate_fixed025:.3f}"
        mc = "n/a" if r.qc_metric_basis in ("tpm", "log2_tpm10") else f"{r.median_counts:.0f}"
        pm = "n/a" if not bool(r.mt_genes_present) else f"{r.median_pct_mito:.2f}"
        L.append(f"| {r.sample_id} | {r.study} | {r.assay} | {'atlas' if r.include_in_atlas else 'out-group'} | {r.cells_before_qc:,} | {r.cells_after_qc:,} | {r.median_genes:.0f} | {mc} | {pm} | {dr} | {df_} | {r.n_fail_low_genes} | {r.n_fail_high_genes} | {r.n_fail_high_mito} | {'' if bool(r.sample_qc_pass) else 'FLAG'} | {'yes' if bool(r.mt_genes_present) else 'no'} |")
    L += ["", "## Figures", ""] + [f"- results/qc/qc_violin_{s}.png — per-sample violins (all cells before QC) of log10(genes detected+1)" +
                                   ("" if any(L_ in ",".join(layer_info[(s, '')][1]) for L_ in ("tpm", "log2_tpm10")) else ", log10(counts+1) and % mitochondrial") +
                                   "; dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; " +
                                   ("two rows of samples" if obs_all[(s, '')].sample_id.nunique() > 24 else "one row") +
                                   ("; TPM/log-scale cohort: counts and % mito panels omitted (no MT- genes, values not counts)" if any(L_ in ",".join(layer_info[(s, '')][1]) for L_ in ("tpm", "log2_tpm10")) else "")
                                   for s in STUDIES if (QC / f"qc_violin_{s}.png").exists() and (s, '') in obs_all]
    (QC / "qc_report.md").write_text("\n".join(L) + "\n")
    summary = {"per_study": tot.to_dict(orient="records"), "assay_contradictions": contradictions, "shared_umi_genes": len(shared_umi), "shared_all_genes": len(shared_all),
               "peak_rss_gb_by_study": peak, "wall_min_by_study": wall, "jessa": jessa_df.to_dict(orient="records")}
    (QC / "qc_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    log.write(f"report written in {time.time() - t0:.0f}s; contradictions: {len(contradictions)}")


if __name__ == "__main__":
    main()
