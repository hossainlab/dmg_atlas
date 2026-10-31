#!/usr/bin/env python
"""Phase 3 aggregation: merge per-study call tables, per-sample summary, validation vs author labels, hallmark CNV
frequencies, CNV heatmaps for representative samples and the Jessa cnv_score vs marker-score scatter.

Usage: python scripts/06_cnv_aggregate.py
Outputs: data/processed/malignant_calls.csv
         results/tables/malignant_calls_per_sample.csv
         results/tables/malignant_call_validation.csv            pooled + sample-level (median/IQR/bootstrap CI) agreement,
                                                                  tuning vs held-out split, event-backed vs score-only
         results/tables/malignant_call_validation_summary.csv    one row per study (and all studies) with the headline numbers
         results/tables/tier3_excluded_author_composition.csv   author-label composition of cells excluded from tier-3 agreement
         results/tables/cnv_event_frequencies.csv                fraction of CNV-malignant cells (eq(True)) + n_evaluable + frac_samples
         results/tables/liu_author_vs_cnv_call_crosstab.csv      tier-2 crosstab (replaces the by-construction consensus crosstab)
         results/tables/jessa_author_celltype_vs_consensus_coarse.csv, jessa_author_celltype_vs_cnv_call.csv
         results/tables/consensus_ambiguous_composition.csv      CNV-05: why cells are ambiguous, by author label
         results/tables/cnv_audit_call_changes.csv, cnv_audit_call_changes_by_author_label.csv   vs malignant_calls.pre_audit.csv.gz
         results/tables/cnv_score_floor_sensitivity.csv          CNV-09: cluster calls re-derived at alternative score floors
         results/qc/cnv_cluster_stats_vs_author.csv               CNV-20: per-cluster evidence vs author labels (scripted)
         results/qc/cnv/cnv_heatmap_<sample>.png, results/qc/cnv/jessa_cnvscore_vs_markerscore.png
Audit corrections 2026-10-19 (docs/audit/cnv_audit.md, docs/audit/fixes_cnv.csv): event fractions use eq(True) and report
evaluable-cell counts (CNV-01); validation is reported at sample level with CIs (CNV-08), split by tuning set (CNV-06/07),
with tier-3 coverage and excluded-cell composition (CNV-04) and by call confidence (CNV-09).
"""
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas import utils  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from dmg_atlas import INTERIM, PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas import cnv as C  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402

STUDIES = ["jessa_scrna", "jessa_multiome", "collot", "zang", "desisto", "scpca", "okada", "liu", "filbin", "labelle"]
EVENTS = list(C.HALLMARK_EVENTS)
HEATMAP_SAMPLES = ["jessa_scrna__P-6337_S-8821", "jessa_multiome__P-1709_S-1709", "okada__1277_pons", "okada__1006_thalamus", "collot__PMC-DMG-6", "liu__BT869"]
# Samples processed in the "benchmark subset" iterations of logs/phase3_cnv.log (17:03-17:16, 2026-10-17) that set the final
# thresholds (audit CNV-06): 9 jessa_scrna + 7 okada + 4 jessa_multiome + 2 collot = 22 samples.
TUNING_SAMPLES = [
    "jessa_scrna__P-1713_S-1713", "jessa_scrna__P-3407_S-3447", "jessa_scrna__P-4198_S-4459", "jessa_scrna__P-4504_S-4916",
    "jessa_scrna__P-6254_S-8499", "jessa_scrna__P-6328_S-8672", "jessa_scrna__P-6337_S-8821", "jessa_scrna__P-6519_S-9084", "jessa_scrna__P-6640_S-9581",
    "okada__1006_thalamus", "okada__1428_pons", "okada__1734_pons", "okada__1890_pons", "okada__1963_pons", "okada__2054_thalamus", "okada__586_pons",
    "jessa_multiome__P-1694_S-1694", "jessa_multiome__P-1701_S-1701", "jessa_multiome__P-1779_S-1781", "jessa_multiome__P-6253_S-8498",
    "collot__DMG-PMC-9_edge", "collot__PMC-DMG-7",
]
PRE_AUDIT = PROCESSED / "malignant_calls.pre_audit.csv.gz"
N_BOOT = 2000


def prf(pred: pd.Series, truth: pd.Series) -> dict:
    """agreement, precision/recall/F1 for 'malignant' on cells where both are defined (and pred != unassigned)."""
    m = truth.notna() & pred.notna() & (pred != "unassigned")
    p, t = pred[m].astype(str), truth[m].astype(str)
    tp = int(((p == "malignant") & (t == "malignant")).sum()); fp = int(((p == "malignant") & (t != "malignant")).sum())
    fn = int(((p != "malignant") & (t == "malignant")).sum())
    out = dict(n_compared=int(m.sum()), agreement=float((p == t).mean()) if m.sum() else np.nan,
               precision_malignant=tp / (tp + fp) if tp + fp else np.nan, recall_malignant=tp / (tp + fn) if tp + fn else np.nan)
    out["f1_malignant"] = (2 * out["precision_malignant"] * out["recall_malignant"] / (out["precision_malignant"] + out["recall_malignant"])
                           if (tp + fp) and (tp + fn) and tp else np.nan)
    return out


def sample_level(vals: np.ndarray, seed: int = 0) -> dict:
    """Sample-level summary of a per-sample metric: n, median, IQR, unweighted mean, bootstrap 95 % CI of the mean (over samples)."""
    v = np.asarray(vals, float); v = v[~np.isnan(v)]
    if len(v) == 0:
        return dict(n_samples=0)
    out = dict(n_samples=int(len(v)), sample_median=float(np.median(v)), sample_q25=float(np.percentile(v, 25)), sample_q75=float(np.percentile(v, 75)),
               sample_mean=float(v.mean()), sample_min=float(v.min()))
    if len(v) >= 3:
        rng = np.random.default_rng(seed)
        boots = rng.choice(v, size=(N_BOOT, len(v)), replace=True).mean(axis=1)
        out["sample_mean_ci95_lo"], out["sample_mean_ci95_hi"] = float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))
        bmed = np.median(rng.choice(v, size=(N_BOOT, len(v)), replace=True), axis=1)
        out["sample_median_ci95_lo"], out["sample_median_ci95_hi"] = float(np.percentile(bmed, 2.5)), float(np.percentile(bmed, 97.5))
    return out


def tier3_pred(gg: pd.DataFrame) -> pd.Series:
    """author-independent tier 3: cnv_call where marker and CNV agree (astrocyte-like cells excluded), else NaN."""
    mk_bin = np.where(gg.marker_call == "normal", "normal", np.where(gg.marker_call == "malignant_candidate", "malignant", None))
    return gg.cnv_call.where(mk_bin == gg.cnv_call.values)


def validation_block(study: str, g: pd.DataFrame, subset: str) -> list[dict]:
    """Pooled + sample-level rows for one study (or 'all') and one subset of samples."""
    rows = []
    preds = {"tier2_cnv_call": lambda d: d.cnv_call,
             "tier2_cnv_call__event_backed": lambda d: d.cnv_call.where((d.cnv_call != "malignant") | (d.cnv_confidence == "event_backed")),
             "tier2_cnv_call__score_only": lambda d: d.cnv_call.where((d.cnv_call != "malignant") | (d.cnv_confidence == "score_only")),
             "tier3_cnv_marker_agree": tier3_pred,
             "tier3_consensus_call_BY_CONSTRUCTION": lambda d: d.consensus_call.where(d.consensus_call != "ambiguous"),
             "consensus_call_eventtrust_BY_CONSTRUCTION": lambda d: d.consensus_call_eventtrust.where(d.consensus_call_eventtrust != "ambiguous")}
    for tier, fn in preds.items():
        r = dict(study=study, level="study", subset=subset, sample_id=study, tier=tier, n_author_labelled=len(g),
                 n_ambiguous=int((g.consensus_call == "ambiguous").sum()))
        r.update(prf(fn(g), g.author_call)); r["coverage_of_labelled"] = r["n_compared"] / len(g) if len(g) else np.nan
        per = []
        for sid, gg in g.groupby("sample_id"):
            rr = dict(study=study, level="sample", subset=subset, sample_id=sid, tier=tier, n_author_labelled=len(gg),
                      n_ambiguous=int((gg.consensus_call == "ambiguous").sum()), tuning_sample=sid in TUNING_SAMPLES,
                      low_confidence_sample=bool(gg.low_confidence_sample.fillna(False).astype(bool).iloc[0]))
            rr.update(prf(fn(gg), gg.author_call)); rr["coverage_of_labelled"] = rr["n_compared"] / len(gg)
            per.append(rr)
        r.update(sample_level([p["agreement"] for p in per]))
        rows.append(r)
        if subset == "all_samples" and study != "all":
            rows.extend(per)
    return rows


def main():
    log = utils.StepLogger(LOGS / "phase3_cnv.log", context="aggregate")
    tabs = []
    for s in STUDIES:
        f = INTERIM / "cnv" / f"{s}.calls.csv"
        if f.exists():
            tabs.append(pd.read_csv(f, low_memory=False))
        else:
            log.write(f"missing {f}")
    calls = pd.concat(tabs, ignore_index=True)
    # author label columns requested in the output
    author_cols = {"jessa": ["author_Malignant_normal_consensus_Jessa2022", "author_Cell_type_consensus_Jessa2022"], "liu": ["author_annotation"],
                   "okada": ["author_cellType"], "collot": ["author_celltype"], "labelle": ["author_broad_annot"]}
    import h5py
    from anndata.io import read_elem
    extra = []
    for s in calls.study.unique():
        key = "jessa" if s.startswith("jessa") else s
        cols = author_cols.get(key)
        if not cols:
            continue
        with h5py.File(INTERIM / f"{s}.qc.h5ad", "r") as fh:
            obs = read_elem(fh["obs"])
        extra.append(obs.loc[calls.loc[calls.study == s, "cell_id"], cols])
    if extra:
        extra = pd.concat(extra)
        for c in extra.columns:
            calls[c] = calls["cell_id"].map(extra[c]).values
    # fine author state label (Jessa cell type where present, otherwise the raw author label)
    calls["author_state"] = calls["author_label"].astype(object)
    if "author_Cell_type_consensus_Jessa2022" in calls:
        jm = calls.study.str.startswith("jessa") & calls["author_Cell_type_consensus_Jessa2022"].notna()
        calls.loc[jm, "author_state"] = calls.loc[jm, "author_Cell_type_consensus_Jessa2022"].astype(str)
    front = ["cell_id", "sample_id", "study", "cnv_score", "cnv_call", "cnv_confidence", "marker_call", "consensus_call", "consensus_call_eventtrust",
             "consensus_celltype_coarse", "author_label", "author_call", "marker_type", "marker_normal_score", "marker_tumour_score", "cnv_cluster",
             "cnv_best_ref", "is_reference", "is_reference_final", "ref_type_final", "ref_source", "high_mito_nuclei", "low_confidence_sample"]
    front = [c for c in front if c in calls.columns]
    calls = calls[front + [c for c in calls.columns if c not in front]]
    calls.to_csv(PROCESSED / "malignant_calls.csv", index=False)
    log.write(f"malignant_calls.csv: {len(calls)} cells, {calls.sample_id.nunique()} samples")
    calls["cnv_confidence"] = calls["cnv_confidence"].astype(object)
    calls["low_confidence_sample"] = calls["low_confidence_sample"].fillna(False).astype(bool)

    # ---- per-sample summary
    rows = []
    for (study, sid), g in calls.groupby(["study", "sample_id"], sort=False):
        r = dict(study=study, sample_id=sid, n_cells=len(g), n_high_mito_relaxed=int(g.high_mito_nuclei.sum()) if g.high_mito_nuclei.any() else 0,
                 ref_source=g.ref_source.iloc[0], n_reference_final=int(g.is_reference_final.fillna(False).astype(bool).sum()),
                 tuning_sample=sid in TUNING_SAMPLES, low_confidence_sample=bool(g.low_confidence_sample.iloc[0]))
        for col in ["cnv_call", "cnv_confidence", "marker_call", "consensus_call", "consensus_call_eventtrust"]:
            vc = g[col].value_counts()
            for k in sorted(vc.index):
                r[f"{col}__{k}"] = int(vc[k])
        for k, v in g.consensus_celltype_coarse.value_counts().items():
            r[f"coarse__{k}"] = int(v)
        mal = g[g.cnv_call == "malignant"]
        r["frac_malignant_consensus"] = float((g.consensus_call == "malignant").mean())
        r["frac_cnv_malignant_score_only"] = float((mal.cnv_confidence == "score_only").mean()) if len(mal) else np.nan
        for e in EVENTS:
            # audit correction CNV-01: eq(True) (NaN = event not evaluable in this sample) + evaluable-cell count
            n_ev = int(mal[e].notna().sum()) if e in mal else 0
            r[f"n_evaluable_{e}"] = n_ev
            r[f"frac_malignant_{e}"] = float(mal[e].eq(True).sum() / n_ev) if n_ev else np.nan
        if g.author_call.notna().any():
            r.update({f"cnv_vs_author__{k}": v for k, v in prf(g.cnv_call, g.author_call).items()})
            r.update({f"consensus_vs_author__{k}": v for k, v in prf(g.consensus_call.where(g.consensus_call != "ambiguous"), g.author_call).items()})
        rows.append(r)
    per_sample = pd.DataFrame(rows)
    samp_meta = pd.concat([pd.read_csv(INTERIM / "cnv" / f"{s}.samples.csv") for s in STUDIES if (INTERIM / "cnv" / f"{s}.samples.csv").exists()])
    keep = ["sample_id", "n_ref_candidates_own", "n_ref_immune_own", "n_pooled_immune", "n_pooled_lineage", "pooled_types", "n_ref_final", "n_genes", "n_windows",
            "ref_cats", "ref_cnv_mean", "ref_cnv_sd", "thr_cnv_score", "n_clusters", "n_clusters_malignant", "n_clusters_event_backed", "n_clusters_score_only",
            "arms_evaluable", "hallmark_regions_missing", "status", "wall_s", "peak_rss_gb"]
    per_sample = per_sample.merge(samp_meta[[c for c in keep if c in samp_meta]], on="sample_id", how="left")
    per_sample.to_csv(RESULTS / "tables" / "malignant_calls_per_sample.csv", index=False)

    # ---- validation (author-labelled studies): pooled, sample-level, tuning split, confidence split
    labelled = calls[calls.author_call.notna()]
    val = []
    for study in [s for s in STUDIES if (labelled.study == s).any()]:
        g = labelled[labelled.study == study]
        val.extend(validation_block(study, g, "all_samples"))
        for subset, m in [("tuning_samples", g.sample_id.isin(TUNING_SAMPLES)), ("held_out_samples", ~g.sample_id.isin(TUNING_SAMPLES)),
                          ("low_confidence_samples", g.low_confidence_sample), ("adequate_samples", ~g.low_confidence_sample)]:
            if m.any():
                val.extend(validation_block(study, g[m], subset))
    val.extend(validation_block("all", labelled, "all_samples"))
    for subset, m in [("tuning_samples", labelled.sample_id.isin(TUNING_SAMPLES)), ("held_out_samples", ~labelled.sample_id.isin(TUNING_SAMPLES))]:
        val.extend(validation_block("all", labelled[m], subset))
    val = pd.DataFrame(val)
    val["flag_low_agreement"] = val.agreement < 0.90
    val.to_csv(RESULTS / "tables" / "malignant_call_validation.csv", index=False)

    # headline summary: one row per study (tier 2, all samples) incl. study-equal-weight mean across studies
    st = val[(val.level == "study") & (val.subset == "all_samples") & (val.tier == "tier2_cnv_call")].set_index("study")
    t3 = val[(val.level == "study") & (val.subset == "all_samples") & (val.tier == "tier3_cnv_marker_agree")].set_index("study")
    tun = val[(val.level == "study") & (val.subset == "tuning_samples") & (val.tier == "tier2_cnv_call")].set_index("study")
    held = val[(val.level == "study") & (val.subset == "held_out_samples") & (val.tier == "tier2_cnv_call")].set_index("study")
    evb = val[(val.level == "study") & (val.subset == "all_samples") & (val.tier == "tier2_cnv_call__event_backed")].set_index("study")
    sco = val[(val.level == "study") & (val.subset == "all_samples") & (val.tier == "tier2_cnv_call__score_only")].set_index("study")
    summ = []
    for s in st.index:
        r = dict(study=s, n_labelled=int(st.loc[s, "n_author_labelled"]), n_samples=int(st.loc[s, "n_samples"]),
                 pooled_agreement=st.loc[s, "agreement"], pooled_precision=st.loc[s, "precision_malignant"], pooled_recall=st.loc[s, "recall_malignant"],
                 sample_median=st.loc[s, "sample_median"], sample_q25=st.loc[s, "sample_q25"], sample_q75=st.loc[s, "sample_q75"],
                 sample_mean=st.loc[s, "sample_mean"], sample_mean_ci95_lo=st.loc[s].get("sample_mean_ci95_lo", np.nan), sample_mean_ci95_hi=st.loc[s].get("sample_mean_ci95_hi", np.nan),
                 tier3_coverage=t3.loc[s, "coverage_of_labelled"], tier3_agreement=t3.loc[s, "agreement"],
                 n_tuning_samples=int(tun.loc[s, "n_samples"]) if s in tun.index else 0, tuning_pooled_agreement=tun.loc[s, "agreement"] if s in tun.index else np.nan,
                 tuning_sample_median=tun.loc[s, "sample_median"] if s in tun.index else np.nan,
                 heldout_pooled_agreement=held.loc[s, "agreement"] if s in held.index else np.nan, heldout_sample_median=held.loc[s, "sample_median"] if s in held.index else np.nan,
                 event_backed_precision=evb.loc[s, "precision_malignant"], score_only_precision=sco.loc[s, "precision_malignant"])
        summ.append(r)
    summ = pd.DataFrame(summ)
    per_study = summ[summ.study != "all"]
    summ.loc[summ.study == "all", "study_equal_weight_mean_agreement"] = per_study.pooled_agreement.mean()
    summ.to_csv(RESULTS / "tables" / "malignant_call_validation_summary.csv", index=False)

    # tier-3 exclusion composition (CNV-04): labelled cells where marker and CNV disagree (or astrocyte-like), by author state
    excl_rows = []
    for study, g in labelled.groupby("study"):
        t3p = tier3_pred(g); ex = g[t3p.isna()]
        # author_call is NOT a function of author_state within a study (Jessa multiome has author_state values that
        # occur under both calls), so it must be part of the grouping key: taking gg.author_call.iloc[0] mislabelled
        # 666 multiome cells as malignant (re-audit CNV-04).
        for (st_lab, acall, cnvc, mk), gg in ex.groupby(["author_state", "author_call", "cnv_call", "marker_call"]):
            excl_rows.append(dict(study=study, author_state=st_lab, author_call=acall, cnv_call=cnvc, marker_call=mk, n_cells=len(gg),
                                  frac_of_labelled_in_study=len(gg) / len(g)))
    pd.DataFrame(excl_rows).sort_values(["study", "n_cells"], ascending=[True, False]).to_csv(RESULTS / "tables" / "tier3_excluded_author_composition.csv", index=False)

    # Liu: tier-2 crosstab (the consensus crosstab is tautological: consensus requires author agreement; CNV-04)
    liu = labelled[labelled.study == "liu"]
    if len(liu):
        ct = pd.crosstab(liu.author_label, liu.cnv_call); ct["n_ambiguous_consensus"] = liu.groupby("author_label").consensus_call.apply(lambda s: int((s == "ambiguous").sum()))
        ct.to_csv(RESULTS / "tables" / "liu_author_vs_cnv_call_crosstab.csv")
    jes = calls[calls.study.str.startswith("jessa") & calls.author_call.notna()]
    pd.crosstab(jes["author_Cell_type_consensus_Jessa2022"], jes.consensus_celltype_coarse).to_csv(RESULTS / "tables" / "jessa_author_celltype_vs_consensus_coarse.csv")
    pd.crosstab(jes["author_Cell_type_consensus_Jessa2022"], jes.cnv_call).to_csv(RESULTS / "tables" / "jessa_author_celltype_vs_cnv_call.csv")

    # ---- ambiguous-class composition (CNV-05)
    amb = calls[calls.consensus_call == "ambiguous"].copy()
    mk_bin = np.where(amb.marker_call == "normal", "normal", np.where(amb.marker_call == "malignant_candidate", "malignant", "astrocyte_like"))
    au = amb.author_call.astype(object)
    reason = np.where(amb.cnv_call == "unassigned", "cnv_unassigned",
             np.where(mk_bin != amb.cnv_call.values, "marker_vs_cnv_disagree",
             np.where(au.notna() & (au.astype(str) != amb.cnv_call.astype(str)), "author_vs_cnv_disagree", "other")))
    amb["reason"] = reason
    amb["pattern"] = "cnv=" + amb.cnv_call.astype(str) + "|marker=" + amb.marker_call.astype(str)
    amb_rows = []
    for (study, reason_, pattern, st_lab), gg in amb.groupby(["study", "reason", "pattern", amb.author_state.fillna("no_author_label")]):
        amb_rows.append(dict(study=study, reason=reason_, pattern=pattern, author_state=st_lab, n_cells=len(gg),
                             n_author_malignant=int((gg.author_call == "malignant").sum()), n_author_normal=int((gg.author_call == "normal").sum()),
                             n_eventtrust_malignant=int((gg.consensus_call_eventtrust == "malignant").sum())))
    pd.DataFrame(amb_rows).sort_values(["study", "n_cells"], ascending=[True, False]).to_csv(RESULTS / "tables" / "consensus_ambiguous_composition.csv", index=False)

    # ---- hallmark event frequencies per study (fraction of CNV-malignant cells in clusters carrying the event; CNV-01 eq(True))
    ev_rows = []
    for study, g in calls.groupby("study"):
        mal = g[g.cnv_call == "malignant"]
        r = dict(study=study, n_malignant=len(mal), n_samples=g.sample_id.nunique(), n_samples_with_malignant=mal.sample_id.nunique(),
                 n_malignant_event_backed=int((mal.cnv_confidence == "event_backed").sum()), n_malignant_score_only=int((mal.cnv_confidence == "score_only").sum()))
        for e in EVENTS:
            ok = mal[e].notna() if len(mal) else pd.Series(dtype=bool)
            n_ev = int(ok.sum()); r[f"n_evaluable_{e}"] = n_ev
            r[e] = float(mal[e].eq(True).sum() / n_ev) if n_ev else np.nan
            # fraction of samples (with >= 20 evaluable malignant cells) where >= 20 % of malignant cells carry the event (CNV-23)
            fs = mal[ok].groupby("sample_id")[e].agg(["mean", "size"]) if n_ev else pd.DataFrame(columns=["mean", "size"])
            fs = fs[fs["size"] >= 20]
            r[f"{e}__n_samples_evaluable"] = int(len(fs)); r[f"{e}__frac_samples"] = float((fs["mean"] >= 0.2).mean()) if len(fs) else np.nan
        ev_rows.append(r)
    ev = pd.DataFrame(ev_rows); ev.to_csv(RESULTS / "tables" / "cnv_event_frequencies.csv", index=False)

    # ---- changes vs pre-audit per-cell file (for downstream phases)
    changes = {}
    if PRE_AUDIT.exists():
        old = pd.read_csv(PRE_AUDIT, usecols=["cell_id", "cnv_call", "consensus_call", "chr17p_loss", "chr4q12_gain"], low_memory=False).set_index("cell_id")
        new = calls.set_index("cell_id")
        common = new.index.intersection(old.index)
        o, n = old.loc[common], new.loc[common]
        ch = pd.DataFrame({"study": n.study, "sample_id": n.sample_id, "author_state": n.author_state, "author_call": n.author_call,
                           "old_cnv": o.cnv_call.astype(str), "new_cnv": n.cnv_call.astype(str), "old_cons": o.consensus_call.astype(str), "new_cons": n.consensus_call.astype(str),
                           "new_conf": n.cnv_confidence})
        ch["cnv_changed"] = ch.old_cnv != ch.new_cnv; ch["cons_changed"] = ch.old_cons != ch.new_cons
        rows = []
        for study, g in ch.groupby("study"):
            r = dict(study=study, n_cells=len(g), n_cells_missing_in_old=int((calls.study == study).sum() - len(g)),
                     n_cnv_call_changed=int(g.cnv_changed.sum()), frac_cnv_call_changed=float(g.cnv_changed.mean()),
                     n_consensus_changed=int(g.cons_changed.sum()), frac_consensus_changed=float(g.cons_changed.mean()))
            for (a, b), v in g[g.cnv_changed].groupby(["old_cnv", "new_cnv"]).size().items():
                r[f"cnv_{a}_to_{b}"] = int(v)
            for (a, b), v in g[g.cons_changed].groupby(["old_cons", "new_cons"]).size().items():
                r[f"consensus_{a}_to_{b}"] = int(v)
            lab = g[g.cons_changed & g.author_call.notna()]
            r["consensus_changed_with_author_label"] = len(lab)
            r["consensus_changed_frac_author_malignant"] = float((lab.author_call == "malignant").mean()) if len(lab) else np.nan
            rows.append(r)
        tot = dict(study="all", n_cells=len(ch), n_cnv_call_changed=int(ch.cnv_changed.sum()), frac_cnv_call_changed=float(ch.cnv_changed.mean()),
                   n_consensus_changed=int(ch.cons_changed.sum()), frac_consensus_changed=float(ch.cons_changed.mean()))
        rows.append(tot)
        pd.DataFrame(rows).to_csv(RESULTS / "tables" / "cnv_audit_call_changes.csv", index=False)
        comp = (ch[ch.cons_changed].groupby(["study", "old_cons", "new_cons", ch.author_state.fillna("no_author_label")]).size().rename("n_cells").reset_index()
                .rename(columns={"author_state": "author_state"}).sort_values(["study", "n_cells"], ascending=[True, False]))
        comp.to_csv(RESULTS / "tables" / "cnv_audit_call_changes_by_author_label.csv", index=False)
        changes = dict(n_cnv_call_changed=tot["n_cnv_call_changed"], n_consensus_changed=tot["n_consensus_changed"], n_cells_compared=len(ch))
        log.write(f"vs pre-audit: cnv_call changed {tot['n_cnv_call_changed']}, consensus changed {tot['n_consensus_changed']} of {len(ch)}")

    # ---- per-cluster evidence vs author labels (CNV-20: scripted) and score-floor sensitivity (CNV-09)
    cl_rows = []
    smeta = samp_meta.set_index("sample_id")
    for s in STUDIES:
        d = INTERIM / "cnv" / s
        if not d.exists():
            continue
        for f in sorted(d.glob("*.clusters.csv")):
            sid = f.name.replace(".clusters.csv", "")
            cl = pd.read_csv(f); cl["cnv_cluster"] = pd.to_numeric(cl.cnv_cluster, errors="coerce").astype("Int64").astype(str)
            g = calls[calls.sample_id == sid]
            gc = g.groupby(pd.to_numeric(g.cnv_cluster, errors="coerce").astype("Int64").astype(str))
            au_frac = gc.author_call.apply(lambda x: float((x == "malignant").mean()) if x.notna().any() else np.nan)
            au_n = gc.author_call.apply(lambda x: int(x.notna().sum()))
            cl["sample_id"] = sid; cl["study"] = s; cl["ref_cnv_mean"] = smeta.loc[sid, "ref_cnv_mean"] if sid in smeta.index else np.nan
            cl["score_rel"] = cl.mean_cnv_score - cl.ref_cnv_mean
            cl["n_own_cells"] = cl.cnv_cluster.map(gc.size()).fillna(0).astype(int)
            cl["n_author_labelled"] = cl.cnv_cluster.map(au_n).fillna(0).astype(int); cl["author_malignant_frac"] = cl.cnv_cluster.map(au_frac)
            cl_rows.append(cl)
    clusters = pd.concat(cl_rows, ignore_index=True)
    clusters.to_csv(RESULTS / "qc" / "cnv_cluster_stats_vs_author.csv", index=False)
    lab_cl = clusters[clusters.n_author_labelled >= 20]
    pct = dict(normal_cluster_score_rel_p90=float(lab_cl[lab_cl.author_malignant_frac <= 0.1].score_rel.quantile(0.9)),
               tumour_cluster_score_rel_p10=float(lab_cl[lab_cl.author_malignant_frac >= 0.9].score_rel.quantile(0.1)), n_clusters_labelled=int(len(lab_cl)))
    sens = []
    has_ev = clusters.events.fillna("").astype(str).str.len().gt(0) | clusters.any_arm_flag.fillna(False).astype(bool)
    for floor in [0.004, 0.006, 0.008, 0.010, 0.012, 0.016, 0.020]:
        elev = clusters.score_rel > floor
        call = np.where(clusters.frac_reference >= C.CALL_PARAMS["ref_frac_guard"], elev & has_ev, elev | has_ev)
        for s, idx in clusters.groupby("study").groups.items():
            c = call[idx]; sub = clusters.loc[idx]
            labd = sub.n_author_labelled > 0
            tp = float((sub.author_malignant_frac[labd & c] * sub.n_author_labelled[labd & c]).sum()); npos = float(sub.n_author_labelled[labd & c].sum())
            fn_ = float((sub.author_malignant_frac[labd & ~c] * sub.n_author_labelled[labd & ~c]).sum())
            sens.append(dict(study=s, score_floor=floor, n_cells_cnv_malignant=int(sub.n_own_cells[c].sum()),
                             frac_cells_cnv_malignant=float(sub.n_own_cells[c].sum() / sub.n_own_cells.sum()),
                             precision_malignant_vs_author=tp / npos if npos else np.nan, recall_malignant_vs_author=tp / (tp + fn_) if (tp + fn_) else np.nan))
    pd.DataFrame(sens).to_csv(RESULTS / "tables" / "cnv_score_floor_sensitivity.csv", index=False)

    # ---- figures (house style: no titles in panels, legends outside, 300 dpi; captions in docs/phases/phase3_cnv_notes.md)
    figdir = RESULTS / "qc" / "cnv"; figdir.mkdir(parents=True, exist_ok=True)
    apply_figure_style()
    cmap = {"reference": "#1b9e77", "pooled-ref": "#66c2a5", "normal": "#7570b3", "malignant": "#d95f02"}
    for sid in HEATMAP_SAMPLES:
        study = sid.split("__")[0]
        f = INTERIM / "cnv" / study / f"{sid}.cnv.npz"
        if not f.exists():
            log.write(f"heatmap: {f} missing"); continue
        z = np.load(f, allow_pickle=True)
        X = z["X_cnv"].astype(np.float32); pooled = z["pooled"]; cl = z["cnv_cluster"].astype(str); ref = z["ref_type"].astype(str)
        cc = calls.set_index("cell_id").loc[z["cell_id"][~pooled], "cnv_call"]
        call = np.array(["pooled-ref"] * len(cl), dtype=object); call[~pooled] = cc.values
        call[(ref != "none") & ~pooled] = "reference"
        order = np.lexsort((cl, pd.Categorical(call, list(cmap)).codes))
        wchrom = z["window_chrom"].astype(str)
        bounds = np.flatnonzero(wchrom[1:] != wchrom[:-1]) + 1
        starts = np.r_[0, bounds]; ends = np.r_[bounds, len(wchrom)]
        fig, (ax0, ax) = plt.subplots(1, 2, figsize=(7.2, 4.2), gridspec_kw=dict(width_ratios=[0.04, 1], wspace=0.02))
        vmax = float(np.percentile(np.abs(X), 99)) or 0.05
        im = ax.imshow(X[order], aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
        for b in bounds:
            ax.axvline(b - 0.5, color="k", lw=0.3)
        ticks = (starts + ends) / 2
        labels = [wchrom[s_].replace("chr", "") for s_ in starts]
        # thin crowded labels (CNV-18): keep a label only if its tick is >= 2.5 % of the windows away from the last kept one
        kept, last = [], -np.inf
        for t, l in zip(ticks, labels):
            if t - last >= 0.025 * len(wchrom):
                kept.append(l); last = t
            else:
                kept.append("")
        ax.set_xticks(ticks); ax.set_xticklabels(kept, fontsize=5); ax.set_yticks([])
        ax.set_xlabel("chromosome (vertical lines = chromosome boundaries; unlabelled blocks are the small chromosomes in order)")
        ax0.set_ylabel(f"{(~pooled).sum()} cells + {pooled.sum()} pooled reference")
        ax0.imshow(np.array([[list(cmap).index(c)] for c in call[order]]), aspect="auto", cmap=matplotlib.colors.ListedColormap(list(cmap.values())), vmin=0, vmax=3)
        ax0.set_xticks([]); ax0.set_yticks([]); ax0.set_xlabel("call", fontsize=6)
        for sp in ax0.spines.values():
            sp.set_visible(False)
        cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label(f"smoothed deviation (clipped at 99th pct |x| = {vmax:.3f})")
        handles = [matplotlib.patches.Patch(color=v, label=k) for k, v in cmap.items()]
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=4, frameon=False, fontsize=6, title="call strip", title_fontsize=6)
        fig.savefig(figdir / f"cnv_heatmap_{sid}.png", dpi=300, bbox_inches="tight"); plt.close(fig)
    # Jessa scatter: larger class first, subsample stated in the legend (CNV-19)
    if len(jes):
        fig, ax = plt.subplots(figsize=(4.2, 3.6))
        n_plot = min(60000, len(jes)); sub = jes.sample(n_plot, random_state=0)
        order_lab = sub.author_label.value_counts().index.tolist()   # draw the larger class first
        cols = {"Malignant": "#d95f02", "Normal": "#7570b3"}
        for lab in order_lab:
            m = sub.author_label == lab
            ax.scatter(sub.loc[m, "marker_tumour_score"], sub.loc[m, "cnv_score"], s=2, alpha=0.3, c=cols.get(lab, "grey"),
                       label=f"author {lab}: {int(m.sum()):,} plotted of {int((jes.author_label == lab).sum()):,}", rasterized=True)
        ax.set_xlabel("marker tumour-program score (mean OPC/AC/MES)"); ax.set_ylabel("cnv_score (mean |CNV|)")
        ax.legend(frameon=False, markerscale=5, fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1.0))
        fig.savefig(figdir / "jessa_cnvscore_vs_markerscore.png", dpi=300, bbox_inches="tight"); plt.close(fig)

    summary = dict(n_cells=int(len(calls)), per_study=calls.groupby("study").consensus_call.value_counts().unstack(fill_value=0).to_dict("index"),
                   per_study_cnv_confidence=calls.groupby("study").cnv_confidence.value_counts(dropna=False).unstack(fill_value=0).rename(columns=str).to_dict("index"),
                   tuning_samples=TUNING_SAMPLES, n_low_confidence_samples=int(per_sample.low_confidence_sample.sum()),
                   low_confidence_samples=per_sample.loc[per_sample.low_confidence_sample, "sample_id"].tolist(),
                   cluster_stats_percentiles=pct, changes_vs_pre_audit=changes)
    (RESULTS / "tables" / "phase3_summary.json").write_text(json.dumps(summary, indent=1, default=int))
    log.write("aggregate done")


if __name__ == "__main__":
    main()
