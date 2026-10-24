#!/usr/bin/env python
"""Phase 8 / step 31 — supplementary figures (house style) and generated captions.

Companion of 30_figures.py: same style and helpers (dmg_atlas.figstyle) - bold UPPERCASE panel letters outside the axes and a short
descriptive title centred above every panel (Cell Reports house style, venue model 43:114804) - every caption number computed from
results/tables/. Drawn from tables only; the exploratory per-phase PNGs in results/figures/ are not reused.

Inputs (results/tables/; md5 recorded in docs/supp_figures.md):
  S1 malignant calls : malignant_call_validation.csv, malignant_calls_per_sample.csv (05 / 06)
  S2 integration     : integration_metrics_summary.csv (08)
  S3 mapping         : mapping_uncertainty.csv (11), malignant_nearest_fetal_class.csv (11b)
  S4 metaprograms    : metaprogram_genes.csv, metaprogram_reference_overlap.csv, metaprogram_cut_sensitivity.csv (16),
                       jessa_cnmf_recovery.csv (19)
  S5 masked MPs      : metaprograms_assay_masked.csv, metaprogram_masking_summary.csv (20)
  S6 TME markers     : tme_subcluster_markers_clean.csv (22)
  S7 spatial         : spatial_niches.csv, spatial_colocalisation.csv, spatial_colocalisation_per_section.csv,
                       spatial_colocalisation_sensitivity.csv (15b, audit SP-01 … SP-12). brainstem_mp_vs_fallback.csv (21) is no longer
                       read: the score_genes-argmax fallback reference it compares against was dropped by the audit (SP-06), so the
                       comparison is superseded; panel c now shows the leave-brainstem-out sensitivity instead.
  S8 bulk survival   : openpedcan_deconvolution.csv, survival_results.csv (25; NNLS variable sets, pre-audit — not re-run),
                       survival_rankscore_audit.csv (25b; the corrected rank-score fits)
Outputs (never overwritten if present, unless --force):
  results/pub_figs/supp/figS<N>_<name>.png (300 dpi) and .pdf for N = 1..8
  docs/supp_figures.md   captions + input md5 table
Run: nice -n 10 python scripts/31_supp_figures.py [--force] [--figs 1 .. 8]

Choices fixed before plotting
  S1: sample-level agreement uses tier2 (own CNV call) and tier3 (consensus) vs author labels; samples with >= 100 compared cells.
  S3: fetal cell-of-origin by the Jessa 2022 author state of malignant cells (>= 100 cells) and the overall malignant set.
  S4: top 8 genes per MP; reference overlap shown as -log10 BH q of the hypergeometric test, capped at 20.
  S6: top 4 marker genes per cleaned subcluster by Wilcoxon score; dot size = fraction of cells expressing, colour = log2 fold change.
  S7: per-section rho for the six state x microenvironment pairs with q < 0.05 and the largest pooled rho in fig 5d.
  S8: every adjusted Cox model (rank scores, NNLS proportions and shares, both signature sets); NNLS rows are shown to document
      why they were not trusted, not as findings.
"""
import sys, time, hashlib, argparse
import re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from dmg_atlas import ROOT, RESULTS, LOGS  # noqa: E402
from dmg_atlas.figstyle import (apply_figure_style, QUERY_COLOURS, QUERY_NAMES, W, FAMILY_COLOURS, FAMILY_NAMES,  # noqa: E402
                                STUDY_NAMES, STUDY_COLOURS, TME_NAME, COMPARTMENT_COLOURS, SPATIAL_CLASS, letter, legend_dots, heat)

log = U.StepLogger(LOGS / "phase8_release.log", context="31_supp_figures")
TAB, FIG, DOC = RESULTS / "tables", RESULTS / "pub_figs/supp", ROOT / "docs/supp_figures.md"
TIER_NAMES = {"tier2_cnv_call": "own CNV call", "tier3_cnv_marker_agree": "CNV + marker agree"}  # consensus uses author labels: not a test
used, captions = {}, {}


def T(name: str, **kw) -> pd.DataFrame:
    used[name] = TAB / name
    return pd.read_csv(TAB / name, **kw)


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)
    log.write(f"wrote {(FIG / name).relative_to(ROOT)}.png/.pdf")


def strip(ax, groups, values, colours, k=0):
    for i, (g, v) in enumerate(zip(groups, values)):
        ax.scatter(np.full(len(v), i) + np.random.default_rng(k + i).uniform(-0.2, 0.2, len(v)), v, s=6, c=colours[i], linewidths=0, zorder=3)
        if len(v):
            ax.hlines(np.median(v), i - 0.3, i + 0.3, color="black", lw=0.8, zorder=4)


# ---------------------------------------------------------------- S1
def s1():
    v = T("malignant_call_validation.csv"); cs = T("malignant_calls_per_sample.csv")
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.8), gridspec_kw=dict(width_ratios=[1.1, 1, 1.2], wspace=0.45))
    # Audit 2026-10-20: the CNV track added a 'subset' column (all_samples / tuning_samples / held_out_samples) to
    # malignant_call_validation.csv, which made (tier, study) non-unique and broke the reindex below. The panel keeps its
    # original meaning - all samples - and the tuning/held-out split stays available in the table and the CNV notes.
    st = v[(v.level == "study") & (v.subset == "all_samples") & v.tier.isin(TIER_NAMES) & v.f1_malignant.notna()]
    assert not st.duplicated(["tier", "study"]).any(), "validation table still has duplicate (tier, study) rows"
    ax = axs[0]; tiers = list(TIER_NAMES); studies = sorted(st.study.unique(), key=lambda s: (s not in STUDY_NAMES, s)); x = np.arange(len(studies))
    for i, t in enumerate(tiers):
        r = st[st.tier == t].set_index("study").reindex(studies)
        ax.bar(x + (i - 0.5) * 0.38, r.f1_malignant, 0.38, color=["#969696", "#252525"][i], label=TIER_NAMES[t])
    ax.set_xticks(x); ax.set_xticklabels([STUDY_NAMES.get(s, QUERY_NAMES.get(s, s).split(" (")[0]) for s in studies], rotation=45, ha="right")
    ax.set_ylabel("F1, malignant vs author label"); ax.set_ylim(0, 1.0)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5, columnspacing=0.6, handlelength=1)
    letter(ax, "a", dx=-0.25, title="Malignant-call F1 per study", tpad=11)
    ax = axs[1]
    sm = v[(v.level == "sample") & (v.n_compared >= 100)]
    vals = [sm.loc[sm.tier == t, "agreement"].dropna().values for t in tiers]
    strip(ax, [TIER_NAMES[t] for t in tiers], vals, ["#969696", "#252525"])
    ax.set_xticks([0, 1]); ax.set_xticklabels([TIER_NAMES[t].replace(" + ", "\n+ ") for t in tiers]); ax.set_ylabel("per-sample agreement with author label")
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(None, 1.01)
    letter(ax, "b", dx=-0.25, title="Agreement with author label")
    ax = axs[2]
    cs = cs[cs.study.isin(STUDY_NAMES)]
    order = sorted(cs.study.unique(), key=lambda s: -cs.loc[cs.study == s, "frac_malignant_consensus"].median())
    strip(ax, order, [cs.loc[cs.study == s, "frac_malignant_consensus"].values for s in order], [STUDY_COLOURS[s] for s in order], k=10)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([STUDY_NAMES[s] for s in order], rotation=45, ha="right")
    ax.set_ylabel("consensus-malignant fraction per sample"); ax.set_ylim(0, 1)
    letter(ax, "c", dx=-0.25, title="Consensus-malignant fraction")
    save(fig, "figS1_malignant_calls")
    c3 = st[st.tier == "tier3_cnv_marker_agree"]; c2 = st[st.tier == "tier2_cnv_call"]
    # audit CNV-08: the headline agreement must be reported at the sample level with a CI, not only cell-pooled.
    allr = v[(v.study == "all") & (v.level == "study") & (v.subset == "all_samples") & (v.tier == "tier2_cnv_call")].iloc[0]
    captions[1] = (
        "**Supplementary Figure 1 | Malignant-cell calls.** "
        "**a**, F1 of the author-independent malignant calls against author labels per study: own inferCNV call, and cells where the CNV "
        f"and marker calls agree. Own CNV call F1 {c2.f1_malignant.min():.2f}–{c2.f1_malignant.max():.2f} across {len(c2)} cohorts with "
        f"author labels; CNV + marker agree {c3.f1_malignant.min():.2f}–{c3.f1_malignant.max():.2f}. The consensus call used downstream is "
        "not shown: it requires agreement with the author label where one exists, so it matches author labels by construction. "
        f"**b**, Per-sample agreement with author labels (samples with ≥ 100 compared cells; n = {len(vals[0])} / {len(vals[1])}; bars: medians "
        f"{np.median(vals[0]):.3f} / {np.median(vals[1]):.3f}). The sample is the unit to quote: over all "
        f"{int(allr.n_samples)} labelled samples the own-CNV call agrees with author labels at a median "
        f"{allr.sample_median:.3f} (IQR {allr.sample_q25:.3f}–{allr.sample_q75:.3f}; unweighted mean {allr.sample_mean:.3f}, "
        f"95 % CI {allr.sample_mean_ci95_lo:.3f}–{allr.sample_mean_ci95_hi:.3f} over 2,000 sample bootstraps), against a cell-pooled "
        f"{allr.agreement:.3f} (audit CNV-08). The CNV thresholds were tuned on 22 of these samples; tuning and held-out samples agree to "
        "within 0.003 pooled (audit CNV-06/CNV-07). "
        f"**c**, Fraction of consensus-malignant cells per reference sample, by study (n = {len(cs)} samples; bars: medians).")


# ---------------------------------------------------------------- S2
def s2():
    im = T("integration_metrics_summary.csv")
    cols = [c for c in im.columns if c not in ("embedding", "chosen", "rank_total", "elbo_validation_best", "batch_score_raw_mean",
                                                "bio_score_raw_mean") and im[c].dtype.kind == "f"]
    m = im.set_index("embedding")[cols].T
    m = m.loc[sorted(m.index, key=lambda c: (not c.endswith("_score"), c))]
    fig, ax = plt.subplots(figsize=(W * 0.6, 5.6))
    heat(ax, m, cmap="viridis", vmax=1, fmt=lambda v: f"{v:.2f}", cbar_label="metric value (scaled 0–1 by scIB)", fig=fig)
    ax.set_title("scIB metrics per embedding", loc="center", fontsize=7.5, pad=3)  # single-panel figure: title, no panel letter
    for t in ax.get_xticklabels():
        if t.get_text() == im.loc[im.chosen, "embedding"].iloc[0]:
            t.set_fontweight("bold")
    save(fig, "figS2_integration_metrics")
    ch = im[im.chosen].iloc[0]
    captions[2] = (
        "**Supplementary Figure 2 | Integration benchmark (scIB).** All batch-correction metrics (batch keys sample, study, assay) and "
        f"biological-conservation metrics (coarse lineage and fine labels) for the {len(im)} candidate embeddings; summary scores on top. "
        f"Bold: chosen reference ({ch.embedding}, overall {ch.total_score:.3f}, rank {int(ch.rank_total)}). Assay remains the hardest batch key "
        f"(kBET assay {ch.kbet__assay:.2f}, iLISI assay {ch.ilisi__assay:.3f} for the chosen model).")


# ---------------------------------------------------------------- S3
def s3():
    mu = T("mapping_uncertainty.csv"); nf = T("malignant_nearest_fetal_class.csv")
    fig, axs = plt.subplots(1, 2, figsize=(W, 3.2), gridspec_kw=dict(width_ratios=[0.8, 1.4], wspace=0.6))
    ax = axs[0]; qs = ["filbin", "labelle", "liu", "braun", "siletti"]
    strip(ax, qs, [mu.loc[mu["query"] == q, "frac_classifier_knn_agree"].values for q in qs], [QUERY_COLOURS[q] for q in qs], k=40)
    ax.set_xticks(range(len(qs))); ax.set_xticklabels([QUERY_NAMES[q].split(" (")[0] for q in qs], rotation=45, ha="right")
    ax.set_ylabel("per sample: classifier label = kNN majority label"); ax.set_ylim(0, 1.02)
    letter(ax, "a", dx=-0.25, title="Classifier vs kNN agreement")
    ax = axs[1]
    d = nf[((nf.grouping == "jessa_author_state") & (nf.n_malignant_cells >= 100)) | (nf.grouping == "all")].copy()
    fc = [c for c in nf.columns if c.startswith("frac_") and c not in ("frac_top_fetal_CellClass", "frac_second_fetal_CellClass",
                                                                      "frac_top_fetal_cell_type")]
    fc = [c for c in fc if d[c].max() >= 0.02]
    m = d.set_index("group")[fc]; m.columns = [c[5:].replace("_", " ") for c in fc]
    m.index = ["all malignant" if g == "all" else f"{g} ({n:,})" for g, n in zip(m.index, d.n_malignant_cells)]
    heat(ax, m, cmap="Purples", vmax=1, fmt=lambda v: f"{v:.2f}" if v >= 0.05 else "", cbar_label="fraction of malignant cells", fig=fig)
    ax.set_xlabel("nearest Braun 2023 fetal class (kNN in scANVI space)")
    letter(ax, "b", dx=-0.45, title="Nearest fetal hindbrain class")
    save(fig, "figS3_mapping_origin")
    al = nf[nf.grouping == "all"].iloc[0]
    captions[3] = (
        "**Supplementary Figure 3 | Mapping confidence and fetal cell-of-origin.** "
        "**a**, Per-sample fraction of query cells whose scANVI classifier label equals the majority label of their 15 nearest reference "
        f"cells ({len(mu)} samples; bars: medians " + ", ".join(f"{QUERY_NAMES[q].split(' (')[0]} {mu.loc[mu['query'] == q, 'frac_classifier_knn_agree'].median():.2f}"
                                                              for q in ["filbin", "labelle", "liu", "braun", "siletti"]) + "). "
        f"Median classifier uncertainty (1 − max probability) is ≤ {mu.median_uncertainty_prob.max():.3f} in every sample, so the classifier "
        "is confident even where the neighbourhood disagrees (fetal hindbrain). "
        "**b**, Nearest fetal hindbrain class of malignant cells (majority of 15 Braun 2023 neighbours), overall and by Jessa 2022 author "
        f"state (≥ 100 cells). Overall {al.frac_top_fetal_CellClass:.1%} of malignant cells sit nearest {al.top_fetal_CellClass.lower()} and "
        f"{al.frac_second_fetal_CellClass:.1%} nearest {al.second_fetal_CellClass.lower()}.")


# ---------------------------------------------------------------- S4
def s4():
    g = T("metaprogram_genes.csv"); ro = T("metaprogram_reference_overlap.csv"); cut = T("metaprogram_cut_sensitivity.csv")
    jr = T("jessa_cnmf_recovery.csv"); fam = T("metaprogram_families.csv").set_index("metaprogram")
    fig = plt.figure(figsize=(W, 7.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1], height_ratios=[1.5, 1], hspace=0.5, wspace=0.45)
    mps = sorted(g.metaprogram.unique(), key=lambda m: int(m[2:]))

    # audit (reporting pass 2026-10-21): the fallback label was the bare word "excluded"; it now names the technical class the
    # programme was excluded for (MP2 ribosomal, MP12 mitochondrial, MP13 IEG/dissociation stress; audit MC-07 / MC-08).
    def fname(x):
        f = fam.loc[x, "family"]
        if isinstance(f, str) and f in FAMILY_NAMES:
            return FAMILY_NAMES[f]
        tc = fam.loc[x, "technical_class"]
        # kept short: panel a prints this next to the gene list, so a long suffix collides with it.
        return {"technical:ribosomal": "ribosomal*", "technical:mitochondrial": "mitochondrial*",
                "technical:IEG_dissociation_stress": "IEG stress*"}.get(tc, str(tc).replace("technical:", "")) \
            if isinstance(tc, str) else "unassigned"
    ax = fig.add_subplot(gs[0, 0]); ax.axis("off")
    top = g[g["rank"] <= 8].sort_values(["metaprogram", "rank"]).groupby("metaprogram").gene.agg(", ".join)
    for i, m in enumerate(mps):
        f = fam.loc[m, "family"]
        ax.add_patch(plt.Rectangle((0, 1 - (i + 1) / len(mps)), 0.012, 1 / len(mps) * 0.85, color=FAMILY_COLOURS.get(f, "#000"),
                                   transform=ax.transAxes, clip_on=False))
        ax.text(0.02, 1 - (i + 0.55) / len(mps), f"{m} {fname(m)}", fontsize=5.5, fontweight="bold",
                transform=ax.transAxes, va="center")
        # gene list starts clear of the longest row label ("MP12 mitochondrial*"), which collided with it at x=0.33
        ax.text(0.38, 1 - (i + 0.55) / len(mps), top[m], fontsize=4.5, transform=ax.transAxes, va="center", style="italic")
    letter(ax, "a", dx=0.0, title="Top genes per metaprogram")
    ax = fig.add_subplot(gs[0, 1])
    q = ro.pivot_table(index="metaprogram", columns="reference_set", values="q_bh").reindex(mps)
    q = q.loc[:, (q < 0.05).any()]
    # The full reference-set keys ran off the bottom of the canvas and were clipped mid-word; abbreviate the source
    # prefixes (the caption spells each one out) so every label fits.
    PREFIX = [("inhouseDE:Liu2022_DE_author_states_rederived_", "Liu DE "), ("published:Liu2022_TableS2a_", "Liu S2a "),
              ("filbin_untraceable:", "Filbin untrac. "), ("curated:", "curated "), ("Liu2022_", "Liu "),
              ("Filbin2018_", "Filbin "), ("state:", "")]
    def refname(c):
        for a, b in PREFIX:
            c = c.replace(a, b)
        return c.replace("_", " ")
    m = -np.log10(q.clip(lower=1e-20)); m.columns = [refname(c) for c in m.columns]
    heat(ax, m.where(m > -np.log10(0.05)), cmap="Greens", vmax=20, cbar_label="−log10 q (overlap)", fig=fig)
    ax.set_xticklabels(m.columns, rotation=90, ha="center", fontsize=5)
    letter(ax, "b", dx=-0.2, title="Overlap with published gene sets")
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(cut.cut, cut.n_metaprograms, "o-", color="#525252", ms=3, label="metaprograms")
    ax.plot(cut.cut, cut.n_replicated, "o-", color="#d95f02", ms=3, label="replicated")
    ax2 = ax.twinx(); ax2.plot(cut.cut, cut.n_programmes_assigned, "s--", color="#9e9ac8", ms=3, lw=0.8)
    ax2.set_ylabel("programmes assigned", color="#756bb1"); ax2.spines["right"].set_visible(True); ax2.tick_params(axis="y", colors="#756bb1")
    ax.axvline(0.85, color="#bdbdbd", lw=0.6, ls="--", zorder=0); ax.set_xlabel("dendrogram cut"); ax.set_ylabel("count")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5, columnspacing=0.6, handlelength=1)  # audit P7-17
    letter(ax, "c", dx=-0.15, title="Metaprograms across dendrogram cuts", tpad=11)
    ax = fig.add_subplot(gs[1, 1])
    fams = [f for f in FAMILY_COLOURS if f in set(jr.best_mp_family)]
    strip(ax, fams, [jr.loc[jr.best_mp_family == f, "best_mp_rho"].values for f in fams], [FAMILY_COLOURS[f] for f in fams], k=20)
    ax.axhline(0.5, color="#bdbdbd", lw=0.6, ls="--", zorder=0)
    ax.set_xticks(range(len(fams))); ax.set_xticklabels([FAMILY_NAMES[f] for f in fams], rotation=45, ha="right")
    ax.set_ylabel("best usage ρ with an MP")
    letter(ax, "d", dx=-0.2, title="Recovery of Jessa 2022 programmes")
    save(fig, "figS4_metaprogram_detail")
    rec = (jr.best_mp_rho >= 0.5).mean()
    captions[4] = (
        "**Supplementary Figure 4 | Metaprogram content and robustness.** "
        f"**a**, Top 8 genes of each of the {len(mps)} metaprograms (rank by consensus across member programmes); bar: family. "
        f"A trailing * on the row label marks a technical programme: {int(fam.technical_class.notna().sum())} carry a technical class and "
        "no family — "
        + ", ".join(f"{m} ({str(fam.loc[m, 'technical_class']).replace('technical:', '').replace('_', ' ')})"
                    for m in fam.index[fam.technical_class.notna()])
        + " — and are excluded from the lineage count (audit MC-07, MC-08); the pre-audit immune-contamination programme MP16 "
        "(DOCK8/APBB1IP) no longer forms with the corrected malignant selection. "
        f"**b**, Overlap of MP top-50 genes with published DMG/glioma state gene sets (hypergeometric test, BH q; blank q ≥ 0.05; "
        f"{m.shape[1]} sets with any significant overlap). "
        "**c**, Number of metaprograms, replicated metaprograms and programmes assigned across dendrogram cuts (dashed: primary cut 0.85). "
        f"**d**, Jessa 2022 author cNMF programmes ({len(jr)} programmes, {jr.sample_id.nunique()} samples): best within-sample Spearman ρ of "
        f"per-cell usage with an atlas MP, by the family of the best MP; {rec:.0%} reach ρ ≥ 0.5 (dashed line).")


# ---------------------------------------------------------------- S5
def s5():
    mm = T("metaprograms_assay_masked.csv"); mk = T("metaprogram_masking_summary.csv")
    fig, axs = plt.subplots(1, 2, figsize=(W, 3.4), gridspec_kw=dict(width_ratios=[1.6, 1], wspace=0.5))
    ax = axs[0]; y = np.arange(len(mm))[::-1]
    ax.barh(y, -mm.n_samples_cells, color="#d95f02", height=0.7, label="cell samples")
    ax.barh(y, mm.n_samples_nuclei, color="#1b9e77", height=0.7, label="nucleus samples")
    lab = [f"{r.metaprogram} [{r.original_mps.split(';')[0].split(':')[0]}]" + (" *" if r.replicated else "") for r in mm.itertuples()]
    ax.set_yticks(y); ax.set_yticklabels(lab); ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("samples (cells ← | → nuclei)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5, columnspacing=0.6, handlelength=1)  # audit P7-17
    xt = ax.get_xticks(); ax.set_xticks(xt); ax.set_xticklabels([f"{abs(t):.0f}" for t in xt])
    letter(ax, "a", dx=-0.2, title="Cell and nucleus samples per mMP", tpad=11)
    ax = axs[1]
    w = mk.pivot_table(index="cut", columns="genes", values="n_mixed_assay")
    for gname, col in [("unmasked", "#525252"), ("masked", "#d95f02")]:
        ax.plot(w.index, w[gname], "o-", color=col, ms=3, label=f"{gname} genes")
    ax.axvline(0.85, color="#bdbdbd", lw=0.6, ls="--", zorder=0)
    ax.set_xlabel("dendrogram cut"); ax.set_ylabel("MPs with both cell and nucleus members")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5, columnspacing=0.6, handlelength=1)  # audit P7-17
    letter(ax, "b", dx=-0.25, title="Mixed-assay metaprograms per cut", tpad=11)
    save(fig, "figS5_masked_metaprograms")
    p = mk[mk.cut == 0.85].set_index("genes")
    captions[5] = (
        "**Supplementary Figure 5 | Metaprograms after masking assay-biased genes.** "
        f"**a**, The {len(mm)} masked metaprograms (mMPs) at the primary cut: number of contributing cell (left) and nucleus (right) samples; "
        f"bracket: most frequent original MP of its members; * replicated ({int(mm.replicated.sum())}: "
        + ", ".join(f"{r.metaprogram} [{r.original_mps.split(';')[0].split(':')[0]}]" for r in mm[mm.replicated].itertuples()) + "). "
        f"**b**, MPs containing both cell- and nucleus-derived programmes across cuts, with all genes and after masking "
        f"({int(p.loc['unmasked', 'n_mixed_assay'])} → {int(p.loc['masked', 'n_mixed_assay'])} at 0.85).")


# ---------------------------------------------------------------- S6
def s6():
    mk = T("tme_subcluster_markers_clean.csv")
    top = mk.sort_values("scores", ascending=False).groupby("subcluster").head(4)
    subs = [s for s in TME_NAME if s in set(top.subcluster)]
    genes = list(dict.fromkeys(top.set_index("subcluster").loc[subs, "names"]))
    sub = mk[mk.names.isin(genes) & mk.subcluster.isin(subs)]
    fig, ax = plt.subplots(figsize=(W, 4.2))
    xi = {g: i for i, g in enumerate(genes)}; yi = {s: i for i, s in enumerate(subs)}
    sc = ax.scatter(sub.names.map(xi), sub.subcluster.map(yi), s=sub.pct_nz_group * 30, c=sub.logfoldchanges.clip(-3, 3), cmap="RdBu_r",
                    vmin=-3, vmax=3, linewidths=0)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, style="italic", fontsize=5)
    ax.set_yticks(range(len(subs))); ax.set_yticklabels([TME_NAME[s] for s in subs]); ax.invert_yaxis()
    for i, s in enumerate(subs):
        ax.add_patch(plt.Rectangle((-1.6, i - 0.5), 0.5, 1, color=COMPARTMENT_COLOURS[s.rsplit("_c", 1)[0]], clip_on=False, lw=0))
    ax.set_xlim(-1.7, len(genes) - 0.5)
    cb = fig.colorbar(sc, ax=ax, fraction=0.02, pad=0.01); cb.set_label("log2 fold change (vs compartment)"); cb.outline.set_visible(False)
    for f in (0.25, 0.5, 1.0):
        ax.scatter([], [], s=f * 30, c="#737373", label=f"{f:.0%}")
    ax.legend(title="cells expressing", loc="upper left", bbox_to_anchor=(1.08, 0.3), fontsize=5, title_fontsize=5)
    ax.set_title("Subcluster marker genes", loc="center", fontsize=7.5, pad=3)  # single-panel figure: title, no panel letter
    save(fig, "figS6_tme_markers")
    captions[6] = (
        "**Supplementary Figure 6 | Marker genes of the cleaned microenvironment subclusters.** Top 4 genes per subcluster (Wilcoxon "
        f"score within compartment; {len(genes)} genes, {len(subs)} subclusters); dot size = fraction of the subcluster's cells expressing "
        "the gene, colour = log2 fold change against the other subclusters of the compartment (clipped at ±3); left bar: compartment. "
        "Names follow docs/phase6_notes.md §11.")


# ---------------------------------------------------------------- S7
def s7():
    ni = T("spatial_niches.csv"); ps = T("spatial_colocalisation_per_section.csv")
    co = T("spatial_colocalisation.csv"); sens = T("spatial_colocalisation_sensitivity.csv")
    pmap = T("spatial_patient_map.csv")
    n_bs = int((pmap.cohort == "brainstem").sum())
    fig = plt.figure(figsize=(W, 6.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.4, 1], hspace=0.55, wspace=0.5)
    ax = fig.add_subplot(gs[0, :])
    cls = [c for c in ni.columns if c.startswith("malignant_") or c in ("myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular")]
    m = ni.sort_values("frac_spots", ascending=False).set_index("niche")[cls]
    m.index = [f"{n} ({f:.0%})" for n, f in zip(m.index, ni.sort_values("frac_spots", ascending=False).frac_spots)]
    m.columns = [SPATIAL_CLASS.get(c, c) for c in cls]
    heat(ax, m.T, cmap="Reds", vmax=m.values.max(), fmt=lambda v: f"{v:.2f}" if v >= 0.15 else "", cbar_label="mean spot proportion", fig=fig)
    ax.set_xticklabels(m.index, rotation=0, ha="center"); ax.set_xlabel("niche (share of DMG spots)")
    letter(ax, "a", dx=-0.08, title="Spot niche composition")
    ax = fig.add_subplot(gs[1, 0])
    # audit SP-01: spatial_colocalisation.csv carries both the primary proportion metric and the CLR sensitivity; select the primary.
    d = co[(co.metric == "proportion") & (co.group == "DMG") & (co.pair_type == "state_vs_TME")
           & (co.q_bh_state_vs_TME < 0.05)].nlargest(6, "rho_re")
    pairs = list(zip(d.class_a, d.class_b))
    p2 = ps[(ps.group == "DMG")]
    vals = [p2[((p2.class_a == a) & (p2.class_b == b)) | ((p2.class_a == b) & (p2.class_b == a))].rho.values for a, b in pairs]
    names = [f"{SPATIAL_CLASS.get(a, a)} – {SPATIAL_CLASS.get(b, b)}" if a.startswith("malignant") else f"{SPATIAL_CLASS.get(b, b)} – {SPATIAL_CLASS.get(a, a)}"
             for a, b in pairs]
    strip(ax, names, vals, ["#525252"] * len(names), k=30)
    ax.axhline(0, color="#bdbdbd", lw=0.6, zorder=0)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=35, ha="right"); ax.set_ylabel("per-section Spearman ρ")
    letter(ax, "b", dx=-0.15, title="Per-section co-localisation")
    # audit (reporting pass 2026-10-21): panel c used to compare the metaprogram reference against the score_genes-argmax
    # fallback reference, which SP-06 dropped from the pipeline; that comparison is superseded and its table
    # (brainstem_mp_vs_fallback.csv, 2026-10-18) was never regenerated. The panel now shows the sensitivity that matters for
    # the unit of replication: dropping the 10-section brainstem autopsy, which is one patient carrying 66 % of the DMG spots.
    ax = fig.add_subplot(gs[1, 1])
    prim = co[(co.metric == "proportion") & (co.group == "DMG") & (co.subset == "all")].set_index(["class_a", "class_b"])
    lbo = sens[(sens.subset == "leave_brainstem_out") & (sens.group == "DMG")].set_index(["class_a", "class_b"])
    key = prim.index.intersection(lbo.index)
    pt = prim.loc[key, "pair_type"].values
    ax.scatter(prim.loc[key, "rho_re"], lbo.loc[key, "rho_re"], s=14, linewidths=0,
               c=np.where(pt == "state_vs_TME", "#b2182b", "#525252"))
    legend_dots(ax, {"state_vs_TME": "#b2182b", "state_vs_state": "#525252"},
                {"state_vs_TME": "state × microenvironment", "state_vs_state": "state × state"},
                loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1, fontsize=5)
    lim = [min(prim.loc[key, "rho_re"].min(), lbo.loc[key, "rho_re"].min()) - 0.05,
           max(prim.loc[key, "rho_re"].max(), lbo.loc[key, "rho_re"].max()) + 0.05]
    ax.plot(lim, lim, color="#bdbdbd", lw=0.6, ls="--", zorder=0); ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.axhline(0, color="#e0e0e0", lw=0.5, zorder=0); ax.axvline(0, color="#e0e0e0", lw=0.5, zorder=0)
    ax.set_xlabel(f"pooled ρ, all {int(prim.n_patients.max())} DMG patients")
    ax.set_ylabel(f"pooled ρ, brainstem autopsy excluded\n({int(lbo.n_patients.max())} patients)")
    letter(ax, "c", dx=-0.42, title="Brainstem-autopsy sensitivity", tpad=18)
    save(fig, "figS7_spatial_detail")
    rho = pd.Series(prim.loc[key, "rho_re"].values).corr(pd.Series(lbo.loc[key, "rho_re"].values), method="spearman")
    dmax = float(np.abs(prim.loc[key, "rho_re"].values - lbo.loc[key, "rho_re"].values).max())
    captions[7] = (
        "**Supplementary Figure 7 | Spatial niches, per-section co-localisation and the brainstem-autopsy sensitivity.** Computed from the "
        "audit-corrected spatial tables (track 'spatial', SP-01 … SP-12, regenerated 2026-10-20): proportions rather than q05 abundances, a "
        "single metaprogram-family reference for all four cohorts, and patient-level pooling. "
        f"**a**, Mean class proportions of the {len(ni)} spot niches (Leiden on cell2location proportions of DMG spots; label: share of DMG "
        "spots; sections containing each niche: " + ", ".join(f"{n} {k}" for n, k in zip(ni.sort_values("frac_spots", ascending=False).niche,
                                                                ni.sort_values("frac_spots", ascending=False).n_sections_present)) + "). "
        f"**b**, Per-section spot-level Spearman ρ for the six strongest significant state × microenvironment pairs of Fig. 5d "
        f"({p2.section_id.nunique()} DMG sections; bars: medians). Every dot is a section, not a patient: {n_bs} of these sections are serial "
        "sections of the single brainstem autopsy, which is why the pooled estimates in Fig. 5d are taken over patients. "
        f"**c**, The same pooled ρ with and without that autopsy, over all {len(key)} pairs "
        f"(Spearman {rho:.2f} between the two; largest shift {dmax:.3f}; dashed: identity). Dropping the autopsy removes "
        f"{int(prim.n_patients.max()) - int(lbo.n_patients.max())} patient and {int(prim.n_sections.max()) - int(lbo.n_sections.max())} of the "
        f"{int(prim.n_sections.max())} sections, and leaves the ranking of the pairs intact. The earlier panel c — section tumour fraction "
        "under the score_genes-argmax fallback reference vs the metaprogram reference — was dropped because SP-06 removed the fallback "
        "reference from the pipeline, so `brainstem_mp_vs_fallback.csv` (2026-10-18) is superseded and was not regenerated.")


# ---------------------------------------------------------------- S8
def s8():
    # audit P7-07 / P7-12: the rank-score fits were re-derived by 25b against the corrected metaprogram families (the pre-audit
    # table still scores the superseded novel_MP15 programme) and carry a Schoenfeld PH test. The NNLS / deconvolution variable
    # sets were NOT re-run, so they are taken from the pre-audit table and labelled as such in the caption.
    dc = T("openpedcan_deconvolution.csv")
    sv_old = T("survival_results.csv"); sv_new = T("survival_rankscore_audit.csv")
    sv = pd.concat([sv_new.assign(provenance="corrected (25b)"),
                    sv_old[sv_old.variable_set != "rank_score"].assign(provenance="pre-audit (25)")], ignore_index=True)
    fig, axs = plt.subplots(1, 2, figsize=(W, 4.6), gridspec_kw=dict(width_ratios=[1, 1.3], wspace=0.75))
    ax = axs[0]
    cols = [c for c in dc.columns if c.startswith("mp:") and ":share" not in c and not c.endswith("malignant_total")]
    cols = sorted(cols, key=lambda c: -dc[c].median())
    ax.boxplot([dc[c].values for c in cols], orientation="horizontal", widths=0.6, showfliers=False, medianprops=dict(color="black", lw=0.8),
               boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    ax.set_yticks(range(1, len(cols) + 1)); ax.set_yticklabels([SPATIAL_CLASS.get(c[3:], c[3:]) for c in cols]); ax.invert_yaxis()
    ax.set_xlabel("NNLS proportion per patient (metaprogram signatures)")
    letter(ax, "a", dx=-0.4, title="NNLS proportions per patient")
    ax = axs[1]
    s = sv[sv.model == "adjusted"].copy()
    order = ["rank_score", "mp_proportion", "mp_share", "mpM_proportion", "mpM_share", "fallback_proportion", "tme_proportion", "tmeM_proportion"]
    s["o"] = s.variable_set.map({k: i for i, k in enumerate(order)}); s = s.sort_values(["o", "variable"]).reset_index(drop=True)
    y = np.arange(len(s))[::-1]
    col = np.where(s.variable_set == "rank_score", "black", "#969696")
    lo, hi = 0.5, 2.0; out = (s.hr_per_sd < lo) | (s.hr_per_sd > hi)
    ok = ~out
    ax.errorbar(s.hr_per_sd[ok], y[ok], xerr=[(s.hr_per_sd - s.ci_lo)[ok], (s.ci_hi - s.hr_per_sd)[ok]], fmt="none", ecolor="#bdbdbd", elinewidth=0.5)
    ax.scatter(s.hr_per_sd[ok], y[ok], c=col[ok], s=np.where(s.q_bh[ok] < 0.05, 14, 6), linewidths=0, zorder=3)
    ax.scatter(np.full(out.sum(), lo * 1.02), y[out], marker="<", c=col[out], s=10, linewidths=0, zorder=3)
    ax.axvline(1, color="#bdbdbd", lw=0.6, zorder=0); ax.set_xscale("log"); ax.minorticks_off(); ax.set_xlim(lo, hi)
    ax.set_xticks([0.5, 0.7, 1, 1.4, 2]); ax.set_xticklabels(["0.5", "0.7", "1", "1.4", "2"])
    ax.set_yticks(y); ax.set_yticklabels([f"{v.split(':', 1)[1]}  [{k}]" for v, k in zip(s.variable, s.variable_set)], fontsize=4.5)
    ax.set_xlabel("hazard ratio per SD (adjusted for age + pons)")
    letter(ax, "b", dx=-0.5, title="Hazard ratio per variable set")
    save(fig, "figS8_bulk_survival")
    sig = s[s.q_bh < 0.05]
    tn = dc["mp:T_NK"].median()
    captions[8] = (
        "**Supplementary Figure 8 | OpenPedCan bulk deconvolution and all survival models.** *Exploratory: the bulk-survival layer (code "
        "Phase 7, build-plan Phase 8) was not part of the approved plan (audit P7-16).* "
        f"**a**, NNLS proportions of the metaprogram-family cell2location signatures in {len(dc)} OpenPedCan H3 K27-altered DMG primary "
        f"tumours (unit: patient). The fit is implausible: median T/NK proportion {tn:.2f}, and several malignant families are zero in most "
        "patients. The NNLS layer was not re-run after the audit, so it still rests on the pre-audit signatures and family definition "
        "(audit P7-07, open). "
        f"**b**, Adjusted Cox hazard ratio per SD for every variable set (unit: patient; n = {int(s.n.iloc[0])}, {int(s.events.iloc[0])} "
        "deaths; rank = rank-based family scores [black, used in Fig. 6c — these rows are the corrected 25b fits against the re-derived "
        "families]; mp / mpM = NNLS with all genes / marker genes, proportion or malignant share; fallback = Filbin/Liu-state signatures; "
        "tme / tmeM = microenvironment proportions — all of the latter are pre-audit fits, retained to document why the NNLS readout was "
        "distrusted rather than as findings). Large dots: BH q < 0.05 within the set "
        f"({len(sig)}: " + ", ".join(f"{v} [{k}] HR {h:.2f}" for v, k, h in zip(sig.variable, sig.variable_set, sig.hr_per_sd)) +
        "); none is supported by the rank scores. Arrows: degenerate fits outside the axis (" +
        ", ".join(f"{v} [{k}] HR {h:.2g}" for v, k, h in zip(s.variable[out], s.variable_set[out], s.hr_per_sd[out])) +
        "; the marker-gene MES-like proportion is zero in most patients).")


FIGS = {1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7, 8: s8}
NAMES = {1: "figS1_malignant_calls", 2: "figS2_integration_metrics", 3: "figS3_mapping_origin", 4: "figS4_metaprogram_detail",
         5: "figS5_masked_metaprograms", 6: "figS6_tme_markers", 7: "figS7_spatial_detail", 8: "figS8_bulk_survival"}


def write_doc():
    old = {}
    if DOC.exists():
        for blk in DOC.read_text().split("\n\n"):
            if blk.startswith("**Supplementary Figure "):
                old[int(blk.split()[2])] = blk
    caps = {**old, **captions}
    lines = ["# Supplementary figures — captions", "", f"Generated by `scripts/31_supp_figures.py` on {time.strftime('%Y-%m-%d')}. Do not edit by "
             "hand; every number is computed from `results/tables/`. Files: `results/pub_figs/supp/figS<N>_<name>.png` (300 dpi) and `.pdf`.", ""]
    lines += [caps[k] + "\n" for k in sorted(caps)]
    # audit P7-06: a partial run (--figs ...) used to drop the inputs of figures not re-drawn. Merge the existing Inputs table
    # (md5 recomputed from the current file so the table always describes the tables on disk) with the inputs of this run.
    inputs = {}
    if DOC.exists():
        for m in re.finditer(r"^\| `([^`]+)` \| `([0-9a-f]{32})` \|$", DOC.read_text(), flags=re.M):
            q = ROOT / m.group(1)
            inputs[m.group(1)] = md5(q) if q.exists() else m.group(2) + " (file missing)"
    inputs.update({str(p.relative_to(ROOT)): md5(p) for _, p in used.items()})
    lines += ["## Inputs", "", "| table | md5 |", "|---|---|"] + [f"| `{k}` | `{v}` |" for k, v in sorted(inputs.items())]
    DOC.write_text("\n".join(lines) + "\n")
    log.write(f"wrote {DOC.relative_to(ROOT)} ({len(caps)} captions, {len(inputs)} inputs; {len(used)} read in this run)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--figs", nargs="+", type=int, default=list(FIGS), choices=list(FIGS))
    args = ap.parse_args()
    apply_figure_style()
    FIG.mkdir(parents=True, exist_ok=True)
    for k in args.figs:
        if (FIG / f"{NAMES[k]}.png").exists() and not args.force:
            log.write(f"{NAMES[k]}.png exists, skipping (use --force)"); continue
        with log.step(f"supplementary figure {k} ({NAMES[k]})"):
            FIGS[k]()
    if captions:
        write_doc()
    log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
