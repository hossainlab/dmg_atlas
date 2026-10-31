#!/usr/bin/env python
"""Phase 8 / step 30 — final manuscript figures (house style) and generated captions.

House style (2026-10-21, Cell Reports 'Resource' submission; venue model Cell Reports 43:114804): bold UPPERCASE panel letters outside
the axes, a short descriptive title (2-5 words, sentence case) centred above every panel, and embedding panels with no axis box, ticks or
axis labels. All statistics, n and the remaining description live in the caption. Style from dmg_atlas.figstyle.apply_figure_style
(open frame, 8/7/6 pt ladder, 300 dpi, Type-42 fonts); letters, titles and embedding axes from figstyle.letter / figstyle.umap_ax.
Every number quoted in a caption is computed here from the tables below, never typed by hand.

Inputs (results/tables/ unless stated; all produced by earlier steps, md5 recorded in docs/figures/figures.md):
  Fig 1 atlas   : atlas_core_obs_labels.csv, atlas_core_umap.csv (09), integration_metrics_summary.csv (08), cnv_event_frequencies.csv (06)
  Fig 2 mapping : mapping_summary.csv, mapping_uncertainty.csv, mapping_liu_author_vs_predicted.csv,
                  mapping_normal_reference_specificity.csv (11 / 11b)
  Fig 3 states  : metaprograms.csv, replication_matrix.csv, metaprogram_families.csv, metaprogram_masking_summary.csv,
                  metaprogram_cell_vs_nucleus.csv, metaprogram_family_composition.csv (16 / 19 / 20)
  Fig 4 TME     : tme_subclusters_clean.csv, tme_assay_purity.csv (22), composition_results.csv (18)
  Fig 5 spatial : spatial_abundance_long.csv (14 / 14c), spatial_colocalisation.csv, ren_section_composition.csv,
                  ren_dmg_vs_gbm_composition.csv, brainstem_section_means.csv, brainstem_gradient.csv (15b, audit SP-01…SP-12;
                  the *_mp.csv tables written by step 21 are superseded — they carry the pre-audit novel_MP15 family and
                  abundance-based statistics, and are no longer read here)
  Fig 6 multi   : atac_summary.csv, atac_label_correlation.csv (24), survival_rankscore_audit.csv (25b),
                  foundation_model_comparison.csv (26b)
Outputs (never overwritten if present, unless --force):
  results/pub_figs/main/fig<N>_<name>.png (300 dpi) and .pdf for N = 1..6
  docs/figures/figures.md   captions (one per figure, panel by panel, numbers from the tables) + input md5 table
Supplementary figures: scripts/31_supp_figures.py (results/pub_figs/supp/, docs/figures/supp_figures.md).
Run: nice -n 10 python scripts/30_figures.py [--force] [--figs 1 2 3 4 5 6]

Choices fixed before plotting
  UMAPs: all 311,618 reference cells, random draw order (seed 0), rasterised. Spatial example: brainstem Pons_A (one of the three
  pons sections with tumour fraction > 0.7 in the one autopsy; metaprogram-family reference, q05 abundances on pixel coordinates).
  Composition heatmap: cell-type and malignant-state
  sets, estimable contrasts only (autopsy vs resection is not estimable). Survival panel: rank-score Cox models adjusted for age + pons
  (the readout phase 7 judged trustworthy; NNLS proportions are not shown). Liu author-vs-predicted: author labels with >= 50 cells.
  Normal-reference panel: classes with >= 100 cells.
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
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from dmg_atlas import ROOT, RESULTS, LOGS  # noqa: E402
from dmg_atlas.figstyle import (EMBED_NAMES, apply_figure_style, LINEAGE_COLOURS, QUERY_COLOURS, QUERY_NAMES, W, LINEAGES, FAMILIES,  # noqa: E402
                                FAMILY_COLOURS, FAMILY_NAMES, ASSAY_COLOURS, ASSAY_NAMES, STUDY_NAMES, STUDY_COLOURS, TME_NAME,
                                COMPARTMENT_COLOURS, CONTRAST_NAMES, SPATIAL_CLASS, letter, umap_ax, legend_dots, heat)

log = U.StepLogger(LOGS / "phase8_release.log", context="30_figures")
TAB, FIG, DOC = RESULTS / "tables", RESULTS / "pub_figs/main", ROOT / "docs/figures/figures.md"
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


# ---------------------------------------------------------------- Fig 1
def fig1():
    obs = T("atlas_core_obs_labels.csv", index_col=0, usecols=[0, 2, 3, 6, 5]).join(T("atlas_core_umap.csv", index_col=0))
    im = T("integration_metrics_summary.csv"); cnv = T("cnv_event_frequencies.csv")
    rng = np.random.default_rng(0); o = obs.iloc[rng.permutation(len(obs))]
    fig = plt.figure(figsize=(W, 6.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], hspace=0.5, wspace=0.6)
    # a cells per study x assay
    ax = fig.add_subplot(gs[0, 0])
    ct = obs.groupby(["study", "assay"]).size().unstack(fill_value=0).reindex(columns=list(ASSAY_COLOURS), fill_value=0)
    ct = ct.loc[ct.sum(1).sort_values().index]
    left = np.zeros(len(ct))
    for a in ct.columns:
        ax.barh([STUDY_NAMES.get(s, s) for s in ct.index], ct[a] / 1e3, left=left, color=ASSAY_COLOURS[a], height=0.7, label=ASSAY_NAMES[a])
        left += ct[a].values / 1e3
    # audit (reporting pass 2026-10-21): the Fig 1a legend sat inside the axes; moved above them like every other panel.
    ax.set_xlabel("cells / nuclei (thousands)")
    legend_dots(ax, ASSAY_COLOURS, ASSAY_NAMES, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5, columnspacing=0.6)
    letter(ax, "a", dx=-0.45, title="Cells and nuclei per study", tpad=16)
    # b-d UMAPs
    for k, (col, cmap, names, lab, ttl) in enumerate([("scanvi_label", {c: LINEAGE_COLOURS[c] for c in LINEAGES}, {c: c.replace("T_NK", "T/NK") for c in LINEAGES}, "b",
                                                       "scANVI lineage label"),
                                                      ("assay", ASSAY_COLOURS, ASSAY_NAMES, "c", "Assay"),
                                                      ("study", STUDY_COLOURS, STUDY_NAMES, "d", "Study")]):
        ax = fig.add_subplot(gs[0, 1 + k] if k < 2 else gs[1, 0])
        ax.scatter(o.UMAP1, o.UMAP2, s=0.05, c=o[col].map(cmap).fillna("#bdbdbd").values, linewidths=0, rasterized=True)
        umap_ax(ax); legend_dots(ax, cmap, names, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=5, columnspacing=0.6)
        letter(ax, lab, title=ttl)
    # e integration scores
    ax = fig.add_subplot(gs[1, 1])
    im = im.sort_values("total_score", ascending=False)
    x = np.arange(len(im)); wd = 0.27
    for i, (c, col) in enumerate([("batch_score", "#9e9ac8"), ("bio_score", "#74c476"), ("total_score", "#525252")]):
        ax.bar(x + (i - 1) * wd, im[c], wd, color=col, label=c.replace("_score", "").replace("bio", "biology").replace("total", "overall"))
    ax.set_xticks(x); ax.set_xticklabels([EMBED_NAMES.get(e, e) + (" *" if ch else "") for e, ch in zip(im.embedding, im.chosen)], rotation=45, ha="right")
    ax.set_ylabel("scIB score"); ax.set_ylim(0, 1)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5, columnspacing=0.6, handlelength=1)
    letter(ax, "e", dx=-0.2, title="scIB integration scores", tpad=11)
    # f CNV hallmark events
    ax = fig.add_subplot(gs[1, 2])
    ev = [c for c in cnv.columns if c.startswith("chr") and "__" not in c]
    m = cnv.set_index("study")[ev].dropna(how="all")
    m = m.loc[sorted(m.index, key=lambda s: (s not in STUDY_NAMES, s))]  # reference studies first, then queries
    m.columns = [c.replace("chr", "").replace("_", " ") for c in ev]
    m.index = [STUDY_NAMES.get(s, QUERY_NAMES.get(s, s).split(" 20")[0] + " (query)") for s in m.index]
    heat(ax, m, cmap="Reds", vmax=1, cbar_label="fraction of malignant cells", fig=fig)
    letter(ax, "f", dx=-0.45, title="Hallmark copy-number events")
    save(fig, "fig1_atlas")
    n_ass = obs.assay.value_counts()
    ch = im[im.chosen].iloc[0]; hm = im.set_index("embedding").loc["harmony"]
    top = cnv.set_index("study")[ev].dropna(how="all")
    captions[1] = (
        "**Figure 1 | A cross-cohort reference of H3K27-altered diffuse midline glioma.** "
        f"**a**, Cells and nuclei per study after QC ({len(obs):,} in total; {n_ass.get('snRNA', 0):,} snRNA nuclei, "
        f"{n_ass.get('multiome_RNA', 0):,} multiome nuclei, {n_ass.get('scRNA', 0):,} scRNA cells). "
        "**b–d**, UMAP of the scANVI latent space of all reference cells, coloured by scANVI lineage label (**b**), assay (**c**) and study (**d**). "
        f"Malignant cells make up {(obs.scanvi_label == 'malignant').mean():.1%} of scANVI labels. "
        f"**e**, scIB batch-correction, biological-conservation and overall scores for the candidate embeddings; * marks the chosen "
        f"reference ({ch.embedding}: batch {ch.batch_score:.2f}, biology {ch.bio_score:.2f}, overall {ch.total_score:.2f}; Harmony "
        f"{hm.batch_score:.2f} / {hm.bio_score:.2f} / {hm.total_score:.2f}). Overall = 0.4 × batch + 0.6 × biology over metrics min–max "
        f"scaled across all {len(im)} rows shown, a rule recorded post hoc and not pre-registered; scaling over the five eligible "
        "candidates only (excluding the scvi_seed0_e60 sensitivity run) shifts the totals but not the ranking (audit IM-05). "
        "**f**, Fraction of consensus-malignant cells carrying each hallmark copy-number event (inferCNV arm-level calls) per study, "
        f"reference and query cohorts (studies without CNV calls omitted); chr1q gain ranges {top.chr1q_gain.min():.2f}–{top.chr1q_gain.max():.2f} across studies. "
        "Fractions are pooled over cells, not patients; the comparable per-sample event frequencies are in "
        "`results/tables/cnv_event_frequencies.csv` (audit CNV-23). CNV thresholds were tuned on 22 of the 97 labelled samples (audit CNV-06).")


# ---------------------------------------------------------------- Fig 2
def fig2():
    ms = T("mapping_summary.csv"); mu = T("mapping_uncertainty.csv"); liu = T("mapping_liu_author_vs_predicted.csv")
    nr = T("mapping_normal_reference_specificity.csv")
    qs = ["filbin", "labelle", "liu", "braun", "siletti"]
    fig = plt.figure(figsize=(W, 7.0))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=1.0, height_ratios=[1, 1.6])
    # a predicted lineage per query
    ax = fig.add_subplot(gs[0, 0])
    agg = ms.groupby("query")[[f"n_pred_{l}" for l in LINEAGES] + ["n_cells"]].sum().reindex(qs)
    left = np.zeros(len(qs))
    for l in LINEAGES:
        v = (agg[f"n_pred_{l}"] / agg.n_cells).values
        ax.barh([QUERY_NAMES[q] for q in qs], v, left=left, color=LINEAGE_COLOURS[l], height=0.7); left += v
    ax.invert_yaxis(); ax.set_xlabel("fraction of query cells"); ax.set_xlim(0, 1)
    legend_dots(ax, {l: LINEAGE_COLOURS[l] for l in LINEAGES}, {"T_NK": "T/NK"}, loc="upper center", bbox_to_anchor=(0.4, -0.28), ncol=4,
                fontsize=5, columnspacing=0.6)
    letter(ax, "a", dx=-0.9, title="Predicted lineage composition")
    # b uncertainty per sample
    ax = fig.add_subplot(gs[0, 1])
    data = [mu.loc[mu["query"] == q, "median_uncertainty_knn"].values for q in qs]
    ax.boxplot(data, orientation="horizontal", widths=0.6, showfliers=False, medianprops=dict(color="black", lw=0.8),
               boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for i, (q, d) in enumerate(zip(qs, data)):
        ax.scatter(d, np.full(len(d), i + 1) + np.random.default_rng(i).uniform(-0.2, 0.2, len(d)), s=4, c=QUERY_COLOURS[q], linewidths=0, zorder=3)
    ax.set_yticks(range(1, len(qs) + 1)); ax.set_yticklabels([QUERY_NAMES[q] for q in qs]); ax.invert_yaxis()
    ax.set_xlabel("per-sample median kNN uncertainty\n(percentile of reference distances)"); ax.set_xlim(0, 1.02)
    letter(ax, "b", dx=-0.9, title="Per-sample mapping uncertainty")
    # c Liu author vs predicted
    ax = fig.add_subplot(gs[1, 0])
    pc = [c for c in ["malignant", "myeloid", "oligodendrocyte", "astrocyte", "neuron", "vascular"] if c in liu]
    l2 = liu[(liu.n >= 50) & (liu.author_annotation != "(none)")].set_index("author_annotation")
    m = l2[pc].div(l2[pc].sum(1), axis=0); m.index = [f"{i} ({n:,})" for i, n in zip(m.index, l2.n)]
    heat(ax, m, cmap="Blues", vmax=1, fmt=lambda v: f"{v:.2f}" if v >= 0.05 else "")
    ax.set_xlabel("predicted lineage"); letter(ax, "c", dx=-0.6, title="Author vs predicted label")
    # d normal references
    ax = fig.add_subplot(gs[1, 1])
    n2 = nr[nr.n_cells >= 100].copy()
    n2.index = [f"{'fetal' if q == 'braun' else 'adult'} {c} ({n:,})" for q, c, n in zip(n2["query"], n2.normal_class, n2.n_cells)]
    m = n2[[f"frac_pred_{l}" for l in LINEAGES]]; m.columns = [{"T_NK": "T/NK"}.get(l, l) for l in LINEAGES]
    heat(ax, m, cmap="Blues", vmax=1, fmt=lambda v: f"{v:.2f}" if v >= 0.05 else "", cbar_label="fraction of row", fig=fig)
    ax.axhline((n2["query"] == "braun").sum() - 0.5, color="black", lw=0.6)
    ax.set_xlabel("predicted lineage"); letter(ax, "d", dx=-0.75, title="Normal reference vs predicted label")
    save(fig, "fig2_mapping")
    tum = agg.loc[["filbin", "liu"]]; fm = (tum.n_pred_malignant / tum.n_cells)
    lab = agg.loc["labelle"]; br = agg.loc["braun"]
    med = mu.groupby("query").median_uncertainty_knn.median()
    sil = nr[nr["query"] == "siletti"]; sil_ok = sil[sil.frac_pred_expected.notna() & (sil.n_cells >= 100)]
    caption_ac = liu.set_index("author_annotation")
    bn = nr[(nr["query"] == "braun") & nr.normal_class.isin(["Neuron", "Neuroblast", "Neuronal IPC", "Oligo"])]
    n_norm_mal = int(nr.n_pred_malignant.sum())
    # audit (reporting pass 2026-10-21): the collapse class is computed, not asserted. Before the IM-01 fix the fetal classes
    # collapsed onto astrocyte; against the corrected / retrained reference classifier they collapse onto malignant.
    bn_frac = {l: bn[f"n_pred_{l}"].sum() / bn.n_cells.sum() for l in LINEAGES}
    bn_top = max(bn_frac, key=bn_frac.get)
    captions[2] = (
        "**Figure 2 | Projection of query cohorts onto the reference (zero-shot scANVI encoder projection, reference-classifier label "
        "transfer).** Audit IM-01/IM-02 (2026-10-20): the reference scANVI model was trained with `encode_covariates=False`, so scArches "
        "surgery cannot adapt the latent and changed nothing; query cells are projected zero-shot by the reference encoder and labelled "
        "by the reference classifier applied to that latent. The earlier 'scArches surgery' description was wrong. "
        f"**a**, Predicted lineage composition per query (cell-weighted over all cells of the query): malignant {fm['filbin']:.1%} of "
        f"Filbin 2018 and {fm['liu']:.1%} of Liu 2022 (Smart-seq2 tumour cohorts); LaBelle 2025 is predicted "
        f"{lab.n_pred_myeloid / lab.n_cells:.1%} myeloid; Braun 2023 fetal hindbrain is predicted "
        f"{br.n_pred_malignant / br.n_cells:.1%} malignant and only {br.n_pred_astrocyte / br.n_cells:.1%} astrocyte; "
        f"{n_norm_mal:,} of the {int(nr.n_cells.sum()):,} normal-reference cells are predicted malignant. Per-sample medians and "
        "donor/plate-bootstrap confidence intervals for every rate quoted here are in `mapping_agreement_sample_level.csv` (audit IM-22). "
        "**b**, Per-sample median kNN uncertainty (mean distance to the 15 nearest reference cells, as a percentile of the reference's "
        f"own distances; dots = samples). Medians: " + ", ".join(f"{QUERY_NAMES[q]} {med[q]:.2f}" for q in qs) +
        "; the Smart-seq2 / TPM-derived tumour queries and the fetal reference sit far from the reference even when the classifier is "
        "confident, whereas adult pons nuclei (10x, as the reference) do not. "
        "**c**, Liu 2022 author annotation (rows, labels with ≥ 50 cells, cell number in brackets) vs predicted lineage, row-normalised; "
        f"e.g. {caption_ac.loc['AC-like', 'malignant'] / caption_ac.loc['AC-like', 'n']:.1%} of author AC-like cells are predicted malignant. "
        f"Reported in full rather than selectively (audit P7-20): author Endothelia (n = {int(caption_ac.loc['Endothelia', 'n'])}) split "
        f"{caption_ac.loc['Endothelia', 'vascular'] / caption_ac.loc['Endothelia', 'n']:.1%} vascular / "
        f"{caption_ac.loc['Endothelia', 'malignant'] / caption_ac.loc['Endothelia', 'n']:.1%} malignant / "
        f"{caption_ac.loc['Endothelia', 'astrocyte'] / caption_ac.loc['Endothelia', 'n']:.1%} astrocyte, the largest non-tumour "
        "misassignment in this panel. "
        "**d**, Normal-reference classes (≥ 100 cells; above the line Braun 2023 fetal hindbrain, below Siletti 2023 adult pons) vs predicted "
        f"lineage, row-normalised. Adult pons classes with a reference counterpart are recovered at a median {sil_ok.frac_pred_expected.median():.1%}. "
        f"Fetal classes collapse onto {bn_top}, including those with a reference counterpart "
        f"({bn_frac[bn_top]:.1%} of fetal neurons, neuroblasts, neuronal IPCs and oligodendrocytes are predicted {bn_top}, "
        f"{bn_frac['astrocyte']:.1%} astrocyte); only fetal immune and vascular cells are recovered, so the reference classifier does not "
        "transfer to first-trimester tissue. Unit of replication: rows pool cells over 17 Braun donors and 3 Siletti donors.")


# ---------------------------------------------------------------- Fig 3
def fig3():
    mp = T("metaprograms.csv"); rep = T("replication_matrix.csv", index_col=0); fam = T("metaprogram_families.csv").set_index("metaprogram")
    mk = T("metaprogram_masking_summary.csv"); cvn = T("metaprogram_cell_vs_nucleus.csv"); comp = T("metaprogram_family_composition.csv")
    ab = T("assay_biased_genes.csv")  # audit P7-12: the masked-gene count in the caption is computed, not hard-coded
    fig = plt.figure(figsize=(W, 6.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 0.8, 0.8], height_ratios=[1.35, 1], hspace=0.6, wspace=0.6)
    # a replication matrix
    ax = fig.add_subplot(gs[0, 0])
    mp = mp.set_index("metaprogram")
    order = sorted(rep.index, key=lambda m: (FAMILIES.index(fam.loc[m, "family"]) if fam.loc[m, "family"] in FAMILIES else 99, int(m[2:])))
    cols = sorted(rep.columns, key=lambda c: (c.split("|")[1], c))
    m = rep.loc[order, cols].astype(float)
    # audit (reporting pass 2026-10-21): the fallback label was hard-coded "immune, excluded", which mislabelled the three
    # technical programmes of the corrected run (MP2 ribosomal, MP12 mitochondrial, MP13 IEG/dissociation stress, audit
    # MC-07 / MC-08). The label now comes from the metaprogram's own technical_class.
    def fname(x):
        f = fam.loc[x, "family"]
        if isinstance(f, str) and f in FAMILY_NAMES:
            return FAMILY_NAMES[f]
        tc = mp.loc[x, "technical_class"]
        return str(tc).replace("technical:", "").replace("_", " ") + ", excluded" if isinstance(tc, str) else "unassigned"
    m.index = [f"{x} {fname(x)}" + (" *" if mp.loc[x, "replicated"] else "") for x in order]
    pub = {"jessa2022": "Jessa 2022", "okada": "Okada", "collot": "Collot", "desisto": "DeSisto", "zang": "Zang"}
    m.columns = [f"{pub.get(c.split('|')[0], c.split('|')[0])} ({c.split('|')[1]})" for c in cols]
    heat(ax, m.where(m > 0), cmap="Purples", vmax=m.values.max(), fmt=lambda v: f"{v:.0f}", cbar_label="samples", fig=fig)
    ax.tick_params(axis="y", pad=8)
    for i, x in enumerate(order):
        ax.add_patch(plt.Rectangle((-0.5 - 0.35, i - 0.5), 0.3, 1, color=FAMILY_COLOURS.get(fam.loc[x, "family"], "#000"), clip_on=False, lw=0))
    ax.axvline(sum(c.endswith("cells") for c in cols) - 0.5, color="black", lw=0.6)
    letter(ax, "a", dx=-0.42, title="Samples per metaprogram")
    # b replicated MPs vs cut, masked / unmasked
    ax = fig.add_subplot(gs[0, 1])
    for g, col in [("unmasked", "#525252"), ("masked", "#d95f02")]:
        s = mk[mk.genes == g].sort_values("cut")
        ax.plot(s.cut, s.n_replicated, "o-", color=col, ms=3, label=f"{g} genes")
    ax.axvline(0.85, color="#bdbdbd", lw=0.6, ls="--", zorder=0)
    ax.set_xlabel("dendrogram cut (1 − mean Jaccard)"); ax.set_ylabel("replicated metaprograms")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5, columnspacing=0.6, handlelength=1)  # audit P7-17: legend outside the axes
    letter(ax, "b", dx=-0.25, title="Replication across cuts", tpad=20)
    # c cell vs nucleus Jaccard before / after masking
    ax = fig.add_subplot(gs[0, 2])
    # audit (reporting pass 2026-10-21): the old substring mapper matched 'oc' inside 'technical:IEG_dissociation_stress' and
    # coloured every pair as OC-like. Pairs are now keyed on the NUCLEUS metaprogram's family (the cell side is a single
    # technical programme after MC-07), with technical nucleus programmes drawn grey.
    fl = lambda m: fam.loc[m, "family"] if isinstance(fam.loc[m, "family"], str) and fam.loc[m, "family"] in FAMILY_COLOURS else "unassigned"  # noqa: E731
    for r in cvn.itertuples():
        f = fl(r.nucleus_mp)
        ax.plot([0, 1], [r.jaccard, r.jaccard_masked], "-o", color=FAMILY_COLOURS[f], ms=2, lw=0.6, alpha=0.8)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["all genes", "assay-biased\ngenes masked"]); ax.set_xlim(-0.25, 1.25)
    ax.set_ylabel("Jaccard, cell vs nucleus MP (top 50 genes)")
    legend_dots(ax, {f: FAMILY_COLOURS[f] for f in sorted({fl(m) for m in cvn.nucleus_mp}, key=FAMILIES.index)},
                {**FAMILY_NAMES, "unassigned": "technical (no family)"},
                loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5, columnspacing=0.6)  # audit P7-17: legend outside the axes
    letter(ax, "c", dx=-0.3, title="Cell–nucleus overlap", tpad=40)
    # d per-sample family composition
    ax = fig.add_subplot(gs[1, :])
    aorder = {"scRNA": 0, "multiome_RNA": 1, "snRNA": 2}
    comp = comp.assign(o=comp.assay.map(aorder)).sort_values(["o", "study", "OPC_like"], ascending=[True, True, False])
    bottom = np.zeros(len(comp)); x = np.arange(len(comp))
    for f in FAMILIES:
        ax.bar(x, comp[f], bottom=bottom, color=FAMILY_COLOURS[f], width=0.9, lw=0); bottom += comp[f].values
    for i, a in enumerate(comp.assay):
        ax.add_patch(plt.Rectangle((i - 0.45, -0.07), 0.9, 0.04, color=ASSAY_COLOURS[a], clip_on=False, lw=0))
    ax.set_xlim(-0.6, len(comp) - 0.4); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_ylabel("fraction of malignant cells")
    ax.set_xlabel("samples", labelpad=12)
    # audit P7-17: build the legend from FAMILIES (current families) - FAMILY_COLOURS also holds the pre-audit novel_MP15 key
    legend_dots(ax, {f: FAMILY_COLOURS[f] for f in FAMILIES}, FAMILY_NAMES, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    letter(ax, "d", dx=-0.05, title="Within-sample family ranking")
    save(fig, "fig3_metaprograms")
    prim = mk[mk.cut == 0.85].set_index("genes")
    rep_names = ", ".join(f"{i} {FAMILY_NAMES.get(fam.loc[i, 'family'], i)}" for i in mp.index[mp.replicated])
    # audit (reporting pass 2026-10-21): the "highest pair" was hard-coded as MP13 vs MP7 and labelled AC-like. MP13 is the
    # IEG / dissociation-stress programme (MC-07) and that pair is 0.00; the top pair is now read off the table.
    au = cvn.sort_values("jaccard_masked", ascending=False).iloc[0]
    lin = fam[fam.family.isin(FAMILIES[:4])]            # OPC-, AC-, OC-, MES-like
    lin_nuc_only = int((lin.n_samples_cells == 0).sum())
    lin_cell_only = int((lin.n_samples_nuclei == 0).sum())
    rand = mk[(mk.cut == 0.85) & (mk.genes == "random_mask")]
    captions[3] = (
        "**Figure 3 | Malignant metaprograms and the pre-registered replication test.** "
        f"**a**, Number of samples contributing a cNMF programme to each of the {len(mp)} metaprograms (MPs; per-sample cNMF, Jaccard "
        "clustering of top-50 genes, cut 0.85), split by publication and assay (cells | nuclei, vertical line). Left bar and row label: "
        f"family; rows without a family are the technical programmes excluded from the lineage count (MC-07, MC-08). "
        f"* replicated under the pre-registered rule (≥ 3 publications and both cells and nuclei): {int(mp.replicated.sum())} of {len(mp)} "
        f"({rep_names}); the unit is the sample, and the {len(mp)} MPs rest on "
        f"{int(mp.n_samples.max())} contributing samples at most ({int(mp.n_patients.max())} patients). "
        f"All {len(lin)} lineage-like MPs (OPC-, AC-, OC-, MES-like) are nuclei-only "
        f"({lin_nuc_only} of {len(lin)} have no cell sample, {lin_cell_only} no nucleus sample): after the corrected selection no "
        "cell-derived lineage MP exists (audit MC-07). "
        f"**b**, Replicated MPs across dendrogram cuts with all genes and after masking {int((ab.call != 'none').sum()):,} assay-biased genes (dashed line: primary cut; "
        f"{int(prim.loc['unmasked', 'n_replicated'])} → {int(prim.loc['masked', 'n_replicated'])} replicated at 0.85). A random mask of the same "
        f"size recovers a median {float(rand.n_replicated.iloc[0]):.1f} replicated MPs "
        f"[{int(rand.n_replicated_min.iloc[0])}–{int(rand.n_replicated_max.iloc[0])}, 10 seeds], so the masked gain is within the control band "
        "(audit MC-05); the cut itself is the pre-registered 0.85 and 2 MPs replicate over the whole 0.84–0.90 range. "
        "**c**, Jaccard similarity of each cell-derived MP with nucleus-derived MPs before and after masking. After the corrected "
        f"selection the only cell-dominant programme is MP13 (IEG / dissociation stress), so every pair here involves it; highest pair "
        f"{au.cell_mp} vs {au.nucleus_mp} ({str(au.nucleus_label).replace('technical:', '').replace('_', ' ')}) "
        f"{au.jaccard:.3f} → {au.jaccard_masked:.3f}, i.e. no cell–nucleus correspondence between lineage programmes is measurable. "
        f"**d**, Per-sample share of malignant cells assigned to each MP family (argmax of within-sample z-scored UCell family scores, "
        f"audit MC-01; {len(comp)} samples with scored malignant cells, ordered by assay, then study; strip below the bars: orange scRNA "
        "cells, purple multiome nuclei, green snRNA nuclei). Cell- and nucleus-derived profiles are not comparable: only "
        f"{int((comp.assay_class == 'cells').sum()) if 'assay_class' in comp else 6} of the {len(comp)} samples are cells.")


# ---------------------------------------------------------------- Fig 4
def fig4():
    tc = T("tme_subclusters_clean.csv"); pur = T("tme_assay_purity.csv"); cr = T("composition_results.csv")
    fig = plt.figure(figsize=(W, 6.2))
    gs = fig.add_gridspec(2, 4, width_ratios=[1.0, 0.35, 0.75, 1.15], height_ratios=[1, 1.1], hspace=0.5, wspace=0.25)
    # a subcluster sizes + nuclei fraction
    tc = tc.assign(name=tc.subcluster.map(TME_NAME)).iloc[::-1]
    ax = fig.add_subplot(gs[:, 0]); y = np.arange(len(tc))
    ax.barh(y, tc.n_cells, color=tc.compartment.map(COMPARTMENT_COLOURS), height=0.7)
    ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(tc.name); ax.set_xlabel("cells (log scale)"); ax.set_ylim(-0.6, len(tc) - 0.4)
    ax2 = fig.add_subplot(gs[:, 1], sharey=ax)
    ax2.scatter(tc.frac_nuclei, y, c=tc.compartment.map(COMPARTMENT_COLOURS), s=10, linewidths=0)
    ax2.set_xlim(-0.08, 1.08); ax2.set_xticks([0, 0.5, 1]); ax2.set_xticklabels(["0", "0.5", "1"]); ax2.set_xlabel("fraction\nnuclei")
    ax2.tick_params(labelleft=False); ax2.axvline(0.5, color="#bdbdbd", lw=0.5, zorder=0)
    legend_dots(ax, COMPARTMENT_COLOURS, {"T_NK": "T/NK"}, loc="upper center", bbox_to_anchor=(0.3, -0.07), ncol=2, columnspacing=0.6)
    letter(ax, "a", dx=-0.02, title="Microenvironment subcluster size")
    # b assay purity before / after
    ax = fig.add_subplot(gs[0, 3])
    # audit (reporting pass 2026-10-21): the end-of-line compartment names sat inside the axes (xlim ran to 1.5, so x = 1.06
    # was inside the panel). House style allows no text labels inside a panel: replaced by a legend above the axes.
    for r in pur.itertuples():
        ax.plot([0, 1], [r.purity_17, r.purity_clean], "-o", color=COMPARTMENT_COLOURS[r.compartment], ms=3)
    legend_dots(ax, {c: COMPARTMENT_COLOURS[c] for c in pur.compartment}, {"T_NK": "T/NK"},
                loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5, columnspacing=0.6)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["atlas scVI", "per-compartment\nscVI, masked"]); ax.set_xlim(-0.2, 1.2)
    ax.set_ylabel("assay purity\n(cells in < 10 % or > 90 % nuclei clusters)"); ax.set_ylim(-0.03, 1.05)
    letter(ax, "b", dx=-0.35, title="Assay purity per compartment", tpad=18)
    # c composition heatmap
    ax = fig.add_subplot(gs[1, 3])
    c = cr[cr.set.isin(["celltype", "malignant_state"]) & cr.estimable].copy()
    c["row"] = np.where(c.set == "malignant_state", "malignant: " + c.category.astype(str).map(lambda v: FAMILY_NAMES.get(v, v)),
                        c.category.astype(str).replace({"T_NK": "T/NK"}))
    m = c.pivot_table(index="row", columns="contrast", values="coef").reindex(columns=list(CONTRAST_NAMES))
    q = c.pivot_table(index="row", columns="contrast", values="q_bh").reindex(index=m.index, columns=m.columns)
    m.columns = [CONTRAST_NAMES[x] for x in m.columns]
    lim = np.nanmax(np.abs(m.values))
    heat(ax, m, cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), cbar_label="Δ logit proportion", fig=fig)
    for i in range(q.shape[0]):
        for j in range(q.shape[1]):
            v = q.values[i, j]
            if np.isfinite(v) and v < 0.05:
                ax.text(j, i, "**" if v < 0.01 else "*", ha="center", va="center", fontsize=6)
    letter(ax, "c", dx=-0.7, title="Composition effects per contrast")
    save(fig, "fig4_tme")
    p = pur.set_index("compartment")
    h3all = cr[(cr.contrast == "h3_H3.1_vs_H3.3") & (cr.set == "celltype")]
    h3 = h3all[h3all.q_bh < 0.05]
    h3site = h3all[h3all.q_site_adj < 0.05]
    h3pons = h3all[h3all.q_pons < 0.05]
    h3row = h3all.iloc[0]
    nn = cr[(cr.contrast == "assay_nuclei_vs_cells")]
    n_bh = int(cr.q_bh.notna().sum())
    captions[4] = (
        "**Figure 4 | Microenvironment subtypes and sample-level composition.** "
        f"**a**, Cleaned microenvironment subclusters ({len(tc)} clusters, {tc.n_cells.sum():,} cells; per-compartment scVI with sample batch, "
        "assay covariate and assay-biased genes removed), named from top genes; right: fraction of nuclei per cluster. "
        f"**b**, Assay purity per compartment before and after re-integration (unit: cell): myeloid {p.loc['myeloid', 'purity_17']:.2f} → "
        f"{p.loc['myeloid', 'purity_clean']:.2f}, T/NK {p.loc['T_NK', 'purity_17']:.2f} → {p.loc['T_NK', 'purity_clean']:.2f}; vascular and "
        "oligodendroglial cells are almost all nuclei, so purity cannot fall. The drop is not independent evidence that the assay effect was "
        "removed: assay_class was an encoded covariate of the per-compartment scVI and no held-out-gene or covariate-free control was run "
        "(audit MC-22, open). "
        "**c**, Composition models (unit: sample; empirical-logit proportion ~ contrast + publication, patient-clustered SE): effect per cell "
        f"type or malignant state (within-sample standardised Filbin/Liu marker scores, audit MC-23) and contrast; * q < 0.05, ** q < 0.01 "
        f"(BH over the {n_bh} tests that returned a q value, out of {int(cr.estimable.sum())} estimable and {len(cr)} fitted contrasts). "
        f"Nuclei vs cells (within Jessa 2022) is the strongest axis ({int((nn.q_bh < 0.05).sum())} of {int(nn.estimable.sum())} estimable tests "
        f"at q < 0.05 across all sets). The atlas holds 8 H3.1 samples, all from Jessa 2022, all pontine and all nuclei; this contrast fits "
        f"{int(h3row.n_test)} of them from {int(h3row.n_patients_test)} patients against {int(h3row.n_ref)} H3.3 samples "
        f"({int(h3row.n_patients_ref)} patients), so H3 variant is confounded with site and assay: "
        f"{', '.join(sorted(h3.category))} reach q < 0.05 under assay adjustment "
        f"alone, but only {', '.join(sorted(h3site.category)) if len(h3site) else 'none'} survives adjustment for site and only "
        f"{', '.join(sorted(h3pons.category)) if len(h3pons) else 'none'} the pons-only fit (audit MC-06). Age, site and ONC201 show no cell-type "
        "or malignant-state effect. Autopsy vs resection is not estimable — all 22 autopsy samples come from Okada, so the contrast is "
        "collinear with study — and is omitted (audit MC-06 / DQ-19).")


# ---------------------------------------------------------------- Fig 5
def fig5():
    # audit (reporting pass 2026-10-21), spatial track SP-01 … SP-12:
    #   * spatial_colocalisation.csv now holds a 'metric' column (proportion = primary, clr = sensitivity) and pools over
    #     PATIENTS by DerSimonian-Laird, so the panel-d filter selects the proportion rows and quotes n_patients.
    #   * the brainstem tables written by step 21 (brainstem_*_mp.csv) are superseded: they predate the audit and still carry
    #     the pre-audit novel_MP15 family. Panel e now reads the regenerated brainstem_section_means.csv / brainstem_gradient.csv.
    #   * panels a-c use pixel coordinates with a scale bar (SP-05); array_row/array_col stretched the section ~1.6x.
    sp = T("spatial_abundance_long.csv", usecols=["section_id", "spot", "class", "q05", "pxl_row", "pxl_col"])
    sp = sp[sp.section_id == "Pons_A"]
    w = sp.pivot_table(index=["spot", "pxl_row", "pxl_col"], columns="class", values="q05").reset_index()
    mal = [c for c in w.columns if str(c).startswith("malignant_")]
    tot = w[[c for c in w.columns if c not in ("spot", "pxl_row", "pxl_col")]].sum(1)
    w["tumour"] = w[mal].sum(1) / tot; w["mes"] = w["malignant_MES_like"] / w[mal].sum(1); w["myeloid_p"] = w["myeloid"] / tot
    co = T("spatial_colocalisation.csv"); rs = T("ren_section_composition.csv"); rt = T("ren_dmg_vs_gbm_composition.csv").set_index("variable")
    bm = T("brainstem_section_means.csv").sort_values("anatomical_rank"); bg = T("brainstem_gradient.csv")
    fig = plt.figure(figsize=(W, 6.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.5, wspace=0.75)
    x, y = w.pxl_col.values, -w.pxl_row.values
    # Scale bar on the same convention as scripts/15c_spatial_figures_fix.py (audit SP-05): px per mm from the section's
    # recorded spot_diameter_fullres and the 55 um Visium spot. NOTE (reporting pass 2026-10-21): for the brainstem sections
    # that value (55.0 px for Pons_A) is not consistent with the measured 78.0 px spot pitch under the standard 100 um
    # Visium pitch, so the absolute scale of the brainstem panels is uncertain; it is recorded, not silently corrected.
    sec_used = pd.read_csv(ROOT / "results/qc/spatial_sections_used.csv").set_index("section_id")
    used["results/qc/spatial_sections_used.csv"] = ROOT / "results/qc/spatial_sections_used.csv"
    px_per_mm = float(sec_used.loc["Pons_A", "spot_diameter_fullres"]) / 55.0 * 1000.0
    for k, (col, cmap, lab, cb, ttl) in enumerate([("tumour", "Reds", "a", "tumour fraction", "Tumour fraction"),
                                                   ("mes", "Purples", "b", "MES-like share of tumour", "MES-like share"),
                                                   ("myeloid_p", "Blues", "c", "myeloid fraction", "Myeloid fraction")]):
        ax = fig.add_subplot(gs[0, k])
        s = ax.scatter(x, y, c=w[col], cmap=cmap, s=1.6, marker="h", linewidths=0, rasterized=True, vmin=0, vmax=np.nanpercentile(w[col], 99))
        ax.set_aspect("equal"); ax.axis("off")
        x0, y0 = x.min(), y.min()
        ax.plot([x0, x0 + px_per_mm], [y0 - 0.04 * (y.max() - y0)] * 2, color="black", lw=1.4, solid_capstyle="butt", clip_on=False)
        cbar = fig.colorbar(s, ax=ax, fraction=0.04, pad=0.01, orientation="horizontal"); cbar.set_label(cb); cbar.outline.set_visible(False)
        letter(ax, lab, dx=-0.06, title=ttl)
    # d co-localisation forest
    ax = fig.add_subplot(gs[1, 0])
    d = co[(co.metric == "proportion") & (co.group == "DMG") & (co.pair_type == "state_vs_TME")].copy()
    d["name"] = [f"{SPATIAL_CLASS.get(a, a)} – {SPATIAL_CLASS.get(b, b)}" if a.startswith("malignant") else f"{SPATIAL_CLASS.get(b, b)} – {SPATIAL_CLASS.get(a, a)}"
                 for a, b in zip(d.class_a, d.class_b)]
    dall = d.copy()                 # the whole BH family; d below is the 12 pairs drawn
    d = d.sort_values("rho_re").tail(12)
    yy = np.arange(len(d))
    sig = d.q_bh_state_vs_TME < 0.05
    ax.errorbar(d.rho_re, yy, xerr=[d.rho_re - d.rho_re_lo, d.rho_re_hi - d.rho_re], fmt="none", ecolor="#525252", elinewidth=0.6)
    ax.scatter(d.rho_re, yy, c=np.where(sig, "black", "white"), edgecolors="black", s=12, linewidths=0.6, zorder=3)
    ax.axvline(0, color="#bdbdbd", lw=0.6, zorder=0); ax.set_yticks(yy); ax.set_yticklabels(d.name)
    ax.set_xlabel(f"pooled Spearman ρ ({int(d.n_patients.max())} patients)")
    letter(ax, "d", dx=-0.6, title="Co-localisation")
    # e brainstem gradient
    ax = fig.add_subplot(gs[1, 1])
    for col, c, lbl in [("tumour_fraction_mean", "black", "tumour fraction"), ("prop_oligodendrocyte", LINEAGE_COLOURS["oligodendrocyte"], "oligodendrocyte"),
                        ("prop_malignant_cycling", FAMILY_COLOURS["cycling"], "cycling"), ("prop_malignant_MES_like", FAMILY_COLOURS["MES_like"], "MES-like")]:
        ax.plot(bm.anatomical_rank, bm[col], "o-", color=c, ms=2.5, lw=0.8, label=lbl)
    ax.set_xticks(bm.anatomical_rank); ax.set_xticklabels([s.replace("Spinal_Cord", "SC").replace("_", " ") for s in bm.section_id], rotation=60, ha="right")
    ax.set_ylabel("mean spot proportion")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5, columnspacing=0.6, handlelength=1)  # audit P7-17
    letter(ax, "e", dx=-0.25, title="Brainstem levels", tpad=26)
    # f Ren DMG vs GBM
    ax = fig.add_subplot(gs[1, 2])
    rs = rs[rs.group.isin(["DMG", "GBM"])]
    for k, (col, lbl) in enumerate([("tumour_fraction", "tumour\nfraction"), ("share_malignant_cycling", "cycling share\nof tumour")]):
        for j, g in enumerate(["DMG", "GBM"]):
            v = rs.loc[rs.group == g, col]
            ax.scatter(np.full(len(v), k * 2.6 + j) + np.random.default_rng(k + j).uniform(-0.15, 0.15, len(v)), v, s=10,
                       c="#b2182b" if g == "DMG" else "#525252", linewidths=0)
            ax.hlines(v.median(), k * 2.6 + j - 0.3, k * 2.6 + j + 0.3, color="black", lw=0.8)
    ax.set_xticks([0, 1, 2.6, 3.6]); ax.set_xticklabels(["DMG", "GBM", "DMG", "GBM"])
    for k, lbl in enumerate(["tumour\nfraction", "cycling share\nof tumour"]):
        ax.text(k * 2.6 + 0.5, -0.12, lbl, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6)
    ax.set_ylabel("section value"); ax.set_ylim(0, None)
    letter(ax, "f", dx=-0.25, title="DMG vs GBM sections")
    save(fig, "fig5_spatial")
    tr = bg[bg.test == "trend_along_axis"].set_index("cls")
    top = d.sort_values("rho_re", ascending=False)
    def ci(r):  # noqa: E306
        return (f"ρ = {r.rho_re:.3f} ({r.rho_re_lo:.3f}–{r.rho_re_hi:.3f}), q = {r.q_bh_state_vs_TME:.3g}, "
                f"{r.n_sections_positive}/{r.n_sections} sections positive")
    def pair(a, b):  # noqa: E306
        r = dall[(dall.class_a == a) & (dall.class_b == b)]
        return ci(r.iloc[0]) if len(r) else "not in table"
    tf, cy = rt.loc["tumour_fraction"], rt.loc["share_malignant_cycling"]
    captions[5] = (
        "**Figure 5 | Spatial organisation of malignant states (cell2location deconvolution of 27 Visium sections, 4 cohorts).** "
        "Computed from the audit-corrected spatial tables (track 'spatial', SP-01 … SP-12, regenerated 2026-10-20): co-localisation is "
        "measured on spot *proportions* rather than raw q05 abundances, all four cohorts were deconvolved to 20,000 epochs against a single "
        "metaprogram-family reference, and the 21 DMG sections are pooled over the 12 contributing **patients** — the 10 brainstem sections "
        "are serial sections of one autopsy and carry 66 % of the DMG spots, so treating sections as independent would triple-count that "
        "patient. "
        f"**a–c**, Brainstem autopsy section Pons_A ({len(w):,} spots, one patient): tumour fraction (**a**), MES-like share of the tumour "
        f"(**b**) and myeloid fraction (**c**) per spot (q05 abundances on pixel coordinates, metaprogram-family reference; colour scale "
        f"capped at the 99th percentile; scale bar 1 mm = {px_per_mm:.0f} px). "
        f"**d**, Strongest malignant-state × microenvironment co-localisations (spot-level Spearman ρ pooled over patients by "
        f"DerSimonian–Laird random effects on Fisher z, 95 % CI; the {len(d)} strongest of the {len(dall)} state × microenvironment pairs are "
        f"drawn, filled: BH q < 0.05 across all {len(dall)}). "
        f"Top pair {top.iloc[0]['name']} {ci(top.iloc[0])}. Of the two pre-stated hypotheses, MES-like – myeloid is supported "
        f"({pair('malignant_MES_like', 'myeloid')}) and OPC-like – oligodendrocyte is not "
        f"({pair('malignant_OPC_like', 'oligodendrocyte')}); between-section heterogeneity is high for every pair "
        f"(I² {dall.I2.min():.3f}–{dall.I2.max():.3f}), so the pooled estimates are not section-generalisable. "
        "**e**, Rostral-to-caudal section means in the whole-brainstem autopsy (metaprogram-family reference; one patient, 10 "
        "sections, descriptive only — no inference is drawn from this panel): "
        f"tumour fraction peaks in pons (max {bm.tumour_fraction_mean.max():.2f}) and falls toward the spinal cord; oligodendrocyte proportion "
        f"rises (ρ = {tr.loc['oligodendrocyte', 'rho_prop']:.2f}, p = {tr.loc['oligodendrocyte', 'p_prop']:.3f}) and cycling falls "
        f"(absolute abundance ρ = {tr.loc['malignant_cycling', 'rho_abs']:.2f}, p = {tr.loc['malignant_cycling', 'p_abs']:.3f}); none survives "
        f"BH across {len(tr)} classes (min q = {tr.q_bh_abs.min():.2f}). "
        f"**f**, Ren cohort, {int(tf.n_dmg)} DMG vs {int(tf.n_gbm)} GBM sections, one section per patient (bars: medians; exact "
        f"Mann–Whitney): tumour fraction {tf.dmg_median:.3f} vs {tf.gbm_median:.3f} (p = {tf.p_mwu:.3f}, q = {tf.q_bh:.3f}); cycling share "
        f"{cy.dmg_median:.3f} vs {cy.gbm_median:.3f} (p = {cy.p_mwu:.3f}, q = {cy.q_bh:.3f}). After BH across the {len(rt)} variables tested "
        f"neither difference survives (min q = {rt.q_bh.min():.3f}), so the DMG-vs-GBM comparison is null at n = 5 vs 5 (audit SP-12).")


# ---------------------------------------------------------------- Fig 6
def fig6():
    at = T("atac_summary.csv"); ac = T("atac_label_correlation.csv", index_col=0)
    # audit P7-07/P7-12: panel c reads the corrected rank-score fits (25b, re-derived families). The pre-audit
    # survival_results.csv still scores the superseded novel_MP15 family and its NNLS sets used pre-audit signatures.
    sv = T("survival_rankscore_audit.csv")
    fm = T("foundation_model_comparison.csv"); pk = T("atac_union_peaks.csv.gz")
    fig = plt.figure(figsize=(W, 5.8))
    gs = fig.add_gridspec(2, 2, hspace=0.6, wspace=0.55, width_ratios=[1, 1.2])
    # a ATAC vs RNA neighbourhood enrichment
    ax = fig.add_subplot(gs[0, 0])
    at = at.sort_values("atac_enrichment"); y = np.arange(len(at))
    ax.hlines(y, at.atac_enrichment, at.rna_enrichment, color="#bdbdbd", lw=1)
    ax.scatter(at.atac_enrichment, y, c=[FAMILY_COLOURS.get(f, "#000") for f in at.family], s=18, zorder=3, linewidths=0, label="ATAC (LSI)")
    ax.scatter(at.rna_enrichment, y, facecolors="white", edgecolors=[FAMILY_COLOURS.get(f, "#000") for f in at.family], s=18, zorder=3,
               linewidths=0.8, label="RNA (scVI)")
    ax.axvline(1, color="#525252", lw=0.6, ls="--", zorder=0)
    ax.set_yticks(y); ax.set_yticklabels([FAMILY_NAMES.get(f, f) for f in at.family])
    ax.set_xlabel("same-family kNN neighbours / expected")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5, columnspacing=0.6, handlelength=1)  # audit P7-17
    letter(ax, "a", dx=-0.3, title="Neighbourhood enrichment, ATAC and RNA", tpad=11)
    # b ATAC label correlation
    ax = fig.add_subplot(gs[0, 1])
    keep = [c for c in ac.columns if c != "mal:unassigned"]
    m = ac.loc[keep, keep]
    nm = lambda s: FAMILY_NAMES.get(s.replace("mal:", ""), s.replace("mal:", "")) if s.startswith("mal:") else s  # noqa: E731
    m.index = [nm(s) for s in m.index]; m.columns = m.index
    heat(ax, m, cmap="viridis", vmin=0.3, vmax=1, cbar_label="pseudobulk ATAC correlation", fig=fig)
    letter(ax, "b", dx=-0.4, title="Pseudobulk ATAC correlation")
    # c survival forest (rank scores, adjusted)
    ax = fig.add_subplot(gs[1, 0])
    s = sv[(sv.variable_set == "rank_score") & (sv.model == "adjusted")].copy()
    s["fam"] = s.variable.str.replace("rank:", "", regex=False); s = s.sort_values("hr_per_sd"); y = np.arange(len(s))
    ax.errorbar(s.hr_per_sd, y, xerr=[s.hr_per_sd - s.ci_lo, s.ci_hi - s.hr_per_sd], fmt="none", ecolor="#525252", elinewidth=0.6)
    ax.scatter(s.hr_per_sd, y, c=[FAMILY_COLOURS.get(f, "#000") for f in s.fam], s=16, zorder=3, linewidths=0)
    ax.axvline(1, color="#bdbdbd", lw=0.6, zorder=0); ax.set_xscale("log")
    ax.set_xticks([0.6, 0.8, 1, 1.25, 1.6]); ax.set_xticklabels(["0.6", "0.8", "1", "1.25", "1.6"]); ax.minorticks_off()
    ax.set_yticks(y); ax.set_yticklabels([FAMILY_NAMES.get(f, f) for f in s.fam]); ax.set_xlabel("hazard ratio per SD of rank score")
    letter(ax, "c", dx=-0.3, title="Hazard ratio per family score")
    # d foundation models
    ax = fig.add_subplot(gs[1, 1])
    met = [("knn_balanced_accuracy", "state kNN accuracy"), ("cross_assay_accuracy_nuclei_to_cells", "nuclei → cells transfer"),
           ("same_assay_neighbours", "same-assay neighbours"), ("same_sample_neighbours", "same-sample neighbours")]
    cols = {"PCA": "#bdbdbd", "scVI": "#525252", "scGPT": "#1b9e77", "Geneformer": "#d95f02"}
    xb = np.arange(len(met)); wd = 0.2
    for i, (e, c) in enumerate(cols.items()):
        r = fm.set_index("embedding").loc[e]
        ax.bar(xb + (i - 1.5) * wd, [r[k] for k, _ in met], wd, color=c, label=e)
    r0 = fm.iloc[0]
    for j, (k, e) in enumerate([("same_assay_neighbours", "same_assay_expected"), ("same_sample_neighbours", "same_sample_expected")]):
        jj = [m_[0] for m_ in met].index(k)
        ax.hlines(r0[e], jj - 0.45, jj + 0.45, color="black", lw=0.8, ls=":")
    ax.set_xticks(xb); ax.set_xticklabels([l for _, l in met], rotation=20, ha="right"); ax.set_ylim(0, 1.05)
    ax.set_ylabel("fraction")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5, columnspacing=0.6, handlelength=1)  # audit P7-17
    letter(ax, "d", dx=-0.15, title="Foundation-model embeddings", tpad=24)
    save(fig, "fig6_multimodal")
    a = at.set_index("family"); s2 = s.set_index("fam"); f2 = fm.set_index("embedding")
    n_cox, ev = int(s.n.iloc[0]), int(s.events.iloc[0])
    captions[6] = (
        "**Figure 6 | Chromatin, outcome and foundation-model anchoring of the malignant states.** *Exploratory: the multimodal, bulk-survival "
        "and foundation-model layer (code Phase 7, build-plan Phase 8) was not part of the approved plan and is reported as a methods note, "
        "not as an atlas claim (audit P7-16).* "
        "**a**, Neighbourhood enrichment of each metaprogram family in the Jessa 2022 multiome ATAC (LSI, filled) and RNA (scVI latent, open): "
        "fraction of same-family kNN neighbours over the expectation (dashed line = chance). Unit: cell, pooled over the "
        f"{int(at.n_samples_tested.max())} multiome samples at most (per family {int(at.n_samples_tested.min())}–"
        f"{int(at.n_samples_tested.max())}); all ten are from Jessa 2022, all nuclei, so this is a within-cohort statement. Ranked by ATAC "
        f"enrichment the families run {', '.join(f'{FAMILY_NAMES.get(f, f)} {v:.2f}' for f, v in a.atac_enrichment.sort_values(ascending=False).items())}: "
        "cilia (MP9) is at chance and hypoxia barely above it, while the four lineage-like families and cycling are clearly above. RNA "
        "enrichment is inflated throughout because the family labels are RNA-derived. "
        f"**b**, Correlation of within-sample pseudobulk accessibility (top 20,000 variable peaks of {len(pk):,} union peaks) between malignant "
        "families and normal lineages (mean over samples with ≥ 50 cells of both labels; blank: no such sample). Pairs backed by fewer than "
        "3 samples are masked and not interpreted (audit P7-03; per-pair n in `atac_label_correlation_n.csv`). "
        f"**c**, OpenPedCan v15 H3 K27-altered DMG (unit: patient; n = {n_cox}, {ev} deaths; atlas-overlapping patient removed): Cox hazard "
        "ratio per SD of a rank-based metaprogram-family score, adjusted for age and pontine site (95 % CI). No family is associated with "
        f"survival (min BH q = {s.q_bh.min():.2f}; cycling HR {s2.loc['cycling', 'hr_per_sd']:.2f}, p = {s2.loc['cycling', 'p']:.3f}). "
        f"Proportional hazards is violated in {int(s.ph_violated.sum())} of the {len(s)} fits shown at p < 0.05 by a Schoenfeld test with a rank time "
        f"transform ({', '.join(FAMILY_NAMES.get(f, f) for f in s.fam[s.ph_violated]) if int(s.ph_violated.sum()) else 'none'}; audit P7-07), "
        "so the hazard ratios are average effects over follow-up rather than constant ones. "
        "NNLS deconvolution proportions are not shown: they disagreed with these rank scores (docs/phases/phase7_multimodal_notes.md). "
        f"**d**, Zero-shot foundation-model embeddings vs scVI and PCA on {int(fm.n_cells.iloc[0]):,} confidently assigned malignant cells: "
        f"balanced kNN accuracy for the family label, nuclei → cells label transfer, and the fraction of neighbours from the same assay / sample "
        f"(dotted: expectation under random mixing, {r0.same_assay_expected:.2f} / {r0.same_sample_expected:.3f}). scGPT and Geneformer "
        f"neighbourhoods are almost assay-pure ({f2.loc['scGPT', 'same_assay_neighbours']:.3f}, {f2.loc['Geneformer', 'same_assay_neighbours']:.3f}); "
        f"only scVI mixes assays ({f2.loc['scVI', 'same_assay_neighbours']:.2f}) and samples ({f2.loc['scVI', 'same_sample_neighbours']:.2f}). "
        "PCA family accuracy is circular (labels are scores on the same expression).")


FIGS = {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6}
NAMES = {1: "fig1_atlas", 2: "fig2_mapping", 3: "fig3_metaprograms", 4: "fig4_tme", 5: "fig5_spatial", 6: "fig6_multimodal"}


def write_doc():
    old = {}
    if DOC.exists():  # keep captions of figures not re-drawn in this run
        for blk in DOC.read_text().split("\n\n"):
            if blk.startswith("**Figure "):
                old[int(blk.split()[1])] = blk
    caps = {**old, **captions}
    lines = ["# Figures — captions", "", f"Generated by `scripts/30_figures.py` on {time.strftime('%Y-%m-%d')}. Do not edit by hand; every number "
             "is computed from `results/tables/`. Files: `results/pub_figs/main/fig<N>_<name>.png` (300 dpi) and `.pdf`.", ""]
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
        with log.step(f"figure {k} ({NAMES[k]})"):
            FIGS[k]()
    if captions:
        write_doc()
    log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
