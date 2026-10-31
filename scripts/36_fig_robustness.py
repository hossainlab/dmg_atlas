#!/usr/bin/env python
"""Phase 9 / step 36 — publication figure 7 (preprint Figure 6): robustness of the reference (Phase 9 quality checks), drawn from results/tables/quality/.

  A  leave-one-publication-out classifier: accuracy and balanced accuracy per held-out publication, with the
     held-out-patient baseline of models/scanvi as reference lines (loso_classifier_summary.csv)
  B  per-class recall per held-out publication (loso_classifier_per_class.csv; classes with < 20 true cells greyed)
  C  fate of labelled cells of Okada when Okada is held out: predicted-label fractions per true class (loso_classifier_confusion.csv)
  D  recovery of the 15 primary metaprograms: drop-one-publication (n/5), publication-balanced and size-matched random pools
     (metaprogram_lopo_summary.csv)
  E  cell-type composition coefficients under scANVI vs consensus-only labels (composition_label_sensitivity_compare.csv)
Writes results/pub_figs/main/fig7_robustness.{pdf,png}. Run: python scripts/36_fig_robustness.py
"""
import sys, pathlib
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
R = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R / "src"))
from dmg_atlas.figstyle import apply_figure_style, letter, W, LINEAGE_COLOURS, LINEAGES, CONTRAST_NAMES  # noqa: E402

Q = R / "results/tables/quality"; OUT = R / "results/pub_figs/main"; OUT.mkdir(parents=True, exist_ok=True)
FOLDS = ["okada", "jessa2022", "collot", "zang", "desisto"]
FOLD_NAMES = {"okada": "Okada", "jessa2022": "Jessa 2022", "collot": "Collot", "zang": "Zang", "desisto": "DeSisto"}
CLS_NAMES = {"T_NK": "T/NK", "oligodendrocyte": "oligodendro."}
apply_figure_style()
summ = pd.read_csv(Q / "loso_classifier_summary.csv").set_index("fold")
pc = pd.read_csv(Q / "loso_classifier_per_class.csv")
conf = pd.read_csv(Q / "loso_classifier_confusion.csv")
mp = pd.read_csv(Q / "metaprogram_lopo_summary.csv")
cmp_ = pd.read_csv(Q / "composition_label_sensitivity_compare.csv")

fig = plt.figure(figsize=(W, 6.6))
outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.1], hspace=0.75)
top = outer[0].subgridspec(1, 3, width_ratios=[1, 1.15, 1.2], wspace=1.0)
bot = outer[1].subgridspec(1, 2, width_ratios=[2.1, 1], wspace=0.45)

# A: accuracy vs balanced accuracy per held-out publication
ax = fig.add_subplot(top[0]); y = np.arange(len(FOLDS))
base = summ.loc["baseline_heldout_patients"]
ax.axvline(base.accuracy, color="#999999", lw=0.8, ls="--", label="held-out patients: accuracy")
ax.axvline(base.balanced_accuracy, color="#999999", lw=0.8, ls=":", label="held-out patients: balanced")
for i, f in enumerate(FOLDS):
    a, b = summ.loc[f, "accuracy"], summ.loc[f, "balanced_accuracy"]
    ax.plot([b, a], [i, i], color="#cccccc", lw=1.2, zorder=1)
    ax.scatter(a, i, s=18, color="#252525", zorder=2, label="accuracy" if i == 0 else None)
    ax.scatter(b, i, s=18, facecolor="white", edgecolor="#b2182b", lw=1.0, zorder=2, label="balanced accuracy" if i == 0 else None)
ax.set_yticks(y, [FOLD_NAMES[f] for f in FOLDS]); ax.invert_yaxis(); ax.set_xlim(0.3, 1.02); ax.set_xlabel("held-out publication score")
ax.legend(loc="lower left", bbox_to_anchor=(-0.45, 1.0), ncol=1, handletextpad=0.3, borderaxespad=0.2, fontsize=5.5, handlelength=1.6)
letter(ax, "a", dx=-0.42, dy=1.40, title="Classifier on an unseen study", tpad=40)

# B: per-class recall heatmap
ax = fig.add_subplot(top[1])
m = pc.pivot(index="cls", columns="fold", values="recall").reindex(index=LINEAGES, columns=FOLDS)
small = pc.pivot(index="cls", columns="fold", values="averaged").reindex(index=LINEAGES, columns=FOLDS)
im = ax.imshow(m.values, cmap="Greys", vmin=0, vmax=1, aspect="auto")
for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        v = m.values[i, j]
        if np.isnan(v):
            ax.text(j, i, "–", ha="center", va="center", fontsize=5.5, color="#999999"); continue
        faded = small.values[i, j] is False or small.values[i, j] == False  # noqa: E712
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5,
                color=("#aaaaaa" if faded else ("white" if v > 0.55 else "black")), style="italic" if faded else "normal")
ax.set_xticks(range(len(FOLDS)), [FOLD_NAMES[f] for f in FOLDS], rotation=40, ha="right")
ax.set_yticks(range(len(LINEAGES)), [CLS_NAMES.get(c, c) for c in LINEAGES])
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04, shrink=0.8); cb.set_label("recall"); cb.outline.set_linewidth(0.4)
letter(ax, "b", dx=-0.45, dy=1.40, title="Per-class recall by held-out study", tpad=40)

# C: Okada held out — predicted fate of each true class
ax = fig.add_subplot(top[2])
c = conf[conf.fold == "okada"].pivot(index="true", columns="pred", values="n").fillna(0).reindex(index=LINEAGES).dropna(how="all")
frac = c.div(c.sum(1), axis=0)
left = np.zeros(len(frac))
for p in LINEAGES:
    if p in frac:
        ax.barh(range(len(frac)), frac[p].values, left=left, color=LINEAGE_COLOURS[p], height=0.72, lw=0)
        left += frac[p].values
ax.set_yticks(range(len(frac)), [CLS_NAMES.get(t, t) for t in frac.index]); ax.invert_yaxis()
for i, t_ in enumerate(frac.index):
    ax.text(1.02, i, f"{int(c.loc[t_].sum()):,}", va="center", fontsize=5, color="#555555", transform=ax.get_yaxis_transform())
ax.set_xlim(0, 1); ax.set_xlabel("fraction of labelled cells")
handles = [mpl.patches.Patch(color=LINEAGE_COLOURS[p], label=CLS_NAMES.get(p, p)) for p in LINEAGES if p in frac]
ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.05, 1.0), ncol=3, fontsize=5.5, handlelength=0.8, columnspacing=0.6, title="predicted as", title_fontsize=5.5)
letter(ax, "c", dx=-0.42, dy=1.40, title="Okada held out: true class → prediction", tpad=40)

# D: metaprogram recovery
ax = fig.add_subplot(bot[0])
mp["drop"] = mp.recovered_drop.str.split("/").apply(lambda s: int(s[0]) / int(s[1]))
FAM = {"novel_MP9": "ciliated", "IEG_dissociation_stress": "dissociation stress"}
lab = [f"{r.metaprogram} {FAM.get(n, n).replace('_', '-')}" + ("*" if r.replicated_primary else "")
       for r in mp.itertuples() for n in [r.family if isinstance(r.family, str) else r.technical_class.split(':')[1]]]
x = np.arange(len(mp)); wbar = 0.27
for k, (col, colr, name) in enumerate([("drop", "#252525", "drop one publication (of 5)"),
                                       ("recovered_balanced_frac", "#4575b4", "≤ 5 samples per publication (20 pools)"),
                                       ("recovered_random_frac", "#bdbdbd", "size-matched random pools (20)")]):
    ax.bar(x + (k - 1) * wbar, mp[col], width=wbar, color=colr, label=name, lw=0)
ax.set_xticks(x, lab, rotation=55, ha="right", fontsize=5.5); ax.set_ylim(0, 1.05); ax.set_ylabel("fraction of scenarios recovered")
ax.legend(loc="upper right", fontsize=5.5, ncol=1)
letter(ax, "d", dx=-0.07, dy=1.04, title="Metaprogram recovery without each publication (* replicated)")

# E: composition coefficients, scANVI vs consensus-only labels
ax = fig.add_subplot(bot[1])
d = cmp_.dropna(subset=["coef_scanvi", "coef_consensus_only"])
sig = d.sig_scanvi | d.sig_consensus_only
lim = np.ceil(max(d.coef_scanvi.abs().max(), d.coef_consensus_only.abs().max()))
ax.plot([-lim, lim], [-lim, lim], color="#999999", lw=0.6, ls="--"); ax.axhline(0, color="#dddddd", lw=0.5); ax.axvline(0, color="#dddddd", lw=0.5)
ax.scatter(d.coef_scanvi[~sig], d.coef_consensus_only[~sig], s=9, color="#bdbdbd", lw=0, label="q ≥ 0.05 under both")
MARK = {"h3_H3.1_vs_H3.3": "o", "assay_nuclei_vs_cells": "s"}
for (cat, con), r in d[sig].set_index(["category", "contrast"]).iterrows():
    ax.scatter(r.coef_scanvi, r.coef_consensus_only, s=20, marker=MARK.get(con, "D"), color=LINEAGE_COLOURS.get(cat, "#252525"),
               edgecolor="black" if r.sig_scanvi else "none", lw=0.6, zorder=3)
cats = sorted(set(d[sig].category))
h = [mpl.lines.Line2D([], [], ls="", marker="o", color=LINEAGE_COLOURS.get(c_, "k"), ms=4, label=CLS_NAMES.get(c_, c_)) for c_ in cats]
h += [mpl.lines.Line2D([], [], ls="", marker=MARK[k], mfc="white", mec="black", ms=4, label=CONTRAST_NAMES.get(k, k)) for k in MARK]
h += [mpl.lines.Line2D([], [], ls="", marker="o", mfc="white", mec="black", mew=1.0, ms=4, label="black edge: q < 0.05 (scANVI)")]
ax.legend(handles=h, loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=5, handletextpad=0.2, borderaxespad=0)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_xlabel("Δ logit, scANVI labels"); ax.set_ylabel("Δ logit, consensus-only labels")
letter(ax, "e", dx=-0.32, dy=1.04, title="Composition vs ambiguous-cell labels")

for ext in ("pdf", "png"):
    fig.savefig(OUT / f"fig7_robustness.{ext}")
print("wrote", OUT / "fig7_robustness.pdf")
