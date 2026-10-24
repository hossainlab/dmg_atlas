"""Graphical abstract for the DMG atlas preprint (run from the repo root).

Three data-driven columns: the integrated reference (real UMAP), what the atlas calibrates
(per-sample malignant-call agreement), and what it finds (metaprogram replication, spatial
co-localisation). House style: the graphical abstract carries no panel letters and no panel titles - it reads as a
flow diagram, so the Cell Reports panel-title convention applied to Figures 1-6 and S1-S8 does not
apply here; identity labels only, everything else in the caption. The UMAP uses the stripped
embedding style of figstyle.umap_ax (no axis box, ticks or axis labels; corner arrow pair kept).
Writes results/pub_figs/graphical_abstract.{png,pdf}.
"""
import sys, pathlib
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

R = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R / "src"))
from dmg_atlas import figstyle as FS

T = lambda n: pd.read_csv(R / "results" / "tables" / n)
OUT = R / "results" / "pub_figs"; OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(0)

# ---------------------------------------------------------------- data
um = pd.read_csv(R / "results/tables/atlas_core_umap.csv", index_col=0)
lab = pd.read_csv(R / "results/tables/atlas_core_obs_labels.csv", index_col=0,
                  usecols=["Unnamed: 0", "study", "assay", "scanvi_label", "consensus_call"])
d = um.join(lab, how="inner")
keep = rng.permutation(len(d))[:60000]
d = d.iloc[np.sort(keep)]

val = T("malignant_call_validation.csv")
samp = val[(val.level == "sample") & (val.tier == "tier2_cnv_call") & (val.subset == "all_samples")].dropna(subset=["agreement"])

cut = T("metaprogram_cut_sensitivity.csv")

sp = T("spatial_colocalisation.csv")
sp = sp[(sp.metric == "proportion") & (sp.subset == "all") & (sp.pair_type == "state_vs_TME")].copy()
sp = sp.sort_values("q_bh_state_vs_TME").head(6).iloc[::-1]
nm = lambda c: (c.replace("malignant_", "").replace("novel_MP9", "cilia").replace("_like", "-like").replace("T_NK", "T/NK"))

# ---------------------------------------------------------------- canvas
FS.apply_figure_style(sizes=(8, 7, 6))
fig = plt.figure(figsize=(FS.W, 3.5))
gs = fig.add_gridspec(3, 3, width_ratios=[0.85, 1.25, 1.05], height_ratios=[1, 1, 1],
                      left=0.015, right=0.985, top=0.95, bottom=0.11, wspace=0.42, hspace=0.75)

# ---- column 1: inputs ------------------------------------------------
axin = fig.add_subplot(gs[:, 0]); axin.axis("off")
ref = (d.groupby("study").size() if False else None)
cohorts = [("Okada 2026", 156936, "nuclei"), ("Jessa 2022", 83076, "both"), ("Jessa multiome", 46288, "nuclei"),
           ("Collot 2026", 17562, "nuclei"), ("Zang 2024", 6546, "cells"), ("DeSisto 2024", 1210, "cells")]
col = {"nuclei": FS.ASSAY_COLOURS["snRNA"], "cells": FS.ASSAY_COLOURS["scRNA"], "both": FS.ASSAY_COLOURS["multiome_RNA"]}
y0, h = 0.93, 0.072
mx = max(c[1] for c in cohorts)
for i, (nmc, n, a) in enumerate(cohorts):
    y = y0 - i * h
    axin.add_patch(plt.Rectangle((0.42, y - 0.022), 0.5 * n / mx, 0.036, color=col[a], transform=axin.transAxes, clip_on=False))
    axin.text(0.40, y - 0.004, nmc, ha="right", va="center", fontsize=6, transform=axin.transAxes)
axin.text(0.0, y0 - 6 * h - 0.03, "6 cohorts  ·  57 samples  ·  52 patients\n311,618 cells and nuclei  ·  18,882 shared genes",
          fontsize=6.5, va="top", transform=axin.transAxes)
for i, (t, c) in enumerate([("nuclei", col["nuclei"]), ("cells", col["cells"]), ("both", col["both"])]):
    axin.add_patch(plt.Rectangle((0.00 + i * 0.34, 0.30), 0.045, 0.03, color=c, transform=axin.transAxes, clip_on=False))
    axin.text(0.055 + i * 0.34, 0.315, t, fontsize=6, va="center", transform=axin.transAxes)
axin.text(0.0, 0.22, "5 further cohorts projected\nwithout fine-tuning", fontsize=6.5, va="top", transform=axin.transAxes)
axin.text(0.0, 0.08, "88.4 % nuclei", fontsize=7, va="top", fontweight="bold",
          color=FS.ASSAY_COLOURS["snRNA"], transform=axin.transAxes)

# ---- column 2: the reference UMAP ------------------------------------
axu = fig.add_subplot(gs[:, 1])
order = [l for l in FS.LINEAGES if l in set(d.scanvi_label)]
for l in order:
    s = d[d.scanvi_label == l]
    axu.scatter(s.UMAP1, s.UMAP2, s=0.45, lw=0, color=FS.LINEAGE_COLOURS[l], rasterized=True,
                alpha=0.75 if l == "malignant" else 0.95)
for l in order:
    s = d[d.scanvi_label == l]
    if len(s) < 200: continue
    axu.text(np.median(s.UMAP1), np.median(s.UMAP2), l.replace("T_NK", "T/NK"), fontsize=6.5, ha="center", va="center",
             color="#111111", bbox=dict(fc="white", ec="none", alpha=0.72, pad=0.8))
# stripped embedding axes + corner arrow pair; aspect=False keeps the panel's existing proportions
FS.umap_ax(axu, arrows=True, aspect=False)
axu.text(0.5, -0.045, "scANVI: 0.975 accuracy on 5 held-out patients", fontsize=6.5, ha="center",
         transform=axu.transAxes)

# ---- column 3a: malignant-call agreement -----------------------------
ax1 = fig.add_subplot(gs[0, 2])
x = rng.normal(0, 0.055, len(samp))
ax1.scatter(x, samp.agreement, s=4, lw=0, color="#555555", alpha=0.6, rasterized=True)
med = samp.agreement.median()
ax1.plot([-0.16, 0.16], [med, med], color="#b2182b", lw=1.6, solid_capstyle="butt")
ax1.text(0.21, med, f"median {med:.3f}", fontsize=6, va="center", color="#b2182b")
ax1.set_xlim(-0.33, 0.75); ax1.set_xticks([]); ax1.set_ylim(0.3, 1.02)
ax1.set_ylabel("agreement with\npublished labels", fontsize=6.5, labelpad=1)
ax1.tick_params(labelsize=6)
ax1.text(0.0, 1.05, f"{len(samp)} samples", fontsize=6, transform=ax1.transAxes)

# ---- column 3b: metaprogram replication ------------------------------
ax2 = fig.add_subplot(gs[1, 2])
ax2.plot(cut.cut, cut.n_metaprograms, "-", lw=1.2, color="#555555")
ax2.plot(cut.cut, cut.n_replicated, "-", lw=1.6, color=FS.FAMILY_COLOURS["MES_like"])
ax2.axvline(0.85, color="#bdbdbd", lw=0.6, ls="--", zorder=0)
ax2.text(0.922, cut.n_metaprograms.iloc[-1], "metaprograms", fontsize=6, ha="right", va="bottom", color="#555555")
ax2.text(0.922, cut.n_replicated.iloc[-1] + 0.6, "replicated", fontsize=6, ha="right", va="bottom",
         color=FS.FAMILY_COLOURS["MES_like"])
ax2.set_xlabel("dendrogram cut", fontsize=6.5, labelpad=1); ax2.set_ylabel("count", fontsize=6.5, labelpad=1)
ax2.tick_params(labelsize=6); ax2.set_ylim(0, 24)
ax2.text(0.846, 21.0, "15 vs 2", fontsize=6.5, ha="right", color="#111111")

# ---- column 3c: spatial co-localisation ------------------------------
ax3 = fig.add_subplot(gs[2, 2])
yy = np.arange(len(sp))
for i, (_, r) in enumerate(sp.iterrows()):
    sig = r.q_bh_state_vs_TME < 0.05
    ax3.plot([r.rho_re_lo, r.rho_re_hi], [i, i], lw=0.9, color="#777777", zorder=1)
    ax3.scatter([r.rho_re], [i], s=16, zorder=2, color="#111111" if sig else "white",
                edgecolor="#111111", lw=0.8)
ax3.axvline(0, color="#bdbdbd", lw=0.6, zorder=0)
ax3.set_yticks(yy); ax3.set_yticklabels([f"{nm(r.class_a)} \u00d7 {nm(r.class_b)}" for _, r in sp.iterrows()], fontsize=6)
ax3.set_xlabel("spot-level \u03c1 (12 patients)", fontsize=6.5, labelpad=1)
ax3.tick_params(labelsize=6); ax3.set_ylim(-0.6, len(sp) - 0.4)
ax3.text(1.0, 1.06, "filled: q < 0.05", fontsize=6, ha="right", transform=ax3.transAxes)

# ---- flow arrows between columns -------------------------------------
for x0, x1 in [(0.268, 0.300), (0.606, 0.638)]:
    fig.add_artist(FancyArrowPatch((x0, 0.52), (x1, 0.52), transform=fig.transFigure,
                                   arrowstyle="-|>", mutation_scale=9, lw=0.9, color="#888888"))

fig.savefig(OUT / "graphical_abstract.png", dpi=400)
fig.savefig(OUT / "graphical_abstract.pdf")
print("wrote", OUT / "graphical_abstract.png", "| umap cells drawn:", len(d), "| samples:", len(samp),
      "| spatial pairs:", len(sp))
