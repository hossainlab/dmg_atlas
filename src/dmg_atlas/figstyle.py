"""Matplotlib rcParams for deliverable figures (reproduces the mechanics of the `figure-style` skill's
apply_figure_style(): role-mapped font-size ladder, open frame, outward ticks, frameless legends, 300 dpi, Type-42 fonts).
House style (2026-10-21, Cell Reports 'Resource' submission; venue model Cell Reports 43:114804): bold UPPERCASE panel letters
outside the axes, a short descriptive title (2-5 words, sentence case, no statistics and no n) centred above every panel, and
embedding panels with no axis box, no ticks and no axis labels. Statistics, n and the rest of the description stay in the caption."""
from __future__ import annotations


def apply_figure_style(*, frame: str = "open", sizes=(8, 7, 6), grid: bool = False) -> None:
    import matplotlib as mpl
    base, secondary, tick = sizes
    boxed = frame == "boxed"
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.size": base, "axes.labelsize": base, "axes.titlesize": base,
        "legend.fontsize": secondary, "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.linewidth": 0.6, "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid), "legend.frameon": False, "figure.dpi": 200, "savefig.dpi": 300,
        "savefig.bbox": "tight", "axes.titleweight": "normal", "axes.titlelocation": "left",
        "axes.labelweight": "normal", "lines.linewidth": 1.2, "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


LINEAGE_COLOURS = {
    "malignant": "#b2182b", "myeloid": "#2166ac", "oligodendrocyte": "#1b7837", "astrocyte": "#e08214",
    "neuron": "#762a83", "vascular": "#35978f", "T_NK": "#c51b7d", "Unknown": "#bdbdbd",
}
QUERY_COLOURS = {"filbin": "#1b9e77", "labelle": "#d95f02", "liu": "#7570b3", "braun": "#e7298a", "siletti": "#66a61e"}
QUERY_NAMES = {"filbin": "Filbin 2018 (Smart-seq2)", "labelle": "LaBelle 2025 (SS2 + 10x)", "liu": "Liu 2022 (Smart-seq2)",
               "braun": "Braun 2023 fetal hindbrain", "siletti": "Siletti 2023 adult pons"}

# Shared palettes, names and panel helpers for the publication figures (scripts/30_figures.py, scripts/31_supp_figures.py).
W = 7.2  # full page width, inches (183 mm)
LINEAGES = ["malignant", "myeloid", "oligodendrocyte", "astrocyte", "neuron", "vascular", "T_NK"]
# Audit 2026-10-20: the metaprogram track re-derived the families; the cilia programme is now novel_MP9 (was novel_MP15).
# FAMILIES is iterated against table columns, so it carries the CURRENT name only; the lookup dicts keep both keys so that
# pre-audit tables (e.g. the not-yet-regenerated c2l signatures) still resolve to a colour and a label.
FAMILIES = ["OPC_like", "AC_like", "OC_like", "MES_like", "cycling", "hypoxia", "novel_MP9", "unassigned"]
FAMILY_COLOURS = {"OPC_like": "#4575b4", "AC_like": "#e08214", "OC_like": "#1b7837", "MES_like": "#b2182b", "cycling": "#762a83",
                  "hypoxia": "#8c510a", "novel_MP9": "#35978f", "novel_MP15": "#35978f", "unassigned": "#bdbdbd"}
FAMILY_NAMES = {"OPC_like": "OPC-like", "AC_like": "AC-like", "OC_like": "OC-like", "MES_like": "MES-like", "cycling": "cycling",
                "hypoxia": "hypoxia", "novel_MP9": "ciliated (MP9)", "novel_MP15": "cilia (MP15, pre-audit)", "unassigned": "unassigned"}
ASSAY_COLOURS = {"scRNA": "#d95f02", "snRNA": "#1b9e77", "multiome_RNA": "#7570b3"}
ASSAY_NAMES = {"scRNA": "cells (scRNA)", "snRNA": "nuclei (snRNA)", "multiome_RNA": "nuclei (multiome)"}
STUDY_NAMES = {"collot": "Collot 2026", "desisto": "DeSisto", "jessa_multiome": "Jessa 2022 multiome", "jessa_scrna": "Jessa 2022",
               "okada": "Okada 2026", "zang": "Zang"}
STUDY_COLOURS = {"collot": "#66c2a5", "desisto": "#fc8d62", "jessa_multiome": "#8da0cb", "jessa_scrna": "#e78ac3", "okada": "#a6d854",
                 "zang": "#ffd92f"}
TME_NAME = {  # cleaned subcluster -> curated name (same mapping as 27_package_release.py; docs/phases/phase6_notes.md section 11)
    "myeloid_c0": "homeostatic microglia", "myeloid_c1": "TGFBI/AOAH macrophage", "myeloid_c2": "hypoxic/glycolytic TAM",
    "myeloid_c3": "interferon TAM", "myeloid_c4": "ribosomal residue (cells)", "myeloid_c5": "proliferating myeloid",
    "T_NK_c0": "THEMIS/CAMK4 T", "T_NK_c1": "CD8 T", "T_NK_c2": "CD4 T", "T_NK_c3": "NK", "T_NK_c4": "IL7R/KLRB1 T",
    "vascular_c0": "endothelial (neural genes)", "vascular_c1": "tip/angiogenic EC", "vascular_c2": "capillary EC",
    "vascular_c3": "capillary EC 2", "vascular_c4": "pericyte", "oligodendrocyte_c0": "mature oligodendrocyte 1",
    "oligodendrocyte_c1": "mature oligodendrocyte 2", "oligodendrocyte_c2": "mature oligodendrocyte 3",
    "oligodendrocyte_c3": "mature oligodendrocyte 4", "oligodendrocyte_c4": "committed OPC"}
COMPARTMENT_COLOURS = {"myeloid": LINEAGE_COLOURS["myeloid"], "T_NK": LINEAGE_COLOURS["T_NK"], "vascular": LINEAGE_COLOURS["vascular"],
                       "oligodendrocyte": LINEAGE_COLOURS["oligodendrocyte"]}
CONTRAST_NAMES = {"assay_nuclei_vs_cells": "nuclei vs cells", "h3_H3.1_vs_H3.3": "H3.1 vs H3.3", "age_adult_vs_pediatric": "adult vs paediatric",
                  "site_thalamus_vs_pons": "thalamus vs pons", "treatment_ONC201_vs_standard": "ONC201 vs standard"}
SPATIAL_CLASS = {"malignant_OPC_like": "OPC-like", "malignant_AC_like": "AC-like", "malignant_OC_like": "OC-like",
                 "malignant_MES_like": "MES-like", "malignant_cycling": "cycling", "malignant_hypoxia": "hypoxia",
                 "malignant_novel_MP9": "ciliated (MP9)", "malignant_novel_MP15": "cilia (MP15, pre-audit)", "T_NK": "T/NK"}


def letter(ax, s, dx=-0.02, dy=1.02, title=None, tpad=3.0):
    """Panel letter and panel title, Cell Press house style.

    The letter is upper-cased here so the convention flips in one place: call sites may keep passing 'a', 'b', ...
    `title` is the panel's short descriptive title (2-5 words, sentence case), centred above the axes; `tpad` is its
    distance from the axes top in points - raise it on panels that carry a legend above the axes.
    """
    ax.text(dx, dy, str(s).upper(), transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")
    if title:
        ax.set_title(title, loc="center", fontsize=7.5, pad=tpad)


def umap_ax(ax, *, arrows: bool = False, aspect: bool = True):
    """Embedding panel: no axis box, no ticks, no axis labels (venue model strips them entirely).

    `arrows=True` keeps the small corner arrow pair and its 'UMAP' label (used by the graphical abstract).
    """
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel("")
    for s in ax.spines.values():
        s.set_visible(False)
    if aspect:
        ax.set_aspect("equal", adjustable="datalim")
    if arrows:
        for xy in ((0.17, 0.015), (0.015, 0.17)):
            ax.annotate("", xy=xy, xytext=(0.015, 0.015), xycoords="axes fraction",
                        arrowprops=dict(arrowstyle="->", lw=0.7, color="#444444"))
        ax.text(0.19, 0.012, "UMAP", fontsize=6, transform=ax.transAxes, va="center", color="#444444")


def legend_dots(ax, colours: dict, names=None, **kw):
    import matplotlib.pyplot as plt
    hs = [plt.Line2D([], [], marker="o", ls="", ms=4, mfc=c, mec="none") for c in colours.values()]
    ax.legend(hs, [(names or {}).get(k, k) for k in colours], handletextpad=0.1, borderaxespad=0, **kw)


def heat(ax, m, cmap="Greys", vmin=0, vmax=None, fmt=None, norm=None, cbar_label=None, fig=None):
    import numpy as np
    im = ax.imshow(m.values, aspect="auto", cmap=cmap, vmin=None if norm else vmin, vmax=None if norm else vmax, norm=norm,
                   interpolation="nearest")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=60, ha="right")
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if fmt:
        hi = np.nanmax(np.abs(m.values)) if m.size else 1
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                v = m.values[i, j]
                if np.isfinite(v) and fmt(v):
                    ax.text(j, i, fmt(v), ha="center", va="center", fontsize=5, color="white" if abs(v) > 0.6 * hi else "black")
    if cbar_label and fig is not None:
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label(cbar_label); cb.outline.set_visible(False)
    return im


# Display names for publication figures (no software-internal identifiers on axes)
EMBED_NAMES = {"scvi_seed0": "scVI seed 0", "scvi_seed1": "scVI seed 1", "scvi_seed2": "scVI seed 2",
               "scvi_seed0_e60": "scVI seed 0, 60 epochs", "harmony": "Harmony", "pca": "PCA"}
_METRIC = {"asw_batch": "batch ASW", "asw_label": "label ASW", "ilisi": "iLISI", "clisi": "cLISI", "kbet": "kBET", "pcr": "PCR",
           "ari": "ARI", "nmi": "NMI", "graph_connectivity": "graph connectivity", "batch_score": "batch score",
           "bio_score": "biology score", "fine_bio_score": "fine biology score", "total_score": "overall score"}
_KEY = {"assay": "assay", "sample_id": "sample", "study": "study", "fine": "fine labels"}


def metric_name(m: str) -> str:
    """scIB column name -> readable label, e.g. 'ilisi__sample_id' -> 'iLISI (sample)'."""
    base, _, key = m.partition("__")
    return _METRIC.get(base, base.replace("_", " ")) + (f" ({_KEY.get(key, key)})" if key else "")
