#!/usr/bin/env python
"""Phase 6 / step 15c (audit remediation 2026-10-20; findings SP-04, SP-05, SP-18) — the three Phase 6 spatial figures,
rebuilt from the corrected (metaprogram-family, 20,000-epoch) deconvolution.

  phase6_spatial_abundance_maps.png   SP-05: spots plotted at full-resolution PIXEL coordinates (pxl_col, pxl_row)
                                      instead of (array_col, -array_row), equal aspect, per-panel 1 mm scale bar
                                      derived from spot_diameter_fullres (Visium spot = 55 um), no in-panel text;
                                      panel identities are given in the caption written to docs/phase6_notes.md.
  phase6_brainstem_gradient.png       SP-04: one shared legend OUTSIDE both panels (below), both y-axes labelled,
                                      OC-like drawn in a dark green (#1b7837) instead of #fddbc7.
  phase6_colocalisation_matrix.png    SP-01/SP-18: the PROPORTION-based random-effects matrix (7 families x 6 TME),
                                      cells annotated with rho and the number of patients with positive rho.
House style: fig, ax = plt.subplots; legends outside; no panel titles; 300 dpi (src/dmg_atlas/figstyle.py).
"""
from __future__ import annotations
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style, FAMILY_COLOURS, FAMILY_NAMES, LINEAGE_COLOURS  # noqa: E402
import numpy as np, pandas as pd, matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_spatial.log", context="15c_spatial_figures_fix")
TABLES = ROOT / "results/tables"; FIGS = ROOT / "results/figures"
FAMILIES = ["OPC_like", "AC_like", "OC_like", "MES_like", "cycling", "hypoxia", "novel_MP9"]
STATES = [f"malignant_{f}" for f in FAMILIES]
TME = ["myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular"]
CLASSES = STATES + TME
SHORT = {**{f"malignant_{f}": FAMILY_NAMES[f] for f in FAMILIES},
         "myeloid": "myeloid", "T_NK": "T/NK", "oligodendrocyte": "oligodendrocyte", "astrocyte": "astrocyte",
         "neuron": "neuron", "vascular": "vascular"}
COLOURS = {**{f"malignant_{f}": FAMILY_COLOURS[f] for f in FAMILIES},
           **{k: LINEAGE_COLOURS[k] for k in ("myeloid", "oligodendrocyte", "astrocyte", "neuron", "vascular", "T_NK")}}
BRAINSTEM_ORDER = ["Midbrain", "Pons_A", "Pons_B", "Pons_C", "Pons_D", "Medulla_A", "Medulla_B", "Medulla_C",
                   "Spinal_Cord_A", "Spinal_Cord_B"]
PICKS = [("huchede", "pDMG_Sample-2"), ("ren", "DMG1"), ("collot_st", "BT01_192_core_HE"), ("brainstem", "Pons_A")]
SPOT_UM = 55.0  # Visium spot diameter


def fig_maps():
    apply_figure_style()
    used = pd.read_csv(ROOT / "results/qc/spatial_sections_used.csv").set_index("section_id")
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.5))
    for ax, (cohort, sec) in zip(axes, PICKS):
        df = pd.read_csv(ROOT / f"data/processed/spatial/{cohort}/{sec}_abundance.csv", index_col=0)
        q = df[[f"q05_{c}" for c in CLASSES]]; q.columns = CLASSES
        dom = q.idxmax(1)
        for c in CLASSES:
            k = (dom == c).values
            if k.any():
                ax.scatter(df.loc[k, "pxl_col"], -df.loc[k, "pxl_row"], s=2.2, c=COLOURS[c], linewidths=0,
                           label=SHORT[c], rasterized=True)
        ax.set_aspect("equal")
        px_per_mm = used.loc[sec, "spot_diameter_fullres"] / SPOT_UM * 1000.0
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        xs = x0 + 0.06 * (x1 - x0); ys = y0 + 0.06 * (y1 - y0)
        ax.plot([xs, xs + px_per_mm], [ys, ys], color="black", lw=1.4, solid_capstyle="butt")
        ax.set_xlim(min(x0, xs - 0.02 * (x1 - x0)), max(x1, xs + px_per_mm + 0.02 * (x1 - x0)))
        ax.set_xticks([]); ax.set_yticks([])
        for s_ in ax.spines.values():
            s_.set_visible(False)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=4, color=COLOURS[c], label=SHORT[c]) for c in CLASSES]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.99, 0.5), ncol=1, handletextpad=0.3)
    fig.savefig(FIGS / "phase6_spatial_abundance_maps.png")
    plt.close(fig)
    LOG.write("phase6_spatial_abundance_maps.png: pixel coordinates, equal aspect, 1 mm scale bars "
              + ", ".join(f"{s} {used.loc[s,'spot_diameter_fullres']/SPOT_UM*1000:.0f} px/mm" for _, s in PICKS))


def fig_gradient():
    apply_figure_style()
    sec_tab = pd.read_csv(TABLES / "brainstem_section_means.csv", index_col=0).loc[BRAINSTEM_ORDER]
    x = np.arange(1, 11)
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 4.2), sharex=True)
    for c in STATES:
        axes[0].plot(x, sec_tab[f"prop_{c}"].values, marker="o", ms=3, color=COLOURS[c], label=SHORT[c])
    axes[0].plot(x, sec_tab["tumour_fraction_mean"].values, marker="s", ms=3, color="black", ls="--",
                 label="tumour fraction")
    for c in TME:
        axes[1].plot(x, sec_tab[f"prop_{c}"].values, marker="o", ms=3, color=COLOURS[c], label=SHORT[c])
    for ax in axes:
        ax.set_xticks(x); ax.set_xticklabels([s.replace("_", " ") for s in BRAINSTEM_ORDER], rotation=45, ha="right")
        ax.set_xlabel("anatomical position (rostral → caudal)")
        ax.set_ylabel("mean share of spot cell abundance (q05)")
    handles = ([plt.Line2D([], [], marker="o", ls="-", ms=3, color=COLOURS[c], label=SHORT[c]) for c in STATES]
               + [plt.Line2D([], [], marker="s", ls="--", ms=3, color="black", label="tumour fraction")]
               + [plt.Line2D([], [], marker="o", ls="-", ms=3, color=COLOURS[c], label=SHORT[c]) for c in TME])
    # Reserve the bottom third of the canvas for the rotated tick labels, the x-axis label and the shared
    # legend, so none of the three collide (the legend previously overlapped the x-axis labels).
    fig.subplots_adjust(bottom=0.46, wspace=0.28)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.14), ncol=5, handletextpad=0.3,
               columnspacing=1.2, frameon=False)
    fig.savefig(FIGS / "phase6_brainstem_gradient.png")
    plt.close(fig)
    LOG.write("phase6_brainstem_gradient.png: shared legend below both panels, both y-axes labelled, "
              f"OC-like colour {COLOURS['malignant_OC_like']}")


def fig_coloc():
    apply_figure_style()
    d = pd.read_csv(TABLES / "spatial_colocalisation.csv")
    d = d[(d.group == "DMG") & (d.metric == "proportion") & (d.pair_type == "state_vs_TME")]
    M = d.pivot(index="class_a", columns="class_b", values="rho_re").loc[STATES, TME]
    N = d.pivot(index="class_a", columns="class_b", values="n_sections_positive").loc[STATES, TME]
    nsec = int(d.n_sections.max())
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    v = float(np.nanmax(np.abs(M.values)))
    im = ax.imshow(M.values.astype(float), cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(len(TME))); ax.set_xticklabels([SHORT[c] for c in TME], rotation=45, ha="right")
    ax.set_yticks(range(len(STATES))); ax.set_yticklabels([SHORT[c] for c in STATES])
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M.values[i, j]:.2f}\n{int(N.values[i, j])}/{nsec}", ha="center", va="center", fontsize=5,
                    color="white" if abs(M.values[i, j]) > 0.6 * v else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cb.set_label("random-effects Spearman ρ (spot proportions)"); cb.outline.set_visible(False)
    fig.savefig(FIGS / "phase6_colocalisation_matrix.png")
    plt.close(fig)
    LOG.write(f"phase6_colocalisation_matrix.png: proportion-based, {M.shape[0]}x{M.shape[1]} pairs, "
              f"rho range {np.nanmin(M.values):.3f}..{np.nanmax(M.values):.3f}")


if __name__ == "__main__":
    import shutil
    for f in ("phase6_spatial_abundance_maps.png", "phase6_brainstem_gradient.png", "phase6_colocalisation_matrix.png"):
        p = FIGS / f; pre = FIGS / f.replace(".png", ".pre_audit.png")
        if p.exists() and not pre.exists():
            shutil.copy(p, pre)
    fig_maps(); fig_gradient(); fig_coloc()
