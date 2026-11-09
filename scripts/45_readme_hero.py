"""README hero figure: the integrated reference as a single UMAP.

Inputs : results/tables/atlas_core_obs_labels.csv, results/tables/atlas_core_umap.csv (written by 09_assemble_atlas_core.py)
Output : results/pub_figs/readme_atlas_umap.{png,pdf}

One embedding, coloured by scANVI lineage label, with no panel letter and no in-figure title - the figure carries
the colour key only and everything else belongs to the README caption. House style from src/dmg_atlas/figstyle.py,
so the colours match the publication figures.

Usage: python scripts/45_readme_hero.py [--force]
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from dmg_atlas.figstyle import (  # noqa: E402
    LINEAGE_COLOURS, LINEAGES, W, apply_figure_style, umap_ax,
)

T = ROOT / "results" / "tables"
OUT = ROOT / "results" / "pub_figs"
LINEAGE_NAMES = {"malignant": "malignant", "myeloid": "myeloid", "oligodendrocyte": "oligodendrocyte",
                 "astrocyte": "astrocyte", "neuron": "neuron", "vascular": "vascular", "T_NK": "T / NK"}


def anchor(d, bins=140, sigma=2.0):
    """Label anchor = densest bin of this lineage, so the text sits on the cluster it names."""
    from scipy.ndimage import gaussian_filter
    h, xe, ye = np.histogram2d(d["UMAP1"], d["UMAP2"], bins=bins)
    h = gaussian_filter(h, sigma)
    i, j = np.unravel_index(h.argmax(), h.shape)
    return (xe[i] + xe[i + 1]) / 2, (ye[j] + ye[j + 1]) / 2


def label_lineages(ax, df, key, colours, names, min_gap=2.2, fontsize=10.0):
    """In-plot labels: neutral dark text on a white plate outlined in the lineage colour.

    Colouring the text itself is unreadable here because every label sits on points of its own colour; a plain dark
    label would be legible but would sever the colour-to-name link. The plate carries both. Anchors are pushed apart
    so neighbouring names do not overlap.
    """
    pts = {k: np.array(anchor(df[df[key] == k])) for k in names if (df[key] == k).any()}
    for _ in range(300):  # short repulsion pass between label anchors only
        moved = False
        for a in pts:
            for b in pts:
                if a >= b:
                    continue
                v = pts[a] - pts[b]
                dist = float(np.hypot(*v)) or 1e-6
                if dist < min_gap:
                    step = v / dist * (min_gap - dist) / 2
                    pts[a] += step
                    pts[b] -= step
                    moved = True
        if not moved:
            break
    for k, (x, y) in pts.items():
        ax.text(x, y, names[k], color="#1a1a1a", fontsize=fontsize, fontweight="bold",
                ha="center", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.30", facecolor="white",
                          edgecolor=colours[k], linewidth=1.3, alpha=0.94))


def scatter(ax, df, key, colours, rng):
    order = rng.permutation(len(df))  # no study drawn systematically on top of another
    d = df.iloc[order]
    ax.scatter(d["UMAP1"], d["UMAP2"], s=0.45, lw=0, rasterized=True, alpha=0.55,
               c=[colours.get(v, "#bdbdbd") for v in d[key]])
    umap_ax(ax)


def main() -> None:
    obs = pd.read_csv(T / "atlas_core_obs_labels.csv", index_col=0, usecols=[0, 2, 6])
    df = obs.join(pd.read_csv(T / "atlas_core_umap.csv", index_col=0), how="inner")
    assert df[["UMAP1", "UMAP2"]].notna().all().all(), "missing UMAP coordinates"
    rng = np.random.default_rng(0)

    apply_figure_style()
    fig, ax = plt.subplots(figsize=(W / 1.5, W / 1.5))

    scatter(ax, df, "scanvi_label", LINEAGE_COLOURS, rng)
    label_lineages(ax, df, "scanvi_label", LINEAGE_COLOURS, LINEAGE_NAMES)

    fig.tight_layout(pad=0.2)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"readme_atlas_umap.{ext}")
    plt.close(fig)
    print(f"cells drawn: {len(df):,} | studies: {df['study'].nunique()} | lineages: {df['scanvi_label'].nunique()}")
    print(f"wrote {OUT / 'readme_atlas_umap.png'}")


if __name__ == "__main__":
    main()
