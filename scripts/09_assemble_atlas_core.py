#!/usr/bin/env python
"""Phase 4, step 7 — assemble data/processed/atlas_core.h5ad and the UMAP overview figure.

  obs   = atlas_input obs + scanvi_label, scanvi_maxprob, leiden_scanvi_r1
  obsm  = X_scvi (chosen seed), X_scanvi, X_umap (UMAP on X_scanvi, n_neighbors 30), X_scvi_malignant (NaN for non-malignant)
  X / layers['counts'] = full-gene int32 counts (18,882 genes) when the file stays <= 6 GB, else HVG counts only.
UMAP / Leiden are computed once here (scanpy, igraph flavour) and stored.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path("/Users/jubayer/Projects/dmg_atlas")
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import anndata as ad  # noqa: E402
import scanpy as sc  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase4_integration.log", context="09_atlas_core")
MODELS = ROOT / "models"
OUT = ROOT / "data/processed"
FIG = ROOT / "results/figures"
FIG.mkdir(parents=True, exist_ok=True)
choice = json.load(open(ROOT / "results/tables/integration_model_choice.json"))
CHOSEN = choice["chosen"]


def load_latent(name, obs_names):
    d = MODELS / name
    z = np.load(d / "latent_f16.npy").astype(np.float32)
    names = pd.Index(open(d / "obs_names.txt").read().split("\n"))
    assert len(names) == len(z)
    return pd.DataFrame(z, index=names).reindex(obs_names)


def main():
    with LOG.step("read atlas_input.h5ad (full genes)"):
        import h5py
        from anndata.io import read_elem
        with h5py.File(OUT / "atlas_input.h5ad") as f:  # X and layers['counts'] are identical on disk; read one copy
            adata = ad.AnnData(X=read_elem(f["X"]), obs=read_elem(f["obs"]), var=read_elem(f["var"]), uns=read_elem(f["uns"]))
    obs_names = adata.obs_names

    with LOG.step("attach latents"):
        adata.obsm["X_scvi"] = load_latent(CHOSEN, obs_names).values
        adata.obsm["X_scanvi"] = load_latent("scanvi", obs_names).values
        zm = load_latent("scvi_malignant", obs_names)
        adata.obsm["X_scvi_malignant"] = zm.values  # NaN rows for non-malignant cells
        assert not np.isnan(adata.obsm["X_scvi"]).any() and not np.isnan(adata.obsm["X_scanvi"]).any()
        sd = MODELS / "scanvi"
        names = pd.Index(open(sd / "obs_names.txt").read().split("\n"))
        pred = pd.Series(open(sd / "predicted_labels.txt").read().split("\n"), index=names).reindex(obs_names)
        maxp = pd.Series(np.load(sd / "maxprob.npy"), index=names).reindex(obs_names)
        adata.obs["scanvi_label"] = pred.astype("category").values
        adata.obs["scanvi_maxprob"] = maxp.astype(np.float32).values
        adata.uns["phase4"]["chosen_scvi"] = CHOSEN
        LOG.write(f"malignant latent rows: {int((~zm.isna().any(axis=1)).sum())}; scanvi labels {adata.obs['scanvi_label'].value_counts().to_dict()}")

    with LOG.step("neighbors (X_scanvi, k=30) + UMAP + Leiden r=1.0"):
        sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=30, random_state=0)
        sc.tl.umap(adata, random_state=0)
        sc.tl.leiden(adata, resolution=1.0, key_added="leiden_scanvi_r1", flavor="igraph", n_iterations=2, directed=False, random_state=0)
        LOG.write(f"leiden clusters: {adata.obs['leiden_scanvi_r1'].nunique()}")

    with LOG.step("UMAP overview figure"):
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        xy = adata.obsm["X_umap"]
        order = np.random.default_rng(0).permutation(adata.n_obs)
        for ax, key in zip(axes.ravel(), ["study", "assay", "lineage_label", "consensus_call"]):
            cats = adata.obs[key].astype(str)
            uniq = sorted(cats.unique(), key=lambda s: (s == "Unknown", s))
            cmap = plt.get_cmap("tab20" if len(uniq) > 10 else "tab10")
            colors = {u: ("lightgrey" if u in ("Unknown", "ambiguous") else cmap(i % 20)) for i, u in enumerate(uniq)}
            ax.scatter(xy[order, 0], xy[order, 1], s=0.3, c=[colors[v] for v in cats.values[order]], linewidths=0, rasterized=True)
            handles = [plt.Line2D([], [], marker="o", ls="", color=colors[u], markersize=6, label=f"{u} ({(cats == u).sum():,})") for u in uniq]
            ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8, title=key, title_fontsize=9)
            ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
        fig.tight_layout()
        fig.savefig(FIG / "phase4_umap_overview.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    with LOG.step("write atlas_core.h5ad"):
        # one on-disk copy of the full-gene counts: layers['counts'] (int32 CSR, 18,882 genes), X left empty.
        # Writing X and layers['counts'] both would double the file to ~10 GB (> 6 GB limit).
        counts = adata.X
        del adata.X
        adata.layers["counts"] = counts
        del adata.obsp["connectivities"], adata.obsp["distances"]
        adata.uns.pop("neighbors", None)
        adata.uns["phase4"]["counts_location"] = "layers['counts'] = int32 UMI counts for all 18,882 shared genes; X is empty (set X from the layer before normalising)"
        adata.write_h5ad(OUT / "atlas_core.h5ad")
        LOG.write(f"atlas_core.h5ad {(OUT / 'atlas_core.h5ad').stat().st_size / 1e9:.2f} GB")
    adata.obs[["sample_id", "study", "assay", "lineage_label", "consensus_call", "scanvi_label", "scanvi_maxprob", "leiden_scanvi_r1"]].to_csv(ROOT / "results/tables/atlas_core_obs_labels.csv")
    pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1", "UMAP2"]).to_csv(ROOT / "results/tables/atlas_core_umap.csv")
    LOG.write("done")


if __name__ == "__main__":
    main()
