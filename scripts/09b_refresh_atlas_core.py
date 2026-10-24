#!/usr/bin/env python
"""Audit remediation (reference-labels track, 19 Oct 2026) — refresh data/processed/atlas_core.h5ad after the scANVI retrain.

  core       obs: lineage_label, scanvi_label, scanvi_maxprob, leiden_scanvi_r1 replaced (previous kept as *_pre_audit); CNV-derived call
             columns (consensus_call, consensus_celltype_coarse, cnv_call, marker_call, marker_type, cnv_score, cnv_cluster, chr*/arm_* events,
             score_*) replaced from the corrected malignant_calls.csv (categorical call columns kept as *_pre_audit; chr4q12_gain -> chr4q12_gain_pre_audit,
             the event was removed; new: cnv_confidence, consensus_call_eventtrust, low_confidence_sample);
             obsm: X_scanvi replaced (previous X_scanvi_pre_audit), X_umap recomputed on the new X_scanvi (n_neighbors 30, scanpy UMAP, seed 0;
             previous X_umap_pre_audit), Leiden 1.0 (igraph) recomputed; uns['phase4']['audit_fix_2026_10_19'] provenance.
             The file is rebuilt (layers/counts and var copied byte-wise with h5py) and the previous file renamed atlas_core.pre_audit.h5ad.
             Tables: results/tables/atlas_core_obs_labels.csv, atlas_core_umap.csv (previous *.pre_audit.csv); figure results/figures/phase4_umap_overview.png.
  malignant  obsm['X_scvi_malignant'] replaced in place from models/scvi_malignant (NaN for non-malignant and for okada__1734_pons cells;
             previous kept as X_scvi_malignant_pre_audit).
"""
from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402
from dmg_atlas import figstyle  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import anndata as ad  # noqa: E402
from anndata.io import read_elem, write_elem  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_reflabels.log", context="09b_refresh_atlas_core")
OUT = ROOT / "data/processed"
CORE = OUT / "atlas_core.h5ad"
TAB = ROOT / "results/tables"
FIG = ROOT / "results/figures"
LABELS = TAB / "reflabels/reference_labels_corrected.csv.gz"
MODELS = ROOT / "models"
CALL_CAT = ["consensus_call", "consensus_celltype_coarse", "cnv_call", "marker_call", "marker_type"]
CALL_NUM = ["cnv_score", "cnv_cluster", "chr1q_gain", "chr2_gain", "chr7_gain", "chr13q_loss", "chr14q_loss", "chr10q_loss", "chr17p_loss",
            "arm_1q", "arm_2", "arm_7", "arm_4q12", "arm_13q", "arm_14q", "arm_10q", "arm_17p",
            "score_myeloid", "score_T_NK", "score_oligodendrocyte", "score_vascular", "score_neuron", "score_astrocyte",
            "score_OPC_like", "score_AC_like", "score_MES_like", "score_cycling"]
CALL_NEW = ["cnv_confidence", "consensus_call_eventtrust", "low_confidence_sample", "gene_dev_PDGFRA", "gene_dev_KIT"]


def load_latent(name, obs_names):
    d = MODELS / name
    z = np.load(d / "latent_f16.npy").astype(np.float32)
    names = pd.Index((d / "obs_names.txt").read_text().split("\n"))
    assert len(names) == len(z)
    return pd.DataFrame(z, index=names).reindex(obs_names)


def backup(path: Path):
    if path.exists():
        pre = path.with_name(f"{path.stem}.pre_audit{path.suffix}")
        if not pre.exists():
            shutil.copy2(path, pre)


def core():
    with LOG.step("read atlas_core obs/obsm/uns (counts untouched)"):
        with h5py.File(CORE) as f:
            obs = read_elem(f["obs"])
            uns = read_elem(f["uns"])
            obsm = {k: read_elem(f["obsm"][k]) for k in f["obsm"]}
        obs_names = obs.index
        for c in obs.columns:
            if isinstance(obs[c].dtype, pd.CategoricalDtype):
                obs[c] = obs[c].astype(object)

    with LOG.step("replace labels and CNV call columns"):
        lab = pd.read_csv(LABELS, index_col="cell_id").loc[obs_names]
        assert (obs["lineage_label"].astype(str).values == lab["lineage_label_pre_audit"].astype(str).values).all()
        calls = pd.read_csv(OUT / "malignant_calls.csv", index_col="cell_id", usecols=["cell_id"] + CALL_CAT + CALL_NUM + CALL_NEW, low_memory=False).loc[obs_names]
        for c in ["lineage_label", "scanvi_label", "scanvi_maxprob", "leiden_scanvi_r1"] + CALL_CAT + ["chr4q12_gain"]:
            obs[f"{c}_pre_audit"] = obs[c].values
        obs = obs.drop(columns=["chr4q12_gain"])
        obs["lineage_label"] = lab["lineage_label"].astype(str).values
        for c in CALL_CAT + CALL_NUM + CALL_NEW:
            obs[c] = calls[c].values
        sd = MODELS / "scanvi"
        names = pd.Index((sd / "obs_names.txt").read_text().split("\n"))
        pred = pd.Series((sd / "predicted_labels.txt").read_text().split("\n"), index=names).reindex(obs_names)
        maxp = pd.Series(np.load(sd / "maxprob.npy"), index=names).reindex(obs_names)
        assert pred.notna().all()
        obs["scanvi_label"] = pred.astype(str).values
        obs["scanvi_maxprob"] = maxp.astype(np.float32).values
        split = pd.read_csv(sd / "split_cells.csv.gz", index_col="cell_id")["set"].reindex(obs_names)
        obs["scanvi_validation_patient"] = (split == "validation").values
        obs["low_confidence_sample"] = obs["low_confidence_sample"].astype(bool)
        n_ll = int((obs["lineage_label"] != obs["lineage_label_pre_audit"]).sum())
        n_sl = int((obs["scanvi_label"] != obs["scanvi_label_pre_audit"]).sum())
        LOG.write(f"lineage_label changed {n_ll:,}; scanvi_label changed {n_sl:,}; consensus_call changed "
                  f"{int((obs['consensus_call'].astype(str) != obs['consensus_call_pre_audit'].astype(str)).sum()):,}; "
                  f"scanvi_label counts {obs['scanvi_label'].value_counts().to_dict()}")

    with LOG.step("new X_scanvi, neighbors (k=30), UMAP, Leiden r=1.0"):
        import scanpy as sc
        sc.settings.n_jobs = 8
        Z = load_latent("scanvi", obs_names).values.astype(np.float32)
        assert not np.isnan(Z).any()
        obsm["X_scanvi_pre_audit"] = obsm["X_scanvi"]
        obsm["X_umap_pre_audit"] = obsm["X_umap"]
        obsm["X_scanvi"] = Z
        a = ad.AnnData(obs=pd.DataFrame(index=obs_names))
        a.obsm["X_scanvi"] = Z
        sc.pp.neighbors(a, use_rep="X_scanvi", n_neighbors=30, random_state=0)
        sc.tl.umap(a, random_state=0)
        sc.tl.leiden(a, resolution=1.0, key_added="leiden_scanvi_r1", flavor="igraph", n_iterations=2, directed=False, random_state=0)
        obsm["X_umap"] = a.obsm["X_umap"].astype(np.float32)
        obs["leiden_scanvi_r1"] = a.obs["leiden_scanvi_r1"].astype(str).values
        uns["umap"] = a.uns["umap"]
        uns["leiden_scanvi_r1"] = a.uns["leiden_scanvi_r1"] if "leiden_scanvi_r1" in a.uns else {"params": {"resolution": 1.0, "random_state": 0, "n_iterations": 2}}
        LOG.write(f"leiden clusters: {obs['leiden_scanvi_r1'].nunique()} (pre-audit {obs['leiden_scanvi_r1_pre_audit'].nunique()}); peak RSS {peak_rss_gb():.2f} GB")
        del a

    with LOG.step("provenance + categorise"):
        run = json.loads((MODELS / "scanvi/run.json").read_text())
        metrics = json.loads((MODELS / "scanvi/validation_metrics.json").read_text())
        uns["phase4"]["audit_fix_2026_10_19"] = {
            "lineage_label": "rebuilt from the corrected data/processed/malignant_calls.csv (CNV track: per-arm smoothing, NaN handling, chr4q12_gain removed); "
                             f"{n_ll:,} cells changed; previous in lineage_label_pre_audit (scripts/07c_rebuild_lineage_labels.py)",
            "scanvi": f"models/scanvi retrained from models/scvi_seed1 with corrected labels, patient-grouped validation ({len(json.loads((MODELS / 'scanvi/split.json').read_text())['validation_patients'])} held-out patients), "
                      f"early stopping on validation_classification_loss (best epoch {run['best_epoch']} of {run['epochs_run']}); train acc {metrics['train_accuracy']:.4f}, "
                      f"held-out acc {metrics['validation_accuracy']:.4f}; previous model models/scanvi_pre_audit, previous columns *_pre_audit (modal/phase4b_scanvi_retrain.py)",
            "cnv_columns": "consensus_call, consensus_celltype_coarse, cnv_call, marker_call, marker_type, cnv_score, cnv_cluster, chr*/arm_* and score_* columns replaced from the corrected "
                           "malignant_calls.csv; categorical call columns kept as *_pre_audit; chr4q12_gain removed (kept as chr4q12_gain_pre_audit); cnv_confidence, consensus_call_eventtrust, low_confidence_sample added",
            "okada__1734_pons": "kept in the scVI/scANVI reference (already in the latent) but low_confidence_sample=True; excluded from models/scvi_malignant and to be excluded from Phase 6 analyses",
            "X_umap": "recomputed on the new X_scanvi (scanpy neighbors k=30, UMAP seed 0); previous X_umap_pre_audit; leiden_scanvi_r1 recomputed (previous *_pre_audit)",
            "scanvi_validation_patient": "True for cells of the 5 patients held out from scANVI classifier training (models/scanvi/split.json)",
        }
        for c in obs.columns:
            if obs[c].dtype == object:
                obs[c] = obs[c].astype(str).where(obs[c].notna(), None).astype("category")

    with LOG.step("rebuild atlas_core.h5ad (counts/var copied byte-wise); previous -> atlas_core.pre_audit.h5ad"):
        new = CORE.with_name("atlas_core.new.h5ad")
        if new.exists():
            new.unlink()
        with h5py.File(CORE) as src, h5py.File(new, "w") as dst:
            for k, v in src.attrs.items():
                dst.attrs[k] = v
            for k in src.keys():
                if k in ("obs", "obsm", "uns"):
                    continue
                src.copy(src[k], dst, name=k)
            write_elem(dst, "obs", obs)
            write_elem(dst, "obsm", obsm)
            write_elem(dst, "uns", uns)
        chk = ad.read_h5ad(new, backed="r")
        assert chk.n_obs == len(obs_names) and chk.obs_names.equals(obs_names)
        chk.file.close()
        pre = CORE.with_name("atlas_core.pre_audit.h5ad")
        if not pre.exists():
            os.rename(CORE, pre)
        else:
            CORE.unlink()
        os.rename(new, CORE)
        LOG.write(f"atlas_core.h5ad {CORE.stat().st_size / 1e9:.2f} GB; pre-audit {pre.stat().st_size / 1e9:.2f} GB; peak RSS {peak_rss_gb():.2f} GB")

    with LOG.step("tables + UMAP overview figure"):
        for t in ["atlas_core_obs_labels.csv", "atlas_core_umap.csv"]:
            backup(TAB / t)
        cols = ["sample_id", "study", "assay", "lineage_label", "consensus_call", "scanvi_label", "scanvi_maxprob", "leiden_scanvi_r1",
                "lineage_label_pre_audit", "consensus_call_pre_audit", "scanvi_label_pre_audit", "scanvi_maxprob_pre_audit", "leiden_scanvi_r1_pre_audit",
                "cnv_confidence", "consensus_call_eventtrust", "low_confidence_sample", "scanvi_validation_patient"]
        obs[cols].to_csv(TAB / "atlas_core_obs_labels.csv")
        pd.DataFrame(obsm["X_umap"], index=obs_names, columns=["UMAP1", "UMAP2"]).to_csv(TAB / "atlas_core_umap.csv")
        figure(obs, obsm["X_umap"])
    LOG.write("core done")


def figure(obs, xy):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figstyle.apply_figure_style()
    backup(FIG / "phase4_umap_overview.png")
    panels = [("study", figstyle.STUDY_COLOURS, figstyle.STUDY_NAMES), ("assay", figstyle.ASSAY_COLOURS, figstyle.ASSAY_NAMES),
              ("lineage_label", figstyle.LINEAGE_COLOURS, None),
              ("consensus_call", {"malignant": "#b2182b", "normal": "#2166ac", "ambiguous": "#bdbdbd"}, None)]
    fig, axes = plt.subplots(2, 2, figsize=(figstyle.W, figstyle.W * 0.95))
    order = np.random.default_rng(0).permutation(len(obs))
    for ax, (key, colours, names), letter in zip(axes.ravel(), panels, "abcd"):
        cats = obs[key].astype(str).values
        cols = {k: v for k, v in colours.items() if k in set(cats)}
        for k in set(cats) - set(cols):
            cols[k] = "#bdbdbd"
        ax.scatter(xy[order, 0], xy[order, 1], s=0.15, c=[cols[v] for v in cats[order]], linewidths=0, rasterized=True)
        figstyle.umap_ax(ax)
        figstyle.letter(ax, letter)
        counts = pd.Series(cats).value_counts()
        lab_cols = {k: cols[k] for k in cols if k in counts.index}
        lab_names = {k: f"{(names or {}).get(k, k)} ({counts[k]:,})" for k in lab_cols}
        figstyle.legend_dots(ax, lab_cols, lab_names, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=6, markerscale=1.2)
    fig.subplots_adjust(wspace=0.9, hspace=0.15)
    fig.savefig(FIG / "phase4_umap_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    LOG.write("figure results/figures/phase4_umap_overview.png written (caption in docs/phase4_integration_notes.md Audit corrections)")


def malignant():
    with LOG.step("replace obsm['X_scvi_malignant'] in place from models/scvi_malignant"):
        with h5py.File(CORE, "r+") as f:
            obs_names = pd.Index(read_elem(f["obs/_index"]).astype(str))
            zm = load_latent("scvi_malignant", obs_names)
            n_rows = int((~zm.isna().any(axis=1)).sum())
            if "X_scvi_malignant_pre_audit" not in f["obsm"]:
                f["obsm"].move("X_scvi_malignant", "X_scvi_malignant_pre_audit")
            else:
                del f["obsm/X_scvi_malignant"]
            write_elem(f["obsm"], "X_scvi_malignant", zm.values.astype(np.float32))
            uns = read_elem(f["uns"])
            run = json.loads((MODELS / "scvi_malignant/run.json").read_text())
            uns["phase4"]["audit_fix_2026_10_19"]["X_scvi_malignant"] = (
                f"models/scvi_malignant retrained (seed 1, {run['n_cells']:,} corrected consensus-malignant cells, okada__1734_pons excluded, "
                f"{run['epochs_run']} epochs, elbo_validation {run['elbo_validation_final']:.1f}); rows are NaN for the {len(obs_names) - n_rows:,} cells not in the "
                "malignant-only model (non-malignant and okada__1734_pons); previous latent X_scvi_malignant_pre_audit (scripts/07d_rebuild_malignant_input.py, modal/phase4_train.py)")
            del f["uns"]
            write_elem(f, "uns", uns)
        LOG.write(f"X_scvi_malignant rows with values {n_rows:,} of {len(obs_names):,}")
    LOG.write("malignant done")


if __name__ == "__main__":
    {"core": core, "malignant": malignant}[sys.argv[1]]()
