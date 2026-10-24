#!/usr/bin/env python
"""Phase 4, step 4 — scIB-style integration metrics for scVI (3 seeds) vs Harmony vs unintegrated PCA.

Embeddings: models/<name>/latent_f16.npy (+ obs_names.txt), names = scvi_seed0, scvi_seed1, scvi_seed2, harmony, pca.
Metrics (scib-metrics 0.5.10, jax CPU) on a stratified subsample (by sample_id) of labelled cells
(lineage_label != Unknown; N_SUB cells), kNN graph with 90 neighbours (pynndescent) per embedding:
  batch (per batch key study / assay / sample_id): iLISI, kBET acceptance (per label), batch ASW, PCR comparison vs PCA
  batch (label-based, batch-key independent): graph connectivity
  bio: cLISI, label ASW, NMI / ARI (Leiden, resolution optimised), on lineage_label
  bio (finer, Jessa malignant cells only, author_Cell_type_consensus_Jessa2022): cLISI, label ASW, NMI / ARI
Scores: each metric min-max scaled across the five embeddings (scIB convention), batch score = mean of scaled batch
metrics (all batch keys), bio score = mean of scaled bio metrics (lineage level), total = 0.4 batch + 0.6 bio.
Selection rule: the scVI seed with the highest total; ties broken by lower best validation ELBO.
Writes results/tables/integration_metrics.csv (long: embedding, metric, batch_key, label_key, value) and
results/tables/integration_metrics_summary.csv (wide + scaled scores + rank).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path("/Users/jubayer/Projects/dmg_atlas")
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402
import h5py  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase4_integration.log", context="08_metrics")
TAB = ROOT / "results/tables"
MODELS = ROOT / "models"
EMBEDDINGS = ["scvi_seed0", "scvi_seed1", "scvi_seed2", "harmony", "pca", "scvi_seed0_e60"]  # last = sensitivity run (60 epochs), not eligible for selection
SELECTABLE = ["scvi_seed0", "scvi_seed1", "scvi_seed2"]
BATCH_KEYS = ["study", "assay", "sample_id"]
N_SUB = 60_000
N_SUB_FINE = 30_000
K = 90
SEED = 0
BATCH_METRICS = ["ilisi", "kbet", "asw_batch", "pcr", "graph_connectivity"]
BIO_METRICS = ["clisi", "asw_label", "nmi", "ari"]


def stratified_idx(groups: pd.Series, n_total: int, rng, floor: int = 50) -> np.ndarray:
    counts = groups.value_counts()
    frac = n_total / len(groups)
    take = {g: int(min(c, max(floor, round(c * frac)))) for g, c in counts.items()}
    idx = []
    for g, n in take.items():
        pos = np.flatnonzero(groups.values == g)
        idx.append(rng.choice(pos, n, replace=False))
    return np.sort(np.concatenate(idx))


def load_embedding(name: str, obs_names: pd.Index) -> np.ndarray:
    d = MODELS / name
    z = np.load(d / "latent_f16.npy").astype(np.float32)
    names = pd.Index(open(d / "obs_names.txt").read().split("\n"))
    assert len(names) == z.shape[0], name
    if not names.equals(obs_names):
        z = z[names.get_indexer(obs_names)]
    return z


def main():
    import scib_metrics as sm
    from scib_metrics.nearest_neighbors import pynndescent

    with h5py.File(ROOT / "data/processed/atlas_input_hvg.h5ad") as f:
        obs = read_elem(f["obs"])
    rng = np.random.default_rng(SEED)
    lab_all = obs["lineage_label"].astype(str)
    labelled = np.flatnonzero((lab_all != "Unknown").values)
    sub = labelled[stratified_idx(obs["sample_id"].astype(str).iloc[labelled].reset_index(drop=True), N_SUB, rng)]
    LOG.write(f"metric subsample: {len(sub)} of {len(labelled)} labelled cells, {obs['sample_id'].iloc[sub].nunique()} samples; "
              f"label counts {lab_all.iloc[sub].value_counts().to_dict()}")
    # finer check: Jessa malignant cells with author consensus cell type
    jm = obs["study"].astype(str).str.startswith("jessa") & (obs["consensus_call"].astype(str) == "malignant")
    fine_lab = obs["author_Cell_type_consensus_Jessa2022"].astype(str)
    jm &= ~fine_lab.isin(["nan", "None", "", "Unresolved", "Uncertain", "Normal"])
    jm_idx = np.flatnonzero(jm.values)
    fine_sub = jm_idx[stratified_idx(obs["sample_id"].astype(str).iloc[jm_idx].reset_index(drop=True), N_SUB_FINE, rng)]
    LOG.write(f"fine subsample (Jessa malignant, author cell type): {len(fine_sub)} cells; labels {fine_lab.iloc[fine_sub].value_counts().to_dict()}")

    emb = {n: load_embedding(n, obs.index) for n in EMBEDDINGS}
    rows = []
    knn = {}
    for name in EMBEDDINGS:
        X = emb[name][sub]
        labels = lab_all.iloc[sub].values
        with LOG.step(f"{name}: kNN k={K} on {len(sub)} cells"):
            nn = pynndescent(X, n_neighbors=K, random_state=SEED, n_jobs=8)
            knn[name] = nn
        with LOG.step(f"{name}: bio metrics"):
            rows.append((name, "clisi", "", "lineage_label", sm.clisi_knn(nn, labels)))
            rows.append((name, "asw_label", "", "lineage_label", sm.silhouette_label(X, labels)))
            na = sm.nmi_ari_cluster_labels_leiden(nn, labels, optimize_resolution=True, n_jobs=1, seed=SEED)
            rows.append((name, "nmi", "", "lineage_label", na["nmi"]))
            rows.append((name, "ari", "", "lineage_label", na["ari"]))
            rows.append((name, "graph_connectivity", "", "lineage_label", sm.graph_connectivity(nn, labels)))
        for bk in BATCH_KEYS:
            batches = obs[bk].astype(str).iloc[sub].values
            with LOG.step(f"{name}: batch metrics {bk}"):
                rows.append((name, "ilisi", bk, "", sm.ilisi_knn(nn, batches)))
                rows.append((name, "kbet", bk, "lineage_label", sm.kbet_per_label(nn, batches, labels)))
                rows.append((name, "asw_batch", bk, "lineage_label", sm.silhouette_batch(X, labels, batches)))
                rows.append((name, "pcr", bk, "", sm.pcr_comparison(emb["pca"][sub], X, pd.Categorical(batches).codes.astype(np.int32), categorical=True)))
        # finer bio check
        Xf = emb[name][fine_sub]
        lf = fine_lab.iloc[fine_sub].values
        with LOG.step(f"{name}: fine bio metrics (Jessa malignant author types)"):
            nnf = pynndescent(Xf, n_neighbors=K, random_state=SEED, n_jobs=8)
            rows.append((name, "clisi", "", "jessa_author_celltype_malignant", sm.clisi_knn(nnf, lf)))
            rows.append((name, "asw_label", "", "jessa_author_celltype_malignant", sm.silhouette_label(Xf, lf)))
            na = sm.nmi_ari_cluster_labels_leiden(nnf, lf, optimize_resolution=True, n_jobs=1, seed=SEED)
            rows.append((name, "nmi", "", "jessa_author_celltype_malignant", na["nmi"]))
            rows.append((name, "ari", "", "jessa_author_celltype_malignant", na["ari"]))
        pd.DataFrame(rows, columns=["embedding", "metric", "batch_key", "label_key", "value"]).to_csv(TAB / "integration_metrics.csv", index=False)

    long = pd.DataFrame(rows, columns=["embedding", "metric", "batch_key", "label_key", "value"])
    long["value"] = long["value"].astype(float)
    long.to_csv(TAB / "integration_metrics.csv", index=False)

    # ---- aggregate ----------------------------------------------------------------------
    long["col"] = long.apply(lambda r: f"{r.metric}" + (f"__{r.batch_key}" if r.batch_key else "") + ("__fine" if r.label_key == "jessa_author_celltype_malignant" else ""), axis=1)
    wide = long.pivot(index="embedding", columns="col", values="value").loc[EMBEDDINGS]
    batch_cols = [c for c in wide.columns if c.split("__")[0] in BATCH_METRICS and not c.endswith("__fine")]
    bio_cols = [c for c in wide.columns if c.split("__")[0] in BIO_METRICS and not c.endswith("__fine")]
    fine_cols = [c for c in wide.columns if c.endswith("__fine")]
    scaled = (wide - wide.min()) / (wide.max() - wide.min()).replace(0, np.nan)
    scaled = scaled.fillna(0.5)  # a metric identical across embeddings carries no information
    summ = wide.copy()
    summ["batch_score"] = scaled[batch_cols].mean(axis=1)
    summ["bio_score"] = scaled[bio_cols].mean(axis=1)
    summ["fine_bio_score"] = scaled[fine_cols].mean(axis=1)
    summ["total_score"] = 0.4 * summ["batch_score"] + 0.6 * summ["bio_score"]
    summ["batch_score_raw_mean"] = wide[batch_cols].mean(axis=1)
    summ["bio_score_raw_mean"] = wide[bio_cols].mean(axis=1)
    # validation ELBO for tie-breaks
    elbo = {}
    for n in EMBEDDINGS:
        rj = MODELS / n / "run.json"
        if rj.exists():
            elbo[n] = json.load(open(rj)).get("elbo_validation_best", np.nan)
    summ["elbo_validation_best"] = pd.Series(elbo)
    summ["rank_total"] = summ["total_score"].rank(ascending=False, method="min").astype(int)
    scvi = summ.loc[SELECTABLE]
    chosen = scvi.sort_values(["total_score", "elbo_validation_best"], ascending=[False, True]).index[0]
    summ["chosen"] = summ.index == chosen
    summ.to_csv(TAB / "integration_metrics_summary.csv")
    json.dump({"chosen": chosen, "rule": "highest 0.4*batch+0.6*bio scaled score among scVI seeds 0/1/2 (26-epoch protocol); tie-break lower best validation ELBO; scvi_seed0_e60 is a sensitivity row only",
               "n_sub": int(len(sub)), "n_sub_fine": int(len(fine_sub)), "k": K,
               "scores": summ[["batch_score", "bio_score", "fine_bio_score", "total_score"]].round(4).to_dict(orient="index")},
              open(TAB / "integration_model_choice.json", "w"), indent=1)
    LOG.write("summary:\n" + summ[["batch_score", "bio_score", "fine_bio_score", "total_score", "rank_total"]].round(3).to_string())
    LOG.write(f"CHOSEN scVI model: {chosen}")


if __name__ == "__main__":
    main()
