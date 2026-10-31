#!/usr/bin/env python
"""Phase 9 quality assessment — leave-one-publication-out (LOPO) test of the reference classifier, run on Modal (A10G).

The existing held-out estimate (models/scanvi, 0.975 accuracy) holds out 5 PATIENTS (4 Okada) from scANVI while scVI saw
every cell, so generalisation to an unseen STUDY is unmeasured. Here, for every publication P (jessa_scrna + jessa_multiome =
jessa2022, 4 shared patients), the whole reference pipeline is refitted without P and P is then mapped as a query:
  1. scVI on the cells of the other publications, settings of models/scvi_seed1 (seed 1, n_latent 30, 2 layers, nb,
     batch sample_id, categorical study + assay, continuous pct_mito, max 26 epochs, early stopping on elbo_validation,
     patience 5, train_size 0.9, batch 1024).
  2. scANVI from that scVI with the corrected lineage_label, settings of models/scanvi (modal/phase4b_scanvi_retrain.py:
     patient-grouped inner validation over the training publications only, early stopping on validation_classification_loss,
     best-epoch weights restored, max 25 epochs).
  3. P mapped by scArches surgery (SCANVI.load_query_data, scvi-tools default freezing, weight_decay 0, 100 epochs, no
     early stopping), labels of P hidden ('Unknown'); settings of scripts/11_map_queries.py.
P's corrected labels are used only for scoring.

Run through the sandbox launcher (stages the slim training file and the corrected label table):
  python modal/phase6_cpu_run.py --gpu A10G \
    --stage data/processed/atlas_input_hvg.train.h5ad results/tables/reflabels/reference_labels_corrected.csv.gz \
            modal/phase4b_scanvi_retrain.py modal/phase9_loso_scanvi.py \
    -- "python modal/phase9_loso_scanvi.py"            [--folds okada zang ...]
Outputs (models/loso_scanvi/<publication>/, pulled back by the launcher; no model weights are kept):
  predictions.csv.gz (cell_id, sample_id, patient_id, study, assay, true, pred, maxprob), split.json, run.json
  and models/loso_scanvi/runs.jsonl (one record per fold). Scored locally by scripts/32_loso_classifier.py.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase4b_scanvi_retrain import patient_grouped_split, gpu_info, peak_rss_gb  # noqa: E402

DATA = "data/processed/atlas_input_hvg.train.h5ad"
LABELS = "results/tables/reflabels/reference_labels_corrected.csv.gz"
OUT = Path("models/loso_scanvi")
BATCH, CAT_COV, CONT_COV, UNLAB = "sample_id", ["study", "assay"], ["pct_mito"], "Unknown"
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
FOLDS = ["okada", "jessa2022", "collot", "zang", "desisto"]
ACC, EP = "gpu", dict(scvi=26, scanvi=25, query=100)   # --smoke: cpu, 2 epochs each, 20,000-cell subsample, outputs under models/loso_scanvi_smoke


def load():
    import anndata as ad
    a = ad.read_h5ad(DATA)
    lab = pd.read_csv(LABELS, index_col="cell_id", usecols=["cell_id", "lineage_label"])
    assert set(lab.index) == set(a.obs_names), "label table / training file cell mismatch"
    a.obs["lineage_label"] = lab.loc[a.obs_names, "lineage_label"].astype(str).values
    for c in CAT_COV + [BATCH, "patient_id"]:
        a.obs[c] = a.obs[c].astype(str)
    a.obs["pct_mito"] = a.obs["pct_mito"].astype(np.float32)
    a.obs["publication"] = a.obs.study.map(lambda s: PUBLICATION.get(s, s))
    return a


def cat(a):
    for c in CAT_COV + [BATCH, "patient_id", "lineage_label"]:
        a.obs[c] = a.obs[c].astype(str).astype("category")
    return a


def fold(adata, pub, seed, d):
    import scvi
    import torch
    from scvi.train import SaveCheckpoint
    scvi.settings.seed = seed
    m = (adata.obs.publication == pub).values
    ref, qry = cat(adata[~m].copy()), adata[m].copy()
    qry.obs["lineage_true"] = qry.obs["lineage_label"].values
    qry.obs["lineage_label"] = UNLAB
    qry = cat(qry)
    t0 = time.time()
    # 1. scVI without P
    scvi.model.SCVI.setup_anndata(ref, layer="counts", batch_key=BATCH, categorical_covariate_keys=CAT_COV, continuous_covariate_keys=CONT_COV)
    vae = scvi.model.SCVI(ref, n_latent=30, n_layers=2, gene_likelihood="nb")
    vae.train(max_epochs=EP["scvi"], batch_size=1024, early_stopping=True, early_stopping_monitor="elbo_validation", early_stopping_patience=5,
              early_stopping_min_delta=0.0, check_val_every_n_epoch=1, train_size=0.9, accelerator=ACC, devices=1)
    scvi_epochs = len(vae.history["elbo_train"])
    t1 = time.time()
    # 2. scANVI, patient-grouped inner validation over the training publications
    held, vm, present, tries = patient_grouped_split(ref.obs, 0.10, seed)
    lvae = scvi.model.SCANVI.from_scvi_model(vae, labels_key="lineage_label", unlabeled_category=UNLAB)
    ck = Path(f"ckpt_{pub}")
    lvae.train(max_epochs=EP["scanvi"], n_samples_per_label=None, batch_size=1024, accelerator=ACC, devices=1, check_val_every_n_epoch=1,
               datasplitter_kwargs={"external_indexing": [np.flatnonzero(~vm), np.flatnonzero(vm), np.array([], dtype=int)]},
               early_stopping=True, early_stopping_monitor="validation_classification_loss", early_stopping_mode="min",
               early_stopping_patience=5, early_stopping_min_delta=0.0,
               callbacks=[SaveCheckpoint(dirpath=str(ck), monitor="validation_classification_loss", mode="min", load_best_on_end=True)])
    vcl = lvae.history["validation_classification_loss"].iloc[:, 0].tolist()
    lab_ref = ref.obs.lineage_label.astype(str).values
    pref = lvae.predict(batch_size=4096)
    inner_val = vm & (lab_ref != UNLAB)
    t2 = time.time()
    # 3. scArches surgery of P
    ref_dir = Path(f"ref_{pub}")
    lvae.save(str(ref_dir), overwrite=True, save_anndata=False)
    scvi.model.SCANVI.prepare_query_anndata(qry, str(ref_dir))
    q = scvi.model.SCANVI.load_query_data(qry, str(ref_dir), accelerator=ACC, device="auto")
    trainable = sum(p.numel() for p in q.module.parameters() if p.requires_grad)
    q.train(max_epochs=EP["query"], batch_size=1024, accelerator=ACC, devices=1, check_val_every_n_epoch=10,
            plan_kwargs={"weight_decay": 0.0}, early_stopping=False, n_samples_per_label=None)
    soft = q.predict(soft=True, batch_size=4096)
    t3 = time.time()
    d.mkdir(parents=True, exist_ok=True)
    o = qry.obs
    pd.DataFrame({"cell_id": qry.obs_names, "sample_id": o.sample_id.astype(str).values, "patient_id": o.patient_id.astype(str).values,
                  "study": o.study.astype(str).values, "assay": o.assay.astype(str).values, "true": o.lineage_true.astype(str).values,
                  "pred": soft.idxmax(axis=1).values.astype(str), "maxprob": soft.max(axis=1).values.astype(np.float32).round(4)}
                 ).to_csv(d / "predictions.csv.gz", index=False)
    split = {"heldout_publication": pub, "n_query_cells": int(m.sum()), "n_ref_cells": int((~m).sum()),
             "query_studies": o.study.astype(str).value_counts().to_dict(), "query_patients": sorted(set(o.patient_id.astype(str))),
             "inner_validation_patients": list(map(str, held)), "inner_validation_class_counts": present, "inner_seed_tries": tries,
             "ref_classes": sorted(set(lab_ref) - {UNLAB}), "query_classes": o.lineage_true.astype(str).value_counts().to_dict()}
    json.dump(split, open(d / "split.json", "w"), indent=1, default=str)
    run = {"fold": pub, "seed": seed, "scvi_epochs": scvi_epochs, "scanvi_epochs": len(vcl), "scanvi_best_epoch": int(np.argmin(vcl)) + 1,
           "scanvi_inner_validation_accuracy": float((pref[inner_val] == lab_ref[inner_val]).mean()),
           "query_trainable_params": int(trainable), "query_epochs": EP["query"],
           "wall_scvi_s": round(t1 - t0, 1), "wall_scanvi_s": round(t2 - t1, 1), "wall_query_s": round(t3 - t2, 1),
           "gpu": gpu_info(), "peak_rss_gb": round(peak_rss_gb(), 2), "data": DATA, "labels": LABELS,
           "scvi_tools": scvi.__version__, "torch": torch.__version__}
    json.dump(run, open(d / "run.json", "w"), indent=1, default=str)
    with open(OUT / "runs.jsonl", "a") as fh:
        fh.write(json.dumps(run, default=str) + "\n")
    print("RUN", json.dumps(run, default=str), flush=True)
    shutil.rmtree(ck, ignore_errors=True); shutil.rmtree(ref_dir, ignore_errors=True)
    del vae, lvae, q
    torch.cuda.empty_cache()


def main(args):
    global OUT, ACC
    import scvi
    scvi.settings.num_threads = 8
    adata = load()
    if args.smoke:
        OUT, ACC = Path("models/loso_scanvi_smoke"), "cpu"
        EP.update(scvi=2, scanvi=2, query=2)
        adata = adata[np.random.default_rng(0).choice(adata.n_obs, 20000, replace=False)].copy()
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"loaded {adata.shape}; publications {adata.obs.publication.value_counts().to_dict()}", flush=True)
    for pub in args.folds:
        d = OUT / pub
        if (d / "predictions.csv.gz").exists() and not args.force:
            print(f"{pub}: exists, skip", flush=True); continue
        print(f"=== fold {pub}", flush=True)
        fold(adata, pub, args.seed, d)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--folds", nargs="+", default=FOLDS, choices=FOLDS)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--force", action="store_true")
    p.add_argument("--smoke", action="store_true", help="local CPU check of the code path only")
    main(p.parse_args())
