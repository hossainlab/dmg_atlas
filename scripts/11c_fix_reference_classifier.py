#!/usr/bin/env python
"""Phase 5 / step 11c — audit fix IM-01 / IM-02 (2026-10-19): re-predict query lineages with the REFERENCE classifier.

Background (docs/audit/integration_mapping_audit.md). In scvi-tools 1.4.2, `SCANVI.load_query_data(..., freeze_classifier=True)` tests
`'.0.' in key` in `_archesmixin._set_params_online_update`; the module prefix `classifier.classifier.0.` matches, so every classifier
parameter gets requires_grad=True. Gradient hooks zero the nn.Linear weights but not the BatchNorm affine parameters, which therefore
trained under the unlabelled scANVI loss for 100 epochs per query. The encoder has no covariate inputs (encode_covariates=False) and its
gradient is zeroed entirely, so the saved query latents are the zero-shot output of the reference encoder. The saved predicted_lineage
columns were produced by five unintentionally re-fitted classifiers.

What this script does
  1. Diffs models/scanvi/model.pt against every models/scanvi_query_<q>/model.pt (per-tensor max |difference|)
     -> results/tables/phase5_state_dict_diff.csv
  2. Rebuilds the reference classifier (scvi.module.Classifier, 30 -> 128 -> 128 -> 7, BatchNorm, logits) from the reference state dict and
     checks that (a) it reproduces models/scanvi/predicted_labels.txt on the reference latent, (b) the per-query classifiers rebuilt the
     same way reproduce the saved query labels (so the rebuilt module is the module scvi used).
  3. Re-encodes Filbin with the reference model through load_query_data WITHOUT training and checks that the latent equals the saved
     X_scanvi (float16 rounding only) and that its predictions equal the reference classifier applied to the saved latent.
  4. Applies the reference classifier to obsm['X_scanvi'] of every query cell in data/processed/query_mapping.h5ad and rewrites the file in
     place: predicted_lineage, scanvi_maxprob, uncertainty_prob, classifier_knn_agree and obsm['scanvi_probs'] are replaced; the previous
     values are kept as predicted_lineage_pre_audit, scanvi_maxprob_pre_audit, uncertainty_prob_pre_audit, classifier_knn_agree_pre_audit
     and obsm['scanvi_probs_pre_audit']. Adds obs['liu_prep'] (fresh cell / frozen nucleus / unannotated; IM-24).
  5. Recomputes results/tables/mapping_uncertainty.csv from the new obs (kNN columns unchanged); previous file kept as *.pre_audit.csv.
  6. Writes results/qc/phase5_reference_classifier_verification.json with every check above.
Downstream: run scripts/11b_mapping_summaries.py and scripts/12_package_v1.py summary figures afterwards.
Log: logs/fix_phase5.log
"""
from __future__ import annotations

import os
import sys
import json
import time
import shutil
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402  (thread caps first)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import h5py  # noqa: E402
import torch  # noqa: E402
import anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_phase5.log", context="11c_fix_reference_classifier")
REF_DIR = ROOT / "models/scanvi"
QM = ROOT / "data/processed/query_mapping.h5ad"
QUERIES = ["filbin", "labelle", "liu", "braun", "siletti"]
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "8")))


def tensors(state: dict) -> dict:
    return {k: v for k, v in state.items() if torch.is_tensor(v)}


def build_classifier(state: dict):
    from scvi.module import Classifier
    clf = Classifier(30, n_hidden=128, n_labels=7, n_layers=2, dropout_rate=0.1, logits=True, use_batch_norm=True, use_layer_norm=False)
    clf.load_state_dict({k[len("classifier."):]: v for k, v in state.items() if k.startswith("classifier.")}, strict=True)
    clf.eval()
    return clf


@torch.no_grad()
def predict_probs(clf, Z: np.ndarray, bs: int = 65536) -> np.ndarray:
    return np.vstack([torch.softmax(clf(torch.as_tensor(Z[i:i + bs], dtype=torch.float32)), -1).numpy() for i in range(0, len(Z), bs)])


def main():
    import scvi
    ver = {"scvi_tools": str(scvi.__version__), "torch": str(torch.__version__), "anndata": str(ad.__version__)}
    ref = torch.load(REF_DIR / "model.pt", map_location="cpu", weights_only=False)
    sd = tensors(ref["model_state_dict"])
    cols = (REF_DIR / "soft_prob_columns.txt").read_text().split("\n")
    verification = {"versions": ver, "init_params_reference": ref["attr_dict"]["init_params_"]["non_kwargs"]}

    # ---------------------------------------------------------------------------------------------- 1. state-dict diff
    with LOG.step("state-dict diff reference vs five query models"):
        rows = []
        for q in QUERIES:
            qsd = tensors(torch.load(ROOT / f"models/scanvi_query_{q}/model.pt", map_location="cpu", weights_only=False)["model_state_dict"])
            for k, v in qsd.items():
                same = k in sd and sd[k].shape == v.shape
                rows.append({"query": q, "tensor": k, "block": k.split(".")[0], "shape_ref": tuple(sd[k].shape) if k in sd else None,
                             "shape_query": tuple(v.shape), "same_shape": same,
                             "max_abs_diff": float((v.float() - sd[k].float()).abs().max()) if same and v.numel() else np.nan})
        diff = pd.DataFrame(rows)
        diff.to_csv(ROOT / "results/tables/phase5_state_dict_diff.csv", index=False)
        block_max = diff[diff.same_shape].groupby(["query", "block"]).max_abs_diff.max().unstack().round(4)
        LOG.write("max |diff| per block:\n" + block_max.to_string())
        verification["state_dict_block_max_abs_diff"] = {q: block_max.loc[q].to_dict() for q in QUERIES}
        verification["classifier_tensors_changed"] = sorted(diff[(diff.block == "classifier") & (diff.max_abs_diff > 0)].tensor.unique().tolist())
        assert (block_max["z_encoder"] == 0).all(), "encoder differs from reference — IM-02 premise not met"

    # ---------------------------------------------------------------------------------------------- 2. rebuild classifier
    with LOG.step("rebuild reference classifier and validate on reference / per-query models"):
        ref_clf = build_classifier(sd)
        Zref = np.load(REF_DIR / "latent_f16.npy").astype(np.float32)
        saved_ref = np.array((REF_DIR / "predicted_labels.txt").read_text().split("\n"))
        Pref = predict_probs(ref_clf, Zref)
        ref_agree = float((np.array(cols)[Pref.argmax(1)] == saved_ref).mean())
        ref_probdiff = float(np.abs(Pref - np.load(REF_DIR / "soft_probs_f16.npy").astype(np.float32)).max())
        LOG.write(f"reference classifier on reference latent: label agreement with saved predicted_labels {ref_agree:.6f}, max |prob diff| {ref_probdiff:.4f}")
        verification["reference_classifier_reproduces_reference_predictions"] = {"label_agreement": ref_agree, "max_abs_prob_diff": ref_probdiff}
        del Zref, Pref

    with h5py.File(QM) as f:
        obs = read_elem(f["obs"])
        Zq = f["obsm/X_scanvi"][:]
        qcols = list(f["uns/phase5/scanvi_prob_columns"][:].astype(str))
    assert qcols == cols
    for c in obs.columns:
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            obs[c] = obs[c].astype(object)
    already = "predicted_lineage_pre_audit" in obs.columns
    old_lab_col = "predicted_lineage_pre_audit" if already else "predicted_lineage"

    new_probs = np.zeros((len(obs), 7), dtype=np.float32)
    verification["per_query"] = {}
    with LOG.step("apply reference classifier to saved query latents"):
        for q in QUERIES:
            m = (obs["query"] == q).values
            qsd = tensors(torch.load(ROOT / f"models/scanvi_query_{q}/model.pt", map_location="cpu", weights_only=False)["model_state_dict"])
            Pq_old = predict_probs(build_classifier(qsd), Zq[m])
            repro = float((np.array(cols)[Pq_old.argmax(1)] == obs[old_lab_col].values[m]).mean())
            Pq_new = predict_probs(ref_clf, Zq[m])
            new_probs[m] = Pq_new
            new_counts = pd.Series(np.array(cols)[Pq_new.argmax(1)]).value_counts().to_dict()
            old_counts = pd.Series(obs[old_lab_col].values[m]).value_counts().to_dict()
            LOG.write(f"{q}: n {m.sum():,}; per-query (drifted) classifier reproduces saved labels {repro:.4f}; "
                      f"reference-classifier counts {new_counts}; previous {old_counts}")
            verification["per_query"][q] = {"n": int(m.sum()), "drifted_classifier_reproduces_saved_labels": repro,
                                            "reference_classifier_counts": {k: int(v) for k, v in new_counts.items()},
                                            "pre_audit_counts": {k: int(v) for k, v in old_counts.items()}}

    # ---------------------------------------------------------------------------------------------- 3. zero-shot check (Filbin)
    with LOG.step("zero-shot check: load_query_data(filbin) without training"):
        from scvi.model import SCANVI
        scvi.settings.seed = 1
        scvi.settings.dl_num_workers = 0
        fa = ad.read_h5ad(ROOT / "data/processed/queries/filbin.hvg.h5ad")
        SCANVI.prepare_query_anndata(fa, str(REF_DIR))
        qmodel = SCANVI.load_query_data(fa, str(REF_DIR), accelerator="cpu", device="auto")
        grad_true = [n for n, p in qmodel.module.named_parameters() if n.startswith("classifier") and p.requires_grad]
        msd = qmodel.module.state_dict()
        same = all(torch.equal(v, msd[k]) for k, v in sd.items() if k in msd and v.shape == msd[k].shape)
        qmodel.is_trained_ = True  # no optimiser step is run; the module holds the reference parameters (checked by `same`)
        qmodel.module.eval()
        Z0 = qmodel.get_latent_representation()
        P0 = qmodel.predict(soft=True)[cols].values
        mf = (obs["query"] == "filbin").values
        zs = {"classifier_params_requires_grad_after_load_query_data": grad_true,
              "untrained_query_model_equals_reference_on_shared_tensors": bool(same),
              "max_abs_latent_diff_vs_saved_X_scanvi": float(np.abs(Z0 - Zq[mf]).max()),
              "pearson_r_latent": float(np.corrcoef(Z0.ravel(), Zq[mf].ravel())[0, 1]),
              "label_agreement_with_reference_classifier_on_saved_latent": float((np.array(cols)[P0.argmax(1)] == np.array(cols)[new_probs[mf].argmax(1)]).mean()),
              "max_abs_prob_diff": float(np.abs(P0 - new_probs[mf]).max())}
        LOG.write(f"zero-shot Filbin: {json.dumps(zs)}")
        verification["zero_shot_check_filbin"] = zs
        del qmodel, fa

    # ---------------------------------------------------------------------------------------------- 4. rewrite query_mapping.h5ad
    with LOG.step("rewrite data/processed/query_mapping.h5ad with reference-classifier predictions (old columns kept as *_pre_audit)"):
        A = ad.read_h5ad(QM)
        if not already:
            for c in ["predicted_lineage", "scanvi_maxprob", "uncertainty_prob", "classifier_knn_agree"]:
                A.obs[f"{c}_pre_audit"] = A.obs[c].values
            A.obsm["scanvi_probs_pre_audit"] = A.obsm["scanvi_probs"]
        A.obsm["scanvi_probs"] = new_probs
        A.obs["predicted_lineage"] = pd.Categorical(np.array(cols)[new_probs.argmax(1)])
        A.obs["scanvi_maxprob"] = new_probs.max(1).astype(np.float32)
        A.obs["uncertainty_prob"] = (1.0 - A.obs["scanvi_maxprob"]).astype(np.float32)
        A.obs["classifier_knn_agree"] = A.obs["knn_lineage"].astype(str).values == A.obs["predicted_lineage"].astype(str).values
        prep = A.obs["author_prep"].astype(object).map({"fresh_cells": "fresh cell", "frozen_nuclei": "frozen nucleus"}).fillna("unannotated")
        A.obs["liu_prep"] = pd.Categorical(np.where(A.obs["query"].astype(str) == "liu", prep.astype(str), "n/a"))
        A.uns["phase5"]["audit_fix_2026_10_19"] = {
            "finding": "IM-01/IM-02 (docs/audit/integration_mapping_audit.md)",
            "predicted_lineage": "argmax of the REFERENCE scANVI classifier (models/scanvi) applied to obsm['X_scanvi']; the *_pre_audit columns hold the "
                                 "values from the per-query models/scanvi_query_<q> whose classifier BatchNorm parameters drifted during 'surgery'",
            "X_scanvi": "zero-shot output of the reference encoder (encode_covariates=False; encoder state identical before/after query training)",
            "scanvi_prob_columns": cols, "script": "scripts/11c_fix_reference_classifier.py", "versions": ver}
        A.uns["phase5"]["knn_lineage"] = ("majority scANVI label (atlas_core scanvi_label) among the 15 nearest reference cells in X_scanvi; "
                                         "primary label for off-manifold queries (see docs/phase5_mapping_notes.md)")
        A.write_h5ad(QM, compression="gzip")
        LOG.write(f"wrote {QM}: {A.n_obs} cells; obs {A.obs.shape[1]} cols; obsm {list(A.obsm)} | peak RSS {peak_rss_gb():.2f} GB")
        obs_new = A.obs.copy()
        del A

    # ---------------------------------------------------------------------------------------------- 5. mapping_uncertainty.csv
    with LOG.step("recompute results/tables/mapping_uncertainty.csv"):
        out = ROOT / "results/tables/mapping_uncertainty.csv"
        pre = out.with_name("mapping_uncertainty.pre_audit.csv")
        if out.exists() and not pre.exists():
            shutil.copy2(out, pre)
        rows = []
        for (q, s), g in obs_new.groupby(["query", "sample_id"], observed=True):
            rows.append({"query": q, "sample_id": s, "n_cells": len(g),
                         "median_uncertainty_prob": g.uncertainty_prob.median(), "p90_uncertainty_prob": g.uncertainty_prob.quantile(0.9),
                         "frac_uncertainty_prob_gt_0.5": (g.uncertainty_prob > 0.5).mean(),
                         "median_uncertainty_knn": g.uncertainty_knn.median(), "p90_uncertainty_knn": g.uncertainty_knn.quantile(0.9),
                         "frac_uncertainty_knn_gt_0.5": (g.uncertainty_knn > 0.5).mean(),
                         "median_knn_mean_dist": g.knn_mean_dist.median(), "median_maxprob": g.scanvi_maxprob.median(),
                         "frac_maxprob_ge_0.8": (g.scanvi_maxprob >= 0.8).mean(),
                         "frac_classifier_knn_agree": g.classifier_knn_agree.astype(float).mean()})
        unc = pd.DataFrame(rows)
        unc = unc.set_index("query").loc[QUERIES].reset_index()
        unc.to_csv(out, index=False)
        LOG.write(f"wrote {out} ({len(unc)} sample rows)")

    (ROOT / "results/qc").mkdir(exist_ok=True)
    (ROOT / "results/qc/phase5_reference_classifier_verification.json").write_text(json.dumps(verification, indent=1, default=str))
    LOG.write(f"done | peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
