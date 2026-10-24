#!/usr/bin/env python
"""Phase 7 / step 26b — audit fix P7-01: recompute the embedding comparison on a common cell set with finite latents.

Why: step 26 zeroed NaN rows of the scVI latent (4,281 of 15,000 cells that the malignant scVI model never embedded), so every
scVI number in foundation_model_comparison.csv was wrong (audit P7-01). The scGPT / Geneformer / PCA embeddings in
data/processed/foundation_embeddings.npz were computed on a PRE-AUDIT cell sample and are NOT recomputed here (no GPU re-run):
the comparison is re-intersected instead.

Inputs : data/processed/foundation_embeddings.npz (step 26: cell_ids, pre-audit family, PCA, scVI(pre-audit), scGPT, Geneformer),
         data/processed/atlas_core.h5ad (refreshed; obs sample_id / assay / consensus_call, obsm['X_scvi_malignant'] = retrained
         malignant scVI, 167,395 cells), results/tables/metaprogram_family_reference_classes.csv (step 19, corrected families
         with within-sample z and the confident flag).
Evaluation set: npz cells that (i) are in the corrected malignant selection with a CONFIDENT corrected family (max within-sample
         family z >= 1; all 7 families incl. novel_MP9 cilia), and (ii) have finite values in ALL four embeddings. Assertion: no NaN row.
Outputs: results/tables/foundation_model_comparison.csv (previous -> .pre_audit.csv), results/tables/foundation_eval_cells.csv
         (cell_id, corrected family, pre-audit family, sample, assay class, in_eval), data/processed/foundation_embeddings.npz
         rewritten with scVI = retrained latent (NaN where not covered), scVI_pre_audit, family (corrected), family_pre_audit,
         eval_mask (previous -> foundation_embeddings.pre_audit.npz), results/figures/phase7_foundation_models.png.
Metrics (identical to step 26, on L2-normalised embeddings): kNN (k = 15, cosine) balanced accuracy, stratified 5-fold (seed 0);
         nuclei -> cells transfer (kNN trained on nuclei, scored on scRNA cells, balanced accuracy); silhouette (cosine, 5,000-cell
         subsample, seed 0); mean fraction of the 15 nearest neighbours from the same sample / same assay class, with the expected
         values under label frequencies.
Run: python scripts/26b_foundation_comparison_rebuild.py
"""
import sys, shutil
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors  # noqa: E402
from sklearn.metrics import balanced_accuracy_score, silhouette_score  # noqa: E402
from sklearn.preprocessing import normalize  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CORE = PROCESSED / "atlas_core.h5ad"; TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"
NPZ = PROCESSED / "foundation_embeddings.npz"; OUT = TABLES / "foundation_model_comparison.csv"
log = U.StepLogger(LOGS / "fix_phase7.log", context="26b_foundation_comparison_rebuild")
K, SEED, MIN_Z = 15, 0, 1.0
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}
EMB_ORDER = ["PCA", "scVI", "scGPT", "Geneformer"]

with log.step("intersect npz cells with corrected families and the retrained malignant latent"):
    src = NPZ.with_name("foundation_embeddings.pre_audit.npz")
    if not src.exists():
        shutil.copy2(NPZ, src)
    z = np.load(src, allow_pickle=True)
    cells = pd.Index(z["cell_ids"].astype(str)); fam_old = z["family"].astype(str)
    emb = {k: np.asarray(z[k], dtype=np.float32) for k in EMB_ORDER}
    scvi_pre = emb.pop("scVI")
    with h5py.File(CORE, "r") as f:
        idx_all = pd.Index(read_elem(f["obs"]["_index"]))
        rows = idx_all.get_indexer(cells); assert (rows >= 0).all(), "npz cells missing from atlas_core"
        order = np.argsort(rows)
        Xs = np.full((len(cells), f["obsm"]["X_scvi_malignant"].shape[1]), np.nan, np.float32)
        Xs[order] = f["obsm"]["X_scvi_malignant"][rows[order]]
        obs = pd.DataFrame({c: np.asarray(read_elem(f["obs"][c])).astype(str)[rows] for c in ("sample_id", "assay", "consensus_call")}, index=cells)
    emb["scVI"] = Xs
    ref = pd.read_csv(TABLES / "metaprogram_family_reference_classes.csv").set_index("cell_id")
    new = ref.reindex(cells)
    in_sel = new.family.notna().values; conf = new.confident.fillna(False).astype(bool).values & (new.max_family_z.fillna(-9).values >= MIN_Z)
    finite = np.ones(len(cells), bool)
    for k, E in emb.items():
        finite &= np.isfinite(E).all(1)
    ev = in_sel & conf & finite
    tab = pd.DataFrame({"cell_id": cells, "family": new.family.fillna("not_in_selection").values, "family_pre_audit": fam_old,
                        "max_family_z": new.max_family_z.values, "confident": conf, "consensus_call": obs.consensus_call.values,
                        "sample_id": obs.sample_id.values, "assay_class": obs.assay.map(ASSAY_CLASS).values,
                        "finite_all_embeddings": finite, "in_eval": ev})
    tab.to_csv(TABLES / "foundation_eval_cells.csv", index=False)
    log.write(f"npz cells {len(cells):,}; in corrected selection {in_sel.sum():,}; confident corrected family {(in_sel & conf).sum():,}; "
              f"finite in all four embeddings {finite.sum():,} (retrained scVI finite {np.isfinite(Xs).all(1).sum():,}; pre-audit scVI finite "
              f"{np.isfinite(scvi_pre).all(1).sum():,}); EVALUATION SET {ev.sum():,}")
    log.write("corrected families in the evaluation set: " + tab[ev].family.value_counts().to_dict().__repr__()
              + "; assay " + tab[ev].assay_class.value_counts().to_dict().__repr__())
    log.write("pre-audit family x corrected family (evaluation set):\n" + pd.crosstab(tab[ev].family_pre_audit, tab[ev].family).to_string())

with log.step("metrics on the common evaluation set"):
    lab = tab.family.values[ev]; samp = tab.sample_id.values[ev]; assay = tab.assay_class.values[ev]
    rng = np.random.default_rng(SEED); skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    cells_mask, nuc_mask = assay == "cells", assay == "nuclei"
    sub = rng.choice(ev.sum(), size=min(5000, ev.sum()), replace=False)
    exp_s = np.mean([(samp == s).mean() for s in samp]); exp_a = np.mean([(assay == x).mean() for x in assay])
    rows_m = []
    for name in EMB_ORDER:
        E = emb[name][ev]
        assert np.isfinite(E).all(), f"{name}: NaN rows in the evaluation set"
        Z = normalize(E)
        pred = np.empty(len(lab), dtype=object)
        for tr, te in skf.split(Z, lab):
            pred[te] = KNeighborsClassifier(K, metric="cosine").fit(Z[tr], lab[tr]).predict(Z[te])
        xa = balanced_accuracy_score(lab[cells_mask], KNeighborsClassifier(K, metric="cosine").fit(Z[nuc_mask], lab[nuc_mask]).predict(Z[cells_mask]))
        nn = NearestNeighbors(n_neighbors=K + 1, metric="cosine").fit(Z).kneighbors(Z, return_distance=False)[:, 1:]
        rows_m.append(dict(embedding=name, status="ok (re-intersected, audit P7-01)" if name != "scVI" else "ok (retrained malignant scVI, audit P7-01)",
                           dims=Z.shape[1], n_cells=len(Z), n_families=len(set(lab)), knn_balanced_accuracy=balanced_accuracy_score(lab, pred),
                           cross_assay_accuracy_nuclei_to_cells=xa, n_cells_scrna=int(cells_mask.sum()), n_cells_nuclei=int(nuc_mask.sum()),
                           silhouette=silhouette_score(Z[sub], lab[sub], metric="cosine"),
                           same_sample_neighbours=(samp[nn] == samp[:, None]).mean(), same_sample_expected=exp_s,
                           same_assay_neighbours=(assay[nn] == assay[:, None]).mean(), same_assay_expected=exp_a))
    res = pd.DataFrame(rows_m)
    log.write("\n" + res.drop(columns=["status"]).round(3).to_string(index=False))

with log.step("write table, npz and figure"):
    if OUT.exists() and not OUT.with_suffix(".pre_audit.csv").exists():
        shutil.copy2(OUT, OUT.with_suffix(".pre_audit.csv"))
    res.to_csv(OUT, index=False)
    np.savez_compressed(NPZ, cell_ids=np.asarray(cells, dtype=str), family=tab.family.values.astype(str), family_pre_audit=fam_old,
                        eval_mask=ev, PCA=emb["PCA"], scVI=emb["scVI"], scVI_pre_audit=scvi_pre, scGPT=emb["scGPT"], Geneformer=emb["Geneformer"])
    apply_figure_style()
    mets = [("knn_balanced_accuracy", "kNN family balanced accuracy", None), ("cross_assay_accuracy_nuclei_to_cells", "nuclei → cells transfer accuracy", None),
            ("same_sample_neighbours", "same-sample neighbour fraction", "same_sample_expected"), ("same_assay_neighbours", "same-assay neighbour fraction", "same_assay_expected")]
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.1), layout="constrained")
    for ax, (m_, ylab, expc) in zip(axes, mets):
        ax.bar(range(len(res)), res[m_], color="#4d4d4d", width=0.7)
        ax.set_xticks(range(len(res))); ax.set_xticklabels(res.embedding, rotation=45, ha="right"); ax.set_ylabel(ylab)
        if expc:
            ax.axhline(res[expc].iloc[0], color="grey", lw=0.6, ls=":")
        ax.set_ylim(0, 1.0 if m_ != "same_sample_neighbours" else max(0.05, res[m_].max() * 1.15))
    fig.savefig(FIGS / "phase7_foundation_models.png"); plt.close(fig)
    log.write(f"wrote {OUT.name} (pre-audit copy {OUT.with_suffix('.pre_audit.csv').name}), {NPZ.name} (pre-audit copy {src.name}), phase7_foundation_models.png")
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
