#!/usr/bin/env python
"""Phase 5 / step 11b — label-transfer summaries, specificity check and cell-of-origin view.

Reads data/processed/query_mapping.h5ad (obs only), models/scanvi latents and atlas_core obs.
Writes
  results/tables/mapping_summary.csv                      per query sample: predicted lineage distribution + quality
  results/tables/mapping_<query>_author_vs_predicted.csv  crosstabs (Filbin/LaBelle use Phase-3 consensus calls: no author states)
  results/tables/mapping_liu_by_group.csv                 Liu quality by prep (fresh cells / frozen nuclei), age group, site
  results/tables/mapping_normal_reference_specificity.csv Braun / Siletti: predicted lineage per normal class, false-malignant fraction
  results/tables/mapping_author_agreement.json            headline agreement numbers
  results/tables/malignant_nearest_fetal_class.csv        kNN (k=15) transfer of Braun CellClass / cell_type onto consensus-malignant
                                                          reference cells, summarised per Leiden cluster and per Jessa author state
  results/tables/malignant_nearest_fetal_class_per_cell.csv.gz
"""
from __future__ import annotations

import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from anndata.io import read_elem  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase5_mapping.log", context="11b_mapping_summaries")
T = ROOT / "results/tables"


def backup_pre_audit(path: Path) -> None:
    """Audit rule (2026-10-19): when a deliverable is regenerated, keep the previous version as <name>.pre_audit.<ext> (never overwritten)."""
    import shutil
    if path.exists():
        pre = path.with_name(path.name.replace(".csv.gz", ".pre_audit.csv.gz") if path.name.endswith(".csv.gz") else f"{path.stem}.pre_audit{path.suffix}")
        if not pre.exists():
            shutil.copy2(path, pre)


def sample_level(flag: np.ndarray, sample: np.ndarray, n_boot: int = 1000, seed: int = 0) -> dict:
    """Pooled cell-level rate plus per-sample median / IQR and a 95 % CI from bootstrapping samples (IM-22: cells within a sample are not independent)."""
    flag = np.asarray(flag, dtype=float); sample = np.asarray(sample).astype(str)
    ok = ~np.isnan(flag); flag, sample = flag[ok], sample[ok]
    if len(flag) == 0:
        return {"n_cells": 0}
    per = pd.Series(flag).groupby(sample).agg(["mean", "size"])
    rng = np.random.default_rng(seed)
    sums = per["mean"].values * per["size"].values; sizes = per["size"].values
    boots = []
    for _ in range(n_boot):
        j = rng.integers(0, len(per), len(per))
        boots.append(sums[j].sum() / sizes[j].sum())
    return {"n_cells": int(len(flag)), "n_samples": int(len(per)), "pooled": float(flag.mean()),
            "per_sample_median": float(per["mean"].median()), "per_sample_q25": float(per["mean"].quantile(0.25)), "per_sample_q75": float(per["mean"].quantile(0.75)),
            "boot_ci_low": float(np.percentile(boots, 2.5)), "boot_ci_high": float(np.percentile(boots, 97.5))}
LINEAGES = ["malignant", "myeloid", "oligodendrocyte", "astrocyte", "neuron", "vascular", "T_NK"]

LIU_MALIGNANT = ["OPC-like-1", "OPC-like-2", "OPC-like-3", "AC-like", "MES-like", "OC-like", "S", "G2M", "frozen-OPC-like-b"]
LIU_NORMAL = {"Myeloid": "myeloid", "Microglia": "myeloid", "Oligodendrocyte": "oligodendrocyte", "Endothelia": "vascular", "T cell": "T_NK"}
LABELLE_NORMAL = {"Myeloid": "myeloid", "Tcell": "T_NK"}
BRAUN_EXPECTED = {"Oligo": "oligodendrocyte", "Vascular": "vascular", "Fibroblast": "vascular", "Immune": "myeloid",
                  "Neuron": "neuron", "Neuroblast": "neuron", "Neuronal IPC": "neuron",
                  "Radial glia": None, "Glioblast": None, "Erythrocyte": None, "Placodes": None, "Neural crest": None}
SILETTI_EXPECTED = {"Oligodendrocyte": "oligodendrocyte", "Astrocyte": "astrocyte", "Bergmann glia": "astrocyte", "Ependymal": "astrocyte",
                    "Microglia": "myeloid", "Vascular": "vascular", "Fibroblast": "vascular", "Oligodendrocyte precursor": None,
                    "Committed oligodendrocyte precursor": "oligodendrocyte", "Choroid plexus": None, "Splatter": None}
# every other Siletti supercluster is neuronal


def read_obs():
    with h5py.File(ROOT / "data/processed/query_mapping.h5ad") as f:
        obs = read_elem(f["obs"])
        probs = f["obsm/scanvi_probs"][:]
        cols = list(f["uns/phase5/scanvi_prob_columns"][:].astype(str))
    for c in obs.columns:
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            obs[c] = obs[c].astype(object)
    return obs, pd.DataFrame(probs, index=obs.index, columns=cols)


def lineage_dist(g: pd.DataFrame) -> dict:
    vc = g.predicted_lineage.value_counts()
    d = {f"n_pred_{l}": int(vc.get(l, 0)) for l in LINEAGES}
    d.update({f"frac_pred_{l}": round(vc.get(l, 0) / len(g), 4) for l in LINEAGES})
    d.update({"median_maxprob": g.scanvi_maxprob.median(), "frac_maxprob_ge_0.8": (g.scanvi_maxprob >= 0.8).mean(),
              "median_uncertainty_knn": g.uncertainty_knn.median(), "frac_uncertainty_knn_gt_0.5": (g.uncertainty_knn > 0.5).mean(),
              "frac_pred_malignant_confident": ((g.predicted_lineage == "malignant") & (g.scanvi_maxprob >= 0.8)).mean(),
              "frac_classifier_knn_agree": g.classifier_knn_agree.astype(float).mean()})
    vk = g.knn_lineage.value_counts()
    d.update({f"frac_knn_{l}": round(vk.get(l, 0) / len(g), 4) for l in LINEAGES})
    return d


def main():
    obs, probs = read_obs()
    agreement = {}
    # ------------------------------------------------------------------ 1. per-sample summary
    rows = []
    for (q, s), g in obs.groupby(["query", "sample_id"]):
        meta = {k: g[k].dropna().iloc[0] if g[k].notna().any() else "" for k in ["study", "assay", "site", "age_group", "patient_id", "donor_id"]}
        rows.append({"query": q, "sample_id": s, **meta, "n_cells": len(g), **lineage_dist(g)})
    summ = pd.DataFrame(rows)
    backup_pre_audit(T / "mapping_summary.csv")
    summ.to_csv(T / "mapping_summary.csv", index=False)
    LOG.write(f"mapping_summary.csv: {len(summ)} sample rows")

    # ------------------------------------------------------------------ 2. crosstabs
    def xtab(df, col, name):
        ct = pd.crosstab(df[col].fillna("(none)"), df.predicted_lineage)
        ct["n"] = ct.sum(1)
        ck = pd.crosstab(df[col].fillna("(none)"), df.knn_lineage)
        ck.columns = [f"knn_{c}" for c in ck.columns]
        ct = ct.join(ck)
        backup_pre_audit(T / f"mapping_{name}_author_vs_predicted.csv")
        ct.to_csv(T / f"mapping_{name}_author_vs_predicted.csv")
        return ct

    liu = obs[obs["query"] == "liu"]
    ct_liu = xtab(liu, "author_annotation", "liu")
    lm = liu[liu.author_annotation.isin(LIU_MALIGNANT)]
    ln = liu[liu.author_annotation.isin(LIU_NORMAL)]
    agreement["liu"] = {
        "n_author_malignant": int(len(lm)), "frac_author_malignant_pred_malignant": float((lm.predicted_lineage == "malignant").mean()),
        "frac_author_malignant_knn_malignant": float((lm.knn_lineage == "malignant").mean()),
        "frac_author_normal_knn_matching": float((ln.knn_lineage.values == ln.author_annotation.map(LIU_NORMAL).values).mean()),
        "per_type_frac_knn_matching": {k: float((v.knn_lineage == LIU_NORMAL[k]).mean()) for k, v in ln.groupby("author_annotation")},
        "per_state_frac_pred_malignant": {k: float((v.predicted_lineage == "malignant").mean()) for k, v in lm.groupby("author_annotation")},
        "n_author_normal": int(len(ln)),
        "frac_author_normal_pred_matching": float((ln.predicted_lineage.values == ln.author_annotation.map(LIU_NORMAL).values).mean()),
        "per_type_frac_matching": {k: float((v.predicted_lineage == LIU_NORMAL[k]).mean()) for k, v in ln.groupby("author_annotation")},
        "frac_author_normal_pred_malignant": float((ln.predicted_lineage == "malignant").mean()),
        "n_unannotated": int(liu.author_annotation.isna().sum()),
        "unannotated_frac_pred_malignant": float((liu[liu.author_annotation.isna()].predicted_lineage == "malignant").mean())}

    fil = obs[obs["query"] == "filbin"]
    fil = fil.assign(phase3_call=np.where(fil.consensus_call == "malignant", "consensus malignant",
                                          np.where(fil.consensus_call == "normal", "consensus normal: " + fil.consensus_celltype_coarse.astype(str),
                                                   "consensus " + fil.consensus_call.astype(str))))
    xtab(fil, "phase3_call", "filbin")
    fm = fil[fil.consensus_call == "malignant"]; fn = fil[fil.consensus_call == "normal"]
    agreement["filbin"] = {
        "note": "Filbin deposit has no per-cell author states; Phase-3 consensus_call (CNV + marker) used instead",
        "n_consensus_malignant": int(len(fm)), "frac_consensus_malignant_pred_malignant": float((fm.predicted_lineage == "malignant").mean()),
        "frac_consensus_malignant_knn_malignant": float((fm.knn_lineage == "malignant").mean()),
        "n_consensus_normal": int(len(fn)),
        "frac_consensus_normal_pred_matching": float((fn.predicted_lineage.values == fn.consensus_celltype_coarse.values).mean()) if len(fn) else None,
        "per_type_frac_matching": {k: float((v.predicted_lineage == k).mean()) for k, v in fn.groupby("consensus_celltype_coarse")},
        "n_ambiguous": int((fil.consensus_call == "ambiguous").sum()),
        "ambiguous_frac_pred_malignant": float((fil[fil.consensus_call == "ambiguous"].predicted_lineage == "malignant").mean()),
        "cnv_call_malignant_frac_pred_malignant": float((fil[fil.cnv_call == "malignant"].predicted_lineage == "malignant").mean())}

    lab = obs[obs["query"] == "labelle"]
    lab = lab.assign(author_or_table=np.where(lab.author_broad_annot.notna(), "immune: " + lab.author_broad_annot.astype(str),
                                              "tumour SS2 table, consensus " + lab.consensus_call.astype(str)))
    xtab(lab, "author_or_table", "labelle")
    li = lab[lab.author_broad_annot.notna()]
    lt = lab[(lab.labelle_table == "tumor_ss2") & (lab.consensus_call == "malignant")]
    agreement["labelle"] = {
        "n_author_immune": int(len(li)),
        "frac_author_immune_pred_matching": float((li.predicted_lineage.values == li.author_broad_annot.map(LABELLE_NORMAL).values).mean()),
        "frac_author_immune_knn_matching": float((li.knn_lineage.values == li.author_broad_annot.map(LABELLE_NORMAL).values).mean()),
        "per_type_frac_knn_matching": {k: float((v.knn_lineage == LABELLE_NORMAL[k]).mean()) for k, v in li.groupby("author_broad_annot")},
        "per_type_frac_matching": {k: float((v.predicted_lineage == LABELLE_NORMAL[k]).mean()) for k, v in li.groupby("author_broad_annot")},
        "frac_author_immune_pred_myeloid_or_T_NK": float(li.predicted_lineage.isin(["myeloid", "T_NK"]).mean()),
        "n_tumour_table_consensus_malignant": int(len(lt)), "frac_tumour_table_consensus_malignant_pred_malignant": float((lt.predicted_lineage == "malignant").mean())}

    # ------------------------------------------------------------------ 3. Liu by group
    grp_rows = []
    liu = liu.assign(prep=liu.author_prep.fillna("unannotated"), age=liu.author_age.fillna(liu.age_group), loc=liu.author_location.fillna(liu.site))
    for axis in ["prep", "age", "loc"]:
        for k, g in liu.groupby(axis):
            gm = g[g.author_annotation.isin(LIU_MALIGNANT)]; gn = g[g.author_annotation.isin(LIU_NORMAL)]
            grp_rows.append({"axis": axis, "group": k, "n_cells": len(g), "n_samples": g.sample_id.nunique(), **lineage_dist(g),
                             "n_author_malignant": len(gm), "frac_author_malignant_pred_malignant": (gm.predicted_lineage == "malignant").mean() if len(gm) else np.nan,
                             "n_author_normal": len(gn),
                             "frac_author_normal_pred_matching": (gn.predicted_lineage.values == gn.author_annotation.map(LIU_NORMAL).values).mean() if len(gn) else np.nan})
    backup_pre_audit(T / "mapping_liu_by_group.csv")
    pd.DataFrame(grp_rows).to_csv(T / "mapping_liu_by_group.csv", index=False)

    # ------------------------------------------------------------------ 4. normal-reference specificity
    spec_rows = []
    fm_summary = {}
    for q, col, expected in [("braun", "author_CellClass", BRAUN_EXPECTED), ("siletti", "author_supercluster_term", SILETTI_EXPECTED)]:
        g0 = obs[obs["query"] == q]
        for k, g in g0.groupby(col):
            exp = expected.get(k, "neuron" if q == "siletti" else None)
            spec_rows.append({"query": q, "normal_class": k, "expected_lineage": exp if exp else "(none in reference)", "n_cells": len(g),
                              **lineage_dist(g),
                              "frac_pred_expected": (g.predicted_lineage == exp).mean() if exp else np.nan,
                              "frac_knn_expected": (g.knn_lineage == exp).mean() if exp else np.nan,
                              "frac_false_malignant_maxprob_ge_0.8": ((g.predicted_lineage == "malignant") & (g.scanvi_maxprob >= 0.8)).mean(),
                              "frac_knn_malignant": (g.knn_lineage == "malignant").mean(),
                              "frac_knn_malignant_vote_ge_0.8": ((g.knn_lineage == "malignant") & (g.knn_vote_frac >= 0.8)).mean()})
        is_diff = g0[col].map(lambda k: expected.get(k, "neuron" if q == "siletti" else None)).notna()
        gd = g0[is_diff]
        fm_summary[q] = {"n_cells": int(len(g0)),
                         "frac_false_malignant_maxprob_ge_0.8_all": float(((g0.predicted_lineage == "malignant") & (g0.scanvi_maxprob >= 0.8)).mean()),
                         "frac_pred_malignant_all": float((g0.predicted_lineage == "malignant").mean()),
                         "frac_knn_malignant_all": float((g0.knn_lineage == "malignant").mean()),
                         "frac_knn_malignant_vote_ge_0.8_all": float(((g0.knn_lineage == "malignant") & (g0.knn_vote_frac >= 0.8)).mean()),
                         "frac_knn_malignant_expected_lineage_classes": float((gd.knn_lineage == "malignant").mean()),
                         "frac_knn_expected_lineage": float((gd.knn_lineage.values == gd[col].map(lambda k: expected.get(k, "neuron" if q == "siletti" else None)).values).mean()),
                         "n_cells_with_expected_lineage": int(len(gd)),
                         "frac_false_malignant_maxprob_ge_0.8_expected_lineage_classes": float(((gd.predicted_lineage == "malignant") & (gd.scanvi_maxprob >= 0.8)).mean()),
                         "frac_pred_expected_lineage": float((gd.predicted_lineage.values == gd[col].map(lambda k: expected.get(k, "neuron" if q == "siletti" else None)).values).mean()),
                         "classes_without_reference_counterpart": [k for k, v in expected.items() if v is None]}
    spec = pd.DataFrame(spec_rows)
    backup_pre_audit(T / "mapping_normal_reference_specificity.csv")
    spec.to_csv(T / "mapping_normal_reference_specificity.csv", index=False)
    agreement["normal_reference_false_malignant"] = fm_summary
    agreement["_provenance"] = {"predicted_lineage": "reference scANVI classifier (models/scanvi, retrained 2026-10-19 with corrected labels) applied to the zero-shot query latent; audit fixes IM-01 (11c) and reference-labels wave 2 (scripts/11g_remap_queries_zeroshot.py)",
                                "previous_values": "results/tables/mapping_author_agreement.pre_audit.json"}
    backup_pre_audit(T / "mapping_author_agreement.json")
    (T / "mapping_author_agreement.json").write_text(json.dumps(agreement, indent=1, default=float))

    # ------------------------------------------------------------------ 4b. sample-level uncertainty of the headline rates (IM-22)
    sl_rows = []
    def add_sl(query, statistic, df, flag):
        sl_rows.append({"query": query, "statistic": statistic, **sample_level(flag, df.sample_id.values)})
    add_sl("liu", "author_malignant -> predicted malignant (classifier)", lm, lm.predicted_lineage == "malignant")
    add_sl("liu", "author_malignant -> kNN malignant", lm, lm.knn_lineage == "malignant")
    add_sl("liu", "author_normal -> matching lineage (classifier)", ln, ln.predicted_lineage.values == ln.author_annotation.map(LIU_NORMAL).values)
    add_sl("liu", "author_normal -> matching lineage (kNN)", ln, ln.knn_lineage.values == ln.author_annotation.map(LIU_NORMAL).values)
    add_sl("liu", "author T cell -> T_NK (classifier)", ln[ln.author_annotation == "T cell"], ln[ln.author_annotation == "T cell"].predicted_lineage == "T_NK")
    add_sl("liu", "author T cell -> T_NK (kNN)", ln[ln.author_annotation == "T cell"], ln[ln.author_annotation == "T cell"].knn_lineage == "T_NK")
    add_sl("filbin", "consensus_malignant -> predicted malignant (classifier)", fm, fm.predicted_lineage == "malignant")
    add_sl("filbin", "consensus_malignant -> kNN malignant", fm, fm.knn_lineage == "malignant")
    add_sl("labelle", "author immune -> matching lineage (classifier)", li, li.predicted_lineage.values == li.author_broad_annot.map(LABELLE_NORMAL).values)
    add_sl("labelle", "author immune -> matching lineage (kNN)", li, li.knn_lineage.values == li.author_broad_annot.map(LABELLE_NORMAL).values)
    add_sl("labelle", "author Tcell -> T_NK (classifier)", li[li.author_broad_annot == "Tcell"], li[li.author_broad_annot == "Tcell"].predicted_lineage == "T_NK")
    add_sl("labelle", "author Tcell -> T_NK (kNN)", li[li.author_broad_annot == "Tcell"], li[li.author_broad_annot == "Tcell"].knn_lineage == "T_NK")
    for q, col, expected in [("braun", "author_CellClass", BRAUN_EXPECTED), ("siletti", "author_supercluster_term", SILETTI_EXPECTED)]:
        g0 = obs[obs["query"] == q]
        exp = g0[col].map(lambda k: expected.get(k, "neuron" if q == "siletti" else None))
        gd = g0[exp.notna()]; ed = exp[exp.notna()]
        add_sl(q, "all cells -> predicted malignant (classifier)", g0, g0.predicted_lineage == "malignant")
        add_sl(q, "all cells -> kNN malignant", g0, g0.knn_lineage == "malignant")
        add_sl(q, "counterpart classes -> expected lineage (classifier)", gd, gd.predicted_lineage.values == ed.values)
        add_sl(q, "counterpart classes -> expected lineage (kNN)", gd, gd.knn_lineage.values == ed.values)
        add_sl(q, "counterpart classes -> predicted malignant (classifier)", gd, gd.predicted_lineage == "malignant")
    sl = pd.DataFrame(sl_rows)
    sl["sample_unit"] = sl["query"].map({"liu": "plate", "filbin": "sample", "labelle": "sample", "braun": "donor", "siletti": "donor x dissection"})
    sl.to_csv(T / "mapping_agreement_sample_level.csv", index=False)
    LOG.write(f"mapping_agreement_sample_level.csv: {len(sl)} rows (per-sample median/IQR and sample-bootstrap 95% CI)")
    LOG.write(f"agreement: {json.dumps({k: {kk: vv for kk, vv in v.items() if not isinstance(vv, dict)} for k, v in agreement.items()}, default=float)[:1500]}")

    # ------------------------------------------------------------------ 5. cell-of-origin: Braun kNN onto malignant reference cells
    import pynndescent
    with LOG.step("cell-of-origin kNN (Braun -> consensus-malignant reference cells, k=15)"):
        Zref = np.load(ROOT / "models/scanvi/latent_f16.npy").astype(np.float32)
        ref_names = np.array((ROOT / "models/scanvi/obs_names.txt").read_text().split("\n"))
        with h5py.File(ROOT / "data/processed/atlas_core.h5ad") as f:
            core = pd.DataFrame({c: read_elem(f["obs"][c]) for c in ["consensus_call", "leiden_scanvi_r1", "author_Cell_type_consensus_Jessa2022",
                                                                       "study", "sample_id", "scanvi_label"]}, index=read_elem(f["obs"]["_index"]).astype(str))
        core = core.loc[ref_names]
        mal = np.flatnonzero((core.consensus_call.astype(str) == "malignant").values)
        bq = obs[obs["query"] == "braun"]
        with h5py.File(ROOT / "data/processed/query_mapping.h5ad") as f:
            Zq_all = f["obsm/X_scanvi"][:]
        Zb = Zq_all[(obs["query"] == "braun").values]
        idx = pynndescent.NNDescent(Zb, n_neighbors=30, metric="euclidean", random_state=0, n_jobs=8)
        idx.prepare()
        nn, dist = idx.query(Zref[mal], k=15)
        cls = bq.author_CellClass.values[nn]
        fine = bq.author_cell_type.values[nn]
        reg = bq.author_Region.values[nn]

        def majority(arr):
            out, frac = [], []
            for row in arr:
                u, c = np.unique(row, return_counts=True)
                j = c.argmax(); out.append(u[j]); frac.append(c[j] / len(row))
            return np.array(out), np.array(frac)

        top_cls, frac_cls = majority(cls)
        top_fine, frac_fine = majority(fine)
        top_reg, _ = majority(reg)
        per_cell = pd.DataFrame({"cell_id": ref_names[mal], "study": core.study.values[mal], "sample_id": core.sample_id.values[mal],
                                 "leiden_scanvi_r1": core.leiden_scanvi_r1.astype(str).values[mal],
                                 "author_Cell_type_consensus_Jessa2022": core.author_Cell_type_consensus_Jessa2022.astype(object).values[mal],
                                 "nearest_fetal_CellClass": top_cls, "frac_knn_CellClass": frac_cls,
                                 "nearest_fetal_cell_type": top_fine, "frac_knn_cell_type": frac_fine,
                                 "nearest_fetal_Region": top_reg, "mean_knn_dist": dist.mean(1).astype(np.float32)})
        backup_pre_audit(T / "malignant_nearest_fetal_class_per_cell.csv.gz")
        per_cell.to_csv(T / "malignant_nearest_fetal_class_per_cell.csv.gz", index=False)
        base = bq.author_CellClass.value_counts(normalize=True)  # Braun subsample class frequencies (IM-10 base rates)
        out_rows = []
        for grouping, col in [("leiden_scanvi_r1", "leiden_scanvi_r1"), ("jessa_author_state", "author_Cell_type_consensus_Jessa2022"), ("all", None)]:
            groups = [("all", per_cell)] if col is None else [(k, v) for k, v in per_cell.dropna(subset=[col]).groupby(col) if k not in ("nan", "Uncertain")]
            for k, g in groups:
                vc = g.nearest_fetal_CellClass.value_counts(normalize=True)
                vf = g.nearest_fetal_cell_type.value_counts(normalize=True)
                out_rows.append({"grouping": grouping, "group": k, "n_malignant_cells": len(g),
                                 "top_fetal_CellClass": vc.index[0], "frac_top_fetal_CellClass": round(vc.iloc[0], 4),
                                 "second_fetal_CellClass": vc.index[1] if len(vc) > 1 else "", "frac_second_fetal_CellClass": round(vc.iloc[1], 4) if len(vc) > 1 else np.nan,
                                 "top_fetal_cell_type": vf.index[0], "frac_top_fetal_cell_type": round(vf.iloc[0], 4),
                                 "median_knn_dist": float(g.mean_knn_dist.median()),
                                 **{f"frac_{c.replace(' ', '_')}": round(vc.get(c, 0.0), 4) for c in sorted(bq.author_CellClass.unique())},
                                 **{f"enrichment_{c.replace(' ', '_')}": round(vc.get(c, 0.0) / base[c], 2) for c in sorted(bq.author_CellClass.unique())},
                                 "n_samples": g.sample_id.nunique()})
        # per-patient (sample) dispersion of the three leading classes (IM-10)
        pp = per_cell.groupby("sample_id").nearest_fetal_CellClass.value_counts(normalize=True).unstack(fill_value=0.0)
        pp["n_malignant_cells"] = per_cell.groupby("sample_id").size()
        pp = pp[pp.n_malignant_cells >= 50]
        disp = pd.DataFrame({"class": [c for c in ["Glioblast", "Oligo", "Radial glia", "Neuroblast", "Neuron"] if c in pp.columns]})
        disp["braun_subsample_frequency"] = disp["class"].map(base).round(4)
        disp["pooled_frac"] = disp["class"].map(per_cell.nearest_fetal_CellClass.value_counts(normalize=True)).round(4)
        disp["enrichment_pooled"] = (disp.pooled_frac / disp.braun_subsample_frequency).round(2)
        disp["per_sample_median"] = disp["class"].map(pp.median()).round(4)
        disp["per_sample_q25"] = disp["class"].map(pp.quantile(0.25)).round(4)
        disp["per_sample_q75"] = disp["class"].map(pp.quantile(0.75)).round(4)
        disp["n_samples_ge50_cells"] = len(pp)
        disp.to_csv(T / "malignant_nearest_fetal_class_per_sample_dispersion.csv", index=False)
        backup_pre_audit(T / "malignant_nearest_fetal_class.csv")
        pd.DataFrame(out_rows).to_csv(T / "malignant_nearest_fetal_class.csv", index=False)
        LOG.write("malignant_nearest_fetal_class_per_sample_dispersion.csv:\n" + disp.to_string())
        LOG.write(f"malignant_nearest_fetal_class.csv: {len(out_rows)} rows; overall top class "
                  f"{per_cell.nearest_fetal_CellClass.value_counts(normalize=True).head(4).round(3).to_dict()}")
    LOG.write(f"done | peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
