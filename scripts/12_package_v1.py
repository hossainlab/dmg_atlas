#!/usr/bin/env python
"""Phase 5 / step 12 — atlas v1 packaging.

Writes
  docs/atlas_schema.md                      every obs / obsm / var / layers / uns field of atlas_core.h5ad and
                                            query_mapping.h5ad: type, allowed values, source (generated from the files)
  docs/label_ontology.csv                   label, definition, marker genes, source tier
  results/tables/atlas_v1_summary.csv       per study: samples, patients, cells, assay, site, H3 variants, malignant fraction
  README.md                                 section 'Atlas v1 contents' appended (existing text untouched)
  results/figures/phase5_query_mapping.png  reference UMAP (grey) + query cells per query, coloured by predicted lineage
  results/figures/phase5_mapping_uncertainty.png  per-sample uncertainty distributions (kNN percentile, 1 - max prob)
"""
from __future__ import annotations

import sys
import json
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style, LINEAGE_COLOURS, QUERY_COLOURS, QUERY_NAMES, letter  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from anndata.io import read_elem  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase5_mapping.log", context="12_package_v1")
CORE = ROOT / "data/processed/atlas_core.h5ad"
QMAP = ROOT / "data/processed/query_mapping.h5ad"
FIG = ROOT / "results/figures"
FIG.mkdir(exist_ok=True, parents=True)
QUERIES = ["filbin", "labelle", "liu", "braun", "siletti"]

SOURCE_RULES = [  # (predicate on column name, source)
    (lambda c: c.startswith("author_"), "author metadata (Phase 2 ingest; never modified)"),
    (lambda c: c in {"sample_id", "patient_id", "study", "role", "assay", "chemistry", "site", "age_years", "age_group", "h3_variant",
                     "k27m_status", "sex", "treatment", "tissue_timing", "include_in_atlas", "exclusion_reason", "counts_type", "label_source"},
     "sample_manifest.csv (Phase 1)"),
    (lambda c: c in {"original_cell_id", "n_counts", "n_genes", "pct_mito", "pct_ribo", "pct_hb", "qc_pass", "qc_fail_reason", "qc_assay_class",
                     "qc_min_genes", "qc_max_genes", "qc_pct_mito_max", "doublet_score", "predicted_doublet", "predicted_doublet_fixed"},
     "Phase 2 QC (scripts/03_ingest_qc.py)"),
    (lambda c: c.startswith(("cnv_", "arm_", "chr", "score_", "marker_", "consensus_", "author_call")) or c in {"high_mito_nuclei"},
     "Phase 3 malignant calling (malignant_calls.csv)"),
    (lambda c: c in {"lineage_label", "scanvi_label", "scanvi_maxprob", "leiden_scanvi_r1"}, "Phase 4 integration"),
    (lambda c: c in {"lineage_label_pre_audit", "scanvi_label_pre_audit", "scanvi_maxprob_pre_audit", "leiden_scanvi_r1_pre_audit"} ,
     "Phase 4 integration — pre-audit value (before the 2026-10-19 label correction / scANVI retrain; superseded)"),
    (lambda c: c.endswith("_pre_audit") and c.startswith(("cnv_", "arm_", "chr", "marker_", "consensus_")),
     "Phase 3 malignant calling — pre-audit value from the superseded malignant_calls (CNV audit 2026-10-19)"),
    (lambda c: c.endswith("_pre_audit"), "Phase 5 query mapping — original Phase 5 value (frozen per-query classifier / kNN on the pre-audit latent; superseded, audit IM-01)"),
    (lambda c: c.endswith("_fix1"), "Phase 5 query mapping — wave-1 value (reference classifier on the pre-audit latent, 11c; superseded by the retrained reference, wave 2)"),
    (lambda c: c in {"cnv_confidence", "consensus_call_eventtrust", "low_confidence_sample", "gene_dev_PDGFRA", "gene_dev_KIT"},
     "Phase 3 malignant calling (malignant_calls.csv; CNV audit additions 2026-10-19)"),
    (lambda c: c == "scanvi_validation_patient", "Phase 4 integration — cell belongs to one of the 5 patients held out from scANVI classifier training (models/scanvi/split.json)"),
    (lambda c: c == "liu_prep", "Phase 5 query mapping (audit IM-24: Liu fresh cell / frozen nucleus flag)"),
    (lambda c: c in {"query", "predicted_lineage", "uncertainty_prob", "knn_mean_dist", "uncertainty_knn", "query_n_counts_hvg", "query_n_genes_hvg",
                     "donor_id", "labelle_platform", "labelle_table"}, "Phase 5 query mapping"),
]


def source_of(c: str) -> str:
    for pred, s in SOURCE_RULES:
        if pred(c):
            return s
    return "Phase 2 ingest"


def describe_obs(path: Path, max_cats: int = 12) -> pd.DataFrame:
    rows = []
    with h5py.File(path) as f:
        g = f["obs"]
        for c in g.keys():
            if c == "_index" or c == g.attrs.get("_index"):
                continue
            el = g[c]
            enc = el.attrs.get("encoding-type", "")
            if enc == "categorical":
                cats = read_elem(el["categories"]).astype(str)
                codes = el["codes"][:]
                n_missing = int((codes < 0).sum())
                vals = "; ".join(sorted(cats)) if len(cats) <= max_cats else f"{len(cats)} categories (e.g. {'; '.join(sorted(cats)[:5])}, ...)"
                rows.append({"field": c, "type": "category", "allowed_values": vals, "n_missing": n_missing, "source": source_of(c)})
            elif enc == "nullable-boolean" or (enc == "array" and el.dtype == bool):
                rows.append({"field": c, "type": "bool", "allowed_values": "True; False", "n_missing": 0, "source": source_of(c)})
            elif enc in ("nullable-integer",):
                rows.append({"field": c, "type": "int (nullable)", "allowed_values": "", "n_missing": int(el["mask"][:].sum()), "source": source_of(c)})
            elif enc == "string-array":
                rows.append({"field": c, "type": "string", "allowed_values": "free text", "n_missing": 0, "source": source_of(c)})
            else:
                arr = el[:]
                if np.issubdtype(arr.dtype, np.number):
                    fin = arr[np.isfinite(arr.astype(float))]
                    rng = f"[{fin.min():.4g}, {fin.max():.4g}]" if len(fin) else ""
                    rows.append({"field": c, "type": str(arr.dtype), "allowed_values": f"numeric range {rng}", "n_missing": int(len(arr) - len(fin)), "source": source_of(c)})
                else:
                    rows.append({"field": c, "type": str(arr.dtype), "allowed_values": "", "n_missing": 0, "source": source_of(c)})
    return pd.DataFrame(rows)


def describe_other(path: Path) -> list[str]:
    out = []
    with h5py.File(path) as f:
        n_obs = f["obs"][f["obs"].attrs["_index"]].shape[0]
        n_var = f["var"][f["var"].attrs["_index"]].shape[0]
        out.append(f"Shape: {n_obs:,} obs x {n_var:,} var.")
        if "X" in f:
            x = f["X"]
            shape = tuple(x.attrs["shape"]) if isinstance(x, h5py.Group) else x.shape
            dt = x["data"].dtype if isinstance(x, h5py.Group) else x.dtype
            out.append(f"X: {x.attrs.get('encoding-type', 'dense')} {shape}, dtype {dt}.")
        for grp, desc in [("layers", "layers"), ("obsm", "obsm"), ("varm", "varm"), ("obsp", "obsp")]:
            if grp in f and len(f[grp]):
                for k in f[grp]:
                    el = f[grp][k]
                    shape = tuple(el.attrs["shape"]) if isinstance(el, h5py.Group) and "shape" in el.attrs else getattr(el, "shape", "")
                    dt = el["data"].dtype if isinstance(el, h5py.Group) and "data" in el else getattr(el, "dtype", "")
                    out.append(f"{desc}['{k}']: shape {shape}, dtype {dt}.")
        out.append("var columns: " + ", ".join(k for k in f["var"].keys() if k != "_index") + f"; var index = {f['var'].attrs['_index']}.")
        if "uns" in f:
            out.append("uns keys: " + ", ".join(f["uns"].keys()) + ".")
    return out


OBSM_DESC = {
    "X_scvi": ("float32 (n, 30)", "scVI latent, model models/scvi_seed1 (Phase 4)"),
    "X_scanvi": ("float32 (n, 30)", "scANVI latent, model models/scanvi (Phase 4; queries: zero-shot projection by the same reference encoder, Phase 5 — see phase5_mapping_notes.md §2)"),
    "X_umap": ("float32 (n, 2)", "scanpy UMAP on X_scanvi (n_neighbors 30, random_state 0); no transform model stored"),
    "X_scvi_malignant": ("float32 (n, 30)", "malignant-only scVI latent (models/scvi_malignant, 167,395 corrected consensus-malignant cells); NaN rows for non-malignant cells and for okada__1734_pons (low-confidence sample)"),
    "X_umap_ref": ("float32 (n, 2)", "umap-learn transform of query X_scanvi with a UMAP refitted on a 100,000-cell reference subsample "
                                     "(results/tables/phase5_umap_reference_fit.csv); not the atlas_core X_umap coordinate system"),
    "scanvi_probs": ("float32 (n, 7)", "scANVI soft class probabilities (queries: reference classifier models/scanvi on X_scanvi); column order in uns['phase5']['scanvi_prob_columns']"),
    "scanvi_probs_pre_audit": ("float32 (n, 7)", "superseded query probabilities from the drifted per-query classifiers (audit IM-01); kept for traceability; NaN for cells absent from the previous file"),
    "scanvi_probs_fix1": ("float32 (n, 7)", "wave-1 probabilities (reference classifier on the pre-audit latent, 11c); superseded by the retrained reference (wave 2); NaN for cells absent from the previous file"),
    "X_scanvi_pre_audit": ("float32 (n, 30)", "pre-audit scANVI latent (models/scanvi_pre_audit); superseded 2026-10-19"),
    "X_umap_pre_audit": ("float32 (n, 2)", "pre-audit UMAP on X_scanvi_pre_audit; superseded 2026-10-19"),
    "X_scvi_malignant_pre_audit": ("float32 (n, 30)", "pre-audit malignant-only latent (models/scvi_malignant_pre_audit, 169,571 cells); superseded 2026-10-19"),
}


def write_schema():
    core_obs = describe_obs(CORE)
    q_obs = describe_obs(QMAP)
    lines = ["# Atlas v1 — object schema", "",
             "Generated by `scripts/12_package_v1.py` from the files themselves (field list, dtype, categories, missing counts). "
             "`source` names the pipeline phase that produced the field; author columns are copied verbatim from the deposits.", "",
             "## 1. `data/processed/atlas_core.h5ad` (tumour + microenvironment reference)", ""]
    lines += ["- " + s for s in describe_other(CORE)]
    lines += ["", "Notes: `X` is intentionally empty (counts live in `layers['counts']` to keep the file < 6 GB); set `adata.X = adata.layers['counts']` before normalising. "
              "obs index = `<sample_id>:<original barcode/cell id>`. var index = GRCh38-2020-A canonical symbol (18,882 genes measured in every UMI sample).", ""]
    lines += ["### obsm", "", "| key | type | meaning |", "|---|---|---|"]
    with h5py.File(CORE) as f:
        for k in f["obsm"]:
            t, d = OBSM_DESC.get(k, ("", ""))
            lines.append(f"| {k} | {t} | {d} |")
    lines += ["", "### obs", "", "| field | type | allowed values / range | missing | source |", "|---|---|---|---|---|"]
    for r in core_obs.itertuples():
        lines.append(f"| {r.field} | {r.type} | {r.allowed_values} | {r.n_missing} | {r.source} |")
    lines += ["", "### var", "", "| field | meaning |", "|---|---|",
              "| (index) symbol | GRCh38-2020-A canonical gene symbol |",
              "| ensembl_id | 2020-A Ensembl gene ID (NaN for symbols reached only through an alias) |",
              "| n_samples_measured / measured_in_all_samples | gene present in the deposited feature list of n samples / all samples (Phase 2) |",
              "| highly_variable | one of the 3,000 seurat_v3_paper HVGs used for scVI/scANVI (results/tables/atlas_input_hvg.csv) |",
              "| hvg_seurat_v3_medianrank_first | HVG flag under scanpy's median-rank-first ranking (kept for the record; not used) |",
              "", "## 2. `data/processed/query_mapping.h5ad` (projected queries)", ""]
    lines += ["- " + s for s in describe_other(QMAP)]
    lines += ["", "Notes: `X` = `layers['counts']` = reference-HVG count matrix per query (UMI for Braun/Siletti; round(RSEM) for Liu; TPM-derived pseudo-counts "
              "for Filbin/LaBelle — see docs/phase5_mapping_notes.md). Query cells never enter the reference model; labels are scANVI predictions.", ""]
    lines += ["### obsm", "", "| key | type | meaning |", "|---|---|---|"]
    with h5py.File(QMAP) as f:
        for k in f["obsm"]:
            t, d = OBSM_DESC.get(k, ("", ""))
            lines.append(f"| {k} | {t} | {d} |")
    lines += ["", "### obs", "", "| field | type | allowed values / range | missing | source |", "|---|---|---|---|---|"]
    for r in q_obs.itertuples():
        lines.append(f"| {r.field} | {r.type} | {r.allowed_values} | {r.n_missing} | {r.source} |")
    lines += ["", "## 3. Label vocabularies", "", "See `docs/label_ontology.csv` (lineage labels, consensus tiers, marker genes, source tier).", "",
              "## 4. Models", "", "| directory | content |", "|---|---|",
              "| models/scvi_seed1 | chosen scVI reference (26 epochs; seeds 0/2 kept as scvi_seed0 / scvi_seed2) |",
              "| models/scanvi | scANVI reference (from scvi_seed1; labels_key lineage_label; all labelled cells) |",
              "| models/scvi_malignant | malignant-only scVI (169,571 cells; HVGs results/tables/atlas_malignant_hvg.csv) |",
              "| models/scanvi_query_<filbin,labelle,liu,braun,siletti> | per-query load_query_data models (Phase 5; provenance only — classifier BatchNorm drifted, see audit IM-01; labels in query_mapping.h5ad come from models/scanvi) |",
              "| models/scvi_seed0_e60, scvi_malignant_seed0, scanvi_nspl100_rejected, pca, harmony | sensitivity / rejected / control runs |",
              "", "Each model directory holds model.pt (save_anndata=False), latent_f16.npy, obs_names.txt, history.json, run.json."]
    (ROOT / "docs/atlas_schema.md").write_text("\n".join(lines) + "\n")
    LOG.write(f"docs/atlas_schema.md: {len(core_obs)} atlas_core obs fields, {len(q_obs)} query_mapping obs fields")


def write_ontology():
    rows = [
        # lineage labels (scANVI classes)
        ("lineage_label / scanvi_label / predicted_lineage", "malignant", "H3K27M DMG tumour cell (any state: OPC-, AC-, OC-, MES-like, cycling)",
         "CNV-defined; states: PDGFRA OLIG1 OLIG2 SOX10 CSPG4 (OPC-like), GFAP AQP4 APOE CLU (AC-like), CD44 VIM CHI3L1 ANXA1 (MES-like), MKI67 TOP2A (cycling)",
         "tier 3 consensus (author label + own CNV + marker agreement); Phase 3"),
        ("lineage_label / scanvi_label / predicted_lineage", "myeloid", "microglia / macrophage / monocyte", "PTPRC CSF1R CX3CR1 P2RY12 C1QA C1QB CD74", "tier 3 consensus; Phase 3"),
        ("lineage_label / scanvi_label / predicted_lineage", "T_NK", "T and NK lymphocytes (smallest class; scANVI recall 0.47 — frequently absorbed by myeloid)", "CD3D CD3E CD2 NKG7", "tier 3 consensus; Phase 3"),
        ("lineage_label / scanvi_label / predicted_lineage", "oligodendrocyte", "mature oligodendrocyte (normal, no CNV)", "MBP PLP1 MOG MAG ST18", "tier 3 consensus; Phase 3"),
        ("lineage_label / scanvi_label / predicted_lineage", "astrocyte", "normal astrocyte (no hallmark CNV; AC-like tumour cells are 'malignant')", "AQP4 GFAP", "tier 3 consensus (CNV-negative astrocyte-like); Phase 3"),
        ("lineage_label / scanvi_label / predicted_lineage", "neuron", "neuron", "RBFOX3 SYT1 SNAP25 GAD1 SLC17A7", "tier 3 consensus; Phase 3"),
        ("lineage_label / scanvi_label / predicted_lineage", "vascular", "endothelial cell, pericyte, vascular smooth muscle, perivascular fibroblast", "CLDN5 FLT1 VWF PDGFRB RGS5", "tier 3 consensus; Phase 3"),
        ("lineage_label", "Unknown", "scANVI unlabelled category: consensus ambiguous / unassigned cells (91,309 in the reference); queries enter with this label", "", "Phase 4"),
        # consensus tiers
        ("consensus_call", "malignant", "tier-3 malignant: author malignant where available, else CNV-malignant with marker agreement", "", "Phase 3 tier 3"),
        ("consensus_call", "normal", "tier-3 normal", "", "Phase 3 tier 3"),
        ("consensus_call", "ambiguous", "CNV and marker calls disagree or no reference available", "", "Phase 3 tier 3"),
        ("consensus_call", "unassigned", "sample without usable CNV reference (scpca SCPCS000637)", "", "Phase 3 tier 3"),
        ("cnv_call", "malignant | normal", "tier-2 own inferCNV cluster call (hallmark events chr1q/2/7/4q12 gain, 13q/14q/10q/17p loss; score elevation)", "", "Phase 3 tier 2"),
        ("author_call", "malignant | normal", "tier-1 author label mapped to malignant/normal (Jessa, Okada, Collot, Liu; immune-only for LaBelle)", "", "Phase 3 tier 1"),
        ("marker_call", "normal | astrocyte_like | malignant_candidate", "marker-score call used to pick CNV references and to arbitrate tier 3", "see lineage marker sets", "Phase 3"),
        ("consensus_celltype_coarse", "malignant | myeloid | T_NK | oligodendrocyte | astrocyte | neuron | vascular | other | ambiguous", "coarse type behind lineage_label", "", "Phase 3"),
        # query-side vocabularies kept verbatim
        ("author_annotation (Liu 2022)", "OPC-like-1/2/3, AC-like, MES-like, OC-like, S, G2M, frozen-OPC-like-b | Myeloid, Microglia, Oligodendrocyte, Endothelia, T cell | score-too-low",
         "Liu 2022 malignant states and TME types (author vocabulary, not modified)", "", "author tier 1"),
        ("author_Cell_type_consensus_Jessa2022", "OPC, Proliferating OPC, Astrocytes, Oligodendrocytes, Glial progenitors, Ependymal, Neurons, RGC, Immune, Vascular, Uncertain, ...",
         "Jessa 2022 consensus cell type (mouse-correlation based); malignant states are developmental-lineage names", "", "author tier 1"),
        ("author_CellClass (Braun 2023)", "Radial glia, Glioblast, Neuroblast, Neuron, Neuronal IPC, Oligo, Vascular, Fibroblast, Immune, Erythrocyte, Placodes, Neural crest",
         "first-trimester hindbrain cell classes; Radial glia / Glioblast have no normal counterpart in the tumour reference", "", "author (cellxgene)"),
        ("author_supercluster_term (Siletti 2023)", "Oligodendrocyte, Oligodendrocyte precursor, Committed oligodendrocyte precursor, Astrocyte, Bergmann glia, Ependymal, Microglia, Vascular, Fibroblast, neuronal superclusters, Splatter, ...",
         "adult pons superclusters; OPC has no normal counterpart in the tumour reference", "", "author (cellxgene)"),
    ]
    df = pd.DataFrame(rows, columns=["field", "label", "definition", "marker_genes", "source_tier"])
    df.to_csv(ROOT / "docs/label_ontology.csv", index=False)
    LOG.write(f"docs/label_ontology.csv: {len(df)} rows")


def write_summary():
    with h5py.File(CORE) as f:
        obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in ["study", "sample_id", "patient_id", "assay", "site", "h3_variant", "consensus_call",
                                                                 "scanvi_label", "lineage_label", "age_group", "tissue_timing"]})
    rows = []
    for st, g in obs.groupby("study", observed=True):
        rows.append({"study": st, "role": "reference", "samples": g.sample_id.nunique(), "patients": g.patient_id.nunique(), "cells": len(g),
                     "assay": ";".join(f"{k}:{v}" for k, v in g.assay.value_counts().items() if v),
                     "site": ";".join(f"{k}:{g[g.site == k].sample_id.nunique()}" for k in g.site.unique()),
                     "h3_variants": ";".join(f"{k}:{g[g.h3_variant == k].sample_id.nunique()}" for k in g.h3_variant.unique()),
                     "age_groups": ";".join(f"{k}:{g[g.age_group == k].sample_id.nunique()}" for k in g.age_group.unique()),
                     "malignant_fraction_consensus": round((g.consensus_call.astype(str) == "malignant").mean(), 4),
                     "malignant_fraction_scanvi": round((g.scanvi_label.astype(str) == "malignant").mean(), 4),
                     "ambiguous_fraction_consensus": round((g.consensus_call.astype(str) == "ambiguous").mean(), 4)})
    with h5py.File(QMAP) as f:
        q = pd.DataFrame({c: read_elem(f["obs"][c]) for c in ["query", "study", "sample_id", "patient_id", "donor_id", "assay", "site", "age_group", "h3_variant",
                                                               "predicted_lineage", "scanvi_maxprob", "consensus_call"]})
    for qn, g in q.groupby("query", observed=True):
        pid = g.patient_id if g.patient_id.notna().any() else g.donor_id
        rows.append({"study": str(g.study.iloc[0]), "role": "query (projection only)", "samples": g.sample_id.nunique(), "patients": pid.nunique(), "cells": len(g),
                     "assay": ";".join(f"{k}:{v}" for k, v in g.assay.value_counts().items() if v),
                     "site": ";".join(f"{k}:{g[g.site == k].sample_id.nunique()}" for k in g.site.dropna().unique()) or "normal brain",
                     "h3_variants": ";".join(f"{k}:{g[g.h3_variant == k].sample_id.nunique()}" for k in g.h3_variant.dropna().unique()) or "n/a",
                     "age_groups": ";".join(f"{k}:{g[g.age_group == k].sample_id.nunique()}" for k in g.age_group.dropna().unique()) or "n/a",
                     "malignant_fraction_consensus": round((g.consensus_call.astype(str) == "malignant").mean(), 4) if g.consensus_call.notna().any() else np.nan,
                     "malignant_fraction_scanvi": round((g.predicted_lineage.astype(str) == "malignant").mean(), 4),
                     "ambiguous_fraction_consensus": round((g.consensus_call.astype(str) == "ambiguous").mean(), 4) if g.consensus_call.notna().any() else np.nan})
    df = pd.DataFrame(rows)
    tot = {"study": "TOTAL reference", "role": "reference", "samples": obs.sample_id.nunique(), "patients": obs.patient_id.nunique(), "cells": len(obs),
           "assay": ";".join(f"{k}:{v}" for k, v in obs.assay.value_counts().items() if v), "site": "", "h3_variants": "", "age_groups": "",
           "malignant_fraction_consensus": round((obs.consensus_call.astype(str) == "malignant").mean(), 4),
           "malignant_fraction_scanvi": round((obs.scanvi_label.astype(str) == "malignant").mean(), 4),
           "ambiguous_fraction_consensus": round((obs.consensus_call.astype(str) == "ambiguous").mean(), 4)}
    df = pd.concat([df, pd.DataFrame([tot])])
    df.to_csv(ROOT / "results/tables/atlas_v1_summary.csv", index=False)
    LOG.write(f"results/tables/atlas_v1_summary.csv: {len(df)} rows")
    return df


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def append_readme(summary: pd.DataFrame):
    readme = ROOT / "README.md"
    txt = readme.read_text()
    if "## Atlas v1 contents" in txt:
        LOG.write("README already has an 'Atlas v1 contents' section; appending a dated update block instead of duplicating")
    ref = summary[summary.role == "reference"]
    tot = ref[ref.study == "TOTAL reference"].iloc[0]
    q = summary[summary.role != "reference"]
    lines = ["", "## Atlas v1 contents", "",
             f"Reference: **{int(tot.cells):,} cells / nuclei, {int(tot.samples)} samples, {int(tot.patients)} patients, "
             f"{len(ref) - 1} studies** ({tot.assay}); consensus-malignant fraction {tot.malignant_fraction_consensus:.3f}. "
             f"Queries (projection only): " + "; ".join(f"{r.study} {int(r.cells):,} cells / {int(r.samples)} samples" for r in q.itertuples()) + ".", "",
             "| file | content |", "|---|---|",
             "| `data/processed/atlas_core.h5ad` | reference counts (`layers['counts']`), obs (manifest, QC, author labels, CNV calls, lineage labels, scANVI labels, Leiden), obsm X_scvi / X_scanvi / X_umap / X_scvi_malignant |",
             "| `data/processed/query_mapping.h5ad` | Filbin 2018, LaBelle 2025, Liu 2022, Braun 2023 hindbrain, Siletti 2023 pons mapped onto the scANVI reference (predicted lineage, soft probabilities, two uncertainty measures, X_scanvi, X_umap_ref) |",
             "| `data/processed/sample_manifest.csv` | per-sample metadata (patient, assay, chemistry, site, age, H3 variant, treatment, tissue timing) |",
             "| `data/processed/malignant_calls.csv` | per-cell three-tier malignant calls (author / own CNV / consensus) with hallmark arm flags |",
             "| `models/` | scvi_seed{0,1,2}, scanvi, scvi_malignant, scanvi_query_<query>; sensitivity and control runs |",
             "| `results/tables/` | atlas_v1_summary.csv, integration_metrics.csv, malignant_call_validation.csv, mapping_summary.csv, mapping_uncertainty.csv, mapping_<query>_author_vs_predicted.csv, mapping_normal_reference_specificity.csv, malignant_nearest_fetal_class.csv, query_preparation.csv |",
             "| `results/figures/` | phase4_umap_overview.png, phase5_query_mapping.png, phase5_mapping_uncertainty.png |",
             "| `docs/` | decisions.md, atlas_schema.md, label_ontology.csv, phase2_qc_notes.md, phase3_cnv_notes.md, phase4_integration_notes.md, phase5_mapping_notes.md |",
             "| `logs/` | phase2_qc.log, phase3_cnv.log, phase4_integration.log, phase4_modal_runs_*.jsonl, phase5_mapping.log (wall time + peak RSS per step) |",
             "", "Label vocabulary and obs schema: `docs/label_ontology.csv`, `docs/atlas_schema.md`. Per-study summary: `results/tables/atlas_v1_summary.csv`.",
             f"Checksums: atlas_core.h5ad md5 `{md5(CORE)}`; query_mapping.h5ad md5 `{md5(QMAP)}`.",
             "", "Known limitations of v1: T_NK is under-resolved (764 reference cells, scANVI recall 0.47; query T cells are absorbed by myeloid); "
             "TPM cohorts are mapped with pseudo-counts and are projection-only; normal progenitor classes (fetal radial glia / glioblast, adult OPC) have no "
             "normal counterpart in the reference and are predicted 'malignant' — see docs/phase5_mapping_notes.md.", ""]
    with open(readme, "a") as fh:
        fh.write("\n".join(lines))
    LOG.write("README.md: appended 'Atlas v1 contents'")


# ----------------------------------------------------------------------------------------------------- figures
def figures():
    """Audit revision 2026-10-19 (IM-13/IM-14): no text inside panels (panel letters outside; all identification in the captions in
    docs/phase5_mapping_notes.md §7); colours encode the reference-classifier labels; the flat '1 - max probability' row is replaced by the
    per-sample classifier-vs-kNN agreement. Previous versions kept as *.pre_audit.png."""
    import shutil
    apply_figure_style()
    for name in ["phase5_query_mapping.png", "phase5_mapping_uncertainty.png"]:
        if (FIG / name).exists() and not (FIG / name.replace(".png", ".pre_audit.png")).exists():
            shutil.copy2(FIG / name, FIG / name.replace(".png", ".pre_audit.png"))
    with h5py.File(QMAP) as f:
        obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in ["query", "sample_id", "predicted_lineage", "knn_lineage", "uncertainty_knn",
                                                                 "uncertainty_prob", "scanvi_maxprob", "classifier_knn_agree"]})
        um = f["obsm/X_umap_ref"][:]
    for c in ["query", "sample_id", "predicted_lineage", "knn_lineage"]:
        obs[c] = obs[c].astype(str)
    ref = pd.read_csv(ROOT / "results/tables/phase5_umap_reference_fit.csv")
    rng = np.random.default_rng(0)
    letters = "abcdefghij"

    # Figure 1 — reference UMAP grey; row 1 query cells by reference-classifier label, row 2 by kNN label; one column per query
    fig, axes = plt.subplots(2, 5, figsize=(13, 5.6))
    xlim = (min(ref.umap_1.min(), um[:, 0].min()), max(ref.umap_1.max(), um[:, 0].max()))
    ylim = (min(ref.umap_2.min(), um[:, 1].min()), max(ref.umap_2.max(), um[:, 1].max()))
    n_shown = {}
    for j, q in enumerate(QUERIES):
        m = (obs["query"] == q).values
        idx = np.flatnonzero(m)
        if len(idx) > 30_000:
            idx = rng.choice(idx, 30_000, replace=False)
        idx = rng.permutation(idx)
        n_shown[q] = (int(m.sum()), int(len(idx)))
        for i, col in enumerate(["predicted_lineage", "knn_lineage"]):
            ax = axes[i, j]
            ax.scatter(ref.umap_1, ref.umap_2, s=0.5, c="#d9d9d9", linewidths=0, rasterized=True)
            cols = obs[col].values[idx]
            ax.scatter(um[idx, 0], um[idx, 1], s=1.2, c=[LINEAGE_COLOURS.get(c, "#000000") for c in cols], linewidths=0, rasterized=True)
            ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            letter(ax, letters[i * 5 + j], dx=0.0, dy=1.0)
    ax0 = axes[1, 0]
    ax0.annotate("", xy=(0.12, 0.02), xytext=(0.02, 0.02), xycoords="axes fraction", arrowprops=dict(arrowstyle="->", lw=0.6))
    ax0.annotate("", xy=(0.02, 0.12), xytext=(0.02, 0.02), xycoords="axes fraction", arrowprops=dict(arrowstyle="->", lw=0.6))
    ax0.set_xlabel("UMAP 1", fontsize=6, loc="left", labelpad=1); ax0.set_ylabel("UMAP 2", fontsize=6, loc="bottom", labelpad=1)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=5, color=LINEAGE_COLOURS[l], label=l) for l in
               ["malignant", "myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular"]]
    handles.append(plt.Line2D([], [], marker="o", ls="", ms=5, color="#d9d9d9", label="reference (100k subsample)"))
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.905, 0.5), title="lineage label")
    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    fig.savefig(FIG / "phase5_query_mapping.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Figure 2 — per query sample: kNN uncertainty distribution (row 1) and classifier-vs-kNN label agreement (row 2)
    sizes = [obs[obs["query"] == q].sample_id.nunique() for q in QUERIES]
    fig, axes = plt.subplots(2, 5, figsize=(16, 4.8), sharey="row",
                             gridspec_kw={"width_ratios": [max(s, 12) * (1.6 if q == "siletti" else 1.0) for s, q in zip(sizes, QUERIES)], "wspace": 0.08, "hspace": 0.15, "height_ratios": [2, 1]})
    order = {}
    for j, q in enumerate(QUERIES):
        g = obs[obs["query"] == q]
        samples = g.groupby("sample_id").uncertainty_knn.median().sort_values().index.tolist()
        order[q] = samples
        ax = axes[0, j]
        data = [g.loc[g.sample_id == s, "uncertainty_knn"].values for s in samples]
        vp = ax.violinplot(data, positions=range(len(samples)), showmedians=False, showextrema=False, widths=0.7)
        for b in vp["bodies"]:
            b.set_facecolor(QUERY_COLOURS[q]); b.set_edgecolor("none"); b.set_alpha(0.45)
        ax.scatter(range(len(samples)), [np.median(d) for d in data], s=6, c="black", zorder=3, linewidths=0)
        ax.axhline(0.5, color="#636363", lw=0.5, ls=":")
        ax.set_xlim(-0.8, len(samples) - 0.2); ax.set_ylim(-0.03, 1.03)
        ax.set_xticks(range(len(samples))); ax.set_xticklabels([])
        ax2 = axes[1, j]
        agree = [g.loc[g.sample_id == s, "classifier_knn_agree"].astype(float).mean() for s in samples]
        ax2.bar(range(len(samples)), agree, width=0.7, color=QUERY_COLOURS[q], alpha=0.7, linewidth=0)
        ax2.set_xlim(-0.8, len(samples) - 0.2); ax2.set_ylim(0, 1.02)
        ax2.set_xticks(range(len(samples)))
        ax2.set_xticklabels([s.split("__", 1)[-1] for s in samples], rotation=90, fontsize=6)
        for i, a in enumerate([ax, ax2]):
            if j == 0:
                a.set_ylabel(["kNN uncertainty\n(reference-distance percentile)", "classifier = kNN\nlabel (fraction)"][i])
            else:
                a.tick_params(axis="y", length=0)
            a.spines["left"].set_visible(j == 0)
        letter(ax, letters[j], dx=0.0, dy=1.0)
    axes[1, 2].set_xlabel("query sample")
    fig.savefig(FIG / "phase5_mapping_uncertainty.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    json.dump({"figure1_cells_total_shown": n_shown, "figure2_sample_order": order},
              open(ROOT / "results/figures/phase5_figure_metadata.json", "w"), indent=1)
    LOG.write("figures: results/figures/phase5_query_mapping.png, phase5_mapping_uncertainty.png (audit revision; metadata in phase5_figure_metadata.json)")


if __name__ == "__main__":
    which = sys.argv[1:] or ["schema", "ontology", "summary", "readme", "figures"]
    summary = None
    for w in which:
        with LOG.step(f"package {w}"):
            if w == "schema":
                write_schema()
            elif w == "ontology":
                write_ontology()
            elif w == "summary":
                summary = write_summary()
            elif w == "readme":
                summary = summary if summary is not None else pd.read_csv(ROOT / "results/tables/atlas_v1_summary.csv")
                append_readme(summary)
            elif w == "figures":
                figures()
    LOG.write(f"done | peak RSS {peak_rss_gb():.2f} GB")
