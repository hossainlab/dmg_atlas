"""Assay-specific per-sample QC. Thresholds come from configs/qc.yaml. Cells are flagged (obs['qc_pass']), not removed."""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import yaml

from . import CONFIGS

NUCLEI_ASSAYS = {"snRNA", "multiome_RNA"}
CELL_ASSAYS = {"scRNA"}
SS2_ASSAYS = {"smartseq2", "smartseq2+scRNA"}


def load_config(path: Path | None = None) -> dict:
    with open(path or CONFIGS / "qc.yaml") as fh:
        return yaml.safe_load(fh)


def assay_class(assay: str, platform: str | None = None) -> str:
    if platform == "10x":  # LaBelle 10x cells inside a smartseq2-labelled sample
        return "cells"
    if assay in NUCLEI_ASSAYS:
        return "nuclei"
    if assay in CELL_ASSAYS:
        return "cells"
    if assay in SS2_ASSAYS:
        return "smartseq2"
    raise ValueError(f"unknown assay {assay}")


EXPRESSION_LAYERS = ("counts", "tpm", "log2_tpm10")  # first present wins; 'log2_tpm10' = LaBelle 2025 log2(TPM/10+1) tables (audit DQ-01)


def _matrix(adata: ad.AnnData):
    for key in EXPRESSION_LAYERS:
        if key in adata.layers:
            return adata.layers[key], key
    raise KeyError(f"no expression layer among {EXPRESSION_LAYERS}: {list(adata.layers)}")


def compute_metrics(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    """n_counts, n_genes, pct_mito, pct_ribo, pct_hb from the counts (or TPM) layer. Sparse throughout."""
    M, kind = _matrix(adata)
    M = sp.csr_matrix(M)
    genes = adata.var_names.astype(str)
    mito = np.char.startswith(genes.values.astype(str), cfg["mito_gene_prefix"])
    ribo = np.zeros(len(genes), dtype=bool)
    for p in cfg["ribo_gene_prefixes"]:
        ribo |= np.char.startswith(genes.values.astype(str), p)
    hb = np.array([bool(re.match(cfg["hb_gene_regex"], g)) for g in genes])
    tot = np.asarray(M.sum(axis=1)).ravel().astype(np.float64)
    adata.obs["n_counts"] = tot
    adata.obs["n_genes"] = np.asarray((M > 0).sum(axis=1)).ravel().astype(np.int32)
    with np.errstate(divide="ignore", invalid="ignore"):
        for name, mask in (("pct_mito", mito), ("pct_ribo", ribo), ("pct_hb", hb)):
            s = np.asarray(M[:, mask].sum(axis=1)).ravel() if mask.any() else np.zeros(M.shape[0])
            adata.obs[name] = np.where(tot > 0, 100.0 * s / tot, np.nan)
    adata.uns["qc_metric_basis"] = kind
    adata.var["is_mito"] = mito; adata.var["is_ribo"] = ribo; adata.var["is_hb"] = hb
    # audit DQ-04: record whether the sample's feature annotation carries MT- genes with any signal (per sample; constant within a sample)
    mt_present = bool(mito.any()) and float(M[:, mito].sum()) > 0
    adata.obs["mt_genes_present"] = mt_present
    adata.uns["mt_genes_present"] = mt_present
    return adata


def apply_thresholds(adata: ad.AnnData, cfg: dict, assay: str, chemistry: str | None = None) -> ad.AnnData:
    """Per-sample thresholds; adata must be a single sample. Adds qc_pass, qc_fail_reason, qc_assay_class, qc_max_genes.

    chemistry (manifest column) enables the author-prep-aware rule (configs/qc.yaml author_prep_aware, audit DQ-04): a Smart-seq2
    sample of frozen nuclei uses the nuclei pct_mito_max when MT- genes are present in its annotation; min_genes stays smartseq2."""
    platform = adata.obs["labelle_platform"].iloc[0] if "labelle_platform" in adata.obs else None
    cls = assay_class(assay, platform)
    th = dict(cfg["thresholds"][cls])
    mt_present = bool(adata.uns.get("mt_genes_present", adata.obs["mt_genes_present"].iloc[0] if "mt_genes_present" in adata.obs else True))
    mito_rule = cls
    if cfg.get("author_prep_aware", False) and cls == "smartseq2" and chemistry and "frozen nuclei" in str(chemistry).lower() and mt_present:
        th["pct_mito_max"] = cfg["thresholds"]["nuclei"]["pct_mito_max"]
        mito_rule = "nuclei (frozen nuclei, author_prep_aware)"
    if not mt_present:
        mito_rule = f"{cls} (inoperative: no MT- genes in annotation)"
    q = float(cfg["max_genes_quantile"])
    max_genes = float(np.quantile(adata.obs["n_genes"].values, q)) if adata.n_obs > 0 else np.inf
    ng = adata.obs["n_genes"].values
    pm = adata.obs["pct_mito"].values
    fail_low = ng < th["min_genes"]
    fail_high = ng > max_genes
    fail_mito = ~(pm < th["pct_mito_max"])  # NaN (zero-count cells) fails
    reasons = []
    for a_, b_, c_ in zip(fail_low, fail_high, fail_mito):
        r = []
        if a_: r.append("low_genes")
        if b_: r.append("high_genes")
        if c_: r.append("high_mito")
        reasons.append(";".join(r))
    adata.obs["qc_fail_reason"] = reasons
    adata.obs["qc_pass"] = ~(fail_low | fail_high | fail_mito)
    adata.obs["qc_assay_class"] = cls
    adata.obs["qc_max_genes"] = max_genes
    adata.obs["qc_min_genes"] = th["min_genes"]
    adata.obs["qc_pct_mito_max"] = th["pct_mito_max"]
    adata.obs["qc_mito_rule"] = mito_rule
    return adata


def sample_gate(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    """Sample-level QC flag (audit DQ-03): median n_genes of qc_pass cells >= cfg['sample_gate']['min_median_genes_pass'].
    Adds obs['sample_median_genes_pass'] and obs['sample_qc_pass'] (constant within the sample). Flag only - nothing is dropped;
    Phase 4 (07_assemble_atlas_input.py) decides how to use it."""
    thr = float(cfg.get("sample_gate", {}).get("min_median_genes_pass", 500))
    p = adata.obs.loc[adata.obs["qc_pass"].astype(bool), "n_genes"]
    med = float(np.median(p)) if len(p) else np.nan
    adata.obs["sample_median_genes_pass"] = med
    adata.obs["sample_qc_pass"] = bool(med >= thr) if np.isfinite(med) else False
    adata.obs["sample_gate_min_median_genes"] = thr
    return adata


def run_scrublet(adata: ad.AnnData, cfg: dict, log=print) -> ad.AnnData:
    """Scrublet on the counts layer of ONE 10x sample. Stores doublet_score / predicted_doublet; flag only."""
    import scanpy as sc
    scfg = cfg["scrublet"]
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = pd.array([pd.NA] * adata.n_obs, dtype="boolean")
    adata.uns["scrublet_threshold"] = np.nan
    if "counts" not in adata.layers:
        return adata
    # run on cells passing the minimal gene filter so empty droplets do not distort the simulated doublets
    use = adata.obs["n_genes"].values >= cfg["thresholds"]["cells"]["min_genes"]
    if use.sum() < scfg["min_cells"]:
        log(f"scrublet skipped ({use.sum()} usable cells < {scfg['min_cells']})")
        return adata
    tmp = ad.AnnData(X=sp.csr_matrix(adata.layers["counts"][use]).astype(np.float32))
    # memory guard: scanpy's scrublet densifies HVGs for observed + 2x simulated cells; pre-select 3,000 HVGs
    # (seurat_v3 on counts) so RSS stays ~2 GB for a 10k-cell sample instead of ~5 GB on the full gene set
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.filter_genes(tmp, min_cells=3)
        if tmp.n_vars > scfg.get("n_hvg", 3000):
            try:
                sc.pp.highly_variable_genes(tmp, flavor="seurat_v3", n_top_genes=scfg.get("n_hvg", 3000), subset=True)
            except Exception as e:  # loess can fail on very small samples; fall back to the full gene set
                log(f"scrublet HVG prefilter skipped ({type(e).__name__}); using all {tmp.n_vars} genes")
        thr_src = "auto"
        try:
            sc.pp.scrublet(tmp, expected_doublet_rate=scfg["expected_doublet_rate"], random_state=scfg["random_state"], verbose=False)
            if not np.isfinite(tmp.uns["scrublet"]["threshold"]):
                raise ValueError("non-finite threshold")
        except Exception as e:  # automatic threshold detection can fail on small/odd samples
            log(f"scrublet auto-threshold failed ({type(e).__name__}: {e}); using fixed threshold {scfg['fixed_threshold']}")
            thr_src = "fixed_fallback"
            sc.pp.scrublet(tmp, expected_doublet_rate=scfg["expected_doublet_rate"], random_state=scfg["random_state"],
                           threshold=scfg["fixed_threshold"], verbose=False)
    adata.obs.loc[use, "doublet_score"] = tmp.obs["doublet_score"].values
    pdbl = pd.array([pd.NA] * adata.n_obs, dtype="boolean")
    pdbl[np.where(use)[0]] = tmp.obs["predicted_doublet"].values.astype(bool)
    adata.obs["predicted_doublet"] = pdbl
    fixed = pd.array([pd.NA] * adata.n_obs, dtype="boolean")
    fixed[np.where(use)[0]] = tmp.obs["doublet_score"].values > scfg["fixed_threshold"]
    adata.obs["predicted_doublet_fixed"] = fixed  # score > fixed_threshold (0.25): alternative call independent of the auto threshold
    adata.uns["scrublet_threshold"] = float(tmp.uns["scrublet"]["threshold"])
    adata.uns["scrublet_threshold_source"] = thr_src
    return adata


def per_sample_summary(adata: ad.AnnData) -> dict:
    o = adata.obs
    p = o[o.qc_pass]
    d = {
        "sample_id": o["sample_id"].iloc[0], "study": o["study"].iloc[0], "assay": o["assay"].iloc[0],
        "qc_assay_class": o["qc_assay_class"].iloc[0], "include_in_atlas": bool(o["include_in_atlas"].iloc[0]),
        "cells_before_qc": int(adata.n_obs), "cells_after_qc": int(len(p)),
        "median_genes_before": float(np.median(o["n_genes"])) if len(o) else np.nan,
        "median_genes": float(np.median(p["n_genes"])) if len(p) else np.nan,
        "median_counts": float(np.median(p["n_counts"])) if len(p) else np.nan,
        "median_pct_mito_before": float(np.nanmedian(o["pct_mito"])) if len(o) else np.nan,
        "median_pct_mito": float(np.nanmedian(p["pct_mito"])) if len(p) else np.nan,
        "median_pct_ribo": float(np.nanmedian(p["pct_ribo"])) if len(p) else np.nan,
        "median_pct_hb": float(np.nanmedian(p["pct_hb"])) if len(p) else np.nan,
        "max_genes_threshold": float(o["qc_max_genes"].iloc[0]),
        "n_fail_low_genes": int(o["qc_fail_reason"].str.contains("low_genes").sum()),
        "n_fail_high_genes": int(o["qc_fail_reason"].str.contains("high_genes").sum()),
        "n_fail_high_mito": int(o["qc_fail_reason"].str.contains("high_mito").sum()),
        "scrublet_threshold": adata.uns.get("scrublet_threshold", np.nan),
        "scrublet_threshold_source": adata.uns.get("scrublet_threshold_source", ""),
        "doublet_rate_fixed025": float(p["predicted_doublet_fixed"].astype("boolean").mean()) if "predicted_doublet_fixed" in p and p["predicted_doublet_fixed"].notna().any() else np.nan,
        "doublet_rate_flagged": float(p["predicted_doublet"].astype("boolean").mean()) if p["predicted_doublet"].notna().any() else np.nan,
        "n_doublets_flagged": int(p["predicted_doublet"].astype("boolean").fillna(False).sum()),
        "qc_metric_basis": adata.uns.get("qc_metric_basis", ""),
        "mt_genes_present": bool(o["mt_genes_present"].iloc[0]) if "mt_genes_present" in o else np.nan,
        "qc_mito_rule": o["qc_mito_rule"].iloc[0] if "qc_mito_rule" in o else "",
        "sample_qc_pass": bool(o["sample_qc_pass"].iloc[0]) if "sample_qc_pass" in o else np.nan,
        "sample_gate_min_median_genes": float(o["sample_gate_min_median_genes"].iloc[0]) if "sample_gate_min_median_genes" in o else np.nan,
    }
    if "author_qc_pass" in o:
        d["author_qc_cells"] = int(o["author_qc_pass"].sum())
        d["author_qc_cells_passing_ours"] = int((o["author_qc_pass"] & o["qc_pass"]).sum())
    if "author_Malignant_normal_consensus_Jessa2022" in o:
        mal = o["author_Malignant_normal_consensus_Jessa2022"].astype(str) == "Malignant"
        d["author_malignant_cells"] = int(mal.sum()); d["author_malignant_passing_qc"] = int((mal & o["qc_pass"]).sum())
    return d
