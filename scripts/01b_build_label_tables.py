#!/usr/bin/env python
"""Phase 1 (audit remediation DQ-02): build every per-study label table in data/interim/labels/ from the raw deposits and
supplementary files already on disk, so that scripts/02_build_sample_manifest.py has a scripted upstream.

Inputs (all under data/raw/, see data/raw/MANIFEST_tumour.csv for URL + sha256)
  jessa_scrna/supp/41588_2022_1205_MOESM3_ESM.xlsx   sheets '1-Sample_summary', '3-scRNAseq', '5-scMultiome'
  jessa_scrna/samples/<GSM>_<ID>/*                   barcodes/features/metadata.NG2022 per sample (GSE210567)
  jessa_multiome/<GSM>/*                              idem (GSE210566)
  collot/supp/mmc2.xlsx                              sheet 'A. Patient samples snRNA Seq '
  collot/scRNA/scRNA/Human/scRNA_metadata.csv.gz     per-cell Study / SampleID
  labelle/supp/mmc2.xlsx                             Table S1 (header row 2)
  labelle/GSE227983_{tumor_ss2,immune_ss2,immune_tenx}_meta.csv.gz
  desisto/supp/SCPCP000001_api_project.json          ScPCA portal API project record (samples[].additional_metadata)
  desisto/samples/<GSM>_<name>/                      GEO sample directories (name == ScPCA submitter_id, 905 == 905_3)
  liu/GSE184357_metadata_{fresh,frozen}_postFilter.txt.gz, liu/raw_counts/cm_counts_*.csv (header line only)
  okada/snrna-seq-data/meta.tsv

Outputs (data/interim/labels/ unless --out):
  jessa_sample_labels.csv, jessa_scrna_geo_sample_summary.csv, jessa_multiome_geo_sample_summary.csv,
  collot_sample_labels.csv, collot_sample_summary.csv, labelle_sample_labels.csv, labelle_geo_metadata_sample_summary.csv,
  desisto_sample_labels.csv, liu_geo_metadata_sample_summary.csv, liu_rawfile_to_sample.csv,
  and data/raw/okada/snrna-seq-data/per_sample_summary.csv (written next to meta.tsv, where 02_build_sample_manifest.py reads it).

Usage
  python scripts/01b_build_label_tables.py --out /tmp/labels --verify     # rebuild elsewhere and diff against data/interim/labels
  python scripts/01b_build_label_tables.py                                # write in place (existing files kept as *.pre_audit.csv when they differ)

Known, documented difference to the pre-audit files (2026-10-19): liu_rawfile_to_sample.csv columns `unmatched` and
`matched_by_normalised_cellid` are now computed with the final ingest regex (dmg_atlas.ingest.liu_harmonise_id); the
pre-audit table had been produced with an earlier, stricter ID normalisation (dash required before the well) and under-counted
matches (28 files). The per-sample sum of the new column equals the 8,984 postFilter cells recovered at ingest (Phase 2 notes).
`file` and `metadata_sample` - the only columns consumed downstream - are unchanged. All other tables reproduce byte-identically.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.ingest import liu_harmonise_id  # noqa: E402

RAW = ROOT / "data" / "raw"
LAB = ROOT / "data" / "interim" / "labels"

JESSA_XLSX = RAW / "jessa_scrna" / "supp" / "41588_2022_1205_MOESM3_ESM.xlsx"
JESSA_SRC = ("Jessa2022 NatGenet 10.1038/s41588-022-01205-w 41588_2022_1205_MOESM3_ESM.xlsx sheet 1-Sample_summary "
             "(cols Age,Sex,Location,Group,Histone_status,scRNAseq) + sheets 3-scRNAseq/5-scMultiome (ID membership, Kit)")
COLLOT_SRC = ('Collot2026 CancerCell 10.1016/j.ccell.2025.12.020 mmc2.xlsx sheet "A. Patient samples snRNA Seq" '
              "(cols Paper ID,Study,Tumor_subtype,Location,Gender,Age_years,Primary/ recurrent)")
LABELLE_SRC = ("LaBelle2025 CellRepMed 10.1016/j.xcrm.2025.102095 mmc2.xlsx Table S1 "
               "(cols Sample,Subtype,Age (y),Sex,Location,Seq type,Treatment before surgery)")
DESISTO_SRC = ("ScPCA portal API SCPCP000001 sample metadata (same Colorado cohort as GSE231859; submitter_id matches GEO sample name): "
               "diagnosis, tissue_location, age, molecular_characteristics, treatment, disease_timing")
LIU_MAP_BASIS = "file name (cell-ID formats differ between raw files and postFilter metadata; harmonise at ingest)"

# Jessa 2022 Table S1 'Group' values that are diffuse midline gliomas with H3 K27M (atlas scope). H3.3K27M-Cortex is
# hemispheric and PFA-H3.1K27M-PF is an ependymoma: both are out-group.
JESSA_DMG_GROUPS = {"HGG-H3.1/2K27M-Pons", "HGG-H3.3K27M-Pons", "HGG-H3.3K27M-Thal."}


# ------------------------------------------------------------------------------------------------------------ helpers
def count_lines_gz(p: Path) -> int:
    n = 0
    with gzip.open(p, "rb") as fh:
        for _ in fh:
            n += 1
    return n


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------------------ Jessa
def jessa_h3_variant(histone_status, group) -> str:
    hs = str(histone_status)
    if "G34" in str(group) or "G34" in hs:
        return "H3.3_G34"
    if hs.startswith("HIST1H3"):
        return "H3.1"
    if hs.startswith("HIST2H3"):
        return "H3.2"
    if hs.startswith("H3F3A"):
        return "H3.3"
    if hs == "WT":
        return "WT"
    return "unknown"


def jessa_sample_labels() -> pd.DataFrame:
    s1 = pd.read_excel(JESSA_XLSX, sheet_name="1-Sample_summary")
    s3 = pd.read_excel(JESSA_XLSX, sheet_name="3-scRNAseq").set_index("ID")
    s5 = pd.read_excel(JESSA_XLSX, sheet_name="5-scMultiome").set_index("ID")
    ids3, ids5 = set(s3.index.dropna().astype(str)), set(s5.index.dropna().astype(str))
    keep = ["ID", "Material", "Age", "Sex", "Location", "Group", "Histone_status", "ACVR1_status", "TP53_status",
            "scRNAseq", "scMultiome", "scATACseq"]
    out = s1[s1.ID.astype(str).isin(ids3 | ids5)][keep].copy()
    out["in_GSE210567"] = out.ID.isin(ids3)
    out["in_GSE210566"] = out.ID.isin(ids5)
    out["h3_variant"] = [jessa_h3_variant(h, g) for h, g in zip(out.Histone_status, out.Group)]
    out["site"] = out.Location.map(lambda L: {"Pons": "pons", "Thal.": "thalamus"}.get(L, "other"))
    out["assay_scrna"] = out.scRNAseq.map(lambda s: "snRNA" if "Nuclei" in str(s) else ("scRNA" if "Cell" in str(s) else np.nan))
    out["patient_id"] = out.ID.str.split("_").str[0]
    out["kit_version"] = out.ID.map(s3["Kit"])
    out["n_cells_after_filtering_S3"] = out.ID.map(s3["N_cells_after_filtering"])
    out["nuclei_prep_multiome"] = out.ID.map(s5["Nuclei_prep_method"])
    out["is_dmg"] = out.Group.isin(JESSA_DMG_GROUPS)
    out["label_source"] = JESSA_SRC
    return out.reset_index(drop=True)


def _jessa_geo_summary(sample_dirs) -> list[dict]:
    recs = []
    for d in sample_dirs:
        mtx = glob.glob(f"{d}/*matrix.mtx.gz")
        if not mtx:
            continue
        pref = os.path.basename(mtx[0]).replace("_matrix.mtx.gz", "")
        md = pd.read_csv(f"{d}/{pref}_metadata.NG2022.tsv.gz", sep="\t")
        recs.append({"sample": pref, "features": count_lines_gz(Path(f"{d}/{pref}_features.tsv.gz")),
                     "barcodes": count_lines_gz(Path(f"{d}/{pref}_barcodes.tsv.gz")), "qc_cells": len(md),
                     "technology": md.Technology.iloc[0], "id_paper": md.ID_paper.iloc[0],
                     "n_malig": int((md.Malignant_normal_consensus_Jessa2022 == "Malignant").sum()),
                     "med_mito": round(float(md["percent.mito"].median()), 2)})
    return recs


def jessa_scrna_geo_sample_summary() -> pd.DataFrame:
    recs = _jessa_geo_summary(sorted(glob.glob(str(RAW / "jessa_scrna" / "samples" / "*"))))
    df = pd.DataFrame(recs).rename(columns={"features": "genes"})
    return df[["sample", "genes", "barcodes", "qc_cells", "technology", "id_paper", "n_malig", "med_mito"]]


def jessa_multiome_geo_sample_summary() -> pd.DataFrame:
    recs = _jessa_geo_summary(sorted(d for d in glob.glob(str(RAW / "jessa_multiome" / "GSM*")) if os.path.isdir(d)))
    df = pd.DataFrame(recs)
    df["gsm"] = df["sample"].str.split("_").str[0]
    df["sample"] = df["sample"].str.split("_", n=1).str[1]
    return df[["gsm", "sample", "features", "barcodes", "qc_cells", "technology", "n_malig", "med_mito"]]


# ------------------------------------------------------------------------------------------------------------ Collot
def collot_sample_labels() -> pd.DataFrame:
    df = pd.read_excel(RAW / "collot" / "supp" / "mmc2.xlsx", sheet_name="A. Patient samples snRNA Seq ")
    df["label_source"] = COLLOT_SRC
    return df


def collot_sample_summary() -> pd.DataFrame:
    cm = pd.read_csv(RAW / "collot" / "scRNA" / "scRNA" / "Human" / "scRNA_metadata.csv.gz", usecols=["Study", "SampleID"])
    return cm.groupby(["Study", "SampleID"]).size().reset_index(name="n_cells")


# ------------------------------------------------------------------------------------------------------------ LaBelle
LABELLE_RENAME = {"Sample": "Sample", "Broad subtype": "broad_subtype", "Subtype": "subtype", "Age (y)": "age_years", "Sex": "sex",
                  "Location": "location", "Seq type": "seq_type", "Immunotherapy Treatment": "immunotherapy",
                  "Treatment before surgery": "treatment_before_surgery", "Steroid Treatment": "steroid",
                  "Tumor mutational Burden": "tmb", "Tumor cells sorted": "tumor_sorted", "Tcells sorted": "tcells_sorted",
                  "Myeloid cells sorted": "myeloid_sorted", "VDJ sequencing (10X only)": "vdj"}


def labelle_sample_labels() -> pd.DataFrame:
    df = pd.read_excel(RAW / "labelle" / "supp" / "mmc2.xlsx", sheet_name=0, header=1).rename(columns=LABELLE_RENAME)
    df["label_source"] = LABELLE_SRC
    return df


def labelle_geo_metadata_sample_summary() -> pd.DataFrame:
    recs = []
    for k in ["tumor_ss2", "immune_ss2", "immune_tenx"]:
        m = pd.read_csv(RAW / "labelle" / f"GSE227983_{k}_meta.csv.gz", index_col=0)
        g = m.groupby("sampleid").agg(Subtype=("Subtype", "first"), Age=("Age", "first"), Sex=("Sex", "first"),
                                     Location=("Location", "first"), n=("Subtype", "size")).reset_index()
        g["file"] = k
        recs.append(g)
    return pd.concat(recs, ignore_index=True)


# ------------------------------------------------------------------------------------------------------------ DeSisto (ScPCA SCPCP000001)
SCPCA_API_PROJECT = "https://api.scpca.alexslemonade.org/v1/projects/SCPCP000001/"
SCPCA_API_JSON = RAW / "desisto" / "supp" / "SCPCP000001_api_project.json"


def desisto_k27m_status(molecular_characteristics) -> str:
    s = str(molecular_characteristics)
    if re.search(r"(No|Negative for)[^.;]*?(H3K27M|H3F3A)", s):
        return "WT"
    if "K27M" in s:
        return "K27M"
    return "unknown"


def desisto_site(tissue_location) -> str:
    L = str(tissue_location).lower()
    if "thalam" in L:
        return "thalamus"
    if "pons" in L:
        return "pons"
    if "spine" in L or "spinal" in L:
        return "spinal"
    return "other"


def desisto_sample_labels(fetch: bool = False) -> pd.DataFrame:
    if fetch or not SCPCA_API_JSON.exists():
        import requests
        r = requests.get(SCPCA_API_PROJECT, timeout=60)
        r.raise_for_status()
        SCPCA_API_JSON.parent.mkdir(parents=True, exist_ok=True)
        SCPCA_API_JSON.write_bytes(r.content)
        print(f"fetched {SCPCA_API_PROJECT} -> {SCPCA_API_JSON} ({len(r.content)} bytes, sha256 {sha256(SCPCA_API_JSON)})")
    proj = json.loads(SCPCA_API_JSON.read_text())
    geo_names = {os.path.basename(d).split("_", 1)[1] for d in glob.glob(str(RAW / "desisto" / "samples" / "*"))}
    recs = []
    for s in proj["samples"]:
        am = dict(s.get("additional_metadata") or {})
        submitter = am.pop("submitter_id", None)
        participant = am.pop("participant_id", None)
        geo = {"905_3": "905"}.get(submitter, submitter)  # ScPCA 905_3 is GEO sample 905
        if geo not in geo_names:
            continue
        recs.append({"scpca_id": s["scpca_id"], "submitter_id": submitter, "participant_id": participant,
                     "diagnosis": s["diagnosis"], "subdiagnosis": s["subdiagnosis"], "tissue_location": s["tissue_location"],
                     "age": s["age"], "sex": s["sex"], "disease_timing": s["disease_timing"],
                     "seq_units": ";".join(s["seq_units"]), "technologies": ";".join(s["technologies"]),
                     "cells": s["sample_cell_count_estimate"], "treatment": s["treatment"], "extra": json.dumps(am),
                     "molecular_characteristics": am.get("molecular_characteristics"), "outcome": am.get("outcome"),
                     "geo_sample": geo, "k27m_status": desisto_k27m_status(am.get("molecular_characteristics")),
                     "site": desisto_site(s["tissue_location"]), "label_source": DESISTO_SRC})
    # the portal's literal "NA" strings are kept verbatim (02_build_sample_manifest.py reads the table with dtype=str)
    return pd.DataFrame(recs).sort_values("scpca_id").reset_index(drop=True)


# ------------------------------------------------------------------------------------------------------------ Liu 2022
# Author annotations that are not malignant tumour states. 'n_malig' below counts annotations containing 'like'
# (OPC-like-1/2/3, AC-like, MES-like, OC-like, frozen-OPC-like-b); cycling malignant cells ('G2M', 'S') are NOT counted -
# this reproduces the pre-audit table and the column is informational only (it enters only the manifest `notes` string).
LIU_FILE_OVERRIDES = {"MGH66": "MGH6666410", "DIPG50aut1": "SUDIPG50Aut", "DIPG50bio": "SUDIPG50Bio", "AAA076530001": "AAA007653"}


def _liu_meta() -> pd.DataFrame:
    fresh = pd.read_csv(RAW / "liu" / "GSE184357_metadata_fresh_postFilter.txt.gz", sep="\t")
    frozen = pd.read_csv(RAW / "liu" / "GSE184357_metadata_frozen_postFilter.txt.gz", sep="\t")
    fresh["prep"] = "fresh"
    frozen["prep"] = "frozen"
    meta = pd.concat([fresh, frozen], ignore_index=True)
    meta.columns = [c.strip() for c in meta.columns]
    meta["cell"] = meta["cell"].str.strip()
    return meta


def liu_geo_metadata_sample_summary() -> pd.DataFrame:
    meta = _liu_meta()
    return meta.groupby("sample").agg(prep=("prep", "first"), age=("age", "first"), location=("location", "first"),
                                      n_cells=("cell", "size"),
                                      n_malig=("annotation", lambda s: int(s.str.contains("like").sum()))).reset_index()


def liu_file_to_sample(stem: str, samples: set) -> str | None:
    """'cm_counts_BCH836_2' -> 'BT836'; 'cm_counts_MUV017' -> 'MUV17'; 'cm_counts_GBM27' -> 'ICGCGBM27'; 'cm_counts_SUDIPG55biop' -> 'SUDIPG55'."""
    s = stem.replace("cm_counts_", "")
    s = re.sub(r"_2$", "", s)
    s = re.sub(r"^BCH", "BT", s)
    s = re.sub(r"^GBM", "ICGCGBM", s)
    s = re.sub(r"^MUV0+", "MUV", s)
    s = re.sub(r"biop$", "", s)
    s = LIU_FILE_OVERRIDES.get(s, s)
    return s if s in samples else None


def liu_rawfile_to_sample() -> pd.DataFrame:
    meta = _liu_meta()
    samples, cells = set(meta["sample"]), set(meta["cell"])
    recs = []
    for fn in sorted(os.listdir(RAW / "liu" / "raw_counts")):
        with open(RAW / "liu" / "raw_counts" / fn) as fh:
            cols = [c.strip('"') for c in fh.readline().rstrip("\n").split(",")[1:]]
        ms = liu_file_to_sample(fn[:-4], samples)
        matched = len({liu_harmonise_id(c, ms) for c in cols} & cells) if ms else 0
        recs.append({"file": fn, "n_cols": len(cols), "example": cols[0], "unmatched": len(cols) - matched,
                     "metadata_sample": ms or "(no postFilter metadata sample)", "mapping_basis": LIU_MAP_BASIS,
                     "matched_by_normalised_cellid": matched})
    return pd.DataFrame(recs)


# ------------------------------------------------------------------------------------------------------------ Okada (Cell Browser meta.tsv)
def okada_per_sample_summary() -> pd.DataFrame:
    mt = pd.read_csv(RAW / "okada" / "snrna-seq-data" / "meta.tsv", sep="\t",
                     usecols=["orig.ident", "sample", "location", "treatement", "age", "sex", "nCount_RNA"])
    return mt.groupby("sample").agg(orig=("orig.ident", "first"), location=("location", "first"), treatment=("treatement", "first"),
                                    age=("age", "first"), sex=("sex", "first"), n_cells=("nCount_RNA", "size"),
                                    median_umi=("nCount_RNA", "median")).reset_index()


# ------------------------------------------------------------------------------------------------------------ driver
BUILDERS = {
    "jessa_sample_labels.csv": jessa_sample_labels,
    "jessa_scrna_geo_sample_summary.csv": jessa_scrna_geo_sample_summary,
    "jessa_multiome_geo_sample_summary.csv": jessa_multiome_geo_sample_summary,
    "collot_sample_labels.csv": collot_sample_labels,
    "collot_sample_summary.csv": collot_sample_summary,
    "labelle_sample_labels.csv": labelle_sample_labels,
    "labelle_geo_metadata_sample_summary.csv": labelle_geo_metadata_sample_summary,
    "desisto_sample_labels.csv": desisto_sample_labels,
    "liu_geo_metadata_sample_summary.csv": liu_geo_metadata_sample_summary,
    "liu_rawfile_to_sample.csv": liu_rawfile_to_sample,
}
OKADA_OUT = RAW / "okada" / "snrna-seq-data" / "per_sample_summary.csv"


def cell_diff(old: pd.DataFrame, new: pd.DataFrame) -> str:
    if list(old.columns) != list(new.columns):
        return f"columns differ: {sorted(set(old.columns) ^ set(new.columns))}"
    if old.shape != new.shape:
        return f"shape {old.shape} vs {new.shape}"
    parts = []
    for c in old.columns:
        a, b = old[c].astype(str), new[c].astype(str)
        n = int((a.values != b.values).sum())
        if n:
            parts.append(f"{c}: {n} cells")
    return "; ".join(parts) if parts else "identical values"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="output directory (default data/interim/labels)")
    ap.add_argument("--verify", action="store_true", help="compare each rebuilt table with the file in data/interim/labels (and the Okada summary)")
    ap.add_argument("--fetch-scpca", action="store_true", help="re-fetch the ScPCA API project record instead of using the local copy")
    args = ap.parse_args()
    out = Path(args.out) if args.out else LAB
    out.mkdir(parents=True, exist_ok=True)
    report = []
    for name, fn in list(BUILDERS.items()) + [("per_sample_summary.csv", okada_per_sample_summary)]:
        df = fn(fetch=args.fetch_scpca) if name == "desisto_sample_labels.csv" else fn()
        target = (OKADA_OUT if args.out is None else out / name) if name == "per_sample_summary.csv" else out / name
        existing = OKADA_OUT if name == "per_sample_summary.csv" else LAB / name
        tmp = target.with_suffix(".tmp.csv")
        df.to_csv(tmp, index=False)
        status = "new"
        if existing.exists():
            same = sha256(tmp) == sha256(existing)
            status = "byte-identical" if same else "DIFFERS (" + cell_diff(pd.read_csv(existing), pd.read_csv(tmp)) + ")"
            if not same and target.resolve() == existing.resolve():
                shutil.copy2(existing, existing.with_suffix(".pre_audit.csv"))
        os.replace(tmp, target)
        report.append({"table": name, "rows": len(df), "cols": df.shape[1], "written_to": str(target.relative_to(ROOT)) if target.is_relative_to(ROOT) else str(target),
                       "vs_existing": status, "sha256": sha256(target)})
        print(f"{name:45s} {len(df):4d} rows  {status}")
    rep = pd.DataFrame(report)
    if args.verify:
        vpath = ROOT / "results" / "qc" / "label_tables_verification.csv"
        vpath.parent.mkdir(parents=True, exist_ok=True)
        rep.to_csv(vpath, index=False)
        print(f"verification table -> {vpath.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
