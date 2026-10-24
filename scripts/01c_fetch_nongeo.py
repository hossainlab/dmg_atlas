#!/usr/bin/env python
"""Phase 1 (audit remediation DQ-02): the non-GEO downloads of the tumour track - Collot 2026 (Zenodo 17922773 + Europe PMC
supplement), ScPCA SCPCP000021 (portal API + dataset zip), ScPCA SCPCP000001 API record (DeSisto labels) and Okada 2026
(UCSC Cell Browser 'brain-dmg-multiomics') - as a script with URLs, local paths and sha256 checksums.

The Phase 1 fetches were performed interactively on 2026-10-17 and recorded in data/raw/MANIFEST_tumour.csv (url, bytes,
sha256, download_utc). This script makes them reproducible:

  python scripts/01c_fetch_nongeo.py --dry-run     print the plan (url -> local path, expected size + sha256) and exit
  python scripts/01c_fetch_nongeo.py --verify      re-hash every local file and compare with the manifest (default action)
  python scripts/01c_fetch_nongeo.py --fetch       download files that are MISSING locally, then verify (never overwrites)

The ScPCA dataset zip is served through a pre-signed S3 URL that is minted per download; it cannot be re-fetched from the
manifest URL. Procedure used (ScPCA portal API, https://api.scpca.alexslemonade.org/v1/):
  1. POST /v1/tokens/ {"email": <contact email>, "is_activated": true}   (= acceptance of the portal Terms of Use)
  2. with header "API-KEY: <token>", request the dataset 0d57349a-f81c-4fd2-b32a-ba5069814238 (format ANN_DATA, samples
     SCPCS000636 + SCPCS000637); the response carries a time-limited download_url to scpca-portal-deployer-prod.s3.amazonaws.com
  3. download the zip; sha256 must equal the manifest value below.
--fetch therefore skips the zip unless --scpca-download-url <presigned url> is given (step 2 is done by the user).

Zenodo also publishes md5 for each file (record API https://zenodo.org/api/records/17922773): scRNA.tar.gz md5
560424543ad4d855b122f73ba82128fc, README.md md5 039b0d7f032c3fca12296321d19a462a; --verify checks both md5 and sha256.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = RAW / "MANIFEST_tumour.csv"

CB = "https://cells-test.gi.ucsc.edu/brain-dmg-multiomics"
ZEN = "https://zenodo.org/api/records/17922773/files"
SCPCA_API = "https://api.scpca.alexslemonade.org/v1"

# (study, accession, url, local_path relative to repo, zenodo_md5 or None, notes)
PLAN = [
    ("collot", "zenodo.17922773", f"{ZEN}/README.md/content", "data/raw/collot/README.md", "039b0d7f032c3fca12296321d19a462a", "Zenodo record README"),
    ("collot", "zenodo.17922773", f"{ZEN}/scRNA.tar.gz/content", "data/raw/collot/scRNA.tar.gz", "560424543ad4d855b122f73ba82128fc",
     "human + mouse IUE sc/snRNA; 1.00 GB; unpacked to data/raw/collot/scRNA/"),
    ("collot", "PMC12994606", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12994606/supplementaryFiles?includeInlineImage=false",
     "data/raw/collot/supp/PMC12994606_supplementaryFiles.zip", None, "Europe PMC supplementary bundle for doi 10.1016/j.ccell.2025.12.020 (mmc1-5 unpacked to data/raw/collot/supp/)"),
    ("okada", "cellbrowser:brain-dmg-multiomics", f"{CB}/dataset.json", "data/raw/okada/dataset.json", None, "collection record"),
    ("okada", "cellbrowser:brain-dmg-multiomics", f"{CB}/desc.json", "data/raw/okada/desc.json", None, "collection description"),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/dataset.json", "data/raw/okada/snrna-seq-data/dataset.json", None, ""),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/desc.json", "data/raw/okada/snrna-seq-data/desc.json", None, ""),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/files.txt", "data/raw/okada/snrna-seq-data/files.txt", None, ""),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/features.tsv.gz", "data/raw/okada/snrna-seq-data/features.tsv.gz", None, "2,000 HVG of the exprMatrix"),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/barcodes.tsv.gz", "data/raw/okada/snrna-seq-data/barcodes.tsv.gz", None, ""),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/meta.tsv", "data/raw/okada/snrna-seq-data/meta.tsv", None,
     "per-nucleus metadata; source of per_sample_summary.csv (scripts/01b_build_label_tables.py)"),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/cellbrowser.json.bak", "data/raw/okada/snrna-seq-data/cellbrowser.json.bak", None, ""),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/brain-dmg-snRNA-Homo-sapiens.h5ad",
     "data/raw/okada/snrna-seq-data/brain-dmg-snRNA-Homo-sapiens.h5ad", None, "Cell Browser h5ad export (1.51 GB); raw.X = Seurat LogNormalize(1e4); counts recovered at ingest"),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/matrix.mtx.gz", "data/raw/okada/snrna-seq-data/matrix.mtx.gz", None, "exprMatrix (2000 features, Float32); not used"),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/UMAP.coords.tsv.gz", "data/raw/okada/snrna-seq-data/UMAP.coords.tsv.gz", None, "embedding"),
    ("okada", "cellbrowser:brain-dmg-multiomics/snrna-seq-data", f"{CB}/snrna-seq-data/tMinusSNE.coords.tsv.gz", "data/raw/okada/snrna-seq-data/tMinusSNE.coords.tsv.gz", None, "embedding"),
    ("scpca", "SCPCP000021", f"{SCPCA_API}/projects/SCPCP000021/", "data/raw/scpca/SCPCP000021_api_project.json", None, "ScPCA API project record (JSON)"),
    ("scpca", "SCPCS000636", f"{SCPCA_API}/samples/SCPCS000636/", "data/raw/scpca/SCPCS000636_api_sample.json", None, "ScPCA API sample record (JSON)"),
    ("scpca", "SCPCS000637", f"{SCPCA_API}/samples/SCPCS000637/", "data/raw/scpca/SCPCS000637_api_sample.json", None, "ScPCA API sample record (JSON)"),
    ("scpca", "SCPCP000021", "https://scpca-portal-deployer-prod.s3.amazonaws.com/0d57349a-f81c-4fd2-b32a-ba5069814238.zip",
     "data/raw/scpca/SCPCP000021_SCPCS000636_SCPCS000637_ann-data.zip", None,
     "ScPCA dataset 0d57349a-f81c-4fd2-b32a-ba5069814238 (ANN_DATA); pre-signed URL (token required, see module docstring); unpacked to data/raw/scpca/unpacked/"),
    ("desisto", "SCPCP000001", f"{SCPCA_API}/projects/SCPCP000001/", "data/raw/desisto/supp/SCPCP000001_api_project.json", None,
     "ScPCA API project record (JSON) for the Colorado pHGG cohort = GSE231859; source of desisto_sample_labels.csv (fetched 2026-10-19 during audit remediation)"),
]


def sha256_md5(p: Path):
    h, m = hashlib.sha256(), hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
            m.update(chunk)
    return h.hexdigest(), m.hexdigest()


def download(url: str, dest: Path):
    import requests
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(1 << 22):
                fh.write(chunk)
    tmp.rename(dest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--fetch", action="store_true", help="download files missing locally (never overwrites)")
    ap.add_argument("--scpca-download-url", default=None, help="pre-signed S3 URL for the ScPCA dataset zip (see docstring)")
    ap.add_argument("--report", default=str(ROOT / "results" / "qc" / "nongeo_fetch_verification.csv"))
    args = ap.parse_args()
    man = pd.read_csv(MANIFEST)
    by_path = man.set_index("local_path")

    if args.dry_run:
        for study, acc, url, lp, md5, note in PLAN:
            exp = by_path.loc[lp] if lp in by_path.index else None
            print(f"[{study}] {url}\n      -> {lp}  expected {int(exp.bytes):,} B sha256 {exp.sha256[:16]}..." if exp is not None else f"[{study}] {url}\n      -> {lp}  (not in manifest)")
        return

    rows = []
    for study, acc, url, lp, md5, note in PLAN:
        dest = ROOT / lp
        exp = by_path.loc[lp] if lp in by_path.index else None
        if not dest.exists() and args.fetch:
            u = args.scpca_download_url if "s3.amazonaws.com" in url else url
            if u is None:
                print(f"SKIP {lp}: pre-signed URL required (--scpca-download-url)")
            else:
                print(f"fetching {u} -> {lp}")
                download(u, dest)
        if dest.exists():
            t0 = time.time()
            sha, m5 = sha256_md5(dest)
            size = dest.stat().st_size
            ok_sha = (exp is None) or (sha == exp.sha256)
            ok_size = (exp is None) or (size == int(exp.bytes))
            ok_md5 = (md5 is None) or (m5 == md5)
            status = "OK" if (ok_sha and ok_size and ok_md5) else "MISMATCH"
            if exp is None:
                status = "present (not in manifest)"
        else:
            sha = m5 = ""
            size = 0
            status = "MISSING"
            t0 = time.time()
        rows.append({"study": study, "accession": acc, "url": url, "local_path": lp, "bytes": size, "sha256": sha, "md5": m5,
                     "manifest_bytes": int(exp.bytes) if exp is not None else None, "manifest_sha256": exp.sha256 if exp is not None else None,
                     "zenodo_md5": md5, "status": status, "hash_seconds": round(time.time() - t0, 1), "notes": note})
        print(f"{status:12s} {lp} ({size:,} B)")
    rep = pd.DataFrame(rows)
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    rep.to_csv(args.report, index=False)
    print(f"{(rep.status == 'OK').sum()}/{len(rep)} OK; report -> {Path(args.report).relative_to(ROOT)}")
    bad = rep[~rep.status.isin(["OK", "present (not in manifest)"])]
    if len(bad):
        print("NOT OK:\n" + bad[["local_path", "status"]].to_string(index=False))
        sys.exit(1)


if __name__ == "__main__":
    main()
