#!/usr/bin/env python
"""Phase 1 (tumour track): fetch GEO deposits of the DMG tumour cohorts.

Sequential queue of files; each file >20 MB is fetched as 4 parallel byte-range chunks
(<= 4 streams at any time). Final size is checked against Content-Length, sha256 is
recorded in data/raw/MANIFEST_tumour.csv, wall time per file in logs/phase1_tumour.log.
Re-running skips files already marked '<file>.manifested'.

Usage: python scripts/01_fetch_tumour_geo.py [queue_name ...]
  queues: jessa_scrna jessa_multiome desisto zang filbin liu labelle labelle2025
"""
import csv, datetime, hashlib, os, re, sys, threading, time
from concurrent.futures import ThreadPoolExecutor
import requests

import pathlib as _pl
ROOT = str(_pl.Path(__file__).resolve().parents[1])  # audit DQ-14: no hard-coded machine path
RAW = f"{ROOT}/data/raw"
LOG = f"{ROOT}/logs/phase1_tumour.log"
MAN = f"{RAW}/MANIFEST_tumour.csv"
GEO_SERIES = "https://ftp.ncbi.nlm.nih.gov/geo/series/"
GEO_SAMPLES = "https://ftp.ncbi.nlm.nih.gov/geo/samples/"
NCHUNKS = 4

S = requests.Session()
S.headers["Accept-Encoding"] = "identity"
man_lock = threading.Lock()


def log(msg):
    line = f"{datetime.datetime.utcnow().isoformat(timespec='seconds')}Z  {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def manifest_add(study, acc, url, path, notes=""):
    b = os.path.getsize(path)
    h = sha256(path)
    new = not os.path.exists(MAN)
    with man_lock, open(MAN, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["study", "accession", "url", "local_path", "bytes", "sha256", "download_utc", "notes"])
        w.writerow([study, acc, url, os.path.relpath(path, ROOT), b, h,
                    datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z", notes])
    return b, h


def head(url):
    r = S.head(url, allow_redirects=True, timeout=60)
    r.raise_for_status()
    return int(r.headers.get("Content-Length", 0) or 0)


def dl_range(url, dest, start, end, retries=6):
    """Resumable range download of bytes start..end (inclusive) into dest."""
    want = end - start + 1
    for attempt in range(retries):
        have = os.path.getsize(dest) if os.path.exists(dest) else 0
        if have >= want:
            return
        try:
            with S.get(url, headers={"Range": f"bytes={start + have}-{end}"}, stream=True, timeout=180) as r:
                if r.status_code != 206:
                    raise RuntimeError(f"status {r.status_code}")
                with open(dest, "ab") as f:
                    for ch in r.iter_content(1 << 20):
                        f.write(ch)
        except Exception as e:  # noqa
            log(f"  retry {attempt + 1} {os.path.basename(dest)}: {e}")
            time.sleep(10 * (attempt + 1))
    if os.path.getsize(dest) != want:
        raise RuntimeError(f"chunk incomplete {dest}")


def download(url, path, study, acc, notes="", nchunks=NCHUNKS):
    if os.path.exists(path + ".manifested"):
        log(f"SKIP (done) {os.path.basename(path)}")
        return
    t0 = time.time()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    size = head(url)
    if size > 20_000_000 and nchunks > 1:
        bounds = [(i * size // nchunks, (i + 1) * size // nchunks - 1) for i in range(nchunks)]
        parts = [f"{path}.part{i}" for i in range(nchunks)]
        with ThreadPoolExecutor(nchunks) as ex:
            list(ex.map(lambda a: dl_range(url, a[1], *a[0]), zip(bounds, parts)))
        with open(path, "wb") as out:
            for p in parts:
                with open(p, "rb") as f:
                    for ch in iter(lambda: f.read(1 << 22), b""):
                        out.write(ch)
                os.remove(p)
    else:
        with S.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(path, "wb") as f:
                for ch in r.iter_content(1 << 20):
                    f.write(ch)
    got = os.path.getsize(path)
    if size and got != size:
        raise RuntimeError(f"size mismatch {got} != {size} for {url}")
    b, h = manifest_add(study, acc, url, path, notes)
    open(path + ".manifested", "w").write(h)
    log(f"DL {study} {os.path.basename(path)} {b} bytes {time.time() - t0:.0f}s "
        f"{b / 1e6 / max(time.time() - t0, 1):.1f} MB/s sha256={h[:12]}")


def series_url(gse, fn):
    return f"{GEO_SERIES}{gse[:-3]}nnn/{gse}/suppl/{fn}"


def sample_url(gsm, fn):
    return f"{GEO_SAMPLES}{gsm[:-3]}nnn/{gsm}/suppl/{fn}"


def multiome_rna_files():
    """RNA/GEX matrices + metadata for the 10 GSMs of GSE210566 that carry them (ATAC bigwigs skipped)."""
    out = []
    with open(f"{RAW}/jessa_multiome/GSE210566_filelist.txt") as fh:
        for line in fh:
            if not line.startswith("File"):
                continue
            name = line.split("\t")[1]
            if name.endswith(".ATAC.bw"):
                continue
            gsm = name.split("_")[0]
            out.append((gsm, name))
    return out


QUEUES = {
    "jessa_scrna": lambda: [(series_url("GSE210567", "GSE210567_RAW.tar"), f"{RAW}/jessa_scrna/GSE210567_RAW.tar",
                             "jessa_scrna", "GSE210567", "37 samples, per-sample MTX + metadata.NG2022.tsv.gz")],
    "desisto": lambda: [(series_url("GSE231859", "GSE231859_RAW.tar"), f"{RAW}/desisto/GSE231859_RAW.tar",
                         "desisto", "GSE231859", "19 samples, per-sample MTX")],
    "jessa_multiome": lambda: [(sample_url(gsm, fn), f"{RAW}/jessa_multiome/{gsm}/{fn}", "jessa_multiome", gsm,
                                "scMultiome GEX part; ATAC bigwigs skipped") for gsm, fn in multiome_rna_files()],
    "zang": lambda: [(sample_url("GSM7846769", fn), f"{RAW}/zang/{fn}", "zang", "GSM7846769",
                      "DMG sample of GSE245674; GSM7846770 (ependymoma) skipped")
                     for fn in ["GSM7846769_barcodes.tsv.gz", "GSM7846769_features.tsv.gz", "GSM7846769_matrix1.mtx.gz"]],
    "filbin": lambda: [(series_url("GSE102130", "GSE102130_K27Mproject.RSEM.vh20170621.txt.gz"),
                        f"{RAW}/filbin/GSE102130_K27Mproject.RSEM.vh20170621.txt.gz", "filbin", "GSE102130",
                        "Smart-seq2 RSEM TPM, all cells")],
    "liu": lambda: [(series_url("GSE184357", fn), f"{RAW}/liu/{fn}", "liu", "GSE184357",
                     "scATAC tar and pre/postFilter count tables not downloaded")
                    for fn in ["GSE184357_DMG_raw_counts.tar.gz", "GSE184357_metadata_fresh_postFilter.txt.gz",
                               "GSE184357_metadata_frozen_postFilter.txt.gz"]],
    "labelle": lambda: [(series_url("GSE227983", fn), f"{RAW}/labelle/{fn}", "labelle", "GSE227983", "processed TPM only")
                        for fn in ["GSE227983_tumor_ss2_tpm.csv.gz", "GSE227983_tumor_ss2_meta.csv.gz",
                                   "GSE227983_immune_ss2_tpm.csv.gz", "GSE227983_immune_ss2_meta.csv.gz",
                                   "GSE227983_immune_tenx_tpm.csv.gz", "GSE227983_immune_tenx_meta.csv.gz"]],
    # GSE271936 is Xenium in situ (4 regions). Only the H3K27M region's cell x gene MTX is fetched; images/transcripts skipped.
    "labelle2025": lambda: [(sample_url("GSM8389375", fn), f"{RAW}/labelle/GSE271936/{fn}", "labelle", "GSM8389375",
                             "Xenium in situ cell x gene matrix (H3K27M_P6_Region_1); morphology/transcripts skipped")
                            for fn in ["GSM8389375_H3K27M_P6_Region_1-barcodes.tsv.gz",
                                       "GSM8389375_H3K27M_P6_Region_1-features.tsv.gz",
                                       "GSM8389375_H3K27M_P6_Region_1-matrix.mtx.gz"]],
}

if __name__ == "__main__":
    names = sys.argv[1:] or list(QUEUES)
    for q in names:
        log(f"QUEUE {q} start")
        t0 = time.time()
        for url, path, study, acc, notes in QUEUES[q]():
            try:
                download(url, path, study, acc, notes)
            except Exception as e:  # noqa
                log(f"FAIL {study} {os.path.basename(path)}: {e}")
        log(f"QUEUE {q} done in {time.time() - t0:.0f}s")
