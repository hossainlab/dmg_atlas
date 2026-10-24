#!/usr/bin/env python
"""Parallel byte-range downloader with size/sha256 verification, manifest and timing log.

Usage: python scripts/dl_ranged.py URL DEST --study S --accession A [--chunks 4] [--notes '...'] [--md5 expected]
Resumable: existing complete part files are kept; partial parts are resumed by Range.
"""
import argparse, csv, hashlib, os, sys, time, threading, datetime
import requests

import pathlib as _pl
ROOT = str(_pl.Path(__file__).resolve().parents[1])  # audit DQ-14: no hard-coded machine path
MANIFEST = f"{ROOT}/data/raw/MANIFEST_spatial_refs.csv"
LOG = f"{ROOT}/logs/phase1_spatial_refs.log"
COLS = ["study","accession","url","local_path","bytes","sha256","download_utc","notes"]

def log(msg):
    line = f"{datetime.datetime.utcnow().isoformat(timespec='seconds')}Z  {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh: fh.write(line + "\n")

def head_size(url):
    r = requests.head(url, allow_redirects=True, timeout=60)
    r.raise_for_status()
    return int(r.headers["Content-Length"]), r.headers.get("Accept-Ranges","") == "bytes"

def fetch_part(url, start, end, path, progress, idx):
    have = os.path.getsize(path) if os.path.exists(path) else 0
    want = end - start + 1
    while have < want:
        try:
            hdr = {"Range": f"bytes={start+have}-{end}"}
            with requests.get(url, headers=hdr, stream=True, timeout=120) as r:
                if r.status_code not in (200, 206):
                    raise RuntimeError(f"status {r.status_code}")
                if r.status_code == 200 and have > 0:
                    raise RuntimeError("server ignored Range")
                with open(path, "ab") as fh:
                    for blk in r.iter_content(1 << 20):
                        fh.write(blk); have += len(blk); progress[idx] = have
        except Exception as e:
            log(f"  part{idx} retry after error: {e}")
            time.sleep(5)
    progress[idx] = have

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""): h.update(blk)
    return h.hexdigest()

def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""): h.update(blk)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("url"); ap.add_argument("dest")
    ap.add_argument("--study", required=True); ap.add_argument("--accession", required=True)
    ap.add_argument("--chunks", type=int, default=4); ap.add_argument("--notes", default="")
    ap.add_argument("--md5", default=None)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.dest), exist_ok=True)
    t0 = time.time()
    size, ranges = head_size(a.url)
    log(f"START {a.study} {os.path.basename(a.dest)} size={size} ranges={ranges} chunks={a.chunks}")
    if os.path.exists(a.dest) and os.path.getsize(a.dest) == size:
        log(f"  already complete: {a.dest}")
    else:
        n = a.chunks if ranges and size > 8 << 20 else 1
        bounds = [(i*size//n, (i+1)*size//n - 1) for i in range(n)]
        p0 = f"{a.dest}.part0"
        if n > 1 and os.path.exists(p0) and os.path.getsize(p0) > bounds[0][1] + 1:
            # an earlier single-stream attempt left a long prefix: keep it as chunk 0, split the tail over the other chunks
            P = os.path.getsize(p0); m = n - 1
            bounds = [(0, P - 1)] + [(P + i*(size-P)//m, P + (i+1)*(size-P)//m - 1) for i in range(m)]
            log(f"  resuming from single-stream prefix of {P} bytes; tail split into {m} chunks")
        parts = [f"{a.dest}.part{i}" for i in range(n)]
        progress = [0]*n
        th = [threading.Thread(target=fetch_part, args=(a.url, s, e, p, progress, i), daemon=True)
              for i,((s,e),p) in enumerate(zip(bounds, parts))]
        for t in th: t.start()
        last = time.time(); lastb = 0
        while any(t.is_alive() for t in th):
            time.sleep(30)
            done = sum(progress)
            rate = (done-lastb)/max(time.time()-last,1)/1e6
            last, lastb = time.time(), done
            log(f"  {a.study}: {done/1e6:.1f}/{size/1e6:.1f} MB  {rate:.2f} MB/s")
        with open(a.dest, "wb") as out:
            for (s_, e_), p in zip(bounds, parts):
                left = e_ - s_ + 1   # copy exactly the chunk's byte range (a part may have grown past it)
                with open(p, "rb") as fh:
                    while left > 0:
                        blk = fh.read(min(1 << 24, left)); out.write(blk); left -= len(blk)
        for p in parts: os.remove(p)
    got = os.path.getsize(a.dest)
    if got != size:
        log(f"FAIL {a.study}: size mismatch {got} != {size}"); sys.exit(2)
    sha = sha256_file(a.dest)
    note = a.notes
    if a.md5:
        m = md5_file(a.dest)
        ok = m == a.md5
        note += f"; md5 {'OK' if ok else 'MISMATCH'} ({m})"
        if not ok: log(f"FAIL {a.study}: md5 mismatch"); sys.exit(3)
    new = not os.path.exists(MANIFEST)
    with open(MANIFEST, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        if new: w.writeheader()
        w.writerow(dict(study=a.study, accession=a.accession, url=a.url,
                        local_path=os.path.relpath(a.dest, ROOT), bytes=got, sha256=sha,
                        download_utc=datetime.datetime.utcnow().isoformat(timespec="seconds")+"Z", notes=note))
    log(f"DONE {a.study} {os.path.basename(a.dest)} bytes={got} sha256={sha[:12]}.. wall={time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
