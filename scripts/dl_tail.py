#!/usr/bin/env python
"""Download explicit byte ranges of a URL to separate files (used to finish a throttled Zenodo transfer in parallel).
Usage: dl_tail.py URL OUTPREFIX START END NCHUNKS   -> writes OUTPREFIX.tail{i} for i in range(NCHUNKS)."""
import sys, os, threading, time
sys.path.insert(0, os.path.dirname(__file__))
from dl_ranged import fetch_part, log
url, pref, start, end, n = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
bounds = [(start + i*(end-start+1)//n, start + (i+1)*(end-start+1)//n - 1) for i in range(n)]
progress = [0]*n
th = [threading.Thread(target=fetch_part, args=(url, s, e, f"{pref}.tail{i}", progress, i), daemon=True) for i,(s,e) in enumerate(bounds)]
log(f"TAIL start {os.path.basename(pref)} ranges={bounds}")
for t in th: t.start()
while any(t.is_alive() for t in th):
    time.sleep(60); log(f"  tail progress: {[round(p/1e6,1) for p in progress]} MB of {[round((e-s+1)/1e6,1) for s,e in bounds]}")
log(f"TAIL done {os.path.basename(pref)} sizes={[os.path.getsize(f'{pref}.tail{i}') for i in range(n)]}")
