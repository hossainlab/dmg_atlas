"""Thread caps, timing and memory logging shared by all phase scripts.

Import this module FIRST (before numpy/scanpy) so the thread caps take effect.
"""
from __future__ import annotations

import os
import time
import resource
import sys
from contextlib import contextmanager
from pathlib import Path

N_THREADS = os.environ.get("DMG_ATLAS_THREADS", "8")
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, N_THREADS)
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")

try:  # pandas >= 3 (Modal image) infers Arrow-backed strings, which anndata 0.11/0.12 cannot write to h5ad
    import pandas as _pd
    if int(_pd.__version__.split(".")[0]) >= 3:
        _pd.set_option("future.infer_string", False)
except Exception:
    pass


def peak_rss_gb() -> float:
    """Peak resident set size of this process in GB (macOS reports bytes, Linux KB)."""
    v = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return v / 1e9
    return v * 1024 / 1e9


def current_rss_gb() -> float:
    try:
        import psutil
        return psutil.Process().memory_info().rss / 1e9
    except Exception:
        return float("nan")


class StepLogger:
    """Append-only text log: one line per step with wall time and peak RSS."""

    def __init__(self, path: Path | str, context: str = ""):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.context = context

    def write(self, msg: str):
        line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} [{self.context}] {msg}"
        with open(self.path, "a") as fh:
            fh.write(line + "\n")
        print(line, flush=True)

    @contextmanager
    def step(self, name: str):
        t0 = time.time()
        self.write(f"START {name}")
        try:
            yield
        finally:
            self.write(f"END   {name} | wall {time.time() - t0:.1f}s | rss_now {current_rss_gb():.2f} GB | peak_rss {peak_rss_gb():.2f} GB")
