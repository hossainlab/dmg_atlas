#!/bin/bash
# Phase 2 driver: one study at a time (thermal rule), 8 threads, nice 10. Usage: bash scripts/run_phase2_all.sh [studies...]
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8 NUMBA_CACHE_DIR=/tmp/numba_cache
STUDIES=${@:-"desisto jessa_scrna jessa_multiome collot okada"}
for s in $STUDIES; do
  echo "=== $s $(date)" >> logs/phase2_driver.log
  nice -n 10 python scripts/03_ingest_qc.py --study $s >> logs/phase2_driver.log 2>&1 || echo "FAILED $s" >> logs/phase2_driver.log
done
echo "=== ALL DONE $(date)" >> logs/phase2_driver.log
