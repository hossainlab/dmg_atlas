#!/usr/bin/env python
"""Subset Braun 2023 'All' h5ad (cellxgene) to hindbrain/pons/medulla cells in a chunked, backed fashion.
Reads CSR row blocks from the source file and appends them to resizable datasets in the output h5ad;
never holds more than one block (~25k cells) in memory. Output follows the anndata on-disk spec (CSR X)."""
import sys, time, datetime, h5py, numpy as np, pandas as pd
from anndata.io import read_elem, write_elem

import pathlib as _pl
ROOT = str(_pl.Path(__file__).resolve().parents[1])  # audit DQ-14: no hard-coded machine path
SRC = f"{ROOT}/data/external/braun2023/braun2023_all_52ff91e6.h5ad"
DST = f"{ROOT}/data/external/braun2023/braun2023_hindbrain_pons_medulla.h5ad"
KEEP_REGIONS = ["Hindbrain", "Pons", "Medulla"]
BLOCK = 25000
LOG = f"{ROOT}/logs/phase1_spatial_refs.log"
def log(m):
    line = f"{datetime.datetime.utcnow().isoformat(timespec='seconds')}Z  {m}"; print(line, flush=True)
    open(LOG, "a").write(line + "\n")

t0 = time.time()
with h5py.File(SRC, "r") as src, h5py.File(DST, "w") as dst:
    obs = read_elem(src["obs"]); var = read_elem(src["var"])
    keep = np.where(obs["Region"].isin(KEEP_REGIONS).to_numpy())[0]
    log(f"braun subset: keeping {len(keep)} of {len(obs)} cells from regions {KEEP_REGIONS}")
    indptr = src["X"]["indptr"][:]
    n_genes = var.shape[0]
    # contiguous runs of selected rows
    breaks = np.where(np.diff(keep) != 1)[0] + 1
    runs = np.split(keep, breaks)
    log(f"  {len(runs)} contiguous row runs")
    Xg = dst.create_group("X")
    Xg.attrs["encoding-type"] = "csr_matrix"; Xg.attrs["encoding-version"] = "0.1.0"
    Xg.attrs["shape"] = np.array([len(keep), n_genes], dtype=np.int64)
    d_data = Xg.create_dataset("data", shape=(0,), maxshape=(None,), dtype=src["X"]["data"].dtype, chunks=(1<<20,), compression="gzip", compression_opts=4)
    d_ind = Xg.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(1<<20,), compression="gzip", compression_opts=4)
    new_indptr = [0]; total = 0; done = 0
    for run in runs:
        for s in range(0, len(run), BLOCK):
            rows = run[s:s+BLOCK]; a, b = indptr[rows[0]], indptr[rows[-1]+1]
            data = src["X"]["data"][a:b]; ind = src["X"]["indices"][a:b]
            d_data.resize((total + len(data),)); d_data[total:] = data
            d_ind.resize((total + len(ind),)); d_ind[total:] = ind.astype(np.int32)
            ptr = indptr[rows[0]:rows[-1]+2] - a + total
            new_indptr.extend(ptr[1:].tolist()); total += len(data); done += len(rows)
            if done % 100000 < BLOCK: log(f"  {done}/{len(keep)} cells, nnz={total}, {time.time()-t0:.0f}s")
    Xg.create_dataset("indptr", data=np.asarray(new_indptr, dtype=np.int64))
    assert len(new_indptr) == len(keep) + 1
    write_elem(dst, "obs", obs.iloc[keep].copy())
    write_elem(dst, "var", var)
    obsm = dst.create_group("obsm"); obsm.attrs["encoding-type"] = "dict"; obsm.attrs["encoding-version"] = "0.1.0"
    for k in src["obsm"].keys():
        arr = src["obsm"][k]
        if isinstance(arr, h5py.Dataset): write_elem(obsm, k, arr[:][keep])
    write_elem(dst, "uns", read_elem(src["uns"]))
    for g in ["layers", "obsp", "varm", "varp"]:
        grp = dst.create_group(g); grp.attrs["encoding-type"] = "dict"; grp.attrs["encoding-version"] = "0.1.0"
    dst.attrs["encoding-type"] = "anndata"; dst.attrs["encoding-version"] = "0.1.0"
    log(f"braun subset DONE: {len(keep)} cells x {n_genes} genes, nnz={total}, wall={time.time()-t0:.0f}s")
