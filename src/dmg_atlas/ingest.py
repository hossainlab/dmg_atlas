"""Per-study loaders -> AnnData with .layers['counts'] (UMI / RSEM counts), .layers['tpm'] (Filbin linear TPM) or
.layers['log2_tpm10'] (LaBelle: GSE227983 tables are log2(TPM/10+1), audit DQ-01). Every loader records uns['value_scale'].

Every loader is a generator `iter_<study>(rows)` yielding (manifest_row, AnnData) one sample at a time so
that QC/scrublet run per sample and memory stays bounded. Conventions:
  * obs index = f"{sample_id}:{original_cell_id}"
  * manifest columns are joined into obs (MANIFEST_OBS_COLS)
  * author per-cell annotations are copied to obs with prefix 'author_' and are never modified
  * matrices are scipy CSR, int32 for integer counts, float32 for TPM/RSEM expected counts
  * var index = gene identifier as deposited; var['ensembl_id'] where the deposit provides it
"""
from __future__ import annotations

import gzip
import io
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import anndata as ad
import h5py

from . import RAW, INTERIM, PROCESSED

MANIFEST_CSV = PROCESSED / "sample_manifest.csv"
MANIFEST_OBS_COLS = ["sample_id", "patient_id", "study", "role", "assay", "chemistry", "site", "age_years", "age_group",
                     "h3_variant", "k27m_status", "sex", "treatment", "tissue_timing", "include_in_atlas",
                     "exclusion_reason", "counts_type", "label_source"]


def load_manifest() -> pd.DataFrame:
    m = pd.read_csv(MANIFEST_CSV)
    m["include_in_atlas"] = m["include_in_atlas"].astype(str).str.lower().eq("true")
    return m


def select_rows(m: pd.DataFrame, study: str) -> pd.DataFrame:
    """Rows to process for a study: included integrated/projection rows + excluded integrated-role out-group rows.
    Collot's excluded rows are re-processed Jessa copies (duplicates), not an out-group, and are skipped."""
    s = m[m.study == study]
    inc = s[s.role.isin(["integrated", "projection"]) & s.include_in_atlas]
    out = s[(s.role == "integrated") & ~s.include_in_atlas]
    if study == "collot":
        out = out.iloc[0:0]
    return pd.concat([inc, out])


def _attach_manifest(adata: ad.AnnData, row: pd.Series, cell_ids) -> ad.AnnData:
    for c in MANIFEST_OBS_COLS:
        v = row.get(c)
        adata.obs[c] = v if not (isinstance(v, float) and np.isnan(v)) else np.nan
    adata.obs["original_cell_id"] = np.asarray(cell_ids).astype(str)
    adata.obs_names = [f"{row.sample_id}:{c}" for c in adata.obs["original_cell_id"]]
    adata.obs_names_make_unique()
    return adata


def _author(adata: ad.AnnData, df: pd.DataFrame, cols=None, prefix="author_"):
    """Join author per-cell table (indexed like adata.obs['original_cell_id']) with prefix."""
    cols = list(df.columns) if cols is None else cols
    df = df[cols].copy()
    df.columns = [prefix + re.sub(r"[^0-9A-Za-z_.]", "_", c) for c in df.columns]
    df = df[~df.index.duplicated()]
    j = df.reindex(adata.obs["original_cell_id"].values)
    j.index = adata.obs_names
    for c in j.columns:
        adata.obs[c] = j[c].values
    return adata


def _to_csr_int(M, allow_float=False):
    M = sp.csr_matrix(M)
    if not allow_float:
        M.data = np.rint(M.data).astype(np.int32)
    else:
        M.data = M.data.astype(np.float32)
    M.eliminate_zeros(); M.sort_indices()
    return M


def read_10x_mtx_dir(d: Path, prefix: str = "") -> ad.AnnData:
    """Read a cellranger MTX triplet (gz) -> AnnData cells x genes, int32 CSR, var index=symbol, var['ensembl_id']."""
    d = Path(d)
    f_mtx = next(d.glob(f"{prefix}*matrix*.mtx.gz"))
    f_bc = next(d.glob(f"{prefix}*barcodes.tsv.gz"))
    f_ft = next(iter(list(d.glob(f"{prefix}*features.tsv.gz")) + list(d.glob(f"{prefix}*genes.tsv.gz"))))
    with gzip.open(f_mtx, "rb") as fh:
        M = scipy.io.mmread(fh)
    feats = pd.read_csv(f_ft, sep="\t", header=None, dtype=str)
    bcs = pd.read_csv(f_bc, sep="\t", header=None, dtype=str)[0].values
    X = _to_csr_int(M.T)
    var = pd.DataFrame({"ensembl_id": feats[0].values, "symbol": feats[1].values}, index=feats[1].values)
    if feats.shape[1] >= 3:
        var["feature_type"] = feats[2].values
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=bcs), var=var)
    if "feature_type" in a.var:
        a = a[:, a.var.feature_type == "Gene Expression"].copy()
    a.var_names = a.var["symbol"].values
    a.var_names_make_unique()
    a.layers["counts"] = a.X
    return a


# ----------------------------------------------------------------------------- Jessa 2022 scRNA / multiome
def _iter_jessa(rows, series_dir: Path, per_sample_subdir: bool):
    for _, row in rows.iterrows():
        d = RAW.parent.parent / row.file_paths if not Path(row.file_paths).is_absolute() else Path(row.file_paths)
        a = read_10x_mtx_dir(d)
        meta_f = next(d.glob("*metadata.NG2022.tsv.gz"))
        meta = pd.read_csv(meta_f, sep="\t", dtype=str).set_index("Cell_barcode")
        a = _attach_manifest(a, row, a.obs_names)
        a = _author(a, meta)
        a.obs["author_qc_pass"] = a.obs["original_cell_id"].isin(meta.index).values
        for c in a.obs.columns:
            if c.startswith("author_") and (c.startswith("author_cNMF") or c.startswith("author_nCount") or c.startswith("author_nFeature")
                                             or c.startswith("author_percent") or "score" in c.lower() or c.startswith("author_G1") or c.startswith("author_G2")):
                a.obs[c] = pd.to_numeric(a.obs[c], errors="coerce")
        yield row, a


def iter_jessa_scrna(rows):
    yield from _iter_jessa(rows, RAW / "jessa_scrna" / "samples", True)


def iter_jessa_multiome(rows):
    yield from _iter_jessa(rows, RAW / "jessa_multiome", True)


# ----------------------------------------------------------------------------- DeSisto, Zang (plain 10x MTX)
def iter_desisto(rows):
    for _, row in rows.iterrows():
        d = RAW.parent.parent / row.file_paths
        a = read_10x_mtx_dir(d)
        a = _attach_manifest(a, row, a.obs_names)
        yield row, a


def iter_zang(rows):
    for _, row in rows.iterrows():
        a = read_10x_mtx_dir(RAW / "zang", prefix="GSM7846769")
        a = _attach_manifest(a, row, a.obs_names)
        yield row, a


# ----------------------------------------------------------------------------- Collot 2025 (Zenodo; streamed MTX)
COLLOT_CACHE = INTERIM / "collot_this_study_counts.h5ad"


def _collot_stream(keep_cells: np.ndarray, n_genes: int, n_cells: int, mtx_gz: Path, chunksize=20_000_000) -> sp.csr_matrix:
    """Stream the 292M-entry MTX and keep only the requested columns (cells). Returns cells x genes CSR int32."""
    col_map = -np.ones(n_cells + 1, dtype=np.int64)
    col_map[keep_cells + 1] = np.arange(len(keep_cells))  # MTX is 1-based
    parts = []
    with gzip.open(mtx_gz, "rt") as fh:
        # skip header/comment lines and the size line
        pos_line = fh.readline()
        while pos_line.startswith("%"):
            pos_line = fh.readline()
        for chunk in pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], dtype={"g": np.int32, "c": np.int32, "v": np.int32},
                                 chunksize=chunksize, engine="c"):
            cc = col_map[chunk["c"].values]
            msk = cc >= 0
            if msk.any():
                parts.append(sp.csr_matrix((chunk["v"].values[msk], (cc[msk], chunk["g"].values[msk] - 1)),
                                           shape=(len(keep_cells), n_genes), dtype=np.int32))
    M = parts[0]
    for p in parts[1:]:
        M = M + p
    return _to_csr_int(M)


def build_collot_cache(log=print) -> ad.AnnData:
    if COLLOT_CACHE.exists():
        return ad.read_h5ad(COLLOT_CACHE)
    d = RAW / "collot" / "scRNA" / "scRNA" / "Human"
    genes = pd.read_csv(d / "scRNA_genes.tsv.gz", sep="\t", header=None, dtype=str)[0].values
    cells = pd.read_csv(d / "scRNA_cells.tsv.gz", sep="\t", header=None, dtype=str)[0].values
    meta = pd.read_csv(d / "scRNA_metadata.csv.gz", dtype=str).set_index("cell_id")
    with gzip.open(d / "scRNA_counts.mtx.gz", "rt") as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_genes, n_cells, nnz = (int(x) for x in line.split())
    assert n_genes == len(genes) and n_cells == len(cells), (n_genes, len(genes), n_cells, len(cells))
    is_new = meta.reindex(cells)["Study"].eq("This_study").values
    keep = np.where(is_new)[0]
    log(f"collot: streaming {nnz:,} entries, keeping {len(keep):,}/{n_cells:,} cells from Study==This_study")
    X = _collot_stream(keep, n_genes, n_cells, d / "scRNA_counts.mtx.gz")
    a = ad.AnnData(X=X, obs=meta.reindex(cells[keep]), var=pd.DataFrame(index=genes))
    a.layers["counts"] = a.X
    INTERIM.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(COLLOT_CACHE)
    return a


def iter_collot(rows, log=print):
    cache = build_collot_cache(log)
    for _, row in rows.iterrows():
        sid = row.sample_id.replace("collot__", "")
        sub = cache[cache.obs["SampleID"] == sid].copy()
        if sub.n_obs == 0:
            raise ValueError(f"collot sample {sid} not found; available: {cache.obs.SampleID.unique()}")
        meta = sub.obs.copy()
        sub.obs = pd.DataFrame(index=sub.obs_names)
        sub = _attach_manifest(sub, row, meta.index.values)
        sub = _author(sub, meta.set_index(meta.index.astype(str)), cols=["SampleID", "Study", "celltype"])
        yield row, sub


# ----------------------------------------------------------------------------- ScPCA (AnnData, filtered)
def iter_scpca(rows):
    for _, row in rows.iterrows():
        d = RAW.parent.parent / row.file_paths
        f = next(Path(d).glob("*_filtered_rna.h5ad"))
        a0 = ad.read_h5ad(f)
        X = _to_csr_int(a0.X)
        var = pd.DataFrame({"ensembl_id": a0.var_names.values, "symbol": a0.var["gene_symbol"].astype(str).values}, index=a0.var_names.values)
        a = ad.AnnData(X=X, obs=pd.DataFrame(index=a0.obs_names), var=var)
        a.var_names = a.var["symbol"].values
        a.var_names_make_unique()
        a.layers["counts"] = a.X
        a = _attach_manifest(a, row, a0.obs_names)
        cols = [c for c in ["scDblFinder_class", "scDblFinder_score", "miQC_pass", "prob_compromised", "scpca_filter",
                            "subsets_mito_percent", "sum", "detected", "library_id", "diagnosis", "subdiagnosis"] if c in a0.obs]
        a = _author(a, a0.obs, cols=cols)
        yield row, a


# ----------------------------------------------------------------------------- Okada 2026 (UCSC Cell Browser h5ad, chunked)
OKADA_H5AD = RAW / "okada" / "snrna-seq-data" / "brain-dmg-snRNA-Homo-sapiens.h5ad"


def _h5_series(g, key):
    it = g[key]
    if isinstance(it, h5py.Group):
        return pd.Categorical.from_codes(it["codes"][:], it["categories"][:].astype(str))
    v = it[:]
    return v.astype(str) if v.dtype.kind in "OS" else v


def okada_obs() -> pd.DataFrame:
    with h5py.File(OKADA_H5AD) as f:
        g = f["obs"]
        idx = g[g.attrs["_index"]][:].astype(str)
        cols = [c for c in g.attrs["column-order"]]
        df = pd.DataFrame({c: _h5_series(g, c) for c in cols}, index=idx)
    return df


def iter_okada(rows, log=print):
    """Read raw.X (LogNormalize 1e4, CSR) one sample at a time by contiguous row block; recover UMI = round(expm1(x)*nCount/1e4)."""
    obs = okada_obs()
    with h5py.File(OKADA_H5AD) as f:
        genes = f["raw/var/_index"][:].astype(str)
        indptr = f["raw/X/indptr"][:]
        for _, row in rows.iterrows():
            sname = row.sample_id.replace("okada__", "")
            idx = np.where(obs["sample"].astype(str).values == sname)[0]
            if len(idx) == 0:
                raise ValueError(f"okada sample {sname} not in obs")
            assert np.all(np.diff(idx) == 1), f"okada sample {sname} rows are not contiguous"
            r0, r1 = idx[0], idx[-1] + 1
            p0, p1 = indptr[r0], indptr[r1]
            data = f["raw/X/data"][p0:p1]
            indices = f["raw/X/indices"][p0:p1]
            ip = indptr[r0:r1 + 1] - p0
            ncount = obs["nCount_RNA"].values[r0:r1].astype(np.float64)
            rows_rep = np.repeat(np.arange(r1 - r0), np.diff(ip))
            raw = np.expm1(data) * ncount[rows_rep] / 1e4
            cnt = np.rint(raw)
            max_dev = float(np.abs(raw - cnt).max())
            X = sp.csr_matrix((cnt.astype(np.int32), indices.astype(np.int32), ip.astype(np.int64)), shape=(r1 - r0, len(genes)))
            rs = np.asarray(X.sum(axis=1)).ravel()
            rs_dev = float(np.abs(rs - ncount).max())
            log(f"okada {sname}: {r1 - r0} nuclei; max |x-round(x)| = {max_dev:.2e}; max |rowsum-nCount_RNA| = {rs_dev:.3g}")
            if max_dev > 1e-6 or rs_dev > 1:
                raise ValueError(f"okada {sname}: count recovery failed (max_dev={max_dev}, rs_dev={rs_dev})")
            a = ad.AnnData(X=X, obs=pd.DataFrame(index=obs.index[r0:r1]), var=pd.DataFrame(index=genes))
            a.var_names_make_unique()
            a.layers["counts"] = a.X
            a = _attach_manifest(a, row, obs.index[r0:r1])
            sub = obs.iloc[r0:r1].copy()
            sub.index = sub.index.astype(str)
            a = _author(a, sub, cols=[c for c in sub.columns if c not in ("_index",)])
            a.uns["okada_count_recovery"] = {"max_abs_dev_from_integer": max_dev, "max_abs_rowsum_minus_nCount": rs_dev}
            yield row, a


# ----------------------------------------------------------------------------- Filbin 2018 (Smart-seq2 TPM)
def iter_filbin(rows):
    f = RAW / "filbin" / "GSE102130_K27Mproject.RSEM.vh20170621.txt.gz"
    df = pd.read_csv(f, sep="\t", index_col=0)
    prefix = pd.Series([c.split("-")[0] for c in df.columns], index=df.columns)
    for _, row in rows.iterrows():
        cols = prefix.index[prefix == row.patient_id]
        if len(cols) == 0:
            raise ValueError(f"filbin: no columns for {row.patient_id}")
        X = sp.csr_matrix(df[cols].values.T.astype(np.float32))
        a = ad.AnnData(X=X, obs=pd.DataFrame(index=cols), var=pd.DataFrame(index=df.index.astype(str)))
        a.var_names_make_unique()
        a.layers["tpm"] = a.X
        a.uns["value_scale"] = "linear TPM (RSEM; per-cell sum 1e6)"
        a = _attach_manifest(a, row, cols)
        a.obs["author_plate"] = [c.split("-")[1] if "-" in c else np.nan for c in cols]
        yield row, a


# ----------------------------------------------------------------------------- Liu 2022 (Smart-seq2 RSEM expected counts)
LIU_ID_RE = re.compile(r"[-_]P(\d+)[-_]?([A-H]\d{2})")


def liu_harmonise_id(raw_id: str, metadata_sample: str):
    """'BT_1126-P01-A01' + 'BT1126' -> 'BT1126.P1.A01' (format of the postFilter metadata)."""
    m = LIU_ID_RE.search(raw_id)
    return f"{metadata_sample}.P{int(m.group(1))}.{m.group(2)}" if m else None


def iter_liu(rows, log=print):
    fresh = pd.read_csv(RAW / "liu" / "GSE184357_metadata_fresh_postFilter.txt.gz", sep="\t")
    frozen = pd.read_csv(RAW / "liu" / "GSE184357_metadata_frozen_postFilter.txt.gz", sep="\t")
    fresh["prep"] = "fresh_cells"; frozen["prep"] = "frozen_nuclei"
    meta = pd.concat([fresh, frozen]); meta.columns = [c.strip() for c in meta.columns]
    meta["cell"] = meta["cell"].str.strip(); meta = meta.set_index("cell")
    fmap = pd.read_csv(INTERIM / "labels" / "liu_rawfile_to_sample.csv")
    for _, row in rows.iterrows():
        msample = row.sample_id.replace("liu__", "")
        files = fmap[fmap.metadata_sample == msample].file.tolist()
        if not files:
            log(f"liu {msample}: no raw count file (author metadata only) - skipped")
            continue
        mats, ids, raw_ids, srcs = [], [], [], []
        genes = None
        for fn in files:
            df = pd.read_csv(RAW / "liu" / "raw_counts" / fn, index_col=0)
            if genes is None:
                genes = df.index.astype(str)
            else:
                df = df.reindex(genes)
            mats.append(sp.csr_matrix(df.values.T.astype(np.float32)))
            raw_ids += list(df.columns); srcs += [fn] * df.shape[1]
            ids += [liu_harmonise_id(c, msample) for c in df.columns]
        X = sp.vstack(mats).tocsr()
        ids = pd.Series(ids, dtype=object)
        # duplicates of a harmonised ID (e.g. MUV82Nuc-P3 vs MUV82-P3): keep the column with the larger library size
        lib = np.asarray(X.sum(axis=1)).ravel()
        order = np.lexsort((-lib, ids.fillna("NA").values))
        keep_mask = np.ones(len(ids), dtype=bool)
        seen = set()
        for i in order:
            k = ids.iloc[i]
            if k is None:
                continue
            if k in seen:
                keep_mask[i] = False
            seen.add(k)
        n_dup = int((~keep_mask).sum())
        X = X[keep_mask]; ids = ids[keep_mask].reset_index(drop=True)
        raw_ids = list(np.array(raw_ids)[keep_mask]); srcs = list(np.array(srcs)[keep_mask])
        cell_ids = [i if i is not None else r for i, r in zip(ids, raw_ids)]
        a = ad.AnnData(X=_to_csr_int(X, allow_float=True), obs=pd.DataFrame(index=cell_ids), var=pd.DataFrame(index=genes))
        a.var_names_make_unique()
        a.layers["counts"] = a.X
        a = _attach_manifest(a, row, cell_ids)
        a.obs["raw_cell_id"] = raw_ids; a.obs["raw_file"] = srcs
        a = _author(a, meta, cols=["sample", "age", "location", "annotation", "prep"])
        a.obs["author_postfilter"] = a.obs["original_cell_id"].isin(meta.index).values
        n_meta = int((meta["sample"] == msample).sum())
        log(f"liu {msample}: {len(files)} file(s), {a.n_obs} raw cells ({n_dup} duplicate IDs dropped), {int(a.obs.author_postfilter.sum())}/{n_meta} author postFilter cells recovered")
        a.uns["liu_recovery"] = {"n_meta": n_meta, "n_recovered": int(a.obs.author_postfilter.sum()), "n_dup_dropped": n_dup}
        yield row, a


# ----------------------------------------------------------------------------- LaBelle 2025 (Smart-seq2 TPM + 10x CPM; GSE227983)
def _labelle_tables():
    d = RAW / "labelle"
    out = {}
    for key in ["tumor_ss2", "immune_ss2", "immune_tenx"]:
        meta = pd.read_csv(d / f"GSE227983_{key}_meta.csv.gz", index_col=0)
        out[key] = meta
    return out


def _stream_csv_columns(path: Path, wanted: set, dtype=np.float32):
    """Stream a dense genes x cells CSV (gz) keeping only `wanted` columns; returns (genes, cells, csr cells x genes).
    pandas on the 24k-column LaBelle 10x table is pathologically slow; a line parser does it in < 1 min."""
    import operator
    rows_data = []
    genes = []
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split(",")
        idx = [i for i, c in enumerate(hdr) if c in wanted]
        cells = [hdr[i] for i in idx]
        if not idx:
            return [], [], None
        getter = operator.itemgetter(*idx) if len(idx) > 1 else (lambda p: (p[idx[0]],))
        for line in fh:
            parts = line.rstrip("\n").split(",")
            genes.append(parts[0].strip('"'))
            rows_data.append(sp.csr_matrix(np.array(getter(parts), dtype=dtype)))
    M = sp.vstack(rows_data).T.tocsr()  # cells x genes
    return genes, cells, M


def iter_labelle(rows, log=print):
    d = RAW / "labelle"
    metas = _labelle_tables()
    want = set(rows.sample_id.str.replace("labelle__", "", regex=False))
    tabs = {}
    for key, meta in metas.items():
        cells = set(meta.index[meta["sampleid"].isin(want)])
        if not cells:
            continue
        genes, cols, M = _stream_csv_columns(d / f"GSE227983_{key}_tpm.csv.gz", cells)
        tabs[key] = (pd.Index(genes), pd.Index(cols), M)
        log(f"labelle {key}: {len(cols)} cells x {len(genes)} genes read")
    for _, row in rows.iterrows():
        sid = row.sample_id.replace("labelle__", "")
        parts = []
        for key, (genes, cols, M) in tabs.items():
            meta = metas[key]
            sel = np.array([c in meta.index and meta.loc[c, "sampleid"] == sid for c in cols])
            if not sel.any():
                continue
            a = ad.AnnData(X=M[sel], obs=pd.DataFrame(index=cols[sel]), var=pd.DataFrame(index=genes))
            a.var_names_make_unique()
            a.obs["original_cell_id"] = cols[sel].values
            a.obs["labelle_table"] = key
            a.obs["labelle_platform"] = "smartseq2" if key.endswith("ss2") else "10x"
            mcols = [c for c in meta.columns if c not in ("sampleid",)]
            a = _author(a, meta.loc[cols[sel]], cols=mcols)
            parts.append(a)
        if not parts:
            raise ValueError(f"labelle {sid}: no cells found")
        # audit DQ-09: tumour and immune tables of one sample have different gene lists (14,642 vs 19,357); zero-filled union
        # (outer join) instead of the previous inner join, which discarded 4,736 measured immune genes
        table_genes = {p_.obs["labelle_table"].iloc[0]: sorted(p_.var_names) for p_ in parts}
        a = parts[0] if len(parts) == 1 else ad.concat(parts, join="outer", merge="first", fill_value=0)
        a.uns["labelle_table_genes"] = table_genes  # measured genes per source table (03_ingest_qc.py turns them into var flags)
        # audit DQ-01: all three GSE227983 tables (tumor_ss2, immune_ss2, immune_tenx) are log2(TPM/10+1), not linear TPM
        # (per-cell maxima 7.6-9.0, column sums 3.7k-12.6k); the layer is named accordingly and the scale recorded in uns
        a.layers["log2_tpm10"] = sp.csr_matrix(a.X)
        a.uns["value_scale"] = "log2(TPM/10+1) as deposited in GSE227983 (SS2 and 10x tables); linear TPM = (2**x - 1) * 10; UMI counts not recoverable"
        orig = a.obs_names.values
        keep_obs = a.obs.copy()
        a = _attach_manifest(a, row, orig)
        for c in keep_obs.columns:
            a.obs[c] = keep_obs[c].values
        yield row, a


LOADERS = {
    "jessa_scrna": iter_jessa_scrna, "jessa_multiome": iter_jessa_multiome, "collot": iter_collot, "desisto": iter_desisto,
    "zang": iter_zang, "scpca": iter_scpca, "okada": iter_okada, "filbin": iter_filbin, "liu": iter_liu, "labelle": iter_labelle,
}
