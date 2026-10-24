"""Gene-space harmonisation to a single symbol set.

Rule (documented in docs/phase2_qc_notes.md):
  1. A reference table is built from the union of the 10x `features.tsv` files present in data/raw
     plus the ScPCA var table. Each row = (ensembl_id without version, symbol, source, priority).
     Sources in priority order (first listed wins; the Zang GRCh38-2020-A table is the most complete 10x annotation among the deposits and therefore ranks first — this is a completeness rule, not a 'newest' rule):
        1 zang          GRCh38-2020-A  (36,601 genes, Ensembl 98 / GENCODE v32)
        2 jessa_multiome cellranger-arc GRCh38-2020-A-2.0.0 GEX rows (60,658; versioned IDs with "_N" suffix)
        3 scpca         Ensembl 104 (60,319)
        4 desisto       GRCh38-3.0.0 (33,538; Ensembl 93)
        5 jessa_scrna   hg19-1.2.0 (32,738; Ensembl 75)
  2. The ATLAS GENE UNIVERSE is the set of symbols of the GRCh38-2020-A reference (source 1), because
     the two largest UMI studies (Jessa multiome superset, Okada) use it and it is the current
     cellranger default. The canonical symbol of an Ensembl ID is its symbol in the highest-priority
     source that contains it.
  3. Mapping a study's gene identifiers to the universe:
        a. Ensembl IDs available (10x features / ScPCA): map by versionless Ensembl ID -> canonical symbol.
        b. Symbols only (Collot, Okada, Filbin, Liu, LaBelle): exact match to a universe symbol; else
           alias lookup = the symbol appears in an older reference for an Ensembl ID whose canonical
           symbol is in the universe (e.g. hg19 'FAM138A' -> 2020-A 'FAM138A'; 'MIR1302-10' -> 'MIR1302-2HG').
           A symbol whose alias lookup resolves to more than one canonical symbol is ambiguous and dropped.
        c. Unmapped genes are dropped and counted; genes that collapse onto the same canonical symbol
           (e.g. PAR_Y copies, hg19 duplicates) are summed.
  4. Each study's output var is the universe subset it measured (no zero-filling of unmeasured genes),
     in universe order, with var['ensembl_id'] where known and var['n_source_features'] (>1 = collapsed).
Genes lost per study are appended to results/qc/gene_harmonisation.csv.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from . import RAW, INTERIM

REFERENCE_CSV = INTERIM / "gene_reference.csv"
UNIVERSE_CSV = INTERIM / "gene_universe.csv"

_SOURCES = [
    # (name, priority, path, columns)
    ("zang", 1, RAW / "zang" / "GSM7846769_features.tsv.gz"),
    ("jessa_multiome", 2, RAW / "jessa_multiome" / "GSM6432689" / "GSM6432689_P-6253_S-8498_features.tsv.gz"),
    ("scpca", 3, RAW / "scpca" / "unpacked" / "SCPCP000021_single-cell" / "SCPCS000636" / "SCPCL000964_filtered_rna.h5ad"),
    ("desisto", 4, RAW / "desisto" / "samples" / "GSM7305260_834" / "GSM7305260_834_features.tsv.gz"),
    ("jessa_scrna", 5, RAW / "jessa_scrna" / "samples" / "GSM6432726_P-6251_S-8496" / "GSM6432726_P-6251_S-8496_features.tsv.gz"),
]


def strip_version(x: str) -> str:
    """ENSG00000223972.5_2 -> ENSG00000223972"""
    return re.sub(r"\..*$", "", x)


def _read_features(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    df = df.iloc[:, :3] if df.shape[1] >= 3 else df
    df.columns = ["ensembl_id", "symbol", "feature_type"][: df.shape[1]]
    if "feature_type" in df:
        df = df[df.feature_type == "Gene Expression"]
    return df[["ensembl_id", "symbol"]]


def build_reference(force: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (reference, universe). Cached in data/interim."""
    if REFERENCE_CSV.exists() and UNIVERSE_CSV.exists() and not force:
        return pd.read_csv(REFERENCE_CSV), pd.read_csv(UNIVERSE_CSV)
    parts = []
    for name, prio, path in _SOURCES:
        if name == "scpca":
            import h5py
            with h5py.File(path) as f:
                ens = f["var/_index"][:].astype(str)
                sym = f["var/gene_symbol"]
                if isinstance(sym, h5py.Group):  # categorical
                    sym = sym["categories"][:].astype(str)[sym["codes"][:]]
                else:
                    sym = sym[:].astype(str)
            df = pd.DataFrame({"ensembl_id": ens, "symbol": sym})
            df = df[df.symbol.notna() & (df.symbol != "") & (df.symbol != "nan")]
        else:
            df = _read_features(path)
        df = df.assign(ensembl_id=df.ensembl_id.map(strip_version), source=name, priority=prio)
        parts.append(df)
    ref = pd.concat(parts, ignore_index=True)
    ref = ref[ref.ensembl_id.str.startswith("ENSG")]
    # canonical symbol per Ensembl ID = symbol in highest-priority source
    canon = (ref.sort_values("priority").drop_duplicates("ensembl_id")[["ensembl_id", "symbol"]]
             .rename(columns={"symbol": "canonical_symbol"}))
    ref = ref.merge(canon, on="ensembl_id", how="left")
    uni = ref[ref.source == "zang"][["ensembl_id", "symbol"]].rename(columns={"symbol": "canonical_symbol"})
    uni = uni.drop_duplicates("canonical_symbol").reset_index(drop=True)  # 2020-A has a handful of duplicated symbols
    uni["universe_rank"] = np.arange(len(uni))
    INTERIM.mkdir(parents=True, exist_ok=True)
    ref.to_csv(REFERENCE_CSV, index=False)
    uni.to_csv(UNIVERSE_CSV, index=False)
    return ref, uni


class GeneMapper:
    def __init__(self):
        self.ref, self.uni = build_reference()
        self.universe = set(self.uni.canonical_symbol)
        self.uni_order = {s: i for i, s in enumerate(self.uni.canonical_symbol)}
        self.ens2sym = dict(zip(self.uni.ensembl_id, self.uni.canonical_symbol))
        self.sym2ens = dict(zip(self.uni.canonical_symbol, self.uni.ensembl_id))
        # Ensembl IDs outside 2020-A but present in other references -> their canonical symbol if it is in the universe
        extra = self.ref[~self.ref.ensembl_id.isin(self.ens2sym)].drop_duplicates("ensembl_id")
        for e, s in zip(extra.ensembl_id, extra.canonical_symbol):
            if s in self.universe:
                self.ens2sym.setdefault(e, s)
        # alias table: symbol (any source) -> set of canonical symbols in universe
        # symbols that are duplicated WITHIN a source (e.g. hg19 '5S_rRNA' x many, ScPCA placeholder 'snoZ196')
        # denote several loci and cannot serve as an alias
        dup_in_source = set(self.ref[self.ref.duplicated(["source", "symbol"], keep=False)].symbol)
        sub = self.ref[self.ref.canonical_symbol.isin(self.universe) & ~self.ref.symbol.isin(dup_in_source)]
        alias = sub.groupby("symbol").canonical_symbol.agg(lambda x: set(x))
        self.alias = {k: next(iter(v)) for k, v in alias.items() if len(v) == 1 and k not in self.universe}
        self.ambiguous = {k for k, v in alias.items() if len(v) > 1 and k not in self.universe}

    def map_symbols(self, symbols: pd.Index) -> tuple[np.ndarray, np.ndarray]:
        """Return (canonical symbol or None per input, method code: exact/alias/unmapped/ambiguous)."""
        out, how = [], []
        for s in symbols:
            if s in self.universe:
                out.append(s); how.append("exact")
            elif s in self.alias:
                out.append(self.alias[s]); how.append("alias")
            elif s in self.ambiguous:
                out.append(None); how.append("ambiguous")
            else:
                out.append(None); how.append("unmapped")
        return np.array(out, dtype=object), np.array(how)

    def map_ensembl(self, ens: pd.Index) -> tuple[np.ndarray, np.ndarray]:
        out, how = [], []
        for e in ens:
            e0 = strip_version(e)
            if e0 in self.ens2sym:
                out.append(self.ens2sym[e0]); how.append("ensembl")
            else:
                out.append(None); how.append("unmapped")
        return np.array(out, dtype=object), np.array(how)

    def harmonise(self, adata: ad.AnnData, study: str, layer_keys=("counts", "tpm", "log2_tpm10")) -> tuple[ad.AnnData, dict]:
        """Map var to the universe; collapse duplicates by summing; return (new AnnData, loss record)."""
        has_ens = "ensembl_id" in adata.var and adata.var["ensembl_id"].notna().mean() > 0.9
        if has_ens:
            canon, how = self.map_ensembl(adata.var["ensembl_id"].astype(str))
        else:
            canon, how = self.map_symbols(adata.var_names)
        keep = np.array([c is not None for c in canon])
        rec = {
            "study": study, "n_input_genes": adata.n_vars, "mapped_ensembl": int((how == "ensembl").sum()),
            "mapped_exact_symbol": int((how == "exact").sum()), "mapped_alias": int((how == "alias").sum()),
            "ambiguous_dropped": int((how == "ambiguous").sum()), "unmapped_dropped": int((how == "unmapped").sum()),
            "unmapped_examples": ";".join(list(adata.var_names[how == "unmapped"][:8])),
        }
        canon_k = canon[keep].astype(str)
        uniq, inv = np.unique(canon_k, return_inverse=True)
        n_dup = int(len(canon_k) - len(uniq))
        # collapse matrix M (cells x kept genes) -> cells x uniq via sparse indicator
        ind = sp.csr_matrix((np.ones(len(canon_k), dtype=np.float32), (np.arange(len(canon_k)), inv)), shape=(len(canon_k), len(uniq)))
        order = np.argsort([self.uni_order[s] for s in uniq])
        uniq = uniq[order]; ind = ind[:, order]

        def collapse(M):
            M = M[:, keep]
            M = M if sp.issparse(M) else sp.csr_matrix(M)
            out = (M @ ind).tocsr()
            out.data = out.data.astype(M.dtype if M.dtype != np.float64 else np.float32)
            out.sort_indices()
            return out

        layers = {k: collapse(adata.layers[k]) for k in adata.layers.keys() if k in layer_keys}
        X = collapse(adata.X) if adata.X is not None else None
        var = pd.DataFrame(index=pd.Index(uniq, name="symbol"))
        var["ensembl_id"] = [self.sym2ens.get(s) for s in uniq]
        cnt = pd.Series(canon_k).value_counts()
        var["n_source_features"] = cnt.reindex(uniq).values.astype(int)
        new = ad.AnnData(X=X, obs=adata.obs.copy(), var=var, layers=layers, uns=adata.uns.copy())
        rec.update({"duplicate_symbols_collapsed": n_dup, "n_output_genes": new.n_vars})
        return new, rec
