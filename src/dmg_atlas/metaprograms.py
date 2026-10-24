"""Phase 6 — per-sample consensus NMF (cNMF-style) on malignant cells, metaprogram consensus and scoring helpers.

Malignant-cell selection (recorded in obs['mal_sel_source']):
    consensus_call == 'malignant'  OR  (scanvi_label == 'malignant' AND scanvi_maxprob >= 0.9).

Per-sample cNMF (run_cnmf_sample):
    1. counts of the selected cells -> genes in >= 3 cells -> 2,000 seurat_v3 HVGs within the sample
       -> normalize_total(1e4) + log1p -> dense float32 (cells x 2,000).
    2. for k in 4..9: N_RESTARTS sklearn NMF fits (init='random', solver='cd', max_iter=500, tol=1e-4,
       seeds 0..N_RESTARTS-1). Components are L2-normalised, pooled (N_RESTARTS*k rows), clustered by
       average-linkage on cosine distance into k clusters; the consensus programme is the cluster median.
       Stability(k) = silhouette of that clustering (precomputed cosine distances); the 'full' criterion is
       that every consensus cluster contains components from >= MIN_RESTART_FRAC of the restarts.
    3. k selection rule (pre-specified here): k* = the LARGEST k in 4..9 with silhouette >= SIL_MIN and
       'full' satisfied; if no k qualifies, k* = argmax silhouette (ties -> smaller k). Reconstruction error
       is recorded for every k but is not used in the rule (it decreases monotonically with k).
    4. usages W: non-negative least squares of each cell on the fixed consensus H; gene spectra scores:
       ordinary least squares of per-gene z-scored expression on the row-normalised usages (cNMF's
       'gene_spectra_score'), top-50 genes per programme by that score.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import nnls
from scipy.cluster.hierarchy import linkage, fcluster, cophenet
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import normalize

KS = tuple(range(4, 10))
N_RESTARTS = 10
SIL_MIN = 0.90
MIN_RESTART_FRAC = 0.8
TOP_N = 50


def select_malignant(obs: pd.DataFrame) -> pd.Series:
    cons = obs["consensus_call"].astype(str) == "malignant"
    scv = (obs["scanvi_label"].astype(str) == "malignant") & (obs["scanvi_maxprob"] >= 0.9)
    src = np.select([cons & scv, cons, scv], ["both", "consensus_only", "scanvi_only"], "not_selected")
    return pd.Series(src, index=obs.index, name="mal_sel_source")


def prep_matrix(X_counts: sp.csr_matrix, genes: np.ndarray, n_top: int = 2000):
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(X=X_counts.astype(np.float32), var=pd.DataFrame(index=pd.Index(genes)))
    sc.pp.filter_genes(a, min_cells=3)
    sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=n_top, subset=False)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a = a[:, a.var.highly_variable].copy()
    return a.X.toarray().astype(np.float32), a.var_names.to_numpy()


def consensus_nmf(M: np.ndarray, ks=KS, n_restarts: int = N_RESTARTS):
    out = {}
    for k in ks:
        Hs, errs = [], []
        for s in range(n_restarts):
            m = NMF(n_components=k, init="random", solver="cd", max_iter=500, random_state=s, tol=1e-4)
            m.fit(M)
            Hs.append(normalize(m.components_, norm="l2"))
            errs.append(m.reconstruction_err_)
        H = np.vstack(Hs)
        rst = np.repeat(np.arange(n_restarts), k)
        D = pdist(H, metric="cosine")
        Z = linkage(D, "average")
        lab = fcluster(Z, k, "maxclust")
        sil = float(silhouette_score(squareform(D), lab, metric="precomputed")) if len(set(lab)) > 1 else np.nan
        frac = [len(set(rst[lab == c])) / n_restarts for c in np.unique(lab)]
        Hc = normalize(np.vstack([np.median(H[lab == c], axis=0) for c in np.unique(lab)]), "l2")
        out[k] = dict(H=Hc, sil=sil, coph=float(cophenet(Z, D)[0]), full=bool(min(frac) >= MIN_RESTART_FRAC),
                      min_restart_frac=float(min(frac)), n_clusters=int(len(np.unique(lab))),
                      err_mean=float(np.mean(errs)), err_sd=float(np.std(errs)))
    return out


def choose_k(res: dict) -> tuple[int, str]:
    ok = [k for k, v in res.items() if v["full"] and v["n_clusters"] == k and v["sil"] >= SIL_MIN]
    if ok:
        return max(ok), "largest_k_with_sil>=%.2f_and_full" % SIL_MIN
    sils = {k: v["sil"] for k, v in res.items()}
    return min(sils, key=lambda k: (-sils[k], k)), "argmax_silhouette_fallback"


def usages_nnls(M: np.ndarray, H: np.ndarray) -> np.ndarray:
    Ht = H.T.astype(np.float64)
    W = np.empty((M.shape[0], H.shape[0]))
    for i in range(M.shape[0]):
        W[i] = nnls(Ht, M[i].astype(np.float64))[0]
    return W


def gene_spectra_scores(M: np.ndarray, W: np.ndarray) -> np.ndarray:
    """OLS of z-scored genes on row-normalised usages -> (k x genes) coefficient matrix."""
    Wn = W / np.clip(W.sum(1, keepdims=True), 1e-12, None)
    Z = (M - M.mean(0)) / np.clip(M.std(0), 1e-8, None)
    beta, *_ = np.linalg.lstsq(Wn, Z, rcond=None)
    return beta


def ucell_scores(X: sp.csr_matrix, genes: np.ndarray, gene_sets: dict, max_rank: int = 1500, chunk: int = 5000) -> pd.DataFrame:
    """UCell (Andreatta & Carmona 2021): per cell, genes ranked by decreasing expression (average ranks for ties;
    ranks > max_rank set to max_rank + 1); U = 1 - (sum(rank_set) - n(n+1)/2) / (n * max_rank). Rank-based, hence
    invariant to per-cell normalisation. Zero-expression genes share the bottom tie and receive max_rank + 1."""
    from scipy.stats import rankdata
    gidx = {g: i for i, g in enumerate(genes)}
    sets = {k: np.array([gidx[g] for g in v if g in gidx]) for k, v in gene_sets.items()}
    out = np.zeros((X.shape[0], len(sets)), dtype=np.float32)
    for s in range(0, X.shape[0], chunk):
        D = X[s:s + chunk].toarray().astype(np.float32)
        R = rankdata(-D, axis=1, method="average")
        R[D <= 0] = max_rank + 1           # all-zero ties -> bottom rank, as in UCell
        R = np.minimum(R, max_rank + 1)
        for j, (k, idx) in enumerate(sets.items()):
            n = len(idx)
            out[s:s + chunk, j] = 1.0 - (R[:, idx].sum(1) - n * (n + 1) / 2.0) / (n * max_rank)
    return pd.DataFrame(out, columns=list(sets.keys()))


def jaccard_matrix(sets: list[set]) -> np.ndarray:
    n = len(sets)
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            u = len(sets[i] | sets[j])
            J[i, j] = J[j, i] = len(sets[i] & sets[j]) / u if u else 0.0
    return J


def run_cnmf_sample(X_counts: sp.csr_matrix, genes: np.ndarray, cell_ids: np.ndarray, sample_id: str):
    M, hv = prep_matrix(X_counts, genes)
    res = consensus_nmf(M)
    k, rule = choose_k(res)
    H = res[k]["H"]
    W = usages_nnls(M, H)
    S = gene_spectra_scores(M, W)
    prog_names = [f"{sample_id}|k{k}|p{i+1}" for i in range(k)]
    rows = []
    for i, p in enumerate(prog_names):
        order = np.argsort(-S[i])[:TOP_N]
        for r, j in enumerate(order):
            rows.append((p, sample_id, k, i + 1, r + 1, hv[j], float(S[i, j]), float(H[i, j])))
    programs = pd.DataFrame(rows, columns=["program", "sample_id", "k", "program_index", "rank", "gene", "spectra_score", "loading"])
    usages = pd.DataFrame(W, index=pd.Index(cell_ids, name="cell_id"), columns=prog_names)
    ksel = pd.DataFrame([{**{"sample_id": sample_id, "k": kk, "chosen": kk == k, "rule": rule},
                          **{a: v[a] for a in ("sil", "coph", "full", "min_restart_frac", "n_clusters", "err_mean", "err_sd")}}
                         for kk, v in res.items()])
    return programs, usages, ksel
