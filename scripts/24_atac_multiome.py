#!/usr/bin/env python
"""Phase 7 / step 24 — Jessa 2022 multiome ATAC: chromatin accessibility per atlas label and per malignant metaprogram family.

Inputs : data/raw/jessa_multiome/GSM*/ (cellranger-arc MTX triplets; features hold 'Gene Expression' rows with TSS coordinates and
         per-sample 'Peaks' rows), data/processed/atlas_core.h5ad (obs, obsm['X_scvi_malignant'] / ['X_scvi']),
         data/processed/metaprogram_scores.parquet + results/tables/metaprogram_families.csv + metaprogram_genes.csv (19 / 16).
Outputs (results/tables/, never overwritten if present, unless --force):
  atac_summary.csv                one row per malignant metaprogram family: cells, samples, DA peaks up/down, ATAC kNN label agreement
                                  vs permutation null (and the same metric in RNA space), promoter concordance with the MP genes
  atac_label_correlation.csv      label x label Pearson r of pseudobulk accessibility (log CPM, union peaks; mean over samples)
  atac_label_correlation_n.csv    label x label number of samples contributing to each mean r (audit P7-03; pairs with n < 3 are
                                  masked in the figure and not interpreted); atac_label_correlation_long.csv (pair, mean r, n, min, max)
  atac_separability_sample_z.csv  per sample x space: overall kNN same-family fraction, permutation null, z (audit P7-02)
  atac_family_da_peaks.csv        per family x union peak (DA calls only): mean log2FC, samples tested / positive, sign-test p
  atac_pseudobulk_cpm.csv.gz      union peak x (sample | group) CPM; atac_union_peaks.csv.gz (chr, start, end, n_samples, gene_tss)
Figure (results/figures/): phase7_atac_separability.png.
Run: nice -n 10 python scripts/24_atac_multiome.py [--force]

Method (fixed before running)
  Cells: atlas cells of the 10 multiome samples (barcode = obs original_cell_id). Group: malignant cells selected in step 19 ->
  'mal:<family>' (incl. 'mal:unassigned'); qc_pass consensus-normal cells -> their scanvi_label; ambiguous cells dropped.
  Peaks: the authors' cellranger-arc per-sample peak-barcode matrices deposited in GEO (GSE210566; 35k-175k peaks per sample) —
  fragments are NOT re-processed and no peaks are called here (audit P7-10). Union peaks: per-sample peak intervals merged where
  they overlap (bedtools-merge logic); counts binarised per cell.
  Promoter peak: union peak overlapping a gene's TSS interval (features columns 4-6) +- 1 kb.
  DA (EXPLORATORY, audit P7-11: binomial sign test over 5-9 samples, no depth normalisation, no multiple-testing control over
  279k peaks; the table is not a DA resource) (family vs all other malignant cells, within sample, samples with >= 50 cells in both): log2((p_f + 0.005) / (p_r + 0.005)) of
  the fraction of cells with the peak accessible; peaks accessible in >= 2 % of either group. Across samples: DA-up if mean log2FC >= 1,
  >= 80 % of tested samples positive and one-sided sign-test p < 0.05 (needs >= 5 samples); DA-down symmetric.
  Separability: per sample, malignant cells with an assigned family; LSI = TF-IDF on peaks accessible in >= 2 % of those cells,
  TruncatedSVD 31 components, drop component 1 (depth), L2-normalise; fraction of the 15 nearest neighbours sharing the cell's family,
  vs 200 label permutations (z); same metric on the atlas scVI malignant latent for the same cells (RNA reference).
  Promoter concordance: per family, promoter log2FC (mean over samples) of genes in the family's MP gene sets vs all other promoter
  genes; Mann-Whitney U (one-sided, greater).
"""
import re, sys, gzip, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py, scipy.sparse as sp, scipy.io  # noqa: E402
import scipy.stats as ss  # noqa: E402
from sklearn.decomposition import TruncatedSVD  # noqa: E402
from sklearn.neighbors import NearestNeighbors  # noqa: E402
from sklearn.preprocessing import normalize  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import RAW, PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
CORE = PROCESSED / "atlas_core.h5ad"; TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
log = U.StepLogger(LOGS / "fix_phase7.log", context="24_atac_multiome")   # audit wave 3 re-run; original runs in logs/phase7_multimodal.log
MIN_CELLS, MIN_FRAC, EPS, LFC, POS_FRAC, MIN_SAMPLES, K, N_PERM, N_SVD, TSS_PAD, SEED = 50, 0.02, 0.005, 1.0, 0.8, 5, 15, 200, 31, 1000, 0
MIN_PAIR_SAMPLES = 3   # audit P7-03: label-correlation pairs with fewer contributing samples are masked / not interpreted
OUTS = [TABLES / f for f in ("atac_summary.csv", "atac_label_correlation.csv", "atac_family_da_peaks.csv")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)
rng = np.random.default_rng(SEED)

COLS = ["sample_id", "original_cell_id", "scanvi_label", "consensus_call", "qc_pass"]
with log.step("read atlas obs + labels"):
    f = h5py.File(CORE, "r")
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in COLS}, index=read_elem(f["obs"]["_index"]))
    for c in ("sample_id", "original_cell_id", "scanvi_label", "consensus_call"):
        obs[c] = obs[c].astype(str)
    m = obs.sample_id.str.startswith("jessa_multiome").values
    obs = obs[m].copy()
    zkey = "X_scvi_malignant" if "X_scvi_malignant" in f["obsm"] else "X_scvi"
    Zrna = pd.DataFrame(f["obsm"][zkey][np.flatnonzero(m)], index=obs.index)
    mps = pd.read_parquet(PROCESSED / "metaprogram_scores.parquet", columns=["metaprogram"])
    obs["group"] = np.where(obs.index.isin(mps.index), "mal:" + mps.metaprogram.reindex(obs.index).fillna("").astype(str),
                            np.where((obs.consensus_call == "normal") & obs.qc_pass.astype(bool) & (obs.scanvi_label != "malignant"), obs.scanvi_label, ""))
    obs = obs[obs.group != ""]
    log.write(f"{len(obs)} labelled multiome cells; groups {obs.group.value_counts().to_dict()}; RNA latent {zkey}")


def read_sample(d: Path):
    ft = pd.read_csv(next(d.glob("*features.tsv.gz")), sep="\t", header=None, dtype=str)
    bc = pd.read_csv(next(d.glob("*barcodes.tsv.gz")), sep="\t", header=None, dtype=str)[0].values
    with gzip.open(next(d.glob("*matrix.mtx.gz")), "rb") as fh:
        M = scipy.io.mmread(fh).tocsr()                     # features x cells
    pk = np.flatnonzero(ft[2].values == "Peaks"); gn = np.flatnonzero(ft[2].values == "Gene Expression")
    peaks = ft.iloc[pk, [3, 4, 5]].set_axis(["chr", "start", "end"], axis=1).astype({"start": int, "end": int}).reset_index(drop=True)
    genes = ft.iloc[gn, [1, 3, 4, 5]].set_axis(["symbol", "chr", "start", "end"], axis=1)
    return M[pk].T.tocsr(), peaks, bc, genes


with log.step("read peak matrices"):
    per, all_peaks, genes = {}, [], None
    for d in sorted((RAW / "jessa_multiome").glob("GSM*")):
        tag = re.search(r"(P-\d+_S-\d+)", next(d.glob("*matrix.mtx.gz")).name).group(1)
        sid = f"jessa_multiome__{tag}"
        X, peaks, bc, g = read_sample(d)
        genes = g if genes is None else genes
        o = obs[obs.sample_id == sid]
        pos = pd.Index(bc).get_indexer(o.original_cell_id)
        keep = pos >= 0
        per[sid] = dict(X=X[pos[keep]], cells=o.index[keep], peaks=peaks)
        peaks = peaks.assign(sample_id=sid, local=np.arange(len(peaks))); all_peaks.append(peaks)
        log.write(f"  {sid}: {X.shape[1]} peaks, {X.shape[0]} barcodes, {int(keep.sum())} labelled atlas cells (missing {int((~keep).sum())})")

with log.step("union peaks"):
    P = pd.concat(all_peaks, ignore_index=True).sort_values(["chr", "start", "end"]).reset_index(drop=True)
    new = np.ones(len(P), bool)
    cur_end = -1; cur_chr = None; ends = P.end.values; starts = P.start.values; chrs = P.chr.values
    for i in range(len(P)):
        if chrs[i] == cur_chr and starts[i] <= cur_end:
            new[i] = False; cur_end = max(cur_end, ends[i])
        else:
            cur_chr, cur_end = chrs[i], ends[i]
    P["union"] = np.cumsum(new) - 1
    U_ = P.groupby("union").agg(chr=("chr", "first"), start=("start", "min"), end=("end", "max"), n_samples=("sample_id", "nunique"))
    tss = genes.astype({"start": int, "end": int})
    U_["gene_tss"] = ""
    for c, g in tss.groupby("chr"):
        u = U_[U_.chr == c]
        if not len(u):
            continue
        us, ue = u.start.values, u.end.values; order = np.argsort(us)
        for sym, s0, e0 in zip(g.symbol.values, g.start.values - TSS_PAD, g.end.values + TSS_PAD):
            j = np.searchsorted(us[order], e0, side="right")
            for k in order[max(0, j - 5):j]:
                if ue[k] >= s0:
                    U_.loc[u.index[k], "gene_tss"] += (";" if U_.loc[u.index[k], "gene_tss"] else "") + sym
    n_union = len(U_)
    for sid, d in per.items():
        mp_ = P[P.sample_id == sid].set_index("local")["union"].reindex(np.arange(d["X"].shape[1])).values
        A = sp.csr_matrix((np.ones(len(mp_)), (np.arange(len(mp_)), mp_)), shape=(len(mp_), n_union))
        d["B"] = ((d["X"] > 0).astype(np.float32) @ A).tocsr(); d["B"].data[:] = 1; d.pop("X")
    log.write(f"{len(P)} sample peaks -> {n_union} union peaks; {int((U_.gene_tss != '').sum())} promoter peaks; "
              f"present in all 10 samples: {int((U_.n_samples == 10).sum())}")

with log.step("pseudobulk + label correlation"):
    pb, frac = {}, {}
    for sid, d in per.items():
        g = obs.group.reindex(d["cells"]).values
        for grp in np.unique(g):
            idx = np.flatnonzero(g == grp)
            if len(idx) >= MIN_CELLS:
                s = np.asarray(d["B"][idx].sum(0)).ravel()
                pb[f"{sid}|{grp}"] = s / s.sum() * 1e6; frac[(sid, grp)] = (s / len(idx), len(idx))
    PB = pd.DataFrame(pb); L = np.log2(PB + 1)
    groups = sorted({k.split("|")[1] for k in PB.columns})
    C = pd.DataFrame(np.nan, index=groups, columns=groups); var = L.var(1); top = var.sort_values(ascending=False).index[:20000]
    Cn = pd.DataFrame(0, index=groups, columns=groups, dtype=int); long = []
    for a in groups:
        for b in groups:
            r = [np.corrcoef(L.loc[top, f"{s}|{a}"], L.loc[top, f"{s}|{b}"])[0, 1] for s in per if f"{s}|{a}" in L and f"{s}|{b}" in L]
            C.loc[a, b] = np.mean(r) if r else np.nan; Cn.loc[a, b] = len(r)
            long.append(dict(label_a=a, label_b=b, mean_r=np.mean(r) if r else np.nan, n_samples=len(r), min_r=np.min(r) if r else np.nan,
                             max_r=np.max(r) if r else np.nan, interpretable=len(r) >= MIN_PAIR_SAMPLES))
    C.to_csv(OUTS[1]); Cn.to_csv(TABLES / "atac_label_correlation_n.csv"); pd.DataFrame(long).to_csv(TABLES / "atac_label_correlation_long.csv", index=False)
    log.write("samples per label pair (audit P7-03):\n" + Cn.to_string())
    PB.to_csv(TABLES / "atac_pseudobulk_cpm.csv.gz"); U_.to_csv(TABLES / "atac_union_peaks.csv.gz")
    log.write("within-sample accessibility correlation (top 20k variable union peaks):\n" + C.round(3).to_string())

fams = sorted({g for g in obs.group.unique() if g.startswith("mal:") and g != "mal:unassigned"})
fam_genes = {}
fm = pd.read_csv(TABLES / "metaprogram_families.csv").dropna(subset=["family"]); mg = pd.read_csv(TABLES / "metaprogram_genes.csv")
for fam, g in fm.groupby("family"):
    fam_genes["mal:" + fam] = set(mg.gene[mg.metaprogram.isin(g.metaprogram)])
prom = U_.gene_tss.str.split(";").explode(); prom = prom[prom != ""]

with log.step("differential accessibility + promoter concordance"):
    rows_da, summ = [], []
    for fam in fams:
        lfcs = {}
        for sid, d in per.items():
            g = obs.group.reindex(d["cells"]).values
            fi = np.flatnonzero(g == fam); ri = np.flatnonzero((np.char.startswith(g.astype(str), "mal:")) & (g != fam))
            if len(fi) < MIN_CELLS or len(ri) < MIN_CELLS:
                continue
            pf = np.asarray(d["B"][fi].mean(0)).ravel(); pr = np.asarray(d["B"][ri].mean(0)).ravel()
            l = np.log2((pf + EPS) / (pr + EPS)); l[(pf < MIN_FRAC) & (pr < MIN_FRAC)] = np.nan
            lfcs[sid] = l
        row = dict(family=fam.split(":", 1)[1], n_cells=int((obs.group == fam).sum()), n_samples_tested=len(lfcs))
        if len(lfcs) >= MIN_SAMPLES:
            Lm = pd.DataFrame(lfcs); n = Lm.notna().sum(1); mean = Lm.mean(1); npos = (Lm > 0).sum(1); nneg = (Lm < 0).sum(1)
            p_up = ss.binom.sf(npos - 1, n, 0.5); p_dn = ss.binom.sf(nneg - 1, n, 0.5)
            up = (mean >= LFC) & (npos >= POS_FRAC * n) & (p_up < 0.05) & (n >= MIN_SAMPLES)
            dn = (mean <= -LFC) & (nneg >= POS_FRAC * n) & (p_dn < 0.05) & (n >= MIN_SAMPLES)
            for k in np.flatnonzero(up | dn):
                rows_da.append(dict(family=row["family"], union_peak=int(k), direction="up" if up.iloc[k] else "down", mean_log2fc=mean.iloc[k],
                                    n_samples=int(n.iloc[k]), n_positive=int(npos.iloc[k]), p_sign=float(p_up[k] if up.iloc[k] else p_dn[k]),
                                    chr=U_.chr.iloc[k], start=int(U_.start.iloc[k]), end=int(U_.end.iloc[k]), gene_tss=U_.gene_tss.iloc[k]))
            pl = mean.reindex(prom.index).groupby(prom.values).mean().dropna()
            inset = pl.index.isin(fam_genes.get(fam, set()))
            u = ss.mannwhitneyu(pl[inset], pl[~inset], alternative="greater") if inset.sum() >= 5 else None
            row.update(n_da_up=int(up.sum()), n_da_down=int(dn.sum()), n_promoter_genes_in_set=int(inset.sum()),
                       promoter_lfc_set_median=float(pl[inset].median()) if inset.any() else np.nan,
                       promoter_lfc_other_median=float(pl[~inset].median()), promoter_p=float(u.pvalue) if u else np.nan)
        summ.append(row)
    summ = pd.DataFrame(summ)
    log.write("\n" + summ.round(4).to_string(index=False))

with log.step("separability (ATAC LSI vs RNA scVI)"):
    sep = []
    for sid, d in per.items():
        g = obs.group.reindex(d["cells"]).values
        idx = np.flatnonzero(np.isin(g, fams) & np.isfinite(Zrna.loc[d["cells"]].values).all(1))   # cells with both ATAC and RNA latent
        if len(idx) < 200 or len(set(g[idx])) < 2:
            continue
        B = d["B"][idx]; keep = np.asarray(B.mean(0)).ravel() >= MIN_FRAC; B = B[:, keep]
        tf = normalize(B, norm="l1", axis=1); idf = np.log1p(B.shape[0] / (1 + np.asarray(B.sum(0)).ravel()))
        Xl = TruncatedSVD(N_SVD, random_state=SEED).fit_transform(tf.multiply(idf).tocsr() * 1e4)
        spaces = {"atac": normalize(Xl[:, 1:]), "rna": Zrna.loc[d["cells"][idx]].values}
        lab = g[idx]
        for sp_name, Xs in spaces.items():
            nn = NearestNeighbors(n_neighbors=K + 1).fit(Xs).kneighbors(Xs, return_distance=False)[:, 1:]
            same = lambda lb: (lb[nn] == lb[:, None]).mean(1)   # noqa: E731
            obs_s = same(lab); null = np.array([same(rng.permutation(lab)).mean() for _ in range(N_PERM)])
            for fam in sorted(set(lab)):
                mm = lab == fam
                sep.append(dict(sample_id=sid, space=sp_name, family=fam.split(":", 1)[1], n_cells=int(mm.sum()), knn_same=float(obs_s[mm].mean()),
                                expected=float(mm.mean()), overall_knn_same=float(obs_s.mean()), null_mean=float(null.mean()),
                                z_overall=float((obs_s.mean() - null.mean()) / null.std())))
    sep = pd.DataFrame(sep)
    agg = sep.groupby(["family", "space"]).apply(lambda x: pd.Series({"knn_same": np.average(x.knn_same, weights=x.n_cells),
                                                                      "expected": np.average(x.expected, weights=x.n_cells)})).unstack("space")
    agg.columns = [f"{a}_{b}" for a, b in agg.columns]
    summ = summ.merge(agg.reset_index(), on="family", how="left")
    summ["atac_enrichment"] = summ.knn_same_atac / summ.expected_atac; summ["rna_enrichment"] = summ.knn_same_rna / summ.expected_rna
    zs = sep.drop_duplicates(["sample_id", "space"]).groupby("space").z_overall.median()
    sz = (sep.groupby(["sample_id", "space"]).agg(n_cells=("n_cells", "sum"), n_families=("family", "nunique"), overall_knn_same=("overall_knn_same", "first"),
                                                null_mean=("null_mean", "first"), z_overall=("z_overall", "first")).reset_index())
    sz.to_csv(TABLES / "atac_separability_sample_z.csv", index=False)
    za = sz[sz.space == "atac"]
    log.write(f"audit P7-02 per-sample ATAC z: {za.set_index('sample_id').z_overall.round(2).to_dict()}; samples with z > 2: {int((za.z_overall > 2).sum())} of {len(za)}; "
              f"z <= 2 (at chance): {za.sample_id[za.z_overall <= 2].tolist()}; min {za.z_overall.min():.2f}, median {za.z_overall.median():.2f}, max {za.z_overall.max():.2f}")
    log.write(f"kNN family agreement, median per-sample z vs permutation: {zs.round(1).to_dict()}\n" +
              summ[["family", "n_cells", "knn_same_atac", "knn_same_rna", "expected_atac", "atac_enrichment", "rna_enrichment"]].round(3).to_string(index=False))

with log.step("write + figure"):
    summ.to_csv(OUTS[0], index=False); pd.DataFrame(rows_da).to_csv(OUTS[2], index=False)
    sep.to_csv(TABLES / "atac_separability_per_sample.csv", index=False)
    apply_figure_style()   # house style: no titles / text inside panels (captions in docs/phase7_multimodal_notes.md), legends outside
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={"width_ratios": [1, 0.8, 1.3]}, layout="constrained")
    x = np.arange(len(summ)); w = 0.38
    axes[0].bar(x - w / 2, summ.atac_enrichment, w, color="#2166ac", label="ATAC (LSI)")
    axes[0].bar(x + w / 2, summ.rna_enrichment, w, color="#b2182b", label="RNA (scVI)")
    axes[0].axhline(1, color="black", lw=0.5, ls=":"); axes[0].set_xticks(x); axes[0].set_xticklabels(summ.family, rotation=45, ha="right")
    axes[0].set_ylabel("kNN same-family / expected"); axes[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    # per-sample z of the overall kNN family agreement vs the permutation null (audit P7-02)
    for i, (spn, col) in enumerate((("atac", "#2166ac"), ("rna", "#b2182b"))):
        v = sz[sz.space == spn].z_overall.values; jit = rng.uniform(-0.12, 0.12, len(v))
        axes[1].scatter(np.full(len(v), i) + jit, v, s=9, color=col, edgecolor="none")
    axes[1].axhline(2, color="black", lw=0.5, ls=":"); axes[1].axhline(0, color="grey", lw=0.5)
    axes[1].set_yscale("symlog", linthresh=2); axes[1].set_xticks([0, 1]); axes[1].set_xticklabels(["ATAC (LSI)", "RNA (scVI)"], rotation=45, ha="right")
    axes[1].set_ylabel("per-sample z vs permutation"); axes[1].set_xlim(-0.5, 1.5)
    # label correlation: pairs with < MIN_PAIR_SAMPLES contributing samples masked (grey), n printed in every cell (audit P7-03)
    Cm = np.ma.masked_where(Cn.values < MIN_PAIR_SAMPLES, C.values.astype(float))
    cmap = plt.get_cmap("viridis").copy(); cmap.set_bad("#d9d9d9")
    im = axes[2].imshow(Cm, cmap=cmap, aspect="auto", vmin=np.nanmin(C.values), vmax=1)
    for i in range(len(C)):
        for j in range(len(C)):
            axes[2].text(j, i, str(int(Cn.values[i, j])), ha="center", va="center", fontsize=4, color="white" if Cn.values[i, j] >= MIN_PAIR_SAMPLES and Cm[i, j] < 0.85 else "black")
    axes[2].set_xticks(range(len(C))); axes[2].set_xticklabels([c.replace("mal:", "") for c in C.columns], rotation=60, ha="right", fontsize=5)
    axes[2].set_yticks(range(len(C))); axes[2].set_yticklabels([c.replace("mal:", "") for c in C.index], fontsize=5)
    fig.colorbar(im, ax=axes[2], fraction=0.04, pad=0.02).set_label("accessibility r (mean over samples)")
    fig.savefig(FIGS / "phase7_atac_separability.png"); plt.close(fig)
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
