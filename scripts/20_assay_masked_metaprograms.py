#!/usr/bin/env python
"""Phase 6 / step 20 — cell vs nucleus metaprograms: assay-biased genes, assay-masked re-clustering and replication re-test.

Inputs : data/processed/atlas_core.h5ad (malignant cells of Jessa 2022, layers['counts'] via h5py),
         results/tables/cnmf/<sample>_programs.csv + samples_analysed.csv, results/tables/metaprograms.csv + metaprogram_genes.csv.
Outputs (results/tables/, never overwritten if present, unless --force):
  assay_biased_genes.csv               per gene: paired log2FC (nuclei - cells) for P-6337 and P-6640, Jessa-wide Welch test, call
  metaprograms_assay_masked.csv        metaprograms after masking (same columns as metaprograms.csv)
  metaprogram_members_assay_masked.csv programme -> masked metaprogram, with the original (unmasked) metaprogram
  metaprogram_masking_summary.csv      n MPs / replicated MPs at each tree cut, unmasked vs masked
  metaprogram_cell_vs_nucleus.csv      Jaccard between every cell-only and nucleus-only MP gene set, unmasked vs masked
Figure (results/figures/): phase6_metaprogram_masked_replication.png.
Run: nice -n 10 python scripts/20_assay_masked_metaprograms.py [--force]

Audit corrections (2026-10-19; MC-05, MC-15, MC-04). EXPLORATORY analysis: the mask is by construction the set of genes that
differ between the assays in Jessa 2022 (Welch 4 cell vs 23 nuclei samples, two of the four cell samples being the paired patients),
so a cross-assay re-test after removing them is partly tautological, and nuclei programmes keep a median of ~20 of 50 genes.
  * malignant cells = primary selection of 13a_malignant_selection.py (no scANVI rescue);
  * random-mask control: N_RANDOM seeds, each removing from every programme the SAME NUMBER of genes as the real mask but chosen at
    random among that programme's genes; MPs / replicated MPs reported as median [min-max] over seeds in metaprogram_masking_summary.csv
    (genes = 'random_mask');
  * replicated_strict (>= 2 samples per assay class in addition to the pre-registered rule) reported for every cut;
  * cut sensitivity in 0.01 steps 0.75-0.92; ribosomal-protein genes by explicit regex (^RP[SL]\d+[A-Z]?$, RPLP0-2, RPSA) instead of the
    RPS\d/RPL\d prefix (which caught 17 kinases / paralogues).

Method (fixed before running)
  Pseudobulk: per Jessa 2022 sample, sum of raw counts over its selected malignant cells (primary selection of 13a, qc_pass), log2 CPM
  (+1). Genes with CPM >= 10 in at least 2 samples are tested.
  Assay-biased gene: (a) the two patients profiled both ways (P-6337, P-6640: scRNA cells vs multiome nuclei) both show
  |log2FC| >= 1 in the same direction, AND (b) Welch t-test across all Jessa samples (4 cell vs 23 nuclei samples) BH q < 0.05 in
  that direction; plus every mitochondrial (MT-) and ribosomal-protein (RPS*/RPL*) gene regardless of the test.
  Masking: assay-biased genes are removed from each programme's top-50 list, then 16_metaprograms.py's clustering is repeated
  unchanged (Jaccard, average linkage, cut 0.85, >= 5 programmes from >= 3 samples, publications, replication rule).
"""
import re, sys, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py, scipy.sparse as sp  # noqa: E402
import scipy.stats as ss  # noqa: E402
from scipy.cluster.hierarchy import linkage, fcluster  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
CORE = PROCESSED / "atlas_core.h5ad"; CNMF = RESULTS / "tables" / "cnmf"; TABLES = RESULTS / "tables"
FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
log = U.StepLogger(LOGS / "fix_metaprograms.log", context="20_assay_masked")
PAIRS = {"P-6337": ("jessa_scrna__P-6337_S-8821", "jessa_multiome__P-6337_S-8821"),
         "P-6640": ("jessa_scrna__P-6640_S-9581", "jessa_multiome__P-6640_S-9581")}
LFC, Q, MIN_CPM = 1.0, 0.05, 10
CUT, MIN_PROGRAMS, MIN_SAMPLES, MIN_GENE_FRAC, TOP_N = 0.85, 5, 3, 0.25, 50   # = 16_metaprograms.py
SENS_CUTS = [round(c, 2) for c in np.arange(0.75, 0.9201, 0.01)]
N_RANDOM, SEED = 10, 0
RP_RE = r"^RP[SL]\d+[A-Z]?$|^RPLP[0-2]$|^RPSA$"
MIN_SAMPLES_PER_ASSAY = 2
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}
OUTS = [TABLES / f for f in ("assay_biased_genes.csv", "metaprograms_assay_masked.csv", "metaprogram_members_assay_masked.csv",
                             "metaprogram_masking_summary.csv", "metaprogram_cell_vs_nucleus.csv")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)


def read_rows(g, rows: np.ndarray, n_genes: int, block: int = 20000) -> sp.csr_matrix:
    """Sorted row subset of an on-disk CSR group (layers/counts), read in contiguous row blocks via h5py."""
    indptr = g["indptr"][:]; parts = []
    for b0 in range(0, len(indptr) - 1, block):
        r = rows[(rows >= b0) & (rows < b0 + block)]
        if not len(r):
            continue
        b1 = min(b0 + block, len(indptr) - 1); s, e = int(indptr[b0]), int(indptr[b1])
        m = sp.csr_matrix((g["data"][s:e], g["indices"][s:e], indptr[b0:b1 + 1] - s), shape=(b1 - b0, n_genes))
        parts.append(m[r - b0])
    return sp.vstack(parts).tocsr()


COLS = ["sample_id", "study", "assay", "qc_pass"]
with log.step("Jessa malignant pseudobulk"):
    f = h5py.File(CORE, "r")
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in COLS}, index=read_elem(f["obs"]["_index"]))
    for c in ("sample_id", "study", "assay"):
        obs[c] = obs[c].astype(str)
    var = read_elem(f["var"]); genes = var.index.to_numpy()
    selp = pd.read_parquet(PROCESSED / "malignant_selection_audit.parquet", columns=["primary"]).reindex(obs.index)
    sel = (selp.primary.fillna(False).astype(bool) & obs.study.str.startswith("jessa")).values
    rows = np.flatnonzero(sel)
    X = read_rows(f["layers/counts"], rows, len(genes))
    sid = obs.sample_id.values[rows]
    samples = sorted(set(sid))
    G = sp.csr_matrix((np.ones(len(rows)), (pd.Index(samples).get_indexer(sid), np.arange(len(rows)))), shape=(len(samples), len(rows)))
    PB = pd.DataFrame(np.asarray((G @ X).todense()), index=samples, columns=genes); del X
    cpm = PB.div(PB.sum(1), axis=0) * 1e6
    L = np.log2(cpm + 1)
    assay = obs.groupby("sample_id").assay.first().reindex(samples).map(ASSAY_CLASS)
    log.write(f"{len(samples)} Jessa samples ({assay.value_counts().to_dict()}), {len(rows)} malignant cells")

with log.step("assay-biased genes"):
    tested = cpm.columns[(cpm >= MIN_CPM).sum() >= 2]
    t = pd.DataFrame(index=tested)
    for pid, (c, n) in PAIRS.items():
        t[f"lfc_{pid}"] = L.loc[n, tested] - L.loc[c, tested]
    cells, nuc = assay.index[assay == "cells"], assay.index[assay == "nuclei"]
    tt = ss.ttest_ind(L.loc[nuc, tested], L.loc[cells, tested], equal_var=False)
    t["lfc_jessa_mean"] = L.loc[nuc, tested].mean() - L.loc[cells, tested].mean()
    t["t_welch"] = tt.statistic; t["p_welch"] = tt.pvalue
    t["q_welch"] = ss.false_discovery_control(np.nan_to_num(tt.pvalue, nan=1.0), method="bh")
    l1, l2 = t[f"lfc_{list(PAIRS)[0]}"], t[f"lfc_{list(PAIRS)[1]}"]
    pair_up = (l1 >= LFC) & (l2 >= LFC); pair_dn = (l1 <= -LFC) & (l2 <= -LFC)
    test_up = (t.q_welch < Q) & (t.t_welch > 0); test_dn = (t.q_welch < Q) & (t.t_welch < 0)
    t["mito_ribo"] = t.index.str.match(r"^MT-") | t.index.str.match(RP_RE)
    t["call"] = np.select([t.mito_ribo, pair_up & test_up, pair_dn & test_dn], ["mito_ribo", "nuclei_biased", "cells_biased"], "none")
    biased = set(t.index[t.call != "none"])
    allmr = {g for g in genes if re.match(r"^MT-", g) or re.match(RP_RE, g)}
    biased |= allmr
    t.sort_values("t_welch").to_csv(OUTS[0])
    log.write(f"{len(tested)} genes tested; calls {t.call.value_counts().to_dict()}; masked set {len(biased)} genes. "
              f"Top nuclei-biased: {', '.join(t[t.call == 'nuclei_biased'].sort_values('t_welch', ascending=False).index[:12])}; "
              f"top cells-biased: {', '.join(t[t.call == 'cells_biased'].sort_values('t_welch').index[:12])}")

with log.step("masked re-clustering"):
    samp = pd.read_csv(CNMF / "samples_analysed.csv", index_col=0); samp = samp[samp.analysed]
    samp["publication"] = samp.study.map(lambda s: PUBLICATION.get(s, s)); samp["assay_class"] = samp.assay.map(ASSAY_CLASS)
    prog = pd.concat([pd.read_csv(CNMF / f"{s}_programs.csv") for s in samp.index], ignore_index=True)
    names = prog.program.drop_duplicates().to_list()
    pmeta = prog.drop_duplicates("program").set_index("program").loc[names, ["sample_id"]].join(samp[["study", "publication", "assay_class"]], on="sample_id")
    lost = prog.gene.isin(biased).groupby(prog.program).mean().reindex(names)
    log.write(f"fraction of programme top-50 genes masked: median {lost.median():.2f}, by assay class "
              + str(lost.groupby(pmeta.assay_class).median().round(3).to_dict()))

    def cluster(pr):
        sets = [set(pr.gene[pr.program == p]) for p in names]
        J = MP.jaccard_matrix(sets); D = 1 - J; np.fill_diagonal(D, 0)
        return J, linkage(squareform(D, checks=False), "average")

    def call(Z, cut):
        lab = fcluster(Z, cut, "distance"); mps = []
        for c in np.unique(lab):
            m = np.flatnonzero(lab == c)
            if len(m) >= MIN_PROGRAMS and pmeta.sample_id.iloc[m].nunique() >= MIN_SAMPLES:
                mps.append(m)
        return sorted(mps, key=lambda m: -len(m))

    def replicated(m, strict=False):
        sub = pmeta.iloc[m]
        ok = sub.publication.nunique() >= 3 and set(sub.assay_class) >= {"cells", "nuclei"}
        if strict:
            ok = ok and all(sub.loc[sub.assay_class == a, "sample_id"].nunique() >= MIN_SAMPLES_PER_ASSAY for a in ("cells", "nuclei"))
        return bool(ok)

    def summarise(Z, cut, tag, seed=None):
        mps = call(Z, cut)
        return dict(cut=cut, genes=tag, seed=seed, n_metaprograms=len(mps), n_replicated=sum(replicated(m) for m in mps),
                    n_replicated_strict=sum(replicated(m, True) for m in mps),
                    n_mixed_assay=sum(set(pmeta.assay_class.iloc[m]) >= {"cells", "nuclei"} for m in mps),
                    n_programmes_assigned=int(sum(len(m) for m in mps)), primary=cut == CUT)

    pm = prog[~prog.gene.isin(biased)]
    J0, Z0 = cluster(prog); J1, Z1 = cluster(pm)
    rows = [summarise(Z, cut, tag) for cut in SENS_CUTS for tag, Z in (("unmasked", Z0), ("masked", Z1))]
    # random-mask control: same number of genes removed per programme, chosen at random within the programme (MC-05)
    n_masked = prog.gene.isin(biased).groupby(prog.program).sum()
    rng = np.random.default_rng(SEED); rand_rows = []
    for seed in range(N_RANDOM):
        keep_idx = []
        for p, g in prog.groupby("program", sort=False):
            k = int(n_masked[p]); idx = g.index.to_numpy()
            keep_idx.append(rng.choice(idx, size=len(idx) - k, replace=False) if k else idx)
        pr = prog.loc[np.sort(np.concatenate(keep_idx))]
        _, Zr = cluster(pr)
        rand_rows += [summarise(Zr, cut, "random_mask", seed) for cut in SENS_CUTS]
    rand = pd.DataFrame(rand_rows)
    agg = rand.groupby("cut").agg(n_metaprograms=("n_metaprograms", "median"), n_metaprograms_min=("n_metaprograms", "min"), n_metaprograms_max=("n_metaprograms", "max"),
                                 n_replicated=("n_replicated", "median"), n_replicated_min=("n_replicated", "min"), n_replicated_max=("n_replicated", "max"),
                                 n_replicated_strict=("n_replicated_strict", "median"), n_mixed_assay=("n_mixed_assay", "median"),
                                 n_programmes_assigned=("n_programmes_assigned", "median")).reset_index()
    agg["genes"] = "random_mask"; agg["seed"] = f"median of {N_RANDOM}"; agg["primary"] = agg.cut == CUT
    summ = pd.concat([pd.DataFrame(rows), agg], ignore_index=True).sort_values(["cut", "genes"])
    rand.to_csv(TABLES / "metaprogram_masking_random_control.csv", index=False)
    log.write("\n" + summ[summ.cut.isin([0.75, 0.80, 0.84, 0.85, 0.86, 0.87, 0.90])].drop(columns="seed").to_string(index=False))

with log.step("masked metaprograms"):
    orig = pd.read_csv(TABLES / "metaprogram_members.csv").set_index("program").metaprogram.reindex(names)
    mps = call(Z1, CUT); rows_mp = []; mem = pd.Series("unassigned", index=names)
    for i, m in enumerate(mps):
        mp = f"mMP{i+1}"; mm = [names[j] for j in m]; mem[mm] = mp
        sub = pm[pm.program.isin(mm)]
        g = sub.groupby("gene").agg(n_programs=("program", "nunique"), mean_spectra=("spectra_score", "mean"))
        g = g[g.n_programs >= MIN_GENE_FRAC * len(mm)].sort_values(["n_programs", "mean_spectra"], ascending=False).head(TOP_N)
        p_ = pmeta.loc[mm]
        rows_mp.append(dict(metaprogram=mp, n_programs=len(mm), n_samples=p_.sample_id.nunique(), n_publications=p_.publication.nunique(),
                            publications=";".join(sorted(p_.publication.unique())),
                            n_samples_cells=p_.loc[p_.assay_class == "cells", "sample_id"].nunique(),
                            n_samples_nuclei=p_.loc[p_.assay_class == "nuclei", "sample_id"].nunique(),
                            replicated=replicated(m), replicated_strict=replicated(m, True), original_mps=";".join(f"{k}:{v}" for k, v in orig.loc[mm].value_counts().items()),
                            top_genes=";".join(g.index[:15])))
    mpt_m = pd.DataFrame(rows_mp)
    log.write("masked metaprograms:\n" + mpt_m[["metaprogram", "n_programs", "n_publications", "n_samples_cells", "n_samples_nuclei", "replicated", "replicated_strict",
                                                "original_mps", "top_genes"]].to_string(index=False))

with log.step("cell-only vs nucleus-only MP overlap"):
    mpt = pd.read_csv(TABLES / "metaprograms.csv"); mpg = pd.read_csv(TABLES / "metaprogram_genes.csv")
    gs = {m: set(g.gene) for m, g in mpg.groupby("metaprogram")}
    cell_only = mpt.metaprogram[(mpt.n_samples_cells > 0) & (mpt.n_samples_nuclei == 0)].tolist()
    nuc_only = mpt.metaprogram[(mpt.n_samples_nuclei > 0) & (mpt.n_samples_cells == 0)].tolist()
    lab = mpt.set_index("metaprogram").label
    cvn = []
    for c in cell_only:
        for n in nuc_only:
            a, b = gs[c], gs[n]; am, bm = a - biased, b - biased
            cvn.append(dict(cell_mp=c, cell_label=lab[c], nucleus_mp=n, nucleus_label=lab[n], jaccard=len(a & b) / len(a | b),
                            jaccard_masked=len(am & bm) / len(am | bm) if am | bm else np.nan, shared_masked=";".join(sorted(am & bm))))
    cvn = pd.DataFrame(cvn).sort_values("jaccard_masked", ascending=False)
    log.write("\n" + cvn.drop(columns="shared_masked").head(12).round(3).to_string(index=False))

with log.step("write + figure"):
    mpt_m.to_csv(OUTS[1], index=False)
    pd.DataFrame({"program": names, "metaprogram_masked": mem.values, "metaprogram_unmasked": orig.values,
                  "frac_genes_masked": lost.values}).join(pmeta.reset_index(drop=True)).to_csv(OUTS[2], index=False)
    summ.to_csv(OUTS[3], index=False); cvn.to_csv(OUTS[4], index=False)
    apply_figure_style()
    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    for tag, ls in (("unmasked", "-"), ("masked", "--")):
        d = summ[summ.genes == tag]
        ax.plot(d.cut, d.n_metaprograms, ls, marker="o", ms=2.5, color="#888888", label=f"MPs ({tag})")
        ax.plot(d.cut, d.n_replicated, ls, marker="s", ms=2.5, color="#b2182b", label=f"replicated ({tag})")
        ax.plot(d.cut, d.n_replicated_strict, ls, marker="^", ms=2.5, color="#2166ac", label=f"replicated, ≥2 samples/assay ({tag})")
    d = summ[summ.genes == "random_mask"]
    ax.fill_between(d.cut, d.n_replicated_min, d.n_replicated_max, color="#b2182b", alpha=0.15, lw=0, label=f"replicated (random mask, {N_RANDOM} seeds)")
    ax.fill_between(d.cut, d.n_metaprograms_min, d.n_metaprograms_max, color="#888888", alpha=0.15, lw=0, label=f"MPs (random mask, {N_RANDOM} seeds)")
    ax.axvline(CUT, color="black", lw=0.5, ls=":"); ax.set_xlabel("tree cut (1 − mean Jaccard)"); ax.set_ylabel("count")
    ax.legend(fontsize=5, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig.savefig(FIGS / "phase6_metaprogram_masked_replication.png"); plt.close(fig)
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
