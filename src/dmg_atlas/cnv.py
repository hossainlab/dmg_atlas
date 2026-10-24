"""Phase 3 — malignant vs normal calling: marker scoring, inferCNV (infercnvpy), CNV clustering, consensus.

Pipeline per sample (see scripts/05_cnv_calls.py):
  1. working cell set           -> working_mask()
  2. log-normalise              -> normalise()
  3. marker scores + calls      -> score_markers(), marker_calls()
  4. reference cells            -> marker_calls()['is_reference'] (pooled study reference if < MIN_REF cells)
  5. inferCNV + arm means       -> run_infercnv()
  6. CNV clustering + calls     -> cnv_calls()
  7. consensus tiers            -> consensus()

Gene positions: GENCODE v32 (GRCh38, same annotation as the cellranger 2020-A gene universe used in Phase 2),
`data/external/gencode/gencode.v32.annotation.gtf.gz`, matched by Ensembl ID first and gene symbol second.
Chromosome arms: CENTROMERES_GRCH38 holds approximate GRCh38 centromere (acen) intervals typed from memory of the UCSC hg38
cytoBand table (the UCSC download is blocked by the network allowlist in this environment and was NOT retrieved; values are
rounded to 0.1 Mb and only used to split arms; validate_arm_boundaries() confirms each midpoint sits inside a protein-coding
gene-free gap of the GTF, so the arm assignment is insensitive to the exact value; acrocentric chromosomes are all 'q'). The 4q12 band
(PDGFRA/KIT region) is taken as chr4:52.7-59.5 Mb (PDGFRA chr4:54.23 Mb and KIT chr4:54.66 Mb from the GTF fall inside).
"""
from __future__ import annotations

import gzip
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

from dmg_atlas import ROOT, INTERIM

GENCODE_GTF = ROOT / "data" / "external" / "gencode" / "gencode.v32.annotation.gtf.gz"
GENEPOS_CSV = ROOT / "data" / "external" / "gencode" / "gene_positions_gencode_v32.csv"
ANNOTATION_SOURCE = "GENCODE v32 (GRCh38.p13, Ensembl 98) gencode.v32.annotation.gtf.gz from ftp.ebi.ac.uk; gene features; match Ensembl ID then symbol"

# ---------------------------------------------------------------- marker sets
NORMAL_MARKERS = {
    "myeloid": ["PTPRC", "CSF1R", "CX3CR1", "P2RY12", "C1QA", "C1QB", "CD74"],
    "T_NK": ["CD3D", "CD3E", "CD2", "NKG7"],
    "oligodendrocyte": ["MBP", "PLP1", "MOG", "MAG", "ST18"],
    "vascular": ["CLDN5", "FLT1", "VWF", "PDGFRB", "RGS5"],
    "neuron": ["RBFOX3", "SYT1", "SNAP25", "GAD1", "SLC17A7"],
}
ASTRO_MARKERS = {"astrocyte": ["AQP4", "GFAP"]}
TUMOUR_MARKERS = {
    "OPC_like": ["PDGFRA", "OLIG1", "OLIG2", "SOX10", "CSPG4"],
    "AC_like": ["GFAP", "AQP4", "APOE", "CLU"],
    "MES_like": ["CD44", "VIM", "CHI3L1", "ANXA1"],
    "cycling": ["MKI67", "TOP2A"],
}
REF_TYPES = ["myeloid", "T_NK", "oligodendrocyte", "vascular", "neuron"]
MIN_REF = 30           # minimum reference cells to run a sample
POOLED_REF_PER_TYPE = 300

# thresholds (log-normalised score_genes units; tuned on jessa_scrna P-6337_S-8821 vs author labels, see docs/phase3_cnv_notes.md;
# audit correction CNV-06: an earlier comment named P-6251, which is not a Jessa scRNA sample)
MARKER_THR = dict(normal_min=0.25, margin=0.10, ref_min=0.50, ref_ratio=2.0, astro_min=0.25)

# ---------------------------------------------------------------- genome
# approximate GRCh38 centromere (acen) start/end per chromosome (bp), typed from memory of UCSC hg38 cytoBand -- NOT downloaded (allowlist); arm boundary = midpoint
CENTROMERES_GRCH38 = {
    "chr1": (121_700_000, 125_100_000), "chr2": (91_800_000, 96_000_000), "chr3": (87_800_000, 94_000_000),
    "chr4": (48_200_000, 51_800_000), "chr5": (46_100_000, 51_400_000), "chr6": (58_500_000, 62_600_000),
    "chr7": (58_100_000, 62_100_000), "chr8": (43_200_000, 47_200_000), "chr9": (42_200_000, 45_500_000),
    "chr10": (38_000_000, 41_600_000), "chr11": (51_000_000, 55_800_000), "chr12": (33_200_000, 37_800_000),
    "chr13": (16_500_000, 18_900_000), "chr14": (16_100_000, 18_200_000), "chr15": (17_500_000, 20_500_000),
    "chr16": (35_300_000, 38_400_000), "chr17": (22_700_000, 27_400_000), "chr18": (15_400_000, 21_500_000),
    "chr19": (24_200_000, 28_100_000), "chr20": (25_700_000, 30_400_000), "chr21": (10_900_000, 13_000_000),
    "chr22": (13_700_000, 17_400_000), "chrX": (58_100_000, 63_800_000), "chrY": (10_300_000, 10_600_000),
}
BAND_4Q12 = ("chr4", 52_700_000, 59_500_000)
# DMG hallmark events: name -> (region, direction)  region = arm name or whole chromosome
# Audit correction 2026-10-19 (CNV-03): 'chr4q12_gain' removed from the calling set. With 100-gene windows the band is
# covered by 0-2 windows (median 1) that span ~40 Mb and measure PDGFRA/OPC-lineage expression against an immune
# reference (97-100 % of malignant cells in every Filbin patient). An exploratory, non-calling arm_4q12 window mean and
# gene-level PDGFRA/KIT deviations (gene_dev_PDGFRA, gene_dev_KIT) are still exported.
HALLMARK_EVENTS = {
    "chr1q_gain": ("1q", +1), "chr2_gain": ("2", +1), "chr7_gain": ("7", +1),
    "chr13q_loss": ("13q", -1), "chr14q_loss": ("14q", -1), "chr10q_loss": ("10q", -1), "chr17p_loss": ("17p", -1),
}
EXPLORATORY_GENES = ["PDGFRA", "KIT"]   # gene-level (unsmoothed) deviations exported for the dropped 4q12 event
MIN_ARM_GENES = 10      # an arm is evaluable (gets an arm mean) only if >= this many used genes map to it
LOW_CONF_WINDOWS = 100  # per-sample low-confidence flag (CNV-09): < 100 windows or < 3,000 used genes
LOW_CONF_GENES = 3000
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def load_gene_positions(gtf: Path = GENCODE_GTF, cache: Path = GENEPOS_CSV) -> pd.DataFrame:
    """Gene-level table from the GENCODE GTF: ensembl_id (unversioned), symbol, chromosome, start, end, gene_type."""
    if cache.exists():
        return pd.read_csv(cache)
    rx_id = re.compile(r'gene_id "([^"]+)"'); rx_name = re.compile(r'gene_name "([^"]+)"'); rx_type = re.compile(r'gene_type "([^"]+)"')
    rows = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "gene":
                continue
            rows.append((rx_id.search(f[8]).group(1).split(".")[0], rx_name.search(f[8]).group(1), f[0], int(f[3]), int(f[4]), rx_type.search(f[8]).group(1)))
    df = pd.DataFrame(rows, columns=["ensembl_id", "symbol", "chromosome", "start", "end", "gene_type"])
    df = df[~df.ensembl_id.str.endswith("_PAR_Y")]
    cache.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    return df


ACROCENTRIC = {"chr13", "chr14", "chr15", "chr21", "chr22"}   # p arms carry no protein-coding genes: whole chromosome -> q


def validate_arm_boundaries(genepos: pd.DataFrame) -> pd.DataFrame:
    """For each chromosome report the protein-coding gene-free gap that contains the assumed centromere midpoint.
    Phase 3 check (data/external/gencode/arm_boundaries_grch38_validation.csv): every midpoint lies inside a gap of
    >= 2.8 Mb (chr22: 0.16 Mb, irrelevant because acrocentric), so no gene's arm depends on the exact coordinate."""
    pc = genepos[genepos.gene_type == "protein_coding"]
    rows = []
    for c, (s, e) in CENTROMERES_GRCH38.items():
        mid = (s + e) / 2; g = pc[pc.chromosome == c]
        before = g[g.end < mid].end.max(); after = g[g.start > mid].start.min()
        rows.append((c, mid, before, after, (after - before) if pd.notna(before) else np.nan, int(((g.start <= mid) & (g.end >= mid)).sum())))
    return pd.DataFrame(rows, columns=["chromosome", "midpoint", "last_pc_gene_end_before", "first_pc_gene_start_after", "gap_bp", "genes_spanning_midpoint"])


def arm_of(chrom: str, pos: int) -> str:
    c = chrom.replace("chr", "")
    if chrom in ACROCENTRIC:
        return f"{c}q"
    cen = CENTROMERES_GRCH38.get(chrom)
    if cen is None:
        return c
    return f"{c}p" if pos < (cen[0] + cen[1]) / 2 else f"{c}q"


def annotate_var_positions(var: pd.DataFrame, genepos: pd.DataFrame) -> pd.DataFrame:
    """Add chromosome/start/end/arm to var (index = symbol, column ensembl_id). Ensembl ID match first, symbol second."""
    gp_id = genepos.drop_duplicates("ensembl_id").set_index("ensembl_id")
    gp_sym = genepos[genepos.chromosome.isin(AUTOSOMES + ["chrX", "chrY"])].drop_duplicates("symbol", keep=False).set_index("symbol")
    out = var.copy()
    eid = out["ensembl_id"].astype(str) if "ensembl_id" in out else pd.Series("", index=out.index)
    hit_id = eid.isin(gp_id.index)
    chrom = pd.Series(np.nan, index=out.index, dtype=object); start = pd.Series(np.nan, index=out.index); end = pd.Series(np.nan, index=out.index)
    chrom[hit_id] = gp_id.loc[eid[hit_id], "chromosome"].values
    start[hit_id] = gp_id.loc[eid[hit_id], "start"].values
    end[hit_id] = gp_id.loc[eid[hit_id], "end"].values
    hit_sym = ~hit_id & out.index.isin(gp_sym.index)
    chrom[hit_sym] = gp_sym.loc[out.index[hit_sym], "chromosome"].values
    start[hit_sym] = gp_sym.loc[out.index[hit_sym], "start"].values
    end[hit_sym] = gp_sym.loc[out.index[hit_sym], "end"].values
    out["chromosome"] = chrom; out["start"] = start; out["end"] = end
    out["pos_source"] = np.where(hit_id, "ensembl_id", np.where(hit_sym, "symbol", "none"))
    mid = ((start + end) / 2)
    out["arm"] = [arm_of(c, p) if isinstance(c, str) else np.nan for c, p in zip(chrom, mid)]
    return out


# ---------------------------------------------------------------- IO
HIGH_MITO_NUCLEI = {"jessa_scrna__P-6254_S-8499", "jessa_multiome__P-6253_S-8498", "jessa_multiome__P-1709_S-1709", "jessa_multiome__P-1779_S-1781"}
RELAXED_MITO = 20.0


def working_mask(obs: pd.DataFrame, study: str) -> tuple[pd.Series, pd.Series]:
    """Cells entering Phase 3: qc_pass (Jessa: also author_qc_pass); four cytoplasm-contaminated nuclei samples get mito < 20 %.
    Returns (mask, high_mito_nuclei flag)."""
    m = obs["qc_pass"].astype(bool).copy()
    hm = obs["sample_id"].isin(HIGH_MITO_NUCLEI)
    if hm.any():
        only_mito = obs["qc_fail_reason"].astype(str).eq("high_mito") & (obs["pct_mito"] < RELAXED_MITO)
        m = m | (hm & only_mito)
    if "author_qc_pass" in obs:
        m = m & obs["author_qc_pass"].fillna(False).astype(bool)
    return m, hm


def read_rows(path: Path, rows: np.ndarray, layer: str) -> ad.AnnData:
    """Read a row subset of one layer (sparse) + obs + var from a .h5ad without loading X."""
    import h5py
    from anndata.io import read_elem, sparse_dataset
    rows = np.sort(np.asarray(rows))
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"]); var = read_elem(f["var"])
        X = sparse_dataset(f["layers"][layer])[rows]
    X = sp.csr_matrix(X)
    return ad.AnnData(X=X, obs=obs.iloc[rows].copy(), var=var.copy())


def normalise(adata: ad.AnnData, counts_type: str) -> None:
    """In place: UMI/RSEM -> normalize_total(1e4) + log1p; TPM/CPM -> log1p(TPM/10) (Tirosh convention)."""
    adata.X = sp.csr_matrix(adata.X, dtype=np.float32)
    if counts_type == "tpm":
        if adata.X.max() < 30:
            # table already log-transformed. LaBelle 2025 deposits: per-cell max 9-12.4, consistent with log2(TPM/10+1)
            # (Filbin's linear TPM tables give per-cell max TPM 6e4-1.2e5 -> log2(TPM/10+1) = 12.5-13.6); rescale to natural log
            adata.X.data = (adata.X.data * np.log(2.0)).astype(np.float32)
            adata.uns["tpm_scale"] = "input already log2(TPM/10+1); multiplied by ln2"
        else:
            adata.X.data = np.log1p(adata.X.data / 10.0).astype(np.float32)
            adata.uns["tpm_scale"] = "linear TPM -> log1p(TPM/10)"
    else:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)


# ---------------------------------------------------------------- markers
def score_markers(adata: ad.AnnData, seed: int = 0) -> pd.DataFrame:
    """scanpy score_genes for every marker set on log-normalised X; missing genes dropped; a set with < 2 present genes -> NaN."""
    present = set(adata.var_names)
    out = pd.DataFrame(index=adata.obs_names)
    for name, genes in {**NORMAL_MARKERS, **ASTRO_MARKERS, **TUMOUR_MARKERS}.items():
        g = [x for x in genes if x in present]
        if len(g) < 2:
            out[f"score_{name}"] = np.nan
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.score_genes(adata, g, score_name="_tmp", random_state=seed, use_raw=False)
        out[f"score_{name}"] = adata.obs["_tmp"].values
    adata.obs.drop(columns="_tmp", inplace=True, errors="ignore")
    return out


def marker_calls(scores: pd.DataFrame, thr: dict = MARKER_THR) -> pd.DataFrame:
    """Marker-based call per cell.
    marker_call: 'normal' (confident non-astrocyte normal type, low tumour programs), 'astrocyte_like' (astro top, tumour
    programs not dominant), 'malignant_candidate' (everything else).  marker_type: best normal type.
    is_reference: strict subset of 'normal' used as inferCNV reference (score >= ref_min and >= ref_ratio x second-best normal type)."""
    S_n = scores[[f"score_{t}" for t in REF_TYPES]].fillna(-np.inf)
    s_n = S_n.max(axis=1); top = S_n.idxmax(axis=1).str.replace("score_", "")
    second = S_n.apply(lambda r: np.sort(r.values)[-2], axis=1)
    s_t = scores[["score_OPC_like", "score_AC_like", "score_MES_like"]].fillna(-np.inf).max(axis=1)
    s_a = scores["score_astrocyte"].fillna(-np.inf)
    s_opc_mes = scores[["score_OPC_like", "score_MES_like"]].fillna(-np.inf).max(axis=1)
    tumour_score = scores[["score_OPC_like", "score_AC_like", "score_MES_like"]].fillna(0).mean(axis=1)
    # tumour-program competitor per normal type, excluding programs whose markers the normal type shares:
    # myeloid express APOE/CLU (AC-like set); oligodendrocytes express OLIG1/2, SOX10 (OPC-like set);
    # vascular cells express VIM/CD44/ANXA1 (MES-like set); T/NK and neurons share none.
    COMPETITOR = {"myeloid": ["OPC_like", "MES_like"], "oligodendrocyte": ["AC_like", "MES_like"], "vascular": ["OPC_like", "AC_like"],
                  "T_NK": ["OPC_like", "AC_like", "MES_like"], "neuron": ["OPC_like", "AC_like", "MES_like"]}
    comp = np.full(len(scores), -np.inf)
    for t, progs in COMPETITOR.items():
        sel = (top == t).values
        comp[sel] = scores.loc[sel, [f"score_{p}" for p in progs]].fillna(-np.inf).max(axis=1).values
    normal = (s_n >= thr["normal_min"]) & (s_n >= comp + thr["margin"]) & (s_n >= s_a)
    # astrocyte-like: AQP4/GFAP high and not dominated by the OPC-like program (reactive astrocytes share MES-like genes)
    astro = ~normal & (s_a >= thr["astro_min"]) & (s_a > scores["score_OPC_like"].fillna(-np.inf) + thr["margin"])
    call = np.where(normal, "normal", np.where(astro, "astrocyte_like", "malignant_candidate"))
    is_ref = normal & (s_n >= thr["ref_min"]) & (s_n >= thr["ref_ratio"] * np.maximum(second, 0.05))

    return pd.DataFrame({"marker_call": call, "marker_type": np.where(normal, top, np.where(astro, "astrocyte", "none")),
                         "marker_normal_score": s_n.values, "marker_tumour_score": tumour_score.values, "is_reference": is_ref.values},
                        index=scores.index)


# ---------------------------------------------------------------- inferCNV
def _window_positions(var: pd.DataFrame, chr_pos: dict, n_windows: int, window: int, step: int, key_col: str = "chromosome") -> pd.DataFrame:
    """Map each smoothed column to chromosome / arm / genomic midpoint / number of genes (mirrors infercnvpy's running-mean indexing).
    `chr_pos` is keyed by the smoothing unit `var[key_col]`: 'chromosome' (chr17; arm of a window = arm of its midpoint, windows may
    straddle the centromere) or 'arm_key' (chr17p; audit correction CNV-02: every window lies within one arm).
    `in_4q12` marks windows whose midpoint lies in chr4:52.7-59.5 Mb (exploratory only)."""
    keys = sorted(chr_pos, key=lambda c: chr_pos[c]); starts = [chr_pos[c] for c in keys] + [n_windows]
    rows = []
    for key, a, b in zip(keys, starts[:-1], starts[1:]):
        g = var.loc[var[key_col] == key].sort_values("start")
        pos = ((g.start + g.end) / 2).values; G = len(g)
        c = g.chromosome.iloc[0]
        for k in range(b - a):
            if G < window:
                lo, hi = 0, G
            else:
                lo = k * step; hi = min(lo + window, G)
            mid = pos[lo:hi].mean()
            arm = g.arm.iloc[0] if key_col == "arm_key" else arm_of(c, mid)
            rows.append((c, arm, mid, hi - lo, (c == "chr4") and (BAND_4Q12[1] <= mid <= BAND_4Q12[2])))
    return pd.DataFrame(rows, columns=["chromosome", "arm", "pos", "n_genes", "in_4q12"])


def _smooth(d: np.ndarray, var_cnv: pd.DataFrame, window: int, step: int):
    """infercnvpy pyramidal running mean per smoothing unit (var_cnv['chromosome']), then per-cell median centring."""
    from infercnvpy.tl._infercnv import _running_mean_by_chromosome
    chr_pos, sm, _ = _running_mean_by_chromosome(np.clip(d, -3, 3), var_cnv, window, step, False)
    return chr_pos, sm - np.median(sm, axis=1, keepdims=True)


def run_infercnv(adata: ad.AnnData, ref_col: str, genepos: pd.DataFrame, window: int = 100, step: int = 10,
                 min_frac_cells: float = 0.10, noise_sd: float = 1.5, chunksize: int = 4000, min_ref_per_type: int = 10,
                 mode: str = "best_match") -> dict:
    """inferCNV-style smoothed copy-number profile on log-normalised `adata.X`.

    Genes: expressed in >= min_frac_cells of cells, autosomal, with a GENCODE position. Reference categories = values of
    adata.obs[ref_col] other than 'none' with >= min_ref_per_type cells; per-gene mean per category.
    mode='best_match' (default): for each cell the deviation (clipped to +-3, smoothed with infercnvpy's pyramidal running
      mean of `window` genes every `step`, median-centred per cell) is computed against EVERY reference category and the
      category with the smallest mean |deviation| is kept (lineage-matched baseline). Unlike infercnvpy's min/max envelope
      across categories, this does not widen with the per-gene sampling noise of small categories (the envelope zeroed most
      of the CNV signal once a 50-200-cell lineage category was added; Jessa P-6519 recall 0.97 -> 0.57).
    mode='envelope': infercnvpy's bounded difference (kept for comparison).
    Noise filtering: |x| < noise_sd * SD(all values) -> 0 (infercnvpy dynamic_threshold applied once over the whole sample).

    Two smoothed matrices are produced from the same per-cell deviations and reference choice (audit correction CNV-02):
      X_cnv  -- smoothing unit = chromosome (infercnvpy default). Drives cnv_score, PCA/Leiden clustering and the heatmaps;
                identical to the pre-audit pipeline so that the validated score/cluster layer is unchanged.
      X_arm  -- smoothing unit = chromosome ARM ('chr17p'): no window straddles a centromere. Drives arm_means (hallmark
                events). Arms with < MIN_ARM_GENES used genes are not evaluable (NaN). The noise threshold estimated on
                X_cnv is applied to both matrices (one noise level per sample).
    Returns dict(X_cnv dense float32 [cells x windows], windows, X_arm, windows_arm, cnv_score, arm_means, arm_genes, gene_dev,
    best_ref per cell, n_genes_used, ref_cats, chr_pos, chr_pos_arm, var_names)."""
    frac = np.asarray((adata.X > 0).mean(axis=0)).ravel()
    var = annotate_var_positions(adata.var, genepos)
    keep = (frac >= min_frac_cells) & var.chromosome.isin(AUTOSOMES).values
    sub = adata[:, keep].copy(); sub.var = var.loc[keep]
    labels = sub.obs[ref_col].astype(str)
    vc = labels.value_counts()
    ref_cats = [c for c in vc.index if c != "none" and vc[c] >= min_ref_per_type]
    if not ref_cats:
        raise ValueError("no reference categories")
    X = sp.csr_matrix(sub.X, dtype=np.float32)
    means = np.vstack([np.asarray(X[(labels == c).values].mean(axis=0)).ravel() for c in ref_cats]).astype(np.float32)
    sub.var = sub.var.assign(arm_key="chr" + sub.var["arm"].astype(str))
    var_chr = pd.DataFrame({"chromosome": sub.var["chromosome"].values, "start": sub.var["start"].values, "end": sub.var["end"].values}, index=sub.var_names)
    var_arm = pd.DataFrame({"chromosome": sub.var["arm_key"].values, "start": sub.var["start"].values, "end": sub.var["end"].values}, index=sub.var_names)
    out_chunks, arm_chunks, best = [], [], []
    for i in range(0, sub.shape[0], chunksize):
        xd = X[i:i + chunksize].toarray()
        if mode == "envelope":
            lo, hi = means.min(axis=0), means.max(axis=0)
            d = np.zeros_like(xd); above = xd > hi; below = xd < lo
            d[above] = (xd - hi)[above]; d[below] = (xd - lo)[below]
            chr_pos, res = _smooth(d, var_chr, window, step); chr_pos_arm, res_arm = _smooth(d, var_arm, window, step)
            b = np.zeros(len(xd), dtype=int)
        else:
            res = None; score = None; b = None
            for k in range(len(ref_cats)):
                chr_pos, r_k = _smooth(xd - means[k], var_chr, window, step); s_k = np.abs(r_k).mean(axis=1)
                if res is None:
                    res, score, b = r_k, s_k, np.zeros(len(xd), dtype=int)
                else:
                    better = s_k < score
                    res[better] = r_k[better]; score[better] = s_k[better]; b[better] = k
            res_arm = None
            for k in range(len(ref_cats)):
                sel = b == k
                if not sel.any():
                    continue
                chr_pos_arm, r_k = _smooth(xd[sel] - means[k], var_arm, window, step)
                if res_arm is None:
                    res_arm = np.zeros((len(xd), r_k.shape[1]), dtype=np.float32)
                res_arm[sel] = r_k
        out_chunks.append(res.astype(np.float32)); arm_chunks.append(res_arm.astype(np.float32)); best.append(b)
    Xc = np.vstack(out_chunks); Xa = np.vstack(arm_chunks); best = np.concatenate(best)
    if noise_sd:
        thr = noise_sd * Xc.std(); Xc[np.abs(Xc) < thr] = 0; Xa[np.abs(Xa) < thr] = 0
    win = _window_positions(sub.var, chr_pos, Xc.shape[1], window, step, key_col="chromosome")
    win_arm = _window_positions(sub.var, chr_pos_arm, Xa.shape[1], window, step, key_col="arm_key")
    cnv_score = np.abs(Xc).mean(axis=1)
    # arm means from the per-arm matrix; only arms with >= MIN_ARM_GENES used genes are evaluable (others NaN; recorded in arm_genes)
    arm_genes = sub.var.groupby("arm").size()
    arms = {}
    for arm in sorted(win_arm.arm.unique(), key=lambda a: (int(re.sub(r"[pq]", "", a)), a)):
        arms[arm] = Xa[:, (win_arm.arm == arm).values].mean(axis=1) if arm_genes.get(arm, 0) >= MIN_ARM_GENES else np.full(Xa.shape[0], np.nan, np.float32)
    for c in ["2", "7"]:
        arms[c] = Xa[:, (win_arm.chromosome == f"chr{c}").values].mean(axis=1)
    arms["4q12"] = Xa[:, win_arm.in_4q12.values].mean(axis=1) if win_arm.in_4q12.any() else np.full(Xa.shape[0], np.nan, np.float32)
    arm_means = pd.DataFrame(arms, index=sub.obs_names)
    # exploratory gene-level (unsmoothed, unfiltered) deviations from the cell's best-match reference mean (CNV-03)
    gene_dev = {}
    for gsym in EXPLORATORY_GENES:
        if gsym in sub.var_names:
            j = sub.var_names.get_loc(gsym)
            gene_dev[f"gene_dev_{gsym}"] = np.asarray(X[:, j].todense()).ravel() - means[best, j]
        else:
            gene_dev[f"gene_dev_{gsym}"] = np.full(Xc.shape[0], np.nan, np.float32)
    gene_dev = pd.DataFrame(gene_dev, index=sub.obs_names)
    return dict(X_cnv=Xc, windows=win, X_arm=Xa, windows_arm=win_arm, cnv_score=pd.Series(cnv_score, index=sub.obs_names), arm_means=arm_means,
                arm_genes=arm_genes, gene_dev=gene_dev, best_ref=pd.Series(np.array(ref_cats)[best], index=sub.obs_names), n_genes_used=int(keep.sum()),
                ref_cats=ref_cats, chr_pos=chr_pos, chr_pos_arm=chr_pos_arm, var_names=sub.var_names.tolist())


def refine_reference(ref_type: pd.Series, cnv_score: pd.Series, n_mad: float = 3.0) -> pd.Series:
    """Drop reference cells whose cnv_score exceeds median + n_mad*MAD of the reference (likely tumour cells with normal-like markers)."""
    is_ref = ref_type != "none"
    s = cnv_score[is_ref]; med = s.median(); mad = np.median(np.abs(s - med)) * 1.4826
    bad = is_ref & (cnv_score > med + n_mad * max(mad, 1e-6))
    out = ref_type.copy(); out[bad] = "none"
    return out


CALL_PARAMS = dict(score_sd=np.inf, z_single=2.0, z_double=1.0, any_arm_z=np.inf, min_sd=0.01, score_floor=0.008, ref_frac_guard=0.2)
# stage A ("tumour-grade" evidence only): used to purge tumour cells from oligodendrocyte/vascular/neuron reference candidates
CALL_PARAMS_STAGEA = dict(score_sd=np.inf, z_single=5.0, z_double=2.5, any_arm_z=6.0, min_sd=0.01, score_floor=0.02, ref_frac_guard=1.1)



def cnv_calls(X: np.ndarray, cnv_score: pd.Series, arm_means: pd.DataFrame, is_ref: pd.Series, score_sd: float = 2.0,
              z_single: float = 2.0, z_double: float = 1.0, any_arm_z: float = np.inf, min_sd: float = 0.01,
              score_floor: float = 0.01, ref_frac_guard: float = 0.2, resolution: float = 1.0, seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Leiden on CNV-PCA; per cluster, z = (cluster-mean arm value - reference mean) / max(reference per-cell SD, min_sd).
    Cluster = malignant if mean cnv_score > ref_mean + max(score_sd*ref_sd, score_floor), OR one hallmark event with z >= z_single
    in the expected direction, OR >= 2 hallmark events with z >= z_double (optionally any arm |z| >= any_arm_z). Clusters holding
    >= ref_frac_guard reference cells need score elevation AND an event (immune-subtype expression can mimic a weak single event).
    Everything else = normal (reference-like).
    Event flags per cell are inherited from the cluster (subclone-level calls).
    Returns (per-cell DataFrame: cnv_cluster, cnv_score, cnv_call, event flags; per-cluster summary)."""
    tmp = ad.AnnData(X=X, obs=pd.DataFrame(index=cnv_score.index))
    n_comps = int(min(30, X.shape[1] - 1, X.shape[0] - 1))
    sc.pp.pca(tmp, n_comps=n_comps, random_state=seed)
    sc.pp.neighbors(tmp, n_neighbors=15, n_pcs=n_comps, random_state=seed)
    sc.tl.leiden(tmp, resolution=resolution, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    cl = tmp.obs["leiden"].astype(str)
    ref = is_ref.values.astype(bool)
    ref_mu, ref_sd = float(cnv_score[ref].mean()), float(cnv_score[ref].std(ddof=1))
    # absolute floor: after noise filtering a clean reference has SD ~0.001, so SD-based thresholds flag trivial deviations.
    # The floor (0.008) was chosen during the benchmark iterations from author-labelled clusters of jessa_scrna, jessa_multiome,
    # okada, liu and collot (results/qc/cnv_cluster_stats_vs_author.csv, now regenerated by scripts/06_cnv_aggregate.py).
    # Audit correction CNV-20: recomputed percentiles of that file are normal-cluster 90th pct 0.0134 / tumour-cluster 10th pct
    # 0.0114 (all studies), not the 0.017 / 0.010 quoted before. Note (CNV-10): because the noise filter zeroes |x| < 1.5 SD before
    # reference statistics are taken, reference arm SDs collapse to min_sd and 'z >= 2' is in practice an absolute arm deviation
    # >= 0.02 log units.
    thr_score = ref_mu + (max(score_sd * ref_sd, score_floor) if np.isfinite(score_sd) else score_floor)
    ev = {e: (r, d) for e, (r, d) in HALLMARK_EVENTS.items() if r in arm_means}
    ref_arm_mu = arm_means.loc[ref].mean(); ref_arm_sd = arm_means.loc[ref].std(ddof=1).fillna(0).clip(lower=min_sd)
    rows = []; cell_events = {e: np.zeros(len(cl), bool) for e in ev}
    per_cell = pd.DataFrame({"cnv_cluster": cl.values, "cnv_score": cnv_score.values}, index=cnv_score.index)
    for c in sorted(cl.unique(), key=int):
        m = (cl == c).values
        mean_score = float(cnv_score[m].mean())
        z = (arm_means.loc[m].mean() - ref_arm_mu) / ref_arm_sd
        zs = {e: float(d * z[r]) for e, (r, d) in ev.items()}          # positive = expected direction
        single = [e for e, v in zs.items() if v >= z_single]
        double = [e for e, v in zs.items() if v >= z_double]
        events = single if single else (double if len(double) >= 2 else [])
        arm_cols = [c for c in arm_means.columns if c not in ("2", "7", "4q12")]
        zabs = (z[arm_cols]).abs(); other_arm = zabs.max() >= any_arm_z
        if ref[m].mean() >= ref_frac_guard:   # reference-rich cluster: demand both score elevation and an arm event
            call = "malignant" if (mean_score > thr_score and (events or other_arm)) else "normal"
        else:
            call = "malignant" if (mean_score > thr_score or events or other_arm) else "normal"
        for e in events:
            cell_events[e][m] = True
        # audit correction CNV-09: confidence of a malignant call = event_backed (hallmark event or any-arm flag) | score_only
        confidence = ("event_backed" if (events or other_arm) else "score_only") if call == "malignant" else ""
        rows.append(dict(cnv_cluster=c, n_cells=int(m.sum()), frac_reference=float(ref[m].mean()), mean_cnv_score=mean_score,
                         score_elevated=bool(mean_score > thr_score), events=";".join(events), max_other_arm=(zabs.idxmax() if zabs.notna().any() else ""),
                         max_other_arm_z=round(float(zabs.max()), 2) if zabs.notna().any() else np.nan,
                         any_arm_flag=bool(other_arm), cnv_call=call, cnv_confidence=confidence,
                         **{f"z_{e}": round(v, 2) for e, v in zs.items()}))
    clusters = pd.DataFrame(rows)
    for e in ev:
        per_cell[e] = cell_events[e]
    cl_idx = clusters.set_index("cnv_cluster")
    per_cell["cnv_call"] = per_cell["cnv_cluster"].map(cl_idx["cnv_call"])
    per_cell["cnv_confidence"] = per_cell["cnv_cluster"].map(cl_idx["cnv_confidence"]).replace("", np.nan)
    per_cell.attrs = dict(ref_mean=ref_mu, ref_sd=ref_sd, thr_score=thr_score, n_ref=int(ref.sum()))
    return per_cell, clusters


# ---------------------------------------------------------------- author labels & consensus
AUTHOR_MAP = {
    "jessa": {"Malignant": "malignant", "Normal": "normal"},
    "okada": {"Tumor": "malignant", "Monocytic": "normal", "Oligodendrocytes": "normal", "Neuron": "normal",
              "Astrocytes": "normal", "Endothelial-cells": "normal", "T-cells": "normal"},
    "collot": {"OPC_like": "malignant", "MES_like": "malignant", "AC_like": "malignant", "OC_like": "malignant",
               "Immune": "normal", "Macroglia": "normal", "Vascular": "normal", "Neuronal": "normal"},
    "liu": {"Myeloid": "normal", "Microglia": "normal", "Oligodendrocyte": "normal", "Endothelia": "normal", "T cell": "normal"},
    "labelle": {"Myeloid": "normal", "Tcell": "normal"},
}
AUTHOR_COARSE = {
    "okada": {"Monocytic": "myeloid", "Oligodendrocytes": "oligodendrocyte", "Neuron": "neuron", "Astrocytes": "astrocyte",
              "Endothelial-cells": "vascular", "T-cells": "T_NK"},
    "collot": {"Immune": "myeloid", "Macroglia": "other", "Vascular": "vascular", "Neuronal": "neuron"},
    "liu": {"Myeloid": "myeloid", "Microglia": "myeloid", "Oligodendrocyte": "oligodendrocyte", "Endothelia": "vascular", "T cell": "T_NK"},
}
LIU_MALIGNANT_RX = re.compile(r"OPC-like|AC-like|MES-like|OC-like|^S$|^G2M$|frozen-")


def author_binary(obs: pd.DataFrame, study: str) -> tuple[pd.Series, pd.Series]:
    """(author_label raw string, author_call in {malignant, normal, NaN})."""
    if study.startswith("jessa"):
        lab = obs["author_Malignant_normal_consensus_Jessa2022"].astype(object)
        return lab, lab.map(AUTHOR_MAP["jessa"])
    if study == "okada":
        lab = obs["author_cellType"].astype(object); return lab, lab.map(AUTHOR_MAP["okada"])
    if study == "collot":
        lab = obs["author_celltype"].astype(object); return lab, lab.map(AUTHOR_MAP["collot"])
    if study == "liu":
        lab = obs["author_annotation"].astype(object)
        call = lab.map(AUTHOR_MAP["liu"])
        mal = lab.astype(str).str.contains(LIU_MALIGNANT_RX)
        call[mal] = "malignant"
        return lab, call
    if study == "labelle":
        lab = obs["author_broad_annot"].astype(object); return lab, lab.map(AUTHOR_MAP["labelle"])
    lab = pd.Series(np.nan, index=obs.index, dtype=object)
    return lab, lab.copy()


def consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Tier-3 consensus from cnv_call, marker_call (+ author_call where present).
    Rules: marker 'normal' agrees with cnv 'normal'; 'malignant_candidate' agrees with cnv 'malignant';
    'astrocyte_like' + cnv 'normal' -> normal astrocyte; 'astrocyte_like' + cnv 'malignant' -> malignant only if the CNV
    cluster carries a hallmark arm event (normal astrocytes have no CNV reference and can reach the cnv_score threshold
    through lineage expression alone), otherwise ambiguous; author label, when present, must agree too; else 'ambiguous'."""
    cnvc = df["cnv_call"].astype(str); mk = df["marker_call"].astype(str); au = df["author_call"]
    ev_cols = [e for e in HALLMARK_EVENTS if e in df.columns]
    has_event = df[ev_cols].fillna(False).astype(bool).any(axis=1).values if ev_cols else np.zeros(len(df), bool)
    astro_bin = np.where(cnvc.values == "normal", "normal", np.where(has_event, "malignant", "ambiguous"))
    marker_bin = np.where(mk == "normal", "normal", np.where(mk == "malignant_candidate", "malignant", astro_bin))
    agree = (marker_bin == cnvc.values) & (cnvc.values != "unassigned")
    has_author = au.notna().values
    agree = agree & (~has_author | (au.astype(str).values == cnvc.values))
    call = np.where(agree, cnvc.values, "ambiguous")
    coarse = np.where(call == "malignant", "malignant", np.where(call == "ambiguous", "ambiguous", df["marker_type"].replace("none", "other").values))
    out = df.copy(); out["consensus_call"] = call; out["consensus_celltype_coarse"] = coarse
    # Audit correction CNV-05 (sensitivity label, does NOT replace consensus_call): marker-normal / astrocyte-like cells that sit
    # in an EVENT-BACKED CNV-malignant cluster are trusted as malignant when no author label contradicts them. The audit showed
    # 79 % of cnv=malignant & marker=normal cells with an author label are author-malignant (OC-/neuron-like tumour states).
    trust = (cnvc.values == "malignant") & has_event & (marker_bin == "normal") & (~has_author | (au.astype(str).values == "malignant"))
    out["consensus_call_eventtrust"] = np.where(trust, "malignant", call)
    return out


# ---------------------------------------------------------------- per-sample driver
IMMUNE_TYPES = ["myeloid", "T_NK"]
LINEAGE_TYPES = ["oligodendrocyte", "vascular", "neuron"]   # astrocytes excluded: a pooled "astrocyte" category was dominated by
                                                           # AC-like tumour cells and masked AC-like tumour (Jessa snRNA recall 0.88 -> 0.69)
POOL_PER_SAMPLE_PER_TYPE = 50   # pooled lineage reference: at most this many cells per donor sample per type (no donor dominates)
MIN_TYPE_REF = 50               # pooled lineage category needs >= this many cells to enter the bounded reference


def purge_candidates(a: ad.AnnData, genepos: pd.DataFrame, ref_col: str = "ref_type", resolution: float = 1.0, seed: int = 0) -> pd.Series:
    """Stage A for pool construction: inferCNV against the immune reference (own immune candidates + pooled immune anchor,
    MAD-refined), tumour-grade cluster calls (CALL_PARAMS_STAGEA). Returns, for the sample's own cells, True where a
    non-immune reference candidate (oligodendrocyte / vascular / neuron) does NOT sit in a tumour-grade cluster."""
    rt = a.obs[ref_col].astype(str)
    rtA = rt.where(rt.isin(IMMUNE_TYPES), "none")
    if int((rtA != "none").sum()) < MIN_REF:
        raise ValueError("too few immune reference cells for purge")
    a.obs["_refA"] = rtA.values
    res = run_infercnv(a, "_refA", genepos)
    rtA = refine_reference(rtA, res["cnv_score"])
    pc, _ = cnv_calls(res["X_cnv"], res["cnv_score"], res["arm_means"], rtA != "none", resolution=resolution, seed=seed, **CALL_PARAMS_STAGEA)
    a.obs.drop(columns="_refA", inplace=True)
    own = ~a.obs["_pooled"].values
    cand = rt.isin(LINEAGE_TYPES).values & own
    ok = cand & (pc["cnv_call"] != "malignant").values
    return pd.Series(ok[own], index=a.obs_names[own])


def call_sample(a: ad.AnnData, genepos: pd.DataFrame, ref_col: str = "ref_type", resolution: float = 1.0, seed: int = 0,
                log=None) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    """Final inferCNV calling for one sample. `a` holds the sample's own cells (obs['_pooled'] False) plus pooled reference
    cells from OTHER samples of the same study/assay stratum (obs['_pooled'] True): immune anchor (myeloid, T_NK, merged with
    the sample's own immune candidates) and purged lineage categories (oligodendrocyte / vascular / neuron). The sample's own
    non-immune candidates are never used as reference (in near-pure tumours they are tumour cells with glial/neuronal marker
    scores; see docs/phase3_cnv_notes.md benchmark). Reference is MAD-refined on cnv_score; pooled lineage categories with
    < MIN_TYPE_REF cells are dropped. Returns (per_cell, clusters, res, record)."""
    rt = a.obs[ref_col].astype(str)
    own = ~a.obs["_pooled"].values
    rt = rt.where(rt.isin(IMMUNE_TYPES) | ~own, "none")          # own cells: immune only
    for t in LINEAGE_TYPES:
        if int((rt == t).sum()) < MIN_TYPE_REF:
            rt = rt.where(rt != t, "none")
    rec = dict(n_ref_immune_own=int((rt[own].isin(IMMUNE_TYPES)).sum()), n_ref_pooled=int((rt[~own] != "none").sum()),
               pooled_types=rt[~own][rt[~own] != "none"].value_counts().to_dict())
    if int((rt != "none").sum()) < MIN_REF:
        raise ValueError(f"only {int((rt != 'none').sum())} reference cells")
    a.obs["_refB"] = rt.values
    res = run_infercnv(a, "_refB", genepos)
    ref_final = refine_reference(a.obs["_refB"], res["cnv_score"])
    pc, clusters = cnv_calls(res["X_cnv"], res["cnv_score"], res["arm_means"], ref_final != "none", resolution=resolution, seed=seed, **CALL_PARAMS)
    pc["is_reference_final"] = (ref_final != "none").values
    pc["ref_type_final"] = ref_final.values
    pc["cnv_best_ref"] = res["best_ref"].values
    for col in res["gene_dev"].columns:
        pc[col] = res["gene_dev"][col].values
    n_genes, n_windows = res["n_genes_used"], int(res["X_cnv"].shape[1]); n_windows_arm = int(res["X_arm"].shape[1])
    arm_ok = res["arm_means"].columns[res["arm_means"].notna().any()].tolist()
    hall_arms = {r for r, _ in HALLMARK_EVENTS.values()}
    # audit correction CNV-09 / CNV-02: per-sample low-confidence flag and the list of evaluable arms
    rec.update(n_ref_final=int((ref_final != "none").sum()), ref_cats=";".join(res["ref_cats"]), n_genes=n_genes,
               n_windows=n_windows, n_windows_arm=n_windows_arm, ref_cnv_mean=pc.attrs["ref_mean"], ref_cnv_sd=pc.attrs["ref_sd"],
               thr_cnv_score=pc.attrs["thr_score"], n_clusters=len(clusters), n_clusters_malignant=int((clusters.cnv_call == "malignant").sum()),
               n_clusters_event_backed=int((clusters.cnv_confidence == "event_backed").sum()),
               n_clusters_score_only=int((clusters.cnv_confidence == "score_only").sum()),
               low_confidence_sample=bool(n_windows < LOW_CONF_WINDOWS or n_genes < LOW_CONF_GENES),
               arms_evaluable=";".join(a for a in arm_ok if a not in ("2", "7", "4q12")),
               hallmark_regions_missing=";".join(sorted(hall_arms - set(arm_ok))),
               arm_genes=str(res["arm_genes"].to_dict()))
    pc["low_confidence_sample"] = bool(n_windows < LOW_CONF_WINDOWS or n_genes < LOW_CONF_GENES)
    res["ref_final"] = ref_final
    a.obs.drop(columns="_refB", inplace=True)
    return pc, clusters, res, rec
