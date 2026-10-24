#!/usr/bin/env python
"""Phase 7 / step 25 — OpenPedCan (v15) DMG bulk RNA-seq: deconvolution with the atlas signatures and survival models.

Downloads (if missing; public AWS bucket d3b-openaccess-us-east-1-prd-pbta, open-targets/v15/) into data/raw/openpedcan/:
  histologies.tsv, gene-expression-rsem-tpm-collapsed.rds, independent-specimens.rnaseqpanel.primary.tsv (+ md5 from md5sum.txt).
Inputs : results/tables/c2l_signatures.csv (fallback states) and c2l_signatures_mp.csv (metaprogram families) — cell2location
         reference signatures (genes x classes); results/tables/metaprogram_families.csv + metaprogram_genes.csv;
         data/processed/sample_manifest.csv (CBTN participant IDs already in the atlas, to exclude).
Outputs (results/tables/, never overwritten if present, unless --force):
  openpedcan_dmg_cohort.csv       selected biospecimens: participant, age, site, molecular subtype, OS, exclusion reason
  openpedcan_deconvolution.csv    per biospecimen: NNLS proportions (both signature sets) and MP-family rank scores
  survival_results.csv            Cox HR per SD (univariate and adjusted for age + pons), 95 % CI, p, BH q, n, events
Figure (results/figures/): phase7_survival_forest.png.
Run: python scripts/25_openpedcan_survival.py [--force]    (needs pyreadr and lifelines: pip install pyreadr lifelines)

Method (fixed before running)
  Cohort: histologies rows with experimental_strategy == 'RNA-Seq', sample_type == 'Tumor', molecular_subtype containing 'H3 K28'
  (OpenPedCan's DMG H3 K27-altered subtypes), restricted to the independent primary RNA-seq specimen list (one per participant);
  participants whose cohort_participant_id matches a CBTN ID in the atlas manifest are excluded (ScPCA SCPCP000021 cases); survival =
  OS_days with OS_status 'DECEASED' as event; specimens without OS are kept for deconvolution but dropped from Cox models.
  Deconvolution: TPM restricted to signature genes; columns of the signature matrix and each bulk profile scaled to sum 1; non-negative
  least squares; weights rescaled to proportions (mRNA-weighted, not cell-number fractions). Run twice: all shared signature genes, and
  a marker-restricted signature (per class the 50 genes with the highest class / max-other-class ratio among genes above the class median;
  added after the all-gene fit gave implausible proportions and disagreed with the rank scores — collinear malignant signatures). Malignant-compartment share of a family =
  its proportion / sum of malignant proportions.
  MP-family rank score: per sample, genes ranked by TPM; score = mean normalised rank (0-1) of the family's MP genes (singscore-like).
  Cox (lifelines CoxPHFitter): each variable z-scored across the cohort; univariate, and adjusted for age at diagnosis (years) and pons
  location (primary_site mentions pons / brain stem). BH within each variable set x model.
"""
import re, sys, argparse, hashlib, urllib.request
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd  # noqa: E402
import scipy.stats as ss  # noqa: E402
from scipy.optimize import nnls  # noqa: E402
from dmg_atlas import RAW, PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
DL = RAW / "openpedcan"; DL.mkdir(parents=True, exist_ok=True)
BASE = "https://d3b-openaccess-us-east-1-prd-pbta.s3.amazonaws.com/open-targets/v15/"
FILES = ["histologies.tsv", "gene-expression-rsem-tpm-collapsed.rds", "independent-specimens.rnaseqpanel.primary.tsv"]
log = U.StepLogger(LOGS / "phase7_multimodal.log", context="25_openpedcan_survival")
N_MARKERS = 50
OUTS = [TABLES / f for f in ("openpedcan_dmg_cohort.csv", "openpedcan_deconvolution.csv", "survival_results.csv")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)
import pyreadr  # noqa: E402
from lifelines import CoxPHFitter  # noqa: E402


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


with log.step("download"):
    sums = {}
    for name in ("md5sum.txt", "analysis_files_release_md5sum.txt"):
        try:
            for line in urllib.request.urlopen(BASE + name, timeout=60).read().decode().splitlines():
                parts = line.split()
                if len(parts) == 2:
                    sums[Path(parts[1]).name] = parts[0]
        except Exception as e:  # md5 lists are optional
            log.write(f"  {name}: {type(e).__name__}")
    for fn in FILES:
        p = DL / fn
        if not p.exists():
            urllib.request.urlretrieve(BASE + fn, p)
        got = md5(p); exp = sums.get(fn)
        log.write(f"  {fn}: {p.stat().st_size / 1e6:.1f} MB md5 {got} " + ("(matches release md5)" if exp == got else f"(release md5 {exp or 'n/a'})"))

with log.step("cohort"):
    h = pd.read_csv(DL / "histologies.tsv", sep="\t", low_memory=False)
    ind = pd.read_csv(DL / "independent-specimens.rnaseqpanel.primary.tsv", sep="\t")
    sel = h[(h.experimental_strategy == "RNA-Seq") & (h.sample_type == "Tumor") & h.molecular_subtype.astype(str).str.contains("H3 K28")]
    log.write(f"histologies {h.shape}; RNA-Seq tumour H3 K28 rows {len(sel)}; molecular subtypes {sel.molecular_subtype.value_counts().to_dict()}")
    sel = sel[sel.Kids_First_Biospecimen_ID.isin(ind.Kids_First_Biospecimen_ID)].copy()
    man = pd.read_csv(PROCESSED / "sample_manifest.csv")
    atlas_cbtn = {m.replace("CBTN-", "") for m in man.patient_id.astype(str) if m.startswith("CBTN-")}
    sel["excluded"] = np.where(sel.cohort_participant_id.astype(str).isin(atlas_cbtn), "in_atlas_scpca", "")
    sel["age_years"] = pd.to_numeric(sel.age_at_diagnosis_days, errors="coerce") / 365.25
    sel["os_days"] = pd.to_numeric(sel.OS_days, errors="coerce"); sel["os_event"] = (sel.OS_status.astype(str).str.upper() == "DECEASED").astype(int)
    sel["pons"] = sel.primary_site.astype(str).str.contains("Pons|Brain Stem|Brainstem", case=False).astype(int)
    keep_cols = ["Kids_First_Biospecimen_ID", "Kids_First_Participant_ID", "cohort_participant_id", "cohort", "primary_site", "pons", "molecular_subtype",
                 "tumor_descriptor", "age_years", "os_days", "OS_status", "os_event", "excluded"]
    cohort = sel[[c for c in keep_cols if c in sel.columns]].reset_index(drop=True)
    log.write(f"independent primary H3 K28 RNA-seq specimens {len(cohort)}; excluded (in atlas) {int((cohort.excluded != '').sum())}; "
              f"with OS {int(cohort.os_days.notna().sum())}, events {int(cohort.os_event[cohort.os_days.notna()].sum())}; pons {int(cohort.pons.sum())}; "
              f"cohorts {cohort.cohort.value_counts().to_dict() if 'cohort' in cohort else {}}")

with log.step("read TPM"):
    tpm = next(iter(pyreadr.read_r(DL / "gene-expression-rsem-tpm-collapsed.rds").values()))
    ids = [b for b in cohort.Kids_First_Biospecimen_ID if b in tpm.columns]
    T = tpm[ids].astype(float); del tpm
    log.write(f"TPM {T.shape} for the cohort ({len(cohort) - len(ids)} specimens missing from the matrix)")

with log.step("deconvolution + rank scores"):
    out = pd.DataFrame(index=ids)
    for tag, fn, markers in (("fallback", "c2l_signatures.csv", False), ("mp", "c2l_signatures_mp.csv", False),
                             ("fallbackM", "c2l_signatures.csv", True), ("mpM", "c2l_signatures_mp.csv", True)):
        S = pd.read_csv(TABLES / fn, index_col=0)
        S.columns = [re.sub(r"^(means_per_cluster_mu_fg_|mean_)", "", c) for c in S.columns]
        g = S.index.intersection(T.index)
        if markers:  # marker-restricted signature: per class the N_MARKERS most specific genes (class value / max of the other classes)
            Sg = S.loc[g]; mk = set()
            for c in Sg.columns:
                other = Sg.drop(columns=c).max(1); ok = Sg[c] > Sg[c].median()
                mk |= set((Sg.loc[ok, c] / (other[ok] + 1e-3)).sort_values(ascending=False).index[:N_MARKERS])
            g = pd.Index(sorted(mk))
        A = (S.loc[g] / S.loc[g].sum()).values
        W = np.array([nnls(A, (T.loc[g, b] / T.loc[g, b].sum()).values)[0] for b in ids])
        Pp = pd.DataFrame(W / W.sum(1, keepdims=True), index=ids, columns=[f"{tag}:{c}" for c in S.columns])
        mal = [c for c in Pp.columns if ":malignant_" in c]
        Pp[f"{tag}:malignant_total"] = Pp[mal].sum(1)
        for c in mal:
            Pp[c.replace(":malignant_", ":share_")] = Pp[c] / Pp[mal].sum(1).replace(0, np.nan)
        out = out.join(Pp)
        log.write(f"  {tag}: {len(g)} signature genes in TPM; median proportions " + Pp.median().round(3).to_dict().__repr__())
    fm = pd.read_csv(TABLES / "metaprogram_families.csv").dropna(subset=["family"]); mg = pd.read_csv(TABLES / "metaprogram_genes.csv")
    R = T.rank(axis=0) / len(T)
    for fam, g in fm.groupby("family"):
        genes = [x for x in set(mg.gene[mg.metaprogram.isin(g.metaprogram)]) if x in R.index]
        out[f"rank:{fam}"] = R.loc[genes].mean().reindex(ids).values
    for tag in ("mp", "mpM"):
        agree = {fam: ss.spearmanr(out[f"{tag}:malignant_{fam}"], out[f"rank:{fam}"]).statistic for fam in fm.family.unique() if f"{tag}:malignant_{fam}" in out}
        log.write(f"Spearman, NNLS ({tag}) proportion vs rank score per family: " + {k: round(float(v), 2) for k, v in agree.items()}.__repr__())
    out.index.name = "Kids_First_Biospecimen_ID"

with log.step("Cox models"):
    d = cohort.set_index("Kids_First_Biospecimen_ID").join(out, how="inner")
    d = d[(d.excluded == "") & d.os_days.notna() & (d.os_days > 0)]
    sets = {"mp_proportion": [c for c in out.columns if c.startswith("mp:malignant_") and not c.endswith("total")] + ["mp:malignant_total"],
            "mp_share": [c for c in out.columns if c.startswith("mp:share_")],
            "fallback_proportion": [c for c in out.columns if c.startswith("fallback:malignant_")],
            "tme_proportion": [c for c in out.columns if c.startswith("mp:") and ":malignant_" not in c and ":share_" not in c],
            "rank_score": [c for c in out.columns if c.startswith("rank:")],
            "mpM_proportion": [c for c in out.columns if c.startswith("mpM:malignant_")],
            "mpM_share": [c for c in out.columns if c.startswith("mpM:share_")],
            "tmeM_proportion": [c for c in out.columns if c.startswith("mpM:") and ":malignant_" not in c and ":share_" not in c]}
    res = []
    for vset, cols in sets.items():
        for c in cols:
            for model in ("univariate", "adjusted"):
                cov = ["age_years", "pons"] if model == "adjusted" else []
                dd = d[["os_days", "os_event", c] + cov].dropna().copy()
                if dd[c].std() == 0 or len(dd) < 20:
                    continue
                dd[c] = (dd[c] - dd[c].mean()) / dd[c].std()
                try:
                    cph = CoxPHFitter().fit(dd, "os_days", "os_event")
                    s = cph.summary.loc[c]
                    res.append(dict(variable_set=vset, variable=c, model=model, n=len(dd), events=int(dd.os_event.sum()), hr_per_sd=s["exp(coef)"],
                                    ci_lo=s["exp(coef) lower 95%"], ci_hi=s["exp(coef) upper 95%"], p=s["p"], concordance=cph.concordance_index_))
                except Exception as e:
                    res.append(dict(variable_set=vset, variable=c, model=model, n=len(dd), error=type(e).__name__))
    res = pd.DataFrame(res); res["q_bh"] = np.nan
    for _, g in res[res.p.notna()].groupby(["variable_set", "model"]):
        res.loc[g.index, "q_bh"] = ss.false_discovery_control(g.p.values, method="bh")
    log.write(f"Cox cohort n={len(d)}, events {int(d.os_event.sum())}, median OS {d.os_days.median():.0f} d\n" +
              res[res.model == "adjusted"][["variable_set", "variable", "n", "events", "hr_per_sd", "ci_lo", "ci_hi", "p", "q_bh"]].round(3).to_string(index=False))

with log.step("write + figure"):
    cohort.to_csv(OUTS[0], index=False); out.to_csv(OUTS[1]); res.to_csv(OUTS[2], index=False)
    apply_figure_style()
    r = res[(res.model == "adjusted") & res.variable_set.isin(["mpM_proportion", "tmeM_proportion"]) & res.p.notna()].sort_values("hr_per_sd")
    fig, ax = plt.subplots(figsize=(3.6, 0.18 * len(r) + 1.0))
    y = np.arange(len(r))
    ax.errorbar(r.hr_per_sd, y, xerr=[r.hr_per_sd - r.ci_lo, r.ci_hi - r.hr_per_sd], fmt="o", ms=3, color="black", lw=0.8)
    ax.axvline(1, color="grey", lw=0.5, ls=":"); ax.set_xscale("log")
    ax.set_yticks(y); ax.set_yticklabels([v.split(":", 1)[1] + (" *" if q < 0.05 else "") for v, q in zip(r.variable, r.q_bh)], fontsize=5)
    ax.set_xlabel(f"hazard ratio per SD (adjusted for age + pons; n = {int(r.n.max()) if len(r) else 0})")
    fig.savefig(FIGS / "phase7_survival_forest.png"); plt.close(fig)
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
