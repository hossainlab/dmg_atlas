#!/usr/bin/env python
"""Phase 6 / step 21 — brainstem infiltration gradient on metaprogram families (metaprogram-based cell2location reference).

Upstream (Modal): 13_reference_signatures.py prepare --use-metaprograms --tag mp -> modal/phase6_c2l.py signatures (GPU) ->
  13_reference_signatures.py collect --tag mp (results/tables/c2l_signatures_mp.csv) -> modal/phase6_c2l.py deconvolve
  --vis visium_brainstem.h5ad --sig c2l_signatures_mp.csv --cohort brainstem_mp (GPU).
Input : --harvest <dir> = harvested out/ of that deconvolve job (<section>_abundance.csv, brainstem_mp_run.json, history),
        results/tables/brainstem_section_means.csv (fallback-state run, 15_spatial_tests.py) for comparison.
Outputs (never overwritten if present, unless --force):
  data/processed/spatial_mp/brainstem/<section>_abundance.csv (copied), results/tables/c2l_deconvolution_runs_mp.json,
  results/tables/brainstem_gradient_mp.csv        same tests as 15_spatial_tests.py::brainstem, classes = 6 TME + MP families
  results/tables/brainstem_section_means_mp.csv   per section: anatomical rank, tumour fraction, mean proportion/abundance
  results/tables/brainstem_mp_vs_fallback.csv     per section: tumour fraction under both references (agreement check)
Figure (results/figures/): phase6_brainstem_gradient_mp.png.
Run: nice -n 10 python scripts/21_brainstem_gradient_mp.py --harvest <dir> [--force]

Tests (identical to 15_spatial_tests.py, q05 abundances, section = unit): Spearman of anatomical rank (midbrain 1 -> spinal cord B 10)
with the section mean of absolute q05 and of the per-spot proportion; within section, Spearman of spot tumour fraction with each
family's share of the malignant compartment (and with TME absolute abundance), pooled by DerSimonian-Laird on Fisher z.
"""
import sys, json, shutil, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd  # noqa: E402
import scipy.stats as ss  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--harvest", required=True); ap.add_argument("--force", action="store_true")
args = ap.parse_args()
TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
DEST = PROCESSED / "spatial_mp" / "brainstem"
log = U.StepLogger(LOGS / "phase6_spatial.log", context="21_brainstem_gradient_mp")
TME = ["myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular"]
ORDER = ["Midbrain", "Pons_A", "Pons_B", "Pons_C", "Pons_D", "Medulla_A", "Medulla_B", "Medulla_C", "Spinal_Cord_A", "Spinal_Cord_B"]
OUTS = [TABLES / f for f in ("brainstem_gradient_mp.csv", "brainstem_section_means_mp.csv", "brainstem_mp_vs_fallback.csv")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)


def dl_random_effects(z, v):
    """DerSimonian-Laird pooled estimate (same as 15_spatial_tests.py)."""
    w = 1 / v; zf = np.sum(w * z) / np.sum(w); Q = np.sum(w * (z - zf) ** 2); df = len(z) - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w); tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    ws = 1 / (v + tau2); zr = np.sum(ws * z) / np.sum(ws); se = np.sqrt(1 / np.sum(ws))
    return zr, se, tau2


def pooled(rows, test, cls):
    rr = pd.DataFrame(rows)
    z = np.arctanh(np.clip(rr.rho.values, -0.999, 0.999)); v = 1 / (rr.n_spots.values - 3)
    zr, se, tau2 = dl_random_effects(z, v)
    return dict(test=test, cls=cls, n_sections=len(rr), rho_re=np.tanh(zr), rho_re_lo=np.tanh(zr - 1.96 * se), rho_re_hi=np.tanh(zr + 1.96 * se),
                p_re=2 * ss.norm.sf(abs(zr / se)), tau2=tau2, n_sections_positive=int((rr.rho > 0).sum()))


with log.step("collect"):
    h = Path(args.harvest); DEST.mkdir(parents=True, exist_ok=True)
    files = sorted(h.glob("*_abundance.csv"))
    for fp in files:
        shutil.copy2(fp, DEST / fp.name)
    rj = h / "brainstem_mp_run.json"
    if rj.exists():
        json.dump({"brainstem_mp": json.load(open(rj))}, open(TABLES / "c2l_deconvolution_runs_mp.json", "w"), indent=1)
    df = pd.concat([pd.read_csv(fp, index_col=0) for fp in files])
    q = df[[c for c in df.columns if c.startswith("q05_")]].rename(columns=lambda c: c[4:])
    fams = [c for c in q.columns if c.startswith("malignant_")]
    sec = df.section_id
    log.write(f"{len(files)} sections, {len(q)} spots; classes {list(q.columns)}")

with log.step("gradient tests"):
    rank = {s: i + 1 for i, s in enumerate(ORDER)}; x = np.arange(1, 11)
    P = q.div(q.sum(1), axis=0); tf = q[fams].sum(1) / q.sum(1)
    m_abs = q.groupby(sec).mean().loc[ORDER]; m_prop = P.groupby(sec).mean().loc[ORDER]
    tab = m_prop.add_prefix("prop_").join(m_abs.add_prefix("abs_"))
    tab.insert(0, "anatomical_rank", [rank[s] for s in tab.index]); tab.insert(1, "n_spots", sec.value_counts().loc[ORDER].values)
    tab.insert(2, "tumour_fraction_mean", tf.groupby(sec).mean().loc[ORDER].values)
    rows = []
    for c in list(q.columns) + ["tumour_fraction"]:
        ya = tab.tumour_fraction_mean.values if c == "tumour_fraction" else m_abs[c].values
        yp = ya if c == "tumour_fraction" else m_prop[c].values
        ra, pa = ss.spearmanr(x, ya); rp, pp = ss.spearmanr(x, yp)
        rows.append(dict(test="trend_along_axis", cls=c, n_sections=10, rho_abs=ra, p_abs=pa, rho_prop=rp, p_prop=pp,
                         rostral_mean=np.mean(yp[:5]), caudal_mean=np.mean(yp[5:])))
    share = q[fams].div(q[fams].sum(1), axis=0)
    for c, test, Y in [(f, "within_section_tumour_fraction_vs_family_share", share) for f in fams] + \
                      [(t, "within_section_tumour_fraction_vs_TME_abundance", q) for t in TME if t in q.columns]:
        per = []
        for s in ORDER:
            m = (sec == s).values
            r, p = ss.spearmanr(tf[m], Y.loc[m, c]); per.append(dict(test=test, cls=c, section_id=s, n_spots=int(m.sum()), rho=r, p=p))
        rows.append(pooled(per, test + "_RE", c)); rows.extend(per)
    grad = pd.DataFrame(rows)
    tr = grad.test == "trend_along_axis"
    grad.loc[tr, "q_bh_abs"] = ss.false_discovery_control(grad.loc[tr, "p_abs"].values, method="bh")
    log.write("trend along axis:\n" + grad[tr][["cls", "rho_abs", "p_abs", "q_bh_abs", "rho_prop", "p_prop", "rostral_mean", "caudal_mean"]].round(3).to_string(index=False))
    log.write("within-section RE:\n" + grad[grad.test.str.endswith("_RE")][["test", "cls", "rho_re", "rho_re_lo", "rho_re_hi", "p_re", "n_sections_positive"]].round(3).to_string(index=False))

with log.step("compare with fallback reference"):
    fb = pd.read_csv(TABLES / "brainstem_section_means.csv", index_col=0)
    cmp_ = pd.DataFrame({"tumour_fraction_mp": tab.tumour_fraction_mean, "tumour_fraction_fallback": fb.tumour_fraction_mean.reindex(ORDER)})
    r = ss.spearmanr(cmp_.tumour_fraction_mp, cmp_.tumour_fraction_fallback).statistic
    log.write(f"section tumour fraction, MP vs fallback reference: Spearman {r:.2f}\n" + cmp_.round(3).to_string())

with log.step("write + figure"):
    grad.to_csv(OUTS[0], index=False); tab.to_csv(OUTS[1]); cmp_.to_csv(OUTS[2])
    apply_figure_style()
    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(figsize=(4.6, 2.5))
    for i, fm in enumerate(fams):
        ax.plot(x, m_prop[fm].values, marker="o", ms=3, color=cmap(i % 10), label=fm.replace("malignant_", ""))
    ax.plot(x, tab.tumour_fraction_mean.values, marker="s", ms=3, color="black", ls="--", label="tumour fraction")
    ax.set_xticks(x); ax.set_xticklabels([s.replace("_", " ") for s in ORDER], rotation=45, ha="right")
    ax.set_ylabel("mean share of spot cells (q05)"); ax.set_xlabel("anatomical position (rostral → caudal)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=5)
    fig.savefig(FIGS / "phase6_brainstem_gradient_mp.png"); plt.close(fig)
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
