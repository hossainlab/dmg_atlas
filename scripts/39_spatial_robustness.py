#!/usr/bin/env python
"""Phase 9 science check / step 39 — spatial co-localisation: small-sample intervals and signature overlap.

(1) Hartung-Knapp-Sidik-Jonkman (HKSJ) intervals. The primary pooling (15b_spatial_tests_fix.py) averages Fisher z within
    patient and pools 12 patients by DerSimonian-Laird with normal-theory intervals, which are too narrow with few patients and
    high heterogeneity. Here the same patient-level z and variances are re-pooled; the DL estimate is asserted to reproduce
    results/tables/spatial_colocalisation.csv, and HKSJ adds se_HK = sqrt(sum w (z - mu)^2 / ((k - 1) sum w)),
    w = 1 / (v + tau2), with t(k - 1) intervals and P values; BH across the 42 state-microenvironment pairs per metric.
(2) Signature overlap. Deconvolution can create spurious co-localisation between classes with similar reference signatures.
    For every class pair, cosine similarity of the cell2location signatures (log1p of the per-class mean expression,
    centred per gene) is computed and related to the pooled correlation across the 42 pairs (Spearman).
Inputs : results/tables/spatial_colocalisation_per_section.csv, spatial_colocalisation.csv, spatial_patient_map.csv (15b),
         results/tables/c2l_signatures.csv (13-class reference signatures)
Outputs (results/tables/quality/): spatial_hksj.csv, spatial_signature_overlap.csv, spatial_robustness_inputs.json
Run: nice -n 10 python scripts/39_spatial_robustness.py [--force]
Log: logs/phase9_quality.log
"""
import sys, json, hashlib, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd  # noqa: E402
import scipy.stats as ss  # noqa: E402
from dmg_atlas import RESULTS, LOGS  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
log = U.StepLogger(LOGS / "phase9_quality.log", context="39_spatial_robustness")
T = RESULTS / "tables"; OUT = T / "quality"; OUT.mkdir(parents=True, exist_ok=True)
OUTS = [OUT / f for f in ("spatial_hksj.csv", "spatial_signature_overlap.csv", "spatial_robustness_inputs.json")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)
INPUTS = [T / "spatial_colocalisation_per_section.csv", T / "spatial_colocalisation.csv", T / "c2l_signatures.csv"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def patient_level(g):
    """Same within-patient averaging as 15b_spatial_tests_fix.two_level."""
    z = np.arctanh(np.clip(g.rho.values, -0.999, 0.999)); v = 1.0 / (g.n_spots.values - 3)
    d = pd.DataFrame(dict(patient=g.patient_id.values, z=z, v=v)); zp, vp = [], []
    for _, h in d.groupby("patient"):
        k = len(h); s2 = h.z.var(ddof=1) if k > 1 else 0.0
        zp.append(h.z.mean()); vp.append(h.v.sum() / k ** 2 + s2 / k)
    zp, vp = np.asarray(zp), np.asarray(vp)
    return zp, np.where(vp <= 0, np.maximum(vp.max(), 1e-6) * 1e-3 + 1e-9, vp)


def pool(z, v):
    w = 1 / v; mu_f = np.sum(w * z) / np.sum(w); Q = np.sum(w * (z - mu_f) ** 2); k = len(z)
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w); tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    ws = 1 / (v + tau2); mu = np.sum(ws * z) / np.sum(ws); se_dl = np.sqrt(1 / np.sum(ws))
    se_hk = np.sqrt(np.sum(ws * (z - mu) ** 2) / ((k - 1) * np.sum(ws))); tq = ss.t.ppf(0.975, k - 1)
    return dict(k=k, rho_dl=np.tanh(mu), se_dl=se_dl, p_dl=2 * ss.norm.sf(abs(mu / se_dl)), tau2=tau2,
                rho_hk_lo=np.tanh(mu - tq * se_hk), rho_hk_hi=np.tanh(mu + tq * se_hk), se_hk=se_hk,
                p_hk=2 * ss.t.sf(abs(mu / se_hk), k - 1) if se_hk > 0 else np.nan)


with log.step("HKSJ re-pooling"):
    per = pd.read_csv(INPUTS[0]); co = pd.read_csv(INPUTS[1])
    per = per[per.group == "DMG"]
    rows = []
    for (m, a, b), g in per.groupby(["metric", "class_a", "class_b"]):
        zp, vp = patient_level(g); r = pool(zp, vp)
        rows.append(dict(metric=m, class_a=a, class_b=b, **r))
    hk = pd.DataFrame(rows)
    ref = co[(co.subset == "all") & (co.group == "DMG")][["metric", "class_a", "class_b", "rho_re", "pair_type", "q_bh_state_vs_TME"]]
    hk = hk.merge(ref, on=["metric", "class_a", "class_b"], how="inner")
    assert np.allclose(hk.rho_dl, hk.rho_re, atol=1e-6), "re-pooling does not reproduce the published DL estimates"
    hk = hk[hk.pair_type == "state_vs_TME"].copy()
    hk["q_hk"] = np.nan
    for m, g in hk.groupby("metric"):
        hk.loc[g.index, "q_hk"] = ss.false_discovery_control(g.p_hk.values, method="bh")
    hk = hk.rename(columns={"q_bh_state_vs_TME": "q_dl"}).drop(columns=["rho_re"])
    hk.to_csv(OUTS[0], index=False)
    for m in ("proportion", "clr"):
        g = hk[hk.metric == m]
        log.write(f"{m}: pairs q_dl < 0.05: {int((g.q_dl < 0.05).sum())}; q_hk < 0.05: {int((g.q_hk < 0.05).sum())}\n"
                  + g.sort_values("p_hk").head(10)[["class_a", "class_b", "rho_dl", "rho_hk_lo", "rho_hk_hi", "q_dl", "q_hk"]].round(4).to_string(index=False))

with log.step("signature overlap"):
    sig = pd.read_csv(INPUTS[2], index_col=0)
    L = np.log1p(sig); L = L.sub(L.mean(axis=1), axis=0)
    nrm = np.sqrt((L ** 2).sum(axis=0))
    cos = (L.T @ L) / np.outer(nrm, nrm)
    cos = pd.DataFrame(cos.values, index=sig.columns, columns=sig.columns)
    p = hk[hk.metric == "proportion"].copy()
    p["signature_cosine"] = [cos.loc[a, b] for a, b in zip(p.class_a, p.class_b)]
    rho_s, p_s = ss.spearmanr(p.signature_cosine, p.rho_dl)
    p[["class_a", "class_b", "signature_cosine", "rho_dl", "q_dl", "q_hk"]].sort_values("signature_cosine", ascending=False).to_csv(OUTS[1], index=False)
    log.write(f"Spearman(signature cosine, pooled rho) over {len(p)} pairs: rho {rho_s:.3f}, P {p_s:.3g}\n"
              + p.sort_values("signature_cosine", ascending=False).head(10)[["class_a", "class_b", "signature_cosine", "rho_dl", "q_dl"]].round(3).to_string(index=False))
    json.dump(dict(script="scripts/39_spatial_robustness.py", spearman_cosine_vs_rho=dict(rho=rho_s, p=p_s),
                   inputs={str(q.relative_to(RESULTS.parent)): md5(q) for q in INPUTS}), open(OUTS[2], "w"), indent=1)
