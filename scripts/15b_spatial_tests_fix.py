#!/usr/bin/env python
"""Phase 6 / step 15b (audit remediation 2026-10-20; findings SP-01, SP-02, SP-07, SP-08, SP-11, SP-12, SP-13, SP-19).

Re-runs every Phase 6 spatial statistic on the deconvolution produced with the corrected metaprogram-family
reference (13 classes: 7 malignant families + 6 normal lineages), with three method changes:

  SP-01  co-localisation is computed on spot PROPORTIONS (q05 / sum of q05 over the 13 classes) as the primary
         metric, with a centred-log-ratio (CLR) variant; the raw q05 correlations are kept as a secondary table.
         Rationale: every class correlates with the spot total q05, which tracks UMI depth, so raw correlations
         are dominated by a shared size factor.
  SP-11  sections are pooled with a TWO-LEVEL random-effects model whose unit is the PATIENT: the 10 GSE280990
         brainstem sections are one autopsy, Collot BT01_192 core/edge sections are one patient; Fisher z is
         averaged within patient (within-patient sampling variance = sum v_i/k^2 + s^2_within/k) and pooled
         across patients by DerSimonian-Laird. A leave-brainstem-out pooled estimate is reported alongside.
  SP-08  BH correction is applied to BOTH brainstem metrics (absolute and proportion) and to the within-section
         random-effects tests; SP-12 adds BH across the Ren DMG-vs-GBM variables.

Sensitivity (SP-07): the pooled co-localisation is recomputed after dropping spots with < 500 UMIs.
Secondary (SP-09/SP-19): no change to the priors; the compositional and spatial-autocorrelation caveats are
reported in the notes.

Input : results/tables/spatial_abundance_long.csv (rebuilt by 14_cell2location.py collect from the 20,000-epoch run)
Outputs (results/tables/): spatial_colocalisation.csv, spatial_colocalisation_per_section.csv,
  spatial_colocalisation_raw.csv, spatial_colocalisation_sensitivity.csv, spatial_niches.csv,
  spatial_niches_by_section.csv, brainstem_gradient.csv, brainstem_section_means.csv,
  ren_dmg_vs_gbm_composition.csv, ren_section_composition.csv, spatial_patient_map.csv
Figures: results/figures/phase6_colocalisation_matrix.png (proportion-based), phase6_spatial_abundance_maps.png,
  phase6_brainstem_gradient.png  (written by 15c_spatial_figures_fix.py)
"""
from __future__ import annotations
import sys, json, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402
import numpy as np, pandas as pd, scipy.stats as ss, anndata as ad, scanpy as sc  # noqa: E402

LOG = StepLogger(ROOT / "logs/fix_spatial.log", context="15b_spatial_tests_fix")
TABLES = ROOT / "results/tables"
FAMILIES = ["OPC_like", "AC_like", "OC_like", "MES_like", "cycling", "hypoxia", "novel_MP9"]
STATES = [f"malignant_{f}" for f in FAMILIES]
TME = ["myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular"]
CLASSES = STATES + TME
BRAINSTEM_ORDER = ["Midbrain", "Pons_A", "Pons_B", "Pons_C", "Pons_D", "Medulla_A", "Medulla_B", "Medulla_C",
                   "Spinal_Cord_A", "Spinal_Cord_B"]
GBM = ["GBM1", "GBM2", "GBM3", "GBM4", "GBM5_1"]; PERITUMOUR = ["GBM5_2"]
SEED = 0
MIN_UMI = 500  # SP-07 sensitivity threshold


def patient_of(section_id: str, cohort: str) -> str:
    """SP-11: unit of replication. The 10 GSE280990 brainstem sections are serial sections of ONE DIPG autopsy;
    Collot core/edge sections of the same tumour block share a patient; every other section is its own patient."""
    if cohort == "brainstem":
        return "brainstem_autopsy_P1"
    if cohort == "collot_st":
        return "collot_" + section_id.replace("_core_HE", "").replace("_edge_HE", "").replace("_HE", "")
    return f"{cohort}_{section_id}"


def load():
    long = pd.read_csv(TABLES / "spatial_abundance_long.csv")
    q = long.pivot_table(index="spot", columns="class", values="q05")[CLASSES]
    meta = (long.drop_duplicates("spot").set_index("spot")
            [["cohort", "section_id", "region_or_diagnosis", "array_row", "array_col"]].loc[q.index])
    meta["group"] = np.where(meta.section_id.isin(GBM), "GBM",
                             np.where(meta.section_id.isin(PERITUMOUR), "peritumour", "DMG"))
    meta["patient_id"] = [patient_of(s, c) for s, c in zip(meta.section_id, meta.cohort)]
    depth = pd.read_csv(TABLES / "spatial_spot_depth.csv", index_col=0)["n_counts"] \
        if (TABLES / "spatial_spot_depth.csv").exists() else None
    if depth is not None:
        meta["n_counts"] = depth.reindex(meta.index).values
    P = q.div(q.sum(1), axis=0)                                  # proportions (primary)
    gm = np.exp(np.log(q.clip(lower=1e-6)).mean(axis=1))
    C = np.log(q.clip(lower=1e-6)).sub(np.log(gm), axis=0)       # CLR
    return q, P, C, meta


def dl(z, v):
    w = 1 / v; zf = np.sum(w * z) / np.sum(w)
    Q = np.sum(w * (z - zf) ** 2); df = len(z) - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    ws = 1 / (v + tau2); zr = np.sum(ws * z) / np.sum(ws); se = np.sqrt(1 / np.sum(ws))
    i2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return zr, se, tau2, Q, i2


def two_level(per_section: pd.DataFrame):
    """per_section: columns rho, n_spots, patient_id. Fisher-z averaged within patient, DL across patients."""
    z = np.arctanh(np.clip(per_section.rho.values, -0.999, 0.999))
    v = 1.0 / (per_section.n_spots.values - 3)
    d = pd.DataFrame(dict(patient=per_section.patient_id.values, z=z, v=v))
    zp, vp = [], []
    for p, g in d.groupby("patient"):
        k = len(g)
        s2 = g.z.var(ddof=1) if k > 1 else 0.0
        zp.append(g.z.mean()); vp.append(g.v.sum() / k ** 2 + s2 / k)
    zp = np.asarray(zp); vp = np.asarray(vp)
    vp = np.where(vp <= 0, np.maximum(vp.max(), 1e-6) * 1e-3 + 1e-9, vp)
    zr, se, tau2, Q, i2 = dl(zp, vp)
    return dict(n_patients=len(zp), rho_re=np.tanh(zr), rho_re_lo=np.tanh(zr - 1.96 * se),
                rho_re_hi=np.tanh(zr + 1.96 * se), z_re=zr, se_z=se,
                p_re=2 * ss.norm.sf(abs(zr / se)) if se > 0 else np.nan, tau2=tau2, I2=i2)


def per_section_rhos(M, meta, metric):
    rows = []
    for sec, idx in meta.groupby("section_id").groups.items():
        sub = M.loc[idx]
        if len(sub) < 10:
            continue
        rho = sub.corr(method="spearman")
        m0 = meta.loc[idx[0]]
        for a in STATES:
            for b in TME:
                rows.append(dict(metric=metric, section_id=sec, cohort=m0.cohort, group=m0.group,
                                 patient_id=m0.patient_id, class_a=a, class_b=b, rho=rho.loc[a, b], n_spots=len(sub)))
        for i, a in enumerate(STATES):
            for b in STATES[i + 1:]:
                rows.append(dict(metric=metric, section_id=sec, cohort=m0.cohort, group=m0.group,
                                 patient_id=m0.patient_id, class_a=a, class_b=b, rho=rho.loc[a, b],
                                 n_spots=len(sub)))
    return pd.DataFrame(rows)


def pool(per, metric, subset_label="all", sections=None):
    d = per if sections is None else per[per.section_id.isin(sections)]
    out = []
    for grp in ("DMG", "GBM"):
        dd = d[d.group == grp]
        if not len(dd):
            continue
        for (a, b), g in dd.groupby(["class_a", "class_b"]):
            r = two_level(g)
            r.update(metric=metric, subset=subset_label, group=grp, class_a=a, class_b=b, n_sections=len(g),
                     rho_mean=g.rho.mean(), n_sections_positive=int((g.rho > 0).sum()),
                     pair_type="state_vs_TME" if b in TME else "state_vs_state")
            out.append(r)
    out = pd.DataFrame(out)
    # BH over the DMG state-vs-TME family. p_re is NaN when the two-level model has no
    # usable between-patient variance for that pair (e.g. a single contributing patient);
    # those pairs are left with q = NaN rather than being passed to the FDR routine.
    m = (out.group == "DMG") & (out.pair_type == "state_vs_TME") & out.p_re.notna()
    if m.any():
        out.loc[m, "q_bh_state_vs_TME"] = ss.false_discovery_control(
            out.loc[m, "p_re"].clip(lower=0.0, upper=1.0).values, method="bh")
    n_nan = int(((out.group == "DMG") & (out.pair_type == "state_vs_TME") & out.p_re.isna()).sum())
    if n_nan:
        LOG.write(f"{metric}/{subset_label}: {n_nan} DMG state-vs-TME pairs had p_re = NaN and were excluded from BH")
    return out


def colocalisation(q, P, C, meta):
    per_p = per_section_rhos(P, meta, "proportion")
    per_r = per_section_rhos(q, meta, "raw_q05")
    per_c = per_section_rhos(C, meta, "clr")
    per = pd.concat([per_p, per_r, per_c], ignore_index=True)
    per.to_csv(TABLES / "spatial_colocalisation_per_section.csv", index=False)
    main = pd.concat([pool(per_p, "proportion"), pool(per_c, "clr")], ignore_index=True)
    main.to_csv(TABLES / "spatial_colocalisation.csv", index=False)
    pool(per_r, "raw_q05").to_csv(TABLES / "spatial_colocalisation_raw.csv", index=False)
    # SP-11 leave-one-cohort-out and leave-brainstem-out, proportion metric, DMG only
    sens = []
    dmg_sections = sorted(meta.loc[meta.group == "DMG", "section_id"].unique())
    for lab, secs in [("leave_brainstem_out", [s for s in dmg_sections if s not in BRAINSTEM_ORDER]),
                      ("brainstem_only", [s for s in dmg_sections if s in BRAINSTEM_ORDER])]:
        t = pool(per_p[per_p.section_id.isin(secs)], "proportion", subset_label=lab)
        sens.append(t)
    for coh in sorted(meta.loc[meta.group == "DMG", "cohort"].unique()):
        secs = [s for s in dmg_sections if s not in set(meta.loc[meta.cohort == coh, "section_id"])]
        sens.append(pool(per_p[per_p.section_id.isin(secs)], "proportion", subset_label=f"leave_out_{coh}"))
    return per, main, pd.concat(sens, ignore_index=True)


def niches(P, meta):
    dmg = meta.index[meta.group == "DMG"]
    X = P.loc[dmg].fillna(0)
    a = ad.AnnData(X=X.values.astype(np.float32), obs=meta.loc[dmg].copy(), var=pd.DataFrame(index=CLASSES))
    a.obsm["X_prop"] = X.values.astype(np.float32)
    sc.pp.neighbors(a, use_rep="X_prop", n_neighbors=15, random_state=SEED)
    sc.tl.leiden(a, resolution=0.5, flavor="igraph", n_iterations=2, directed=False, random_state=SEED, key_added="niche")
    lab = a.obs["niche"].astype(int)
    sizes = lab.value_counts()
    comp = X.groupby(lab.values).mean().loc[sizes.index]
    comp.insert(0, "n_spots", sizes.values)
    comp.insert(1, "frac_spots", (sizes / sizes.sum()).values)
    comp.insert(2, "dominant_class", X.groupby(lab.values).mean().loc[sizes.index].idxmax(1).values)
    comp.insert(3, "n_sections_present", [lab.index[lab == k].map(meta.section_id).nunique() for k in sizes.index])
    comp.insert(4, "n_patients_present", [lab.index[lab == k].map(meta.patient_id).nunique() for k in sizes.index])
    comp.index.name = "niche"
    comp.to_csv(TABLES / "spatial_niches.csv")
    pd.crosstab(meta.loc[dmg, "section_id"], lab.values, normalize="index").to_csv(TABLES / "spatial_niches_by_section.csv")
    labels = meta.loc[dmg, ["cohort", "section_id", "patient_id"]].assign(niche=lab.values)
    labels.to_csv(ROOT / "data/processed/spatial/spatial_niche_labels.csv.gz")
    return comp


def brainstem(q, P, meta):
    bs = meta.index[meta.cohort == "brainstem"]
    qb, Pb, mb = q.loc[bs], P.loc[bs], meta.loc[bs]
    rank = {s: i + 1 for i, s in enumerate(BRAINSTEM_ORDER)}
    means_abs = qb.groupby(mb.section_id).mean().loc[BRAINSTEM_ORDER]
    means_prop = Pb.groupby(mb.section_id).mean().loc[BRAINSTEM_ORDER]
    tf = qb[STATES].sum(1) / qb.sum(1)
    sec_tab = means_prop.add_prefix("prop_").join(means_abs.add_prefix("abs_"))
    sec_tab.insert(0, "anatomical_rank", [rank[s] for s in sec_tab.index])
    sec_tab.insert(1, "n_spots", mb.section_id.value_counts().loc[BRAINSTEM_ORDER].values)
    sec_tab.insert(2, "tumour_fraction_mean", tf.groupby(mb.section_id).mean().loc[BRAINSTEM_ORDER].values)
    sec_tab.insert(3, "tumour_fraction_median", tf.groupby(mb.section_id).median().loc[BRAINSTEM_ORDER].values)
    sec_tab.to_csv(TABLES / "brainstem_section_means.csv")
    x = np.array([rank[s] for s in BRAINSTEM_ORDER])
    rows = []
    for c in CLASSES + ["tumour_fraction"]:
        ya = sec_tab["tumour_fraction_mean"].values if c == "tumour_fraction" else means_abs[c].values
        yp = ya if c == "tumour_fraction" else means_prop[c].values
        ra, pa = ss.spearmanr(x, ya); rp, pp_ = ss.spearmanr(x, yp)
        rows.append(dict(test="trend_along_axis", cls=c, n_sections=10, n_patients=1, rho_abs=ra, p_abs=pa,
                         rho_prop=rp, p_prop=pp_, rostral_mean=np.mean(yp[:5]), caudal_mean=np.mean(yp[5:])))
    share = qb[STATES].div(qb[STATES].sum(1), axis=0)
    for st in STATES:
        rr = []
        for s in BRAINSTEM_ORDER:
            m = (mb.section_id == s).values
            r, p = ss.spearmanr(tf[m], share.loc[m, st])
            rr.append(dict(test="within_section_tumour_fraction_vs_state_share", cls=st, section_id=s,
                           n_spots=int(m.sum()), rho=r, p=p))
        rr = pd.DataFrame(rr)
        z = np.arctanh(np.clip(rr.rho.values, -0.999, 0.999)); v = 1 / (rr.n_spots.values - 3)
        zr, se, tau2, Q, i2 = dl(z, v)
        rows.append(dict(test="within_section_tumour_fraction_vs_state_share_RE", cls=st, n_sections=10, n_patients=1,
                         rho_re=np.tanh(zr), rho_re_lo=np.tanh(zr - 1.96 * se), rho_re_hi=np.tanh(zr + 1.96 * se),
                         p_re=2 * ss.norm.sf(abs(zr / se)), tau2=tau2, n_sections_positive=int((rr.rho > 0).sum())))
        rows.extend(rr.to_dict("records"))
    for c in TME:
        rr = []
        for s in BRAINSTEM_ORDER:
            m = (mb.section_id == s).values
            r, p = ss.spearmanr(tf[m], qb.loc[m, c])
            rr.append(dict(test="within_section_tumour_fraction_vs_TME_abundance", cls=c, section_id=s,
                           n_spots=int(m.sum()), rho=r, p=p))
        rr = pd.DataFrame(rr)
        z = np.arctanh(np.clip(rr.rho.values, -0.999, 0.999)); v = 1 / (rr.n_spots.values - 3)
        zr, se, tau2, Q, i2 = dl(z, v)
        rows.append(dict(test="within_section_tumour_fraction_vs_TME_abundance_RE", cls=c, n_sections=10, n_patients=1,
                         rho_re=np.tanh(zr), rho_re_lo=np.tanh(zr - 1.96 * se), rho_re_hi=np.tanh(zr + 1.96 * se),
                         p_re=2 * ss.norm.sf(abs(zr / se)), tau2=tau2, n_sections_positive=int((rr.rho > 0).sum())))
        rows.extend(rr.to_dict("records"))
    grad = pd.DataFrame(rows)
    # SP-08: BH within each family of tests, both metrics
    m = grad.test == "trend_along_axis"
    grad.loc[m, "q_bh_abs"] = ss.false_discovery_control(grad.loc[m, "p_abs"].values, method="bh")
    grad.loc[m, "q_bh_prop"] = ss.false_discovery_control(grad.loc[m, "p_prop"].values, method="bh")
    for t in ("within_section_tumour_fraction_vs_state_share_RE", "within_section_tumour_fraction_vs_TME_abundance_RE"):
        m = grad.test == t
        grad.loc[m, "q_bh_re"] = ss.false_discovery_control(grad.loc[m, "p_re"].values, method="bh")
    grad.to_csv(TABLES / "brainstem_gradient.csv", index=False)
    return grad, sec_tab, means_prop


def ren_outgroup(q, meta):
    ren = meta.index[meta.cohort == "ren"]
    qr, mr = q.loc[ren], meta.loc[ren]
    rows = []
    for s, idx in mr.groupby("section_id").groups.items():
        sub = qr.loc[idx]; tot_mal = sub[STATES].sum().sum()
        d = dict(section_id=s, group=mr.loc[idx[0], "group"], diagnosis=mr.loc[idx[0], "region_or_diagnosis"],
                 n_spots=len(sub), tumour_fraction=tot_mal / sub.sum().sum())
        for st in STATES:
            d[f"share_{st}"] = sub[st].sum() / tot_mal
        for c in TME:
            d[f"frac_{c}"] = sub[c].sum() / sub.sum().sum()
        rows.append(d)
    sec = pd.DataFrame(rows).sort_values(["group", "section_id"])
    sec.to_csv(TABLES / "ren_section_composition.csv", index=False)
    out = []
    for col in [f"share_{s}" for s in STATES] + ["tumour_fraction"] + [f"frac_{c}" for c in TME]:
        a = sec.loc[sec.group == "DMG", col].values; b = sec.loc[sec.group == "GBM", col].values
        u, p = ss.mannwhitneyu(a, b, alternative="two-sided", method="exact")
        out.append(dict(variable=col, dmg_median=np.median(a), dmg_min=a.min(), dmg_max=a.max(),
                        gbm_median=np.median(b), gbm_min=b.min(), gbm_max=b.max(), n_dmg=len(a), n_gbm=len(b),
                        U=u, p_mwu=p, peritumour_GBM5_2=sec.loc[sec.section_id == "GBM5_2", col].iloc[0]))
    out = pd.DataFrame(out)
    out["q_bh"] = ss.false_discovery_control(out.p_mwu.values, method="bh")  # SP-12
    out.to_csv(TABLES / "ren_dmg_vs_gbm_composition.csv", index=False)
    return out, sec


def main():
    for f in ("spatial_colocalisation.csv", "spatial_colocalisation_per_section.csv", "spatial_niches.csv",
              "spatial_niches_by_section.csv", "brainstem_gradient.csv", "brainstem_section_means.csv",
              "ren_dmg_vs_gbm_composition.csv", "ren_section_composition.csv"):
        p = TABLES / f; pre = TABLES / f.replace(".csv", ".pre_audit.csv")
        if p.exists() and not pre.exists():
            shutil.copy(p, pre)
    with LOG.step("load"):
        q, P, C, meta = load()
        meta[["cohort", "section_id", "patient_id", "group"]].drop_duplicates("section_id") \
            .to_csv(TABLES / "spatial_patient_map.csv", index=False)
        LOG.write(f"spots {len(q)}, sections {meta.section_id.nunique()}, patients {meta.patient_id.nunique()}, "
                  f"groups {meta.group.value_counts().to_dict()}; DMG patients "
                  f"{meta.loc[meta.group=='DMG','patient_id'].nunique()}")
    with LOG.step("colocalisation (proportion primary, CLR, raw secondary)"):
        per, coloc, sens = colocalisation(q, P, C, meta)
        if "n_counts" in meta:
            keep = meta.index[meta.n_counts >= MIN_UMI]
            pq = per_section_rhos(P.loc[keep], meta.loc[keep], "proportion_min500umi")
            sens = pd.concat([sens, pool(pq, "proportion_min500umi", subset_label="min500umi")], ignore_index=True)
            LOG.write(f"SP-07 sensitivity: {len(meta)-len(keep)} spots (<{MIN_UMI} UMI) dropped")
        sens.to_csv(TABLES / "spatial_colocalisation_sensitivity.csv", index=False)
        k = coloc[(coloc.group == "DMG") & (coloc.metric == "proportion") & (coloc.pair_type == "state_vs_TME")]
        LOG.write("DMG proportion-based state x TME (sorted):\n" + k.sort_values("rho_re", ascending=False)
                  [["class_a", "class_b", "n_sections", "n_patients", "rho_re", "rho_re_lo", "rho_re_hi",
                    "n_sections_positive", "p_re", "q_bh_state_vs_TME"]].round(4).to_string(index=False))
    with LOG.step("niches"):
        comp = niches(P, meta); LOG.write("niches:\n" + comp.round(3).to_string())
    with LOG.step("brainstem gradient"):
        grad, sec_tab, means_prop = brainstem(q, P, meta)
        LOG.write("trend along axis:\n" + grad[grad.test == "trend_along_axis"]
                  [["cls", "rho_abs", "p_abs", "q_bh_abs", "rho_prop", "p_prop", "q_bh_prop"]].round(4).to_string(index=False))
        LOG.write("within-section RE:\n" + grad[grad.test.str.endswith("_RE")]
                  [["test", "cls", "rho_re", "rho_re_lo", "rho_re_hi", "p_re", "q_bh_re", "n_sections_positive"]]
                  .round(4).to_string(index=False))
    with LOG.step("Ren DMG vs GBM"):
        out, sec = ren_outgroup(q, meta)
        LOG.write("DMG vs GBM:\n" + out.round(4).to_string(index=False))
    LOG.write(f"done | peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
