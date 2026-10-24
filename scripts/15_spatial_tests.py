#!/usr/bin/env python
"""Phase 6 / step 15 — spatial tests on the cell2location abundances (section = unit of replication).

Input : results/tables/spatial_abundance_long.csv (q05 and posterior-mean cell abundance per spot x class).
        All statistics use the q05 abundances (cell2location's recommended conservative estimate).
Outputs (results/tables/): spatial_colocalisation.csv (per-pair random-effects summary + per-section rhos in
        spatial_colocalisation_per_section.csv), spatial_niches.csv (+ spatial_niches_by_section.csv, per-spot labels in
        data/processed/spatial/spatial_niche_labels.csv.gz), brainstem_gradient.csv (+ brainstem_section_means.csv),
        ren_dmg_vs_gbm_composition.csv (+ per-section table ren_section_composition.csv).
Figures (results/figures/): phase6_spatial_abundance_maps.png, phase6_colocalisation_matrix.png, phase6_brainstem_gradient.png.

Methods
  (a) co-localisation: per section, Spearman rho between the q05 abundance of each malignant state and each TME class
      across spots. Sections pooled with a DerSimonian-Laird random-effects model on Fisher z (within-section variance
      1/(n_spots-3)); CI back-transformed. A section-level bootstrap (10,000 resamples of sections, simple mean of rho)
      is reported alongside. Primary set = 21 DMG sections; GBM (Ren GBM1-GBM5_1) summarised separately.
  (b) niches: DMG spots jointly; per-spot abundance vector scaled to proportions (q05 / sum of q05); kNN graph
      (k = 15, Euclidean on proportions, 11 dims) and Leiden (igraph flavour, resolution 0.5, seed 0).
  (c) brainstem gradient (GSE280990, one DIPG autopsy): sections ordered midbrain -> pons A-D -> medulla A-C -> spinal
      cord A-B (rank 1..10). Per class, Spearman rho between the anatomical rank and the section mean of (i) absolute
      q05 abundance and (ii) the proportion of the spot's cells. Within section: Spearman between the spot tumour
      fraction (sum malignant / total) and each state's share of the malignant compartment.
  (d) out-group: Ren DMG1-5 vs GBM1-GBM5_1 (GBM5_2 peritumour shown separately): per-section malignant-state
      composition (state q05 / sum malignant q05 over all spots); Mann-Whitney U (two-sided, exact) per state.
"""
from __future__ import annotations

import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.stats as ss  # noqa: E402
import anndata as ad  # noqa: E402
import scanpy as sc  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase6_spatial.log", context="15_spatial_tests")
TABLES = ROOT / "results/tables"; FIGS = ROOT / "results/figures"; FIGS.mkdir(exist_ok=True, parents=True)
STATES = ["malignant_OPC_like", "malignant_AC_like", "malignant_OC_like", "malignant_MES_like", "malignant_cycling"]
TME = ["myeloid", "T_NK", "oligodendrocyte", "astrocyte", "neuron", "vascular"]
CLASSES = STATES + TME
SHORT = {"malignant_OPC_like": "OPC-like", "malignant_AC_like": "AC-like", "malignant_OC_like": "OC-like",
         "malignant_MES_like": "MES-like", "malignant_cycling": "cycling", "myeloid": "myeloid", "T_NK": "T/NK",
         "oligodendrocyte": "oligodendrocyte", "astrocyte": "astrocyte", "neuron": "neuron", "vascular": "vascular"}
COLOURS = {"malignant_OPC_like": "#b2182b", "malignant_AC_like": "#ef8a62", "malignant_OC_like": "#fddbc7",
           "malignant_MES_like": "#7f0000", "malignant_cycling": "#f4a582", "myeloid": "#2166ac", "T_NK": "#c51b7d",
           "oligodendrocyte": "#1b7837", "astrocyte": "#e08214", "neuron": "#762a83", "vascular": "#35978f"}
BRAINSTEM_ORDER = ["Midbrain", "Pons_A", "Pons_B", "Pons_C", "Pons_D", "Medulla_A", "Medulla_B", "Medulla_C", "Spinal_Cord_A", "Spinal_Cord_B"]
GBM = ["GBM1", "GBM2", "GBM3", "GBM4", "GBM5_1"]; PERITUMOUR = ["GBM5_2"]
SEED = 0


def load_wide():
    long = pd.read_csv(TABLES / "spatial_abundance_long.csv")
    q = long.pivot_table(index="spot", columns="class", values="q05")[CLASSES]
    meta = long.drop_duplicates("spot").set_index("spot")[["cohort", "section_id", "region_or_diagnosis", "array_row", "array_col"]].loc[q.index]
    meta["group"] = np.where(meta.section_id.isin(GBM), "GBM", np.where(meta.section_id.isin(PERITUMOUR), "peritumour", "DMG"))
    return q, meta


def dl_random_effects(z: np.ndarray, v: np.ndarray):
    """DerSimonian-Laird pooled estimate on effect sizes z with variances v."""
    w = 1 / v; zf = np.sum(w * z) / np.sum(w)
    Q = np.sum(w * (z - zf) ** 2); df = len(z) - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    ws = 1 / (v + tau2); zr = np.sum(ws * z) / np.sum(ws); se = np.sqrt(1 / np.sum(ws))
    i2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return zr, se, tau2, Q, i2


def colocalisation(q, meta):
    rng = np.random.default_rng(SEED)
    rows = []
    for sec, idx in meta.groupby("section_id").groups.items():
        sub = q.loc[idx]
        if len(sub) < 10:
            continue
        rho = sub.corr(method="spearman")
        for a in CLASSES:
            for b in CLASSES:
                if a < b:
                    rows.append(dict(section_id=sec, cohort=meta.loc[idx[0], "cohort"], group=meta.loc[idx[0], "group"],
                                     class_a=a, class_b=b, rho=rho.loc[a, b], n_spots=len(sub)))
    per = pd.DataFrame(rows)
    per.to_csv(TABLES / "spatial_colocalisation_per_section.csv", index=False)
    out = []
    for grp in ("DMG", "GBM"):
        for (a, b), d in per[per.group == grp].groupby(["class_a", "class_b"]):
            r = np.clip(d.rho.values, -0.999, 0.999); z = np.arctanh(r); v = 1 / (d.n_spots.values - 3)
            zr, se, tau2, Q, i2 = dl_random_effects(z, v)
            boots = np.array([np.mean(rng.choice(r, size=len(r), replace=True)) for _ in range(10000)]) if len(r) > 1 else np.array([r.mean()])
            out.append(dict(group=grp, class_a=a, class_b=b, n_sections=len(r),
                            rho_re=np.tanh(zr), rho_re_lo=np.tanh(zr - 1.96 * se), rho_re_hi=np.tanh(zr + 1.96 * se),
                            z_re=zr, se_z=se, p_re=2 * ss.norm.sf(abs(zr / se)) if se > 0 else np.nan, tau2=tau2, I2=i2,
                            rho_mean=r.mean(), rho_boot_lo=np.quantile(boots, 0.025), rho_boot_hi=np.quantile(boots, 0.975),
                            n_sections_positive=int((r > 0).sum()),
                            pair_type="state_vs_TME" if (a in STATES) != (b in STATES) else ("state_vs_state" if a in STATES else "TME_vs_TME")))
    out = pd.DataFrame(out)
    # BH over the 30 state-vs-TME pairs within DMG
    m = (out.group == "DMG") & (out.pair_type == "state_vs_TME")
    out.loc[m, "q_bh_state_vs_TME"] = ss.false_discovery_control(out.loc[m, "p_re"].values, method="bh")
    out.to_csv(TABLES / "spatial_colocalisation.csv", index=False)
    return out, per


def niches(q, meta):
    dmg = meta.index[meta.group == "DMG"]
    P = q.loc[dmg].div(q.loc[dmg].sum(1), axis=0).fillna(0)
    a = ad.AnnData(X=P.values.astype(np.float32), obs=meta.loc[dmg].copy(), var=pd.DataFrame(index=CLASSES))
    a.obsm["X_prop"] = P.values.astype(np.float32)
    sc.pp.neighbors(a, use_rep="X_prop", n_neighbors=15, random_state=SEED)
    sc.tl.leiden(a, resolution=0.5, flavor="igraph", n_iterations=2, directed=False, random_state=SEED, key_added="niche")
    lab = a.obs["niche"].astype(int)
    order = P.groupby(lab.values).mean()
    # name niches by dominant class, ordered by size
    sizes = lab.value_counts()
    comp = order.loc[sizes.index]
    comp.insert(0, "n_spots", sizes.values)
    comp.insert(1, "frac_spots", (sizes / sizes.sum()).values)
    comp.insert(2, "dominant_class", order.loc[sizes.index].idxmax(1).values)
    comp.insert(3, "n_sections_present", [lab.index[lab == k].map(meta.section_id).nunique() for k in sizes.index])
    comp.index.name = "niche"
    comp.to_csv(TABLES / "spatial_niches.csv")
    bysec = pd.crosstab(meta.loc[dmg, "section_id"], lab.values, normalize="index")
    bysec.to_csv(TABLES / "spatial_niches_by_section.csv")
    labels = meta.loc[dmg, ["cohort", "section_id"]].assign(niche=lab.values)
    (ROOT / "data/processed/spatial").mkdir(parents=True, exist_ok=True)
    labels.to_csv(ROOT / "data/processed/spatial/spatial_niche_labels.csv.gz")
    return comp, bysec, labels


def brainstem(q, meta):
    bs = meta.index[meta.cohort == "brainstem"]
    qb = q.loc[bs]; mb = meta.loc[bs]
    P = qb.div(qb.sum(1), axis=0)
    rank = {s: i + 1 for i, s in enumerate(BRAINSTEM_ORDER)}
    means_abs = qb.groupby(mb.section_id).mean().loc[BRAINSTEM_ORDER]
    means_prop = P.groupby(mb.section_id).mean().loc[BRAINSTEM_ORDER]
    tf = (qb[STATES].sum(1) / qb.sum(1))
    sec_tab = means_prop.add_prefix("prop_").join(means_abs.add_prefix("abs_"))
    sec_tab.insert(0, "anatomical_rank", [rank[s] for s in sec_tab.index])
    sec_tab.insert(1, "n_spots", mb.section_id.value_counts().loc[BRAINSTEM_ORDER].values)
    sec_tab.insert(2, "tumour_fraction_mean", tf.groupby(mb.section_id).mean().loc[BRAINSTEM_ORDER].values)
    sec_tab.insert(3, "tumour_fraction_median", tf.groupby(mb.section_id).median().loc[BRAINSTEM_ORDER].values)
    sec_tab.to_csv(TABLES / "brainstem_section_means.csv")
    rows = []
    x = np.array([rank[s] for s in BRAINSTEM_ORDER])
    for c in CLASSES + ["tumour_fraction"]:
        if c == "tumour_fraction":
            ya = sec_tab["tumour_fraction_mean"].values; yp = ya
        else:
            ya = means_abs[c].values; yp = means_prop[c].values
        ra, pa = ss.spearmanr(x, ya); rp, pp_ = ss.spearmanr(x, yp)
        rows.append(dict(test="trend_along_axis", cls=c, n_sections=10, rho_abs=ra, p_abs=pa, rho_prop=rp, p_prop=pp_,
                         rostral_mean=np.mean(yp[:5]), caudal_mean=np.mean(yp[5:])))
    # within-section tumour fraction vs state share of the malignant compartment
    share = qb[STATES].div(qb[STATES].sum(1), axis=0)
    for st in STATES:
        rhos = []
        for s in BRAINSTEM_ORDER:
            m = mb.section_id == s
            r, p = ss.spearmanr(tf[m], share.loc[m, st])
            rhos.append(dict(test="within_section_tumour_fraction_vs_state_share", cls=st, section_id=s, n_spots=int(m.sum()), rho=r, p=p))
        rr = pd.DataFrame(rhos)
        z = np.arctanh(np.clip(rr.rho.values, -0.999, 0.999)); v = 1 / (rr.n_spots.values - 3)
        zr, se, tau2, Q, i2 = dl_random_effects(z, v)
        rows.append(dict(test="within_section_tumour_fraction_vs_state_share_RE", cls=st, n_sections=10, rho_re=np.tanh(zr),
                         rho_re_lo=np.tanh(zr - 1.96 * se), rho_re_hi=np.tanh(zr + 1.96 * se), p_re=2 * ss.norm.sf(abs(zr / se)),
                         tau2=tau2, n_sections_positive=int((rr.rho > 0).sum())))
        rows.extend(rhos)
    # also TME classes vs tumour fraction (absolute abundance), within section
    for c in TME:
        rhos = []
        for s in BRAINSTEM_ORDER:
            m = mb.section_id == s
            r, p = ss.spearmanr(tf[m], qb.loc[m, c])
            rhos.append(dict(test="within_section_tumour_fraction_vs_TME_abundance", cls=c, section_id=s, n_spots=int(m.sum()), rho=r, p=p))
        rr = pd.DataFrame(rhos)
        z = np.arctanh(np.clip(rr.rho.values, -0.999, 0.999)); v = 1 / (rr.n_spots.values - 3)
        zr, se, tau2, Q, i2 = dl_random_effects(z, v)
        rows.append(dict(test="within_section_tumour_fraction_vs_TME_abundance_RE", cls=c, n_sections=10, rho_re=np.tanh(zr),
                         rho_re_lo=np.tanh(zr - 1.96 * se), rho_re_hi=np.tanh(zr + 1.96 * se), p_re=2 * ss.norm.sf(abs(zr / se)),
                         tau2=tau2, n_sections_positive=int((rr.rho > 0).sum())))
        rows.extend(rhos)
    grad = pd.DataFrame(rows)
    grad.to_csv(TABLES / "brainstem_gradient.csv", index=False)
    return grad, sec_tab, means_prop


def ren_outgroup(q, meta):
    ren = meta.index[meta.cohort == "ren"]
    qr = q.loc[ren]; mr = meta.loc[ren]
    rows = []
    for s, idx in mr.groupby("section_id").groups.items():
        sub = qr.loc[idx]; tot_mal = sub[STATES].sum().sum()
        d = dict(section_id=s, group=mr.loc[idx[0], "group"], diagnosis=mr.loc[idx[0], "region_or_diagnosis"], n_spots=len(sub),
                 tumour_fraction=tot_mal / sub.sum().sum())
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
        out.append(dict(variable=col, dmg_median=np.median(a), dmg_min=a.min(), dmg_max=a.max(), gbm_median=np.median(b),
                        gbm_min=b.min(), gbm_max=b.max(), n_dmg=len(a), n_gbm=len(b), U=u, p_mwu=p,
                        peritumour_GBM5_2=sec.loc[sec.section_id == "GBM5_2", col].iloc[0]))
    out = pd.DataFrame(out)
    out.to_csv(TABLES / "ren_dmg_vs_gbm_composition.csv", index=False)
    return out, sec


def fig_maps(q, meta):
    apply_figure_style()
    picks = {"huchede": "pDMG_Sample-2", "ren": "DMG1", "collot_st": "BT01_192_core_HE", "brainstem": "Pons_A"}
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.4))
    dom_all = q.idxmax(1)
    for ax, (cohort, sec) in zip(axes, picks.items()):
        idx = meta.index[meta.section_id == sec]
        m = meta.loc[idx]; dom = dom_all.loc[idx]
        for c in CLASSES:
            k = dom == c
            if k.any():
                ax.scatter(m.loc[k, "array_col"], -m.loc[k, "array_row"], s=2.2, c=COLOURS[c], linewidths=0, label=SHORT[c], rasterized=True)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp_ in ax.spines.values():
            sp_.set_visible(False)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=4, color=COLOURS[c], label=SHORT[c]) for c in CLASSES]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=1, handletextpad=0.3)
    fig.savefig(FIGS / "phase6_spatial_abundance_maps.png")
    plt.close(fig)
    return picks


def fig_coloc(coloc):
    apply_figure_style()
    d = coloc[(coloc.group == "DMG")]
    M = pd.DataFrame(index=STATES, columns=TME, dtype=float); N = M.copy()
    for _, r in d.iterrows():
        a, b = r.class_a, r.class_b
        if a in STATES and b in TME:
            M.loc[a, b] = r.rho_re; N.loc[a, b] = r.n_sections_positive
        elif b in STATES and a in TME:
            M.loc[b, a] = r.rho_re; N.loc[b, a] = r.n_sections_positive
    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    im = ax.imshow(M.values.astype(float), cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(len(TME))); ax.set_xticklabels([SHORT[c] for c in TME], rotation=45, ha="right")
    ax.set_yticks(range(len(STATES))); ax.set_yticklabels([SHORT[c] for c in STATES])
    for i in range(len(STATES)):
        for j in range(len(TME)):
            ax.text(j, i, f"{M.iloc[i, j]:.2f}\n{int(N.iloc[i, j])}/{int(d.n_sections.max())}", ha="center", va="center", fontsize=5,
                    color="white" if abs(M.iloc[i, j]) > 0.35 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03); cb.set_label("random-effects Spearman ρ")
    fig.savefig(FIGS / "phase6_colocalisation_matrix.png")
    plt.close(fig)


def fig_gradient(means_prop, sec_tab):
    apply_figure_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.5), gridspec_kw={"width_ratios": [1.6, 1]})
    x = np.arange(1, 11)
    ax = axes[0]
    for c in STATES:
        ax.plot(x, means_prop[c].values, marker="o", ms=3, color=COLOURS[c], label=SHORT[c])
    ax.plot(x, sec_tab["tumour_fraction_mean"].values, marker="s", ms=3, color="black", ls="--", label="tumour fraction")
    ax.set_xticks(x); ax.set_xticklabels([s.replace("_", " ") for s in BRAINSTEM_ORDER], rotation=45, ha="right")
    ax.set_ylabel("mean share of spot cells (q05)"); ax.set_xlabel("anatomical position (rostral → caudal)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax = axes[1]
    for c in TME:
        ax.plot(x, means_prop[c].values, marker="o", ms=3, color=COLOURS[c], label=SHORT[c])
    ax.set_xticks(x); ax.set_xticklabels([s.replace("_", " ") for s in BRAINSTEM_ORDER], rotation=45, ha="right")
    ax.set_xlabel("anatomical position (rostral → caudal)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.savefig(FIGS / "phase6_brainstem_gradient.png")
    plt.close(fig)


def main():
    with LOG.step("load"):
        q, meta = load_wide()
        LOG.write(f"spots {len(q)}, sections {meta.section_id.nunique()}, groups {meta.group.value_counts().to_dict()}")
    with LOG.step("colocalisation"):
        coloc, per = colocalisation(q, meta)
        key = coloc[(coloc.group == "DMG") & coloc.pair_type.eq("state_vs_TME")].sort_values("rho_re", ascending=False)
        LOG.write("DMG state-vs-TME (top/bottom):\n" + key[["class_a", "class_b", "n_sections", "rho_re", "rho_re_lo", "rho_re_hi", "rho_boot_lo", "rho_boot_hi", "n_sections_positive", "q_bh_state_vs_TME"]].round(3).to_string(index=False))
    with LOG.step("niches"):
        comp, bysec, labels = niches(q, meta)
        LOG.write("niches:\n" + comp.round(3).to_string())
    with LOG.step("brainstem gradient"):
        grad, sec_tab, means_prop = brainstem(q, meta)
        LOG.write("brainstem trend along axis:\n" + grad[grad.test == "trend_along_axis"][["cls", "rho_abs", "p_abs", "rho_prop", "p_prop", "rostral_mean", "caudal_mean"]].round(3).to_string(index=False))
        LOG.write("brainstem section means:\n" + sec_tab.iloc[:, :4].round(3).to_string())
        LOG.write("within-section RE:\n" + grad[grad.test.str.endswith("_RE")][["test", "cls", "rho_re", "rho_re_lo", "rho_re_hi", "p_re", "n_sections_positive"]].round(3).to_string(index=False))
    with LOG.step("Ren DMG vs GBM"):
        out, sec = ren_outgroup(q, meta)
        LOG.write("DMG vs GBM:\n" + out.round(3).to_string(index=False))
    with LOG.step("figures"):
        picks = fig_maps(q, meta); fig_coloc(coloc); fig_gradient(means_prop, sec_tab)
        LOG.write(f"figures written; map sections {picks}")
    LOG.write(f"done | peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
