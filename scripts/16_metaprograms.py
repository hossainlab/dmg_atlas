#!/usr/bin/env python
"""Phase 6 / step 16 — consensus metaprograms from the per-sample cNMF programmes and the pre-registered replication test.

Inputs : results/tables/cnmf<tag>/<sample>_programs.csv (top-50 genes per programme, from 13_cnmf_per_sample.py --selection ...),
         results/tables/cnmf<tag>/samples_analysed.csv (study, assay per sample), data/processed/malignant_selection_audit.parquet,
         reference gene sets (audit MC-03, all scripted):
           published : results/tables/genesets/reference_gene_sets_scripted.csv, Liu2022_TableS2a_* (Liu 2022 Nat Genet Table S2a)
           inhouseDE : same file, Liu2022_DE_author_states_rederived_* (in-house Wilcoxon DE on the Liu author states; NOT published lists)
           curated   : results/tables/genesets/phase6_curated_gene_sets.csv (dmg_atlas.states.CURATED_PHASE6_SETS; hand-curated)
           filbin_untraceable : dmg_atlas.states.STATE_GENE_SETS (Filbin 2018-derived lists that cannot be traced to a downloaded table)
Outputs (results/tables/, suffix <tag> = '' for the primary selection, '_eventtrust' for the sensitivity selection):
         metaprograms<tag>.csv                  one row per metaprogram (size, studies, assays, patients, replication calls, label, class)
         metaprogram_genes<tag>.csv             MP gene set (rank, gene, n member programmes containing it, mean spectra score)
         metaprogram_members<tag>.csv           every programme -> metaprogram (or 'unassigned') + mean Jaccard to its MP
         replication_matrix<tag>.csv            MP x (publication, assay class): number of samples contributing a member programme
         metaprogram_reference_overlap<tag>.csv MP x reference set: overlap, Jaccard, hypergeometric p, BH q, provenance tier
         metaprogram_cut_sensitivity<tag>.csv   MPs / replicated MPs at tree cuts 0.75-0.92 in 0.01 steps (audit MC-02)
Figures (results/figures/, primary only): phase6_metaprogram_jaccard.png, phase6_replication_matrix.png.
Run: nice -n 10 python scripts/16_metaprograms.py [--selection primary|eventtrust] [--force]

Method (fixed before looking at the clustering; unchanged by the audit except where marked)
  1. Every programme's top-50 genes (by spectra score) is a set; pairwise Jaccard J over all programmes of all analysed
     samples. Average-linkage hierarchical clustering on 1 - J, tree cut at distance CUT (0.85 => mean pairwise
     Jaccard >= 0.15, i.e. ~13 of 50 shared genes). A cluster is a metaprogram if it holds >= MIN_PROGRAMS programmes
     from >= MIN_SAMPLES distinct samples; everything else is 'unassigned'.
  2. MP gene set: genes ranked by the number of member programmes listing them (ties: mean spectra score), restricted
     to genes in >= MIN_GENE_FRAC of members, top 50.
  3. Replication (docs/decisions.md #9, pre-registered): replicated if the MP has member programmes from >= 3 publications
     AND from both cells (scRNA) and nuclei (snRNA, multiome_RNA). jessa_scrna and jessa_multiome are both Jessa 2022 and
     count once. Audit additions: replicated_strict additionally requires >= 2 samples per assay class (MC-05); n_patients
     counts distinct tumours (MC-17); the cut-sensitivity curve runs 0.75-0.92 in 0.01 steps (MC-02).
  4. Class and label (audit MC-03/07/08, P7-04). Technical classes are assigned BEFORE lineage labelling from the MP gene set:
     technical:ribosomal (>= 15 ribosomal-protein genes, explicit RP regex, MC-15), technical:mitochondrial (>= 5 MT- genes or
     MT- + RP >= 15; technical:mitochondrial_ribosomal when both rules fire), technical:IEG_dissociation_stress (>= 10 genes of the curated IEG/heat-shock dissociation set),
     contamination:immune (>= 5 genes of the curated pan-immune set). These MPs get no lineage family and are excluded from
     the lineage counts. Lineage MPs are labelled by the reference set with the largest overlap among sets with BH q < 0.05
     and >= 5 shared genes (universe = all genes in any programme's top-50 list), taking the highest provenance tier that has a
     qualifying set: published (Liu Table S2a) and curated functional sets first, then in-house DE sets, then the untraceable
     Filbin-derived sets; label_provenance records the tier. MPs without any qualifying set are 'novel:' + top-3 genes.
"""
import re, sys, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd  # noqa: E402
import scipy.stats as ss  # noqa: E402
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402
from dmg_atlas import RESULTS, LOGS, PROCESSED  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402
from dmg_atlas.states import STATE_GENE_SETS, CURATED_PHASE6_SETS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true")
ap.add_argument("--selection", default="primary", choices=["primary", "eventtrust"]); args = ap.parse_args()
TAG = "" if args.selection == "primary" else "_eventtrust"
CNMF = RESULTS / "tables" / f"cnmf{TAG}"; TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
log = U.StepLogger(LOGS / "fix_metaprograms.log", context=f"16_metaprograms[{args.selection}]")
CUT, MIN_PROGRAMS, MIN_SAMPLES, MIN_GENE_FRAC, TOP_N = 0.85, 5, 3, 0.25, 50
SENS_CUTS = [round(c, 2) for c in np.arange(0.75, 0.9201, 0.01)]
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}
RP_RE = re.compile(r"^RP[SL]\d+[A-Z]?$|^RPLP[0-2]$|^RPSA$")      # ribosomal proteins only (MC-15: no RPS6K*, no *L paralogues)
MT_RE = re.compile(r"^MT-")
TECH = dict(ribo=15, mito=5, ieg=10, immune=5)
FAMILY_PATTERNS = ((r"G2M|_S$|^S$|[Cc]ycling", "cycling"), (r"OPC", "OPC_like"), (r"(^|_)AC", "AC_like"), (r"(^|_)OC", "OC_like"),
                   (r"MES", "MES_like"), (r"[Hh]ypoxia", "hypoxia"))
OUTS = [TABLES / f"{n}{TAG}.csv" for n in ("metaprograms", "metaprogram_genes", "metaprogram_members", "replication_matrix",
                                          "metaprogram_reference_overlap", "metaprogram_cut_sensitivity")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)


def family_of(label: str) -> str:
    r = label.split(":")[-1]
    for pat, fam in FAMILY_PATTERNS:
        if re.search(pat, r):
            return fam
    return r


def technical_class(genes: set) -> tuple[str | None, dict]:
    n = dict(n_ribo=sum(bool(RP_RE.match(g)) for g in genes), n_mito=sum(bool(MT_RE.match(g)) for g in genes),
             n_ieg=len(genes & set(CURATED_PHASE6_SETS["IEG_dissociation"])), n_immune=len(genes & set(CURATED_PHASE6_SETS["Immune_pan"])))
    if n["n_mito"] >= TECH["mito"] and n["n_ribo"] >= TECH["ribo"]:
        return "technical:mitochondrial_ribosomal", n
    if n["n_ribo"] >= TECH["ribo"]:
        return "technical:ribosomal", n
    if n["n_mito"] >= TECH["mito"] or n["n_mito"] + n["n_ribo"] >= TECH["ribo"]:
        return "technical:mitochondrial", n
    if n["n_ieg"] >= TECH["ieg"]:
        return "technical:IEG_dissociation_stress", n
    if n["n_immune"] >= TECH["immune"]:
        return "contamination:immune", n
    return None, n


with log.step("load programmes"):
    samp = pd.read_csv(CNMF / "samples_analysed.csv", index_col=0)
    samp = samp[samp.analysed]
    samp["publication"] = samp.study.map(lambda s: PUBLICATION.get(s, s))
    samp["assay_class"] = samp.assay.map(ASSAY_CLASS)
    assert samp.assay_class.notna().all(), samp.assay[samp.assay_class.isna()].unique()
    selp = pd.read_parquet(PROCESSED / "malignant_selection_audit.parquet", columns=["sample_id", "patient_id", args.selection, "cnv_confidence"])
    selp = selp[selp[args.selection]]
    samp["patient_id"] = selp.groupby("sample_id").patient_id.first().reindex(samp.index)
    samp["frac_score_only"] = (selp.cnv_confidence == "score_only").groupby(selp.sample_id).mean().reindex(samp.index).fillna(0)
    prog = pd.concat([pd.read_csv(CNMF / f"{s}_programs.csv") for s in samp.index], ignore_index=True)
    names = prog.program.drop_duplicates().to_list()
    sets = [set(prog.gene[prog.program == p]) for p in names]
    pmeta = prog.drop_duplicates("program").set_index("program").loc[names, ["sample_id", "k"]]
    pmeta = pmeta.join(samp[["study", "publication", "assay", "assay_class", "patient_id", "frac_score_only"]], on="sample_id")
    log.write(f"{len(names)} programmes from {len(samp)} samples ({samp.patient_id.nunique()} patients); publications "
              f"{samp.publication.value_counts().to_dict()}; assay classes (samples) {samp.assay_class.value_counts().to_dict()}; "
              f"samples with > 50 % score-only cells: {sorted(samp.index[samp.frac_score_only > 0.5])}")

with log.step("reference sets"):
    refs, tier = {}, {}
    sc = pd.read_csv(TABLES / "genesets" / "reference_gene_sets_scripted.csv")
    for k, g in sc.groupby("gene_set"):
        if k.startswith("Liu2022_TableS2a_"):
            refs[f"published:{k}"] = set(g.gene); tier[f"published:{k}"] = 1
        elif "rederived" in k:
            refs[f"inhouseDE:{k}"] = set(g.gene); tier[f"inhouseDE:{k}"] = 2
    cu = pd.read_csv(TABLES / "genesets" / "phase6_curated_gene_sets.csv")
    for k, g in cu.groupby("gene_set"):
        refs[f"curated:{k}"] = set(g.gene); tier[f"curated:{k}"] = 1
    for k, v in STATE_GENE_SETS.items():
        refs[f"filbin_untraceable:{k}"] = set(v); tier[f"filbin_untraceable:{k}"] = 3
    universe = set(prog.gene); N = len(universe)
    log.write(f"{len(refs)} reference sets (tiers {pd.Series(tier).value_counts().to_dict()}); universe {N} genes")

with log.step("cluster"):
    J = MP.jaccard_matrix(sets)
    D = 1 - J; np.fill_diagonal(D, 0)
    Z = linkage(squareform(D, checks=False), "average")

    def call(cut):
        lab = fcluster(Z, cut, "distance")
        mps = []
        for c in np.unique(lab):
            m = np.flatnonzero(lab == c)
            if len(m) >= MIN_PROGRAMS and pmeta.sample_id.iloc[m].nunique() >= MIN_SAMPLES:
                mps.append(m)
        mps.sort(key=lambda m: -len(m))
        return mps

    def replicated(m, strict=False):
        sub = pmeta.iloc[m]
        ok = sub.publication.nunique() >= 3 and set(sub.assay_class) >= {"cells", "nuclei"}
        if strict:
            ok = ok and sub.loc[sub.assay_class == "cells", "sample_id"].nunique() >= 2 and sub.loc[sub.assay_class == "nuclei", "sample_id"].nunique() >= 2
        return bool(ok)

    def mp_genes(m):
        mem = [names[j] for j in m]
        sub = prog[prog.program.isin(mem)]
        g = sub.groupby("gene").agg(n_programs=("program", "nunique"), mean_spectra=("spectra_score", "mean"))
        g = g[g.n_programs >= MIN_GENE_FRAC * len(mem)].sort_values(["n_programs", "mean_spectra"], ascending=False).head(TOP_N)
        return mem, g

    def best_label(genes):
        rows = []
        for ref, rs in refs.items():
            rs_u = rs & universe; k = len(genes & rs_u)
            p = ss.hypergeom.sf(k - 1, N, len(rs_u), len(genes)) if k else 1.0
            rows.append(dict(reference_set=ref, tier=tier[ref], n_ref_in_universe=len(rs_u), overlap=k,
                             jaccard=k / len(genes | rs_u) if genes | rs_u else 0.0, p_hyper=p, shared_genes=";".join(sorted(genes & rs_u))))
        return pd.DataFrame(rows)

    sens = []
    for cut in SENS_CUTS:
        mps = call(cut)
        reps, reps_s, lin_reps, lin_mps = [], [], [], 0
        for m in mps:
            _, g = mp_genes(m); tech, _ = technical_class(set(g.index))
            r, rs_ = replicated(m), replicated(m, strict=True)
            if tech is None:
                lin_mps += 1
                if r:
                    lin_reps.append(";".join(sorted(set(g.index[:3]))))
            reps.append(r); reps_s.append(rs_)
        sens.append(dict(cut=cut, mean_jaccard_min=round(1 - cut, 2), n_metaprograms=len(mps), n_replicated=int(sum(reps)),
                         n_replicated_strict=int(sum(reps_s)), n_lineage_metaprograms=lin_mps, n_lineage_replicated=len(lin_reps),
                         n_programmes_assigned=int(sum(len(m) for m in mps)), primary=cut == CUT,
                         lineage_replicated_top_genes="|".join(lin_reps)))
    sens = pd.DataFrame(sens)
    mps = call(CUT)
    log.write(f"cut {CUT}: {len(mps)} metaprograms, {sum(len(m) for m in mps)}/{len(names)} programmes assigned\n"
              + sens.drop(columns="lineage_replicated_top_genes").to_string(index=False))

with log.step("gene sets, replication, labels"):
    rows_mp, rows_genes, rows_ov, rows_mem = [], [], [], []
    member_of = pd.Series("unassigned", index=names)
    for i, m in enumerate(mps):
        mp = f"MP{i+1}"
        mem, g = mp_genes(m); member_of[mem] = mp
        genes = set(g.index)
        for r, (gene, x) in enumerate(g.iterrows()):
            rows_genes.append(dict(metaprogram=mp, rank=r + 1, gene=gene, n_programs=int(x.n_programs),
                                   frac_programs=x.n_programs / len(mem), mean_spectra_score=x.mean_spectra))
        ov = best_label(genes); ov.insert(0, "metaprogram", mp); rows_ov.append(ov)
        tech, tn = technical_class(genes)
        pm = pmeta.loc[mem]
        Jm = J[np.ix_(m, m)]
        for j, p in zip(m, mem):
            rows_mem.append(dict(program=p, metaprogram=mp, mean_jaccard_to_mp=(J[j, m].sum() - 1) / max(len(m) - 1, 1)))
        rows_mp.append(dict(metaprogram=mp, n_programs=len(mem), n_samples=pm.sample_id.nunique(), n_patients=pm.patient_id.nunique(),
                            n_publications=pm.publication.nunique(), n_datasets=pm.study.nunique(),
                            publications=";".join(sorted(pm.publication.unique())), datasets=";".join(sorted(pm.study.unique())),
                            n_samples_cells=pm.loc[pm.assay_class == "cells", "sample_id"].nunique(),
                            n_samples_nuclei=pm.loc[pm.assay_class == "nuclei", "sample_id"].nunique(),
                            n_patients_cells=pm.loc[pm.assay_class == "cells", "patient_id"].nunique(),
                            n_patients_nuclei=pm.loc[pm.assay_class == "nuclei", "patient_id"].nunique(),
                            n_samples_score_only_dominant=pm.loc[pm.frac_score_only > 0.5, "sample_id"].nunique(),
                            mean_within_jaccard=(Jm.sum() - len(m)) / (len(m) * (len(m) - 1)),
                            replicated=replicated(m), replicated_strict=replicated(m, strict=True), n_genes=len(genes),
                            technical_class=tech, **tn, top_genes=";".join(list(g.index[:15]))))
    ov = pd.concat(rows_ov, ignore_index=True)
    ov["q_bh"] = ss.false_discovery_control(ov.p_hyper.values, method="bh") if len(ov) else []
    mpt = pd.DataFrame(rows_mp)
    labels, prov, bo, bq = [], [], [], []
    for r in mpt.itertuples():
        o = ov[(ov.metaprogram == r.metaprogram) & (ov.q_bh < 0.05) & (ov.overlap >= 5)]
        if r.technical_class is not None:
            labels.append(r.technical_class); prov.append("technical_rule")
            bo.append(int(o.overlap.max()) if len(o) else 0); bq.append(o.q_bh.min() if len(o) else np.nan); continue
        if len(o):
            t = o.tier.min(); ot = o[o.tier == t].sort_values(["overlap", "q_bh"], ascending=[False, True]).iloc[0]
            labels.append(ot.reference_set); prov.append({1: "published_or_curated", 2: "inhouse_DE", 3: "filbin_untraceable"}[t])
            bo.append(int(ot.overlap)); bq.append(ot.q_bh)
        else:
            labels.append("novel:" + "/".join(r.top_genes.split(";")[:3])); prov.append("none"); bo.append(0); bq.append(np.nan)
    mpt["label"] = labels; mpt["label_provenance"] = prov; mpt["best_overlap"] = bo; mpt["best_q"] = bq
    mpt["best_reference_set"] = [l if p not in ("technical_rule", "none") else None for l, p in zip(labels, prov)]
    mpt["family"] = [None if p in ("technical_rule",) else (family_of(l) if p != "none" else f"novel_{m}")
                     for l, p, m in zip(labels, prov, mpt.metaprogram)]
    # best published-set overlap alongside, so a label from a lower tier can be read against the published lists
    pub = ov[ov.reference_set.str.startswith("published:")].sort_values("overlap", ascending=False).drop_duplicates("metaprogram").set_index("metaprogram")
    mpt["best_published_set"] = mpt.metaprogram.map(pub.reference_set); mpt["best_published_overlap"] = mpt.metaprogram.map(pub.overlap)
    mpt["best_published_q"] = mpt.metaprogram.map(pub.q_bh)
    mem = pmeta.join(member_of.rename("metaprogram")).reset_index(names="program")
    mem = mem.merge(pd.DataFrame(rows_mem)[["program", "mean_jaccard_to_mp"]], on="program", how="left")
    rep = (mem[mem.metaprogram != "unassigned"].groupby(["metaprogram", "publication", "assay_class"]).sample_id.nunique()
           .unstack(["publication", "assay_class"], fill_value=0))
    rep.columns = [f"{a}|{b}" for a, b in rep.columns]
    rep = rep.reindex(mpt.metaprogram).fillna(0).astype(int).reindex(sorted(rep.columns), axis=1)
    log.write("metaprograms:\n" + mpt[["metaprogram", "n_programs", "n_samples", "n_patients", "n_publications", "n_samples_cells", "n_samples_nuclei",
                                      "replicated", "replicated_strict", "label", "label_provenance", "family", "best_overlap",
                                      "best_published_set", "best_published_overlap"]].round(3).to_string(index=False))
    lin = mpt[mpt.technical_class.isna()]
    log.write(f"replicated: {int(mpt.replicated.sum())} of {len(mpt)} (strict {int(mpt.replicated_strict.sum())}); lineage MPs "
              f"{len(lin)}, lineage replicated {int(lin.replicated.sum())}; technical/contamination: "
              + str(mpt.set_index('metaprogram').technical_class.dropna().to_dict()))

with log.step("write"):
    mpt.to_csv(OUTS[0], index=False); pd.DataFrame(rows_genes).to_csv(OUTS[1], index=False); mem.to_csv(OUTS[2], index=False)
    rep.to_csv(OUTS[3]); ov.to_csv(OUTS[4], index=False); sens.to_csv(OUTS[5], index=False)

if TAG == "":
    with log.step("figures"):
        apply_figure_style()
        o = leaves_list(Z)
        fig, ax = plt.subplots(figsize=(4.6, 3.8))
        im = ax.imshow(J[np.ix_(o, o)], cmap="magma", vmin=0, vmax=0.5, interpolation="nearest")
        lab_o = member_of.iloc[o].to_numpy()
        centres = []
        for mp in mpt.metaprogram:
            pos = np.flatnonzero(lab_o == mp)
            if len(pos):
                ax.add_patch(plt.Rectangle((pos.min() - 0.5, pos.min() - 0.5), len(pos), len(pos), fill=False, ec="cyan", lw=0.6))
                centres.append((pos.mean(), mp))
        centres.sort(); ys = [c[0] for c in centres]; min_gap = len(names) / 55  # de-overlap right-margin labels (MC-20)
        for k in range(1, len(ys)):
            ys[k] = max(ys[k], ys[k - 1] + min_gap)
        for k in range(len(ys) - 2, -1, -1):
            ys[k] = min(ys[k], ys[k + 1] - min_gap)
        for (y0, mp), y in zip(centres, ys):
            ax.annotate(mp, xy=(len(names) - 0.5, y0), xytext=(len(names) + 6, y), fontsize=4.5, va="center", annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", lw=0.3, color="0.4"))
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(f"{len(names)} programmes (average linkage on 1 − Jaccard)")
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.16); cb.set_label("Jaccard (top-50 genes)")
        fig.savefig(FIGS / "phase6_metaprogram_jaccard.png"); plt.close(fig)
        fig, ax = plt.subplots(figsize=(max(3.0, 0.45 * rep.shape[1] + 1.8), max(2.0, 0.22 * len(rep) + 1.0)))
        ax.imshow(rep.values > 0, cmap="Greys", vmin=0, vmax=1.6, aspect="auto")
        for i in range(rep.shape[0]):
            for j in range(rep.shape[1]):
                if rep.iat[i, j]:
                    ax.text(j, i, rep.iat[i, j], ha="center", va="center", fontsize=5, color="white")
        ax.set_xticks(range(rep.shape[1])); ax.set_xticklabels([c.replace("|", " | ") for c in rep.columns], rotation=45, ha="right")
        ax.set_yticks(range(len(rep)))
        ax.set_yticklabels([f"{m} {'✓' if r else ''} {l.split(':', 1)[-1]}"[:44] for m, r, l in zip(mpt.metaprogram, mpt.replicated, mpt.label)])
        ax.set_xlabel("publication | assay class"); ax.set_ylabel("metaprogram (✓ replicated)")
        fig.savefig(FIGS / "phase6_replication_matrix.png"); plt.close(fig)
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
