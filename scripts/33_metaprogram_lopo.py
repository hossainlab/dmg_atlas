#!/usr/bin/env python
"""Phase 9 quality assessment / step 33 — publication-dependence of the metaprograms (leave-one-publication-out + balanced).

Question: are the 15 primary metaprograms (16_metaprograms.py, primary selection) and the 2 replicated ones (MP1 cycling,
MP4 hypoxia) driven by one publication? The per-sample cNMF programmes are re-clustered with the UNCHANGED 16_metaprograms.py
rule (top-50 gene sets, Jaccard, average linkage, cut 0.85, MP = >= 5 programmes from >= 3 samples; replicated = members
from >= 3 publications AND from cells and nuclei; jessa_scrna + jessa_multiome = jessa2022) on:
  all          : every analysed sample (must reproduce 16_metaprograms.py: asserted)
  drop_<pub>   : every sample except those of one publication (5 scenarios)
  balanced_s<i>: at most N_BAL samples per publication (seeded random subsample, N_SEEDS seeds), so no publication can
                 dominate the programme pool (by samples Jessa 2022 contributes 25, Okada 18, Collot 5, DeSisto 1, Zang 1)
  random_s<i>  : size-matched control — as many samples as a balanced scenario, drawn at random ignoring publication; a
                 balanced-vs-random gap is the publication effect, a gap between both and 'all' is the pool-size effect
Recovery of a primary MP in a scenario = largest Jaccard between its 50-gene set and any scenario MP's gene set (built
with the 16_metaprograms.py gene rule); recovered if >= RECOVER_J. Dropping a publication leaves 4, so the >= 3-publication
replication rule is harder to meet there; the drop rows report it unchanged and say so.
Inputs : results/tables/cnmf/{samples_analysed.csv,<sample>_programs.csv}, results/tables/{metaprograms,metaprogram_genes}.csv
Outputs (results/tables/quality/, skipped if present unless --force):
         metaprogram_lopo_scenarios.csv   one row per scenario: samples, programmes, MPs, replicated MPs, replicated MP genes
         metaprogram_lopo_recovery.csv    primary MP x scenario: best Jaccard, matched scenario MP, recovered, replicated there
         metaprogram_lopo_summary.csv     per primary MP: recovery rate over drop / balanced scenarios, which drop loses it
         metaprogram_lopo_inputs.json     md5 of every input
Run: nice -n 10 python scripts/33_metaprogram_lopo.py [--force]
Log: logs/phase9_quality.log
"""
import sys, json, hashlib, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd  # noqa: E402
from scipy.cluster.hierarchy import linkage, fcluster  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402
from dmg_atlas import RESULTS, LOGS  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
log = U.StepLogger(LOGS / "phase9_quality.log", context="33_metaprogram_lopo")
TABLES = RESULTS / "tables"; CNMF = TABLES / "cnmf"; OUT = TABLES / "quality"; OUT.mkdir(parents=True, exist_ok=True)
CUT, MIN_PROGRAMS, MIN_SAMPLES, MIN_GENE_FRAC, TOP_N = 0.85, 5, 3, 0.25, 50      # = 16_metaprograms.py
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}
N_BAL, N_SEEDS, RECOVER_J = 5, 20, 0.5
OUTS = [OUT / f"metaprogram_lopo_{n}" for n in ("scenarios.csv", "recovery.csv", "summary.csv", "inputs.json")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def cluster(prog, samples):
    """16_metaprograms.py clustering + replication + gene rule on the programmes of `samples`."""
    p = prog[prog.sample_id.isin(samples)]
    names = p.program.drop_duplicates().to_list()
    sets = [set(p.gene[p.program == n]) for n in names]
    meta = p.drop_duplicates("program").set_index("program").loc[names, ["sample_id"]].join(samp[["publication", "assay_class"]], on="sample_id")
    D = 1 - MP.jaccard_matrix(sets); np.fill_diagonal(D, 0)
    lab = fcluster(linkage(squareform(D, checks=False), "average"), CUT, "distance")
    mps = [m for m in (np.flatnonzero(lab == c) for c in np.unique(lab)) if len(m) >= MIN_PROGRAMS and meta.sample_id.iloc[m].nunique() >= MIN_SAMPLES]
    mps.sort(key=lambda m: -len(m))
    out = []
    for m in mps:
        sub = meta.iloc[m]; mem = [names[j] for j in m]
        g = p[p.program.isin(mem)].groupby("gene").agg(n=("program", "nunique"), s=("spectra_score", "mean"))
        g = g[g.n >= MIN_GENE_FRAC * len(mem)].sort_values(["n", "s"], ascending=False).head(TOP_N)
        out.append(dict(genes=set(g.index), n_programs=len(m), publications=sorted(set(sub.publication)),
                        replicated=bool(sub.publication.nunique() >= 3 and set(sub.assay_class) >= {"cells", "nuclei"})))
    return out, len(names)


with log.step("load programmes"):
    samp = pd.read_csv(CNMF / "samples_analysed.csv", index_col=0)
    samp = samp[samp.analysed].copy()
    samp["publication"] = samp.study.map(lambda s: PUBLICATION.get(s, s))
    samp["assay_class"] = samp.assay.map(ASSAY_CLASS)
    prog = pd.concat([pd.read_csv(CNMF / f"{s}_programs.csv") for s in samp.index], ignore_index=True)
    prim = pd.read_csv(TABLES / "metaprograms.csv").set_index("metaprogram")
    pg = pd.read_csv(TABLES / "metaprogram_genes.csv")
    prim_genes = {m: set(g.gene) for m, g in pg.groupby("metaprogram")}
    inputs = {str(p.relative_to(RESULTS.parent)): md5(p) for p in
              [CNMF / "samples_analysed.csv", TABLES / "metaprograms.csv", TABLES / "metaprogram_genes.csv"] + [CNMF / f"{s}_programs.csv" for s in samp.index]}
    log.write(f"{prog.program.nunique()} programmes, {len(samp)} samples; samples per publication {samp.publication.value_counts().to_dict()}")

with log.step("scenarios"):
    pubs = sorted(samp.publication.unique())
    scen = {"all": list(samp.index)} | {f"drop_{p}": list(samp.index[samp.publication != p]) for p in pubs}
    for i in range(N_SEEDS):
        rng = np.random.default_rng(i)
        scen[f"balanced_s{i}"] = [s for p in pubs for s in rng.permutation(samp.index[samp.publication == p])[:N_BAL]]
    n_bal = len(scen["balanced_s0"])
    for i in range(N_SEEDS):  # size-matched control: same number of samples, drawn ignoring publication
        scen[f"random_s{i}"] = list(np.random.default_rng(1000 + i).choice(samp.index, n_bal, replace=False))
    rows, rec = [], []
    for name, ss in scen.items():
        mps, n_prog = cluster(prog, ss)
        reps = [m for m in mps if m["replicated"]]
        rows.append(dict(scenario=name, kind=name.split("_")[0], n_samples=len(ss), n_programmes=n_prog, n_metaprograms=len(mps),
                         n_replicated=len(reps), samples_per_publication=json.dumps(samp.loc[ss].publication.value_counts().to_dict())))
        for mp, G in prim_genes.items():
            j = [len(G & m["genes"]) / len(G | m["genes"]) for m in mps]
            b = int(np.argmax(j)) if j else -1
            rec.append(dict(metaprogram=mp, scenario=name, best_jaccard=round(j[b], 4) if j else 0.0, matched_rank=b + 1 if j else None,
                            recovered=bool(j and j[b] >= RECOVER_J), replicated_in_scenario=bool(j and j[b] >= RECOVER_J and mps[b]["replicated"]),
                            matched_publications=";".join(mps[b]["publications"]) if j else ""))
        rows[-1]["replicated_matches"] = ";".join(r["metaprogram"] for r in rec[-len(prim_genes):] if r["replicated_in_scenario"])
    sc = pd.DataFrame(rows); rc = pd.DataFrame(rec)
    a = sc.set_index("scenario").loc["all"]
    assert a.n_metaprograms == len(prim) and a.n_replicated == int(prim.replicated.sum()), \
        f"re-clustering does not reproduce 16_metaprograms.py: {a.n_metaprograms} MPs / {a.n_replicated} replicated vs {len(prim)} / {int(prim.replicated.sum())}"
    assert rc[rc.scenario == "all"].best_jaccard.min() == 1.0, "primary MP gene sets not reproduced exactly"
    log.write("scenarios:\n" + sc.to_string(index=False))

with log.step("summary + write"):
    rc["kind"] = rc.scenario.str.split("_").str[0]
    summ = []
    for mp, g in rc.groupby("metaprogram", sort=False):
        d, b, r = g[g.kind == "drop"], g[g.kind == "balanced"], g[g.kind == "random"]
        summ.append(dict(metaprogram=mp, label=prim.loc[mp, "label"], family=prim.loc[mp, "family"], technical_class=prim.loc[mp, "technical_class"],
                         publications=prim.loc[mp, "publications"], replicated_primary=bool(prim.loc[mp, "replicated"]),
                         recovered_drop=f"{int(d.recovered.sum())}/{len(d)}", lost_when_dropping=";".join(d.scenario[~d.recovered].str[5:]),
                         median_best_jaccard_drop=round(d.best_jaccard.median(), 3),
                         recovered_balanced_frac=round(b.recovered.mean(), 3), replicated_balanced_frac=round(b.replicated_in_scenario.mean(), 3),
                         median_best_jaccard_balanced=round(b.best_jaccard.median(), 3),
                         recovered_random_frac=round(r.recovered.mean(), 3), replicated_random_frac=round(r.replicated_in_scenario.mean(), 3)))
    summ = pd.DataFrame(summ)
    summ["metaprogram"] = pd.Categorical(summ.metaprogram, [f"MP{i}" for i in range(1, len(summ) + 1)]); summ = summ.sort_values("metaprogram")
    sc.to_csv(OUTS[0], index=False); rc.drop(columns="kind").to_csv(OUTS[1], index=False); summ.to_csv(OUTS[2], index=False)
    json.dump(dict(script="scripts/33_metaprogram_lopo.py", params=dict(cut=CUT, min_programs=MIN_PROGRAMS, min_samples=MIN_SAMPLES, n_bal=N_BAL,
                                                                     n_seeds=N_SEEDS, recover_j=RECOVER_J), inputs=inputs), open(OUTS[3], "w"), indent=1)
    for k in ("balanced", "random"):
        s = sc[sc.kind == k]
        log.write(f"{k} ({N_SEEDS} seeds, {s.n_samples.iloc[0]} samples): MPs median {s.n_metaprograms.median():.0f} [{s.n_metaprograms.min()}-{s.n_metaprograms.max()}], "
                  f"replicated median {s.n_replicated.median():.0f} [{s.n_replicated.min()}-{s.n_replicated.max()}]")
    log.write(summ.to_string(index=False))
