#!/usr/bin/env python
"""Phase 9 science check / step 41 — do lineage metaprograms appear in whole cells when whole-cell cohorts are added?

The primary metaprograms (16_metaprograms.py) were discovered on six whole-cell samples, all 10x; every lineage-like
metaprogram was found only in nuclei. Filbin 2018 and Liu 2022, the Smart-seq2 cohorts in which OPC-, AC- and OC-like
states were first described, were projection-only. Here their malignant cells are factorised per sample with the unchanged
per-sample cNMF (dmg_atlas.metaprograms.run_cnmf_sample: 2,000 HVGs, k 4-9, 10 restarts, silhouette rule, top-50 genes),
their programmes are pooled with the 400 primary programmes, and the 16_metaprograms.py clustering and replication rule are
re-applied unchanged.
  cells        : consensus_call == 'malignant' and qc_pass (malignant_calls.csv; interim QC files)
  min cells    : MIN_CELLS per sample (100; Smart-seq2 detects ~3x more genes per cell than 10x, and the primary rule of 200
                 would drop 28 of 43 Liu samples) — a documented deviation for these two cohorts only
  input values : Liu RSEM expected counts (layers['counts'], rounded); Filbin TPM (rounded) treated as counts
  assay class  : Liu 'fresh cells' samples = cells, 'frozen nuclei' samples = nuclei (Smart-seq2 nuclei); Filbin = cells
  publications : liu, filbin added to the five primary publications
Recovery of a primary metaprogram = Jaccard >= 0.5 between its 50 genes and an extended-run metaprogram's genes.
Inputs : data/interim/{liu,filbin}.qc.h5ad, data/processed/malignant_calls.csv, results/tables/cnmf/*, metaprograms.csv,
         metaprogram_genes.csv
Outputs: results/tables/cnmf_wholecell/<sample>_programs.csv, results/tables/cnmf_wholecell/samples_analysed.csv,
         results/tables/quality/wholecell_metaprograms.csv (extended-run metaprograms), wholecell_mp_recovery.csv
         (per primary MP: matched extended MP, whole-cell samples and publications, replicated), wholecell_inputs.json
Run: nice -n 10 python scripts/41_wholecell_metaprograms.py [--force]
Log: logs/phase9_quality.log
"""
import sys, json, hashlib, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py, scipy.sparse as sp  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from scipy.cluster.hierarchy import linkage, fcluster  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402
from dmg_atlas import INTERIM, PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
log = U.StepLogger(LOGS / "phase9_quality.log", context="41_wholecell_metaprograms")
T = RESULTS / "tables"; CN = T / "cnmf"; WC = T / "cnmf_wholecell"; WC.mkdir(exist_ok=True); OUT = T / "quality"
OUTS = [OUT / f for f in ("wholecell_metaprograms.csv", "wholecell_mp_recovery.csv", "wholecell_inputs.json")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)
CUT, MIN_PROGRAMS, MIN_SAMPLES, MIN_GENE_FRAC, TOP_N, MIN_CELLS, RECOVER_J = 0.85, 5, 3, 0.25, 50, 100, 0.5
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


with log.step("per-sample cNMF on Liu and Filbin malignant cells"):
    calls = pd.read_csv(PROCESSED / "malignant_calls.csv", usecols=["cell_id", "study", "consensus_call"], low_memory=False)
    mal = set(calls.cell_id[calls.study.isin(["liu", "filbin"]) & calls.consensus_call.eq("malignant")])
    samp_rows, progs = [], []
    for study, layer in [("liu", "counts"), ("filbin", "tpm")]:
        with h5py.File(INTERIM / f"{study}.qc.h5ad", "r") as f:
            obs = read_elem(f["obs"]); genes = read_elem(f["var"]).index.to_numpy()
            X = read_elem(f["layers"][layer]); X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
        keep = obs.index.isin(mal) & obs.qc_pass.astype(bool).values
        for s, idx in pd.Series(np.flatnonzero(keep)).groupby(obs.sample_id.values[keep]):
            rows = np.sort(idx.values)
            ac = "nuclei" if study == "liu" and (obs.chemistry.values[rows] == "Smart-seq2 plates (frozen nuclei)").mean() > 0.5 else "cells"
            rec = dict(sample_id=s, study=study, assay_class=ac, n_malignant=len(rows), analysed=len(rows) >= MIN_CELLS)
            samp_rows.append(rec)
            if not rec["analysed"]:
                continue
            out = WC / f"{s}_programs.csv"
            if out.exists() and not args.force:
                pr = pd.read_csv(out)
            else:
                Xc = X[rows].copy(); Xc.data = np.rint(Xc.data)
                pr, _, _ = MP.run_cnmf_sample(Xc, genes, obs.index.values[rows], s)
                pr.to_csv(out, index=False)
            progs.append(pr)
    samp = pd.DataFrame(samp_rows); samp.to_csv(WC / "samples_analysed.csv", index=False)
    log.write(f"samples: {samp.groupby(['study', 'assay_class']).analysed.agg(['size', 'sum']).to_dict()}")

with log.step("re-cluster with the primary programmes"):
    s0 = pd.read_csv(CN / "samples_analysed.csv", index_col=0); s0 = s0[s0.analysed]
    s0 = pd.DataFrame(dict(study=s0.study, assay_class=s0.assay.map(ASSAY_CLASS)))
    s1 = samp[samp.analysed].set_index("sample_id")[["study", "assay_class"]]
    allsamp = pd.concat([s0, s1]); allsamp["publication"] = allsamp.study.map(lambda x: PUBLICATION.get(x, x))
    prog = pd.concat([pd.read_csv(CN / f"{s}_programs.csv") for s in s0.index] + progs, ignore_index=True)
    names = prog.program.drop_duplicates().to_list(); sets = [set(prog.gene[prog.program == n]) for n in names]
    meta = prog.drop_duplicates("program").set_index("program").loc[names, ["sample_id"]].join(allsamp, on="sample_id")
    D = 1 - MP.jaccard_matrix(sets); np.fill_diagonal(D, 0)
    lab = fcluster(linkage(squareform(D, checks=False), "average"), CUT, "distance")
    mps = [m for m in (np.flatnonzero(lab == c) for c in np.unique(lab)) if len(m) >= MIN_PROGRAMS and meta.sample_id.iloc[m].nunique() >= MIN_SAMPLES]
    mps.sort(key=lambda m: -len(m))
    rows = []
    for i, m in enumerate(mps):
        sub = meta.iloc[m]; mem = [names[j] for j in m]
        g = prog[prog.program.isin(mem)].groupby("gene").agg(n=("program", "nunique"), s=("spectra_score", "mean"))
        g = g[g.n >= MIN_GENE_FRAC * len(mem)].sort_values(["n", "s"], ascending=False).head(TOP_N)
        cells = sub[sub.assay_class == "cells"]
        rows.append(dict(metaprogram=f"xMP{i+1}", n_programs=len(m), n_samples=sub.sample_id.nunique(), publications=";".join(sorted(set(sub.publication))),
                         n_publications=sub.publication.nunique(), n_samples_cells=cells.sample_id.nunique(),
                         cell_publications=";".join(sorted(set(cells.publication))), n_samples_nuclei=sub[sub.assay_class == "nuclei"].sample_id.nunique(),
                         replicated=bool(sub.publication.nunique() >= 3 and set(sub.assay_class) >= {"cells", "nuclei"}),
                         replicated_strict=bool(sub.publication.nunique() >= 3 and cells.sample_id.nunique() >= 2 and (sub.assay_class == "nuclei").sum() >= 2),
                         genes=";".join(g.index)))
    ext = pd.DataFrame(rows); ext.to_csv(OUTS[0], index=False)

with log.step("recovery of primary metaprograms"):
    prim = pd.read_csv(T / "metaprograms.csv").set_index("metaprogram")
    pg = {m: set(g.gene) for m, g in pd.read_csv(T / "metaprogram_genes.csv").groupby("metaprogram")}
    eg = [set(x.split(";")) for x in ext.genes]
    rec = []
    for mp, G in pg.items():
        j = [len(G & e) / len(G | e) for e in eg]; b = int(np.argmax(j)) if j else -1
        r = ext.iloc[b] if j and j[b] >= RECOVER_J else None
        rec.append(dict(metaprogram=mp, family=prim.loc[mp, "family"], technical_class=prim.loc[mp, "technical_class"],
                        replicated_primary=bool(prim.loc[mp, "replicated"]), best_jaccard=round(j[b], 3) if j else 0.0,
                        matched=r.metaprogram if r is not None else "", n_samples_cells=int(r.n_samples_cells) if r is not None else 0,
                        cell_publications=r.cell_publications if r is not None else "",
                        replicated_extended=bool(r.replicated) if r is not None else False,
                        replicated_strict_extended=bool(r.replicated_strict) if r is not None else False))
    rec = pd.DataFrame(rec); rec.to_csv(OUTS[1], index=False)
    new = ext[~ext.metaprogram.isin(set(rec.matched))]
    log.write(f"extended run: {len(names)} programmes, {len(ext)} metaprograms, {int(ext.replicated.sum())} replicated "
              f"({int(ext.replicated_strict.sum())} strict)\n" + rec.to_string(index=False)
              + "\nextended metaprograms without a primary match:\n" + new[["metaprogram", "n_samples", "publications", "n_samples_cells", "replicated", "genes"]].assign(genes=new.genes.str[:80]).to_string(index=False))
    json.dump(dict(script="scripts/41_wholecell_metaprograms.py", params=dict(min_cells=MIN_CELLS, cut=CUT, recover_j=RECOVER_J),
                   inputs={str(p.relative_to(RESULTS.parent)): md5(p) for p in [INTERIM / "liu.qc.h5ad", INTERIM / "filbin.qc.h5ad",
                                                                               PROCESSED / "malignant_calls.csv", T / "metaprograms.csv"]}),
              open(OUTS[2], "w"), indent=1)
