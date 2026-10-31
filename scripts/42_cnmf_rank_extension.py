#!/usr/bin/env python
"""Phase 9 science check / step 42 — was the cNMF rank range (k 4-9) binding?

21 of 50 samples selected the largest rank tested (k = 9). For those samples, consensus NMF is repeated for k = 10-12 with the
unchanged procedure (dmg_atlas.metaprograms.prep_matrix / consensus_nmf / choose_k; 10 restarts; selection = largest k with
silhouette >= 0.90 and every consensus cluster drawn from >= 80% of restarts). If a k in 10-12 qualifies, the sample's programmes
are replaced by those at the new k; otherwise its k = 9 programmes are kept. The 16_metaprograms.py clustering and replication
rule are then re-applied and compared with the primary result.
Inputs : data/processed/atlas_core.h5ad (layers['counts'] rows of selected malignant cells), data/processed/malignant_selection_audit.parquet,
         results/tables/cnmf/{k_selection,samples_analysed,<sample>_programs}.csv, metaprograms.csv, metaprogram_genes.csv
Outputs: results/tables/cnmf_rank12/<sample>_programs.csv and <sample>_kselection.csv (extended samples);
         results/tables/quality/cnmf_rank_extension.csv (per sample: old / new k, silhouettes) and
         cnmf_rank_extension_recovery.csv (primary MP recovery and replication after re-clustering), cnmf_rank_inputs.json
Run: nice -n 10 python scripts/42_cnmf_rank_extension.py [--force]
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
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
log = U.StepLogger(LOGS / "phase9_quality.log", context="42_cnmf_rank_extension")
T = RESULTS / "tables"; CN = T / "cnmf"; R12 = T / "cnmf_rank12"; R12.mkdir(exist_ok=True); OUT = T / "quality"
OUTS = [OUT / f for f in ("cnmf_rank_extension.csv", "cnmf_rank_extension_recovery.csv", "cnmf_rank_inputs.json")]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)
CORE = PROCESSED / "atlas_core.h5ad"; KS_EXT = (10, 11, 12)
CUT, MIN_PROGRAMS, MIN_SAMPLES, MIN_GENE_FRAC, TOP_N, RECOVER_J = 0.85, 5, 3, 0.25, 50, 0.5
PUBLICATION = {"jessa_scrna": "jessa2022", "jessa_multiome": "jessa2022"}
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


with h5py.File(CORE) as f:
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in ["sample_id"]}, index=read_elem(f["obs"]["_index"])); var = read_elem(f["var"])
sel = pd.read_parquet(PROCESSED / "malignant_selection_audit.parquet", columns=["primary"]).reindex(obs.index).primary.fillna(False).astype(bool)


def read_rows(idx):
    idx = np.asarray(idx); runs = np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1); parts = []
    with h5py.File(CORE) as f:
        g = f["layers/counts"]; indptr = g["indptr"]
        for r in runs:
            a, b = int(r[0]), int(r[-1]) + 1
            ip = indptr[a:b + 1]; s, e = int(ip[0]), int(ip[-1])
            parts.append(sp.csr_matrix((g["data"][s:e], g["indices"][s:e], ip - s), shape=(b - a, var.shape[0])))
    return sp.vstack(parts).tocsr()


ksel = pd.read_csv(CN / "k_selection.csv"); chosen = ksel[ksel.chosen.astype(str) == "True"].set_index("sample_id")
todo = sorted(chosen.index[chosen.k == 9])
recs = []
for i, sid in enumerate(todo):
    p_csv, p_k = R12 / f"{sid}_programs.csv", R12 / f"{sid}_kselection.csv"
    if p_k.exists() and not args.force:
        kd = pd.read_csv(p_k); recs.append(kd.iloc[0].to_dict()); continue
    with log.step(f"rank extension {sid} ({i+1}/{len(todo)})"):
        rows = np.flatnonzero(((obs.sample_id == sid) & sel).values)
        M, hv = MP.prep_matrix(read_rows(rows), var.index.to_numpy())
        res = MP.consensus_nmf(M, ks=KS_EXT)
        ok = [k for k, v in res.items() if v["full"] and v["n_clusters"] == k and v["sil"] >= MP.SIL_MIN]
        knew = max(ok) if ok else 9
        rec = dict(sample_id=sid, n_cells=len(rows), k_primary=9, k_new=knew,
                   **{f"sil_k{k}": round(res[k]["sil"], 4) for k in KS_EXT}, **{f"full_k{k}": bool(res[k]["full"]) for k in KS_EXT})
        if knew > 9:
            H = res[knew]["H"]; W = MP.usages_nnls(M, H); S = MP.gene_spectra_scores(M, W); prow = []
            for j in range(knew):
                order = np.argsort(-S[j])[:TOP_N]
                for r, g in enumerate(order):
                    prow.append(dict(program=f"{sid}|k{knew}|p{j+1}", sample_id=sid, k=knew, program_index=j + 1, rank=r + 1,
                                     gene=hv[g], spectra_score=float(S[j, g]), loading=float(H[j, g])))
            pd.DataFrame(prow).to_csv(p_csv, index=False)
        pd.DataFrame([rec]).to_csv(p_k, index=False); recs.append(rec)
        log.write(f"  {sid}: k_new={knew}; sil " + " ".join(f"{k}:{res[k]['sil']:.3f}{'' if res[k]['full'] else '(not full)'}" for k in KS_EXT))
ext = pd.DataFrame(recs); ext.to_csv(OUTS[0], index=False)

with log.step("re-cluster with extended ranks"):
    samp = pd.read_csv(CN / "samples_analysed.csv", index_col=0); samp = samp[samp.analysed].copy()
    samp["publication"] = samp.study.map(lambda x: PUBLICATION.get(x, x)); samp["assay_class"] = samp.assay.map(ASSAY_CLASS)
    replaced = set(ext.sample_id[ext.k_new > 9])
    prog = pd.concat([pd.read_csv((R12 if s in replaced else CN) / f"{s}_programs.csv") for s in samp.index], ignore_index=True)
    names = prog.program.drop_duplicates().to_list(); sets = [set(prog.gene[prog.program == n]) for n in names]
    meta = prog.drop_duplicates("program").set_index("program").loc[names, ["sample_id"]].join(samp[["publication", "assay_class"]], on="sample_id")
    D = 1 - MP.jaccard_matrix(sets); np.fill_diagonal(D, 0)
    lab = fcluster(linkage(squareform(D, checks=False), "average"), CUT, "distance")
    mps = [m for m in (np.flatnonzero(lab == c) for c in np.unique(lab)) if len(m) >= MIN_PROGRAMS and meta.sample_id.iloc[m].nunique() >= MIN_SAMPLES]
    out = []
    for m in mps:
        sub = meta.iloc[m]; mem = [names[j] for j in m]
        g = prog[prog.program.isin(mem)].groupby("gene").agg(n=("program", "nunique"), s=("spectra_score", "mean"))
        g = g[g.n >= MIN_GENE_FRAC * len(mem)].sort_values(["n", "s"], ascending=False).head(TOP_N)
        out.append(dict(genes=set(g.index), replicated=bool(sub.publication.nunique() >= 3 and set(sub.assay_class) >= {"cells", "nuclei"})))
    prim = pd.read_csv(T / "metaprograms.csv").set_index("metaprogram")
    pg = {m: set(g.gene) for m, g in pd.read_csv(T / "metaprogram_genes.csv").groupby("metaprogram")}
    rec = []
    for mp, G in pg.items():
        j = [len(G & o["genes"]) / len(G | o["genes"]) for o in out]; b = int(np.argmax(j))
        rec.append(dict(metaprogram=mp, family=prim.loc[mp, "family"], replicated_primary=bool(prim.loc[mp, "replicated"]),
                        best_jaccard=round(j[b], 3), recovered=j[b] >= RECOVER_J, replicated_extended=bool(j[b] >= RECOVER_J and out[b]["replicated"])))
    rec = pd.DataFrame(rec); rec.to_csv(OUTS[1], index=False)
    log.write(f"{len(replaced)} of {len(todo)} samples moved to k > 9 ({ext.k_new.value_counts().to_dict()}); re-clustered: "
              f"{len(names)} programmes, {len(out)} metaprograms, {sum(o['replicated'] for o in out)} replicated; primary recovered "
              f"{int(rec.recovered.sum())}/{len(rec)}\n" + rec.to_string(index=False))
    json.dump(dict(script="scripts/42_cnmf_rank_extension.py", ks_ext=KS_EXT, n_samples=len(todo), n_replaced=len(replaced),
                   inputs={str(p.relative_to(RESULTS.parent)): md5(p) for p in [CN / "k_selection.csv", T / "metaprograms.csv"]}),
              open(OUTS[2], "w"), indent=1)
