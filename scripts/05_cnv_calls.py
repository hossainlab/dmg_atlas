#!/usr/bin/env python
"""Phase 3: malignant vs normal calls for one study (marker scoring, inferCNV, CNV clustering, consensus tiers).

Usage: python scripts/05_cnv_calls.py --study jessa_scrna [--samples S1 S2] [--resolution 1.0]
Inputs : data/interim/<study>.qc.h5ad (Phase 2), data/external/gencode/gencode.v32.annotation.gtf.gz
Outputs: data/interim/cnv/<study>/<sample>.cnv.npz        smoothed CNV matrix (float16) + window table + reference flags
         data/interim/cnv/<study>/<sample>.clusters.csv   per-CNV-cluster summary (z-scores, events, call)
         data/interim/cnv/<study>.calls.csv               per-cell table (scores, marker/CNV/consensus calls, author labels)
         data/interim/cnv/<study>.samples.csv             per-sample run summary (reference counts, genes, thresholds, wall, RSS)
Reference design: own immune marker cells + pooled study immune anchor + pooled, tumour-purged lineage cells from OTHER samples.
Audit corrections 2026-10-19 (docs/audit/cnv_audit.md): arm means / hallmark events from a per-chromosome-ARM smoothed matrix (CNV-02; cnv_score and clustering unchanged), 4q12 dropped from the hallmark
set (CNV-03), per-cell cnv_confidence event_backed|score_only and per-sample low_confidence_sample flag (CNV-09),
consensus_call_eventtrust sensitivity label (CNV-05). Thresholds unchanged.
         logs/phase3_cnv.log                              wall time + peak RSS per step
Never touches the Phase 2 files.
"""
import sys
import time
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from dmg_atlas import utils  # noqa: E402  (thread caps before numpy import)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import anndata as ad  # noqa: E402

from dmg_atlas import INTERIM, LOGS  # noqa: E402
from dmg_atlas import cnv as C  # noqa: E402

STUDIES = ["jessa_scrna", "jessa_multiome", "collot", "zang", "desisto", "scpca", "okada", "liu", "filbin", "labelle"]
HALLMARK_REGIONS = ["1q", "2", "7", "4q12", "13q", "14q", "10q", "17p"]   # arm_4q12 is exported as exploratory only (not a hallmark since audit CNV-03)


def pool_key(obs: pd.DataFrame) -> pd.Series:
    """Pooled-reference stratum: assay (+ author_prep for Liu fresh cells vs frozen nuclei; + platform for LaBelle)."""
    key = obs["assay"].astype(str)
    if "author_prep" in obs:
        key = key + "|" + obs["author_prep"].astype(str)
    if "labelle_platform" in obs:
        key = key + "|" + obs["labelle_platform"].astype(str)
    return key


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True, choices=STUDIES)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    study = args.study
    log = utils.StepLogger(LOGS / "phase3_cnv.log", context=study)
    t_all = time.time()
    outdir = INTERIM / "cnv" / study; outdir.mkdir(parents=True, exist_ok=True)
    path = INTERIM / f"{study}.qc.h5ad"
    genepos = C.load_gene_positions()

    import h5py
    from anndata.io import read_elem
    with h5py.File(path, "r") as f:
        obs_all = read_elem(f["obs"])
        # audit DQ-01: LaBelle's layer is now 'log2_tpm10' (log2(TPM/10+1)); cnv.normalise() handles both scales via counts_type='tpm'
        layer = next((L for L in ("tpm", "log2_tpm10") if L in f["layers"]), "counts")
    counts_type = "tpm" if layer in ("tpm", "log2_tpm10") else "umi"
    mask, high_mito = C.working_mask(obs_all, study)
    log.write(f"working set: {int(mask.sum())}/{len(obs_all)} cells (qc_pass {int(obs_all.qc_pass.sum())}; relaxed-mito nuclei samples add {int((mask & ~obs_all.qc_pass.astype(bool)).sum())}); layer={layer}")
    samples = sorted(obs_all.loc[mask, "sample_id"].unique())
    if args.samples:
        samples = [s for s in samples if s in set(args.samples)]
    pkey = pool_key(obs_all)

    # ---- pass 1: marker scores, marker calls, reference candidates
    scores_all, calls_all = [], []
    with log.step("pass1 marker scoring"):
        for sid in samples:
            rows = np.where((obs_all.sample_id == sid).values & mask.values)[0]
            if len(rows) < 20:
                log.write(f"{sid}: only {len(rows)} cells, skipped"); continue
            a = C.read_rows(path, rows, layer)
            C.normalise(a, counts_type)
            sc_ = C.score_markers(a, seed=args.seed)
            mc = C.marker_calls(sc_)
            sc_["sample_id"] = sid
            scores_all.append(sc_); calls_all.append(mc)
            del a
    scores = pd.concat(scores_all); mcalls = pd.concat(calls_all)
    ref_type_init = pd.Series(np.where(mcalls.is_reference, mcalls.marker_type, "none"), index=mcalls.index)
    log.write(f"pass1: {len(scores)} cells scored; reference candidates {int((ref_type_init != 'none').sum())} " + str(ref_type_init[ref_type_init != 'none'].value_counts().to_dict()))

    # ---- pooled immune anchor per stratum (<= POOLED_REF_PER_TYPE myeloid + T_NK candidates), read once
    pools_imm = {}
    rng = np.random.default_rng(args.seed)
    strata = pkey[scores.index]
    for key in strata.unique():
        cells = ref_type_init[(ref_type_init.isin(C.IMMUNE_TYPES)) & (strata == key).values]
        pick = []
        for t, grp in cells.groupby(cells.values):
            pick.extend(rng.choice(grp.index.values, min(len(grp), C.POOLED_REF_PER_TYPE), replace=False))
        if pick:
            pool = C.read_rows(path, obs_all.index.get_indexer(pick), layer); C.normalise(pool, counts_type)
            pool.obs["ref_type"] = ref_type_init[pool.obs_names].values; pool.obs["_pooled"] = True
            pools_imm[key] = pool
            log.write(f"pooled immune anchor [{key}]: {pool.n_obs} cells from {pool.obs.sample_id.nunique()} samples " + str(pool.obs.ref_type.value_counts().to_dict()))

    def with_pool(a, pool):
        if pool is None:
            return a, 0
        use = pool[~pool.obs_names.isin(a.obs_names)].copy()
        if use.n_obs == 0:
            return a, 0
        return ad.concat([a, use], join="inner", merge="first"), use.n_obs

    # ---- pass 1b: per-sample purge of lineage candidates (tumour-grade CNV clusters vs immune reference) -> pooled lineage refs
    survivors = {}
    with log.step("pass1b purge lineage candidates"):
        for sid in samples:
            rows = np.where((obs_all.sample_id == sid).values & mask.values)[0]
            if len(rows) < 20:
                continue
            a = C.read_rows(path, rows, layer); C.normalise(a, counts_type)
            a.obs["ref_type"] = ref_type_init[a.obs_names].values; a.obs["_pooled"] = False
            if not a.obs.ref_type.isin(C.LINEAGE_TYPES).any():
                continue
            a, _ = with_pool(a, pools_imm.get(pkey[a.obs_names[0]]))
            try:
                ok = C.purge_candidates(a, genepos, resolution=args.resolution, seed=args.seed)
            except ValueError as e:
                log.write(f"{sid}: purge skipped ({e}); lineage candidates of this sample not pooled"); continue
            keep = ok[ok].index
            survivors[sid] = ref_type_init[keep]
            n_c = int(ref_type_init[a.obs_names[~a.obs['_pooled'].values]].isin(C.LINEAGE_TYPES).sum())
            log.write(f"{sid}: lineage candidates {n_c} -> {len(keep)} survive tumour-grade purge " + str(survivors[sid].value_counts().to_dict()))
            del a
    pools_lin = {}
    for key in strata.unique():
        pick = []
        for sid, sv in survivors.items():
            sv = sv[(strata[sv.index] == key).values]
            for t, grp in sv.groupby(sv.values):
                pick.extend(rng.choice(grp.index.values, min(len(grp), C.POOL_PER_SAMPLE_PER_TYPE), replace=False))
        if not pick:
            continue
        pick = pd.Series(ref_type_init[pick])
        final = []
        for t, grp in pick.groupby(pick.values):
            final.extend(rng.choice(grp.index.values, min(len(grp), C.POOLED_REF_PER_TYPE), replace=False))
        pool = C.read_rows(path, obs_all.index.get_indexer(final), layer); C.normalise(pool, counts_type)
        pool.obs["ref_type"] = ref_type_init[pool.obs_names].values; pool.obs["_pooled"] = True
        pools_lin[key] = pool
        log.write(f"pooled lineage reference [{key}]: {pool.n_obs} cells from {pool.obs.sample_id.nunique()} samples " + str(pool.obs.ref_type.value_counts().to_dict()))

    # ---- pass 2: inferCNV, clustering, calls
    per_cell_out, per_sample = [], []
    for sid in samples:
        t0 = time.time()
        rows = np.where((obs_all.sample_id == sid).values & mask.values)[0]
        if len(rows) < 20:
            continue
        with log.step(f"cnv {sid} ({len(rows)} cells)"):
            a = C.read_rows(path, rows, layer); C.normalise(a, counts_type)
            a.obs["ref_type"] = ref_type_init[a.obs_names].values; a.obs["_pooled"] = False
            n_ref_own = int((a.obs.ref_type != "none").sum()); n_imm_own = int(a.obs.ref_type.isin(C.IMMUNE_TYPES).sum())
            key = pkey[a.obs_names[0]]
            a, n_pi = with_pool(a, pools_imm.get(key)); a, n_pl = with_pool(a, pools_lin.get(key))
            ref_source = f"own_immune[{n_imm_own}]+pooled_immune[{n_pi}]+pooled_lineage[{n_pl}]"
            rec = dict(sample_id=sid, study=study, n_cells=len(rows), n_ref_candidates_own=n_ref_own, n_ref_immune_own=n_imm_own, n_pooled_immune=n_pi, n_pooled_lineage=n_pl, ref_source=ref_source)
            cells = a.obs_names[~a.obs["_pooled"].values]
            try:
                pc, clusters, res, rec2 = C.call_sample(a, genepos, resolution=args.resolution, seed=args.seed, log=log)
            except ValueError as e:
                log.write(f"{sid}: {e} -> cnv_call unassigned")
                pc = pd.DataFrame({"cnv_cluster": "NA", "cnv_score": np.nan, "cnv_call": "unassigned"}, index=cells)
                rec.update(status="no_reference")
            else:
                for r in HALLMARK_REGIONS:
                    pc[f"arm_{r}"] = res["arm_means"][r].values if r in res["arm_means"] else np.nan
                clusters.to_csv(outdir / f"{sid}.clusters.csv", index=False)
                np.savez_compressed(outdir / f"{sid}.cnv.npz", X_cnv=res["X_cnv"].astype(np.float16), cell_id=np.array(a.obs_names, dtype=object),
                                    pooled=a.obs["_pooled"].values, ref_type=res["ref_final"].values.astype(str), cnv_cluster=pc.cnv_cluster.values.astype(str),
                                    window_chrom=res["windows"].chromosome.values.astype(str), window_pos=res["windows"].pos.values, window_arm=res["windows"].arm.values.astype(str),
                                    window_n_genes=res["windows"].n_genes.values.astype(int),
                                    chr_pos_keys=np.array(list(res["chr_pos"].keys())), chr_pos_vals=np.array(list(res["chr_pos"].values())),
                                    # per-ARM smoothed matrix behind arm_means / hallmark events (audit correction CNV-02)
                                    X_arm=res["X_arm"].astype(np.float16), arm_window_arm=res["windows_arm"].arm.values.astype(str),
                                    arm_window_pos=res["windows_arm"].pos.values, arm_window_n_genes=res["windows_arm"].n_genes.values.astype(int))
                rec.update({k: (str(v) if isinstance(v, dict) else v) for k, v in rec2.items()}); rec["status"] = "ok"
                pc = pc.loc[cells]
            df = pd.concat([scores.loc[cells].drop(columns="sample_id"), mcalls.loc[cells], pc], axis=1)
            df["author_label"], df["author_call"] = C.author_binary(obs_all.loc[cells], study)
            df = C.consensus(df)
            df.insert(0, "cell_id", df.index); df.insert(1, "sample_id", sid); df.insert(2, "study", study)
            df["high_mito_nuclei"] = high_mito[cells].values
            df["ref_source"] = ref_source
            per_cell_out.append(df)
            rec.update(n_malignant=int((df.consensus_call == "malignant").sum()), n_normal=int((df.consensus_call == "normal").sum()),
                       n_ambiguous=int((df.consensus_call == "ambiguous").sum()), cnv_malignant=int((df.cnv_call == "malignant").sum()),
                       wall_s=round(time.time() - t0, 1), peak_rss_gb=round(utils.peak_rss_gb(), 2))
            per_sample.append(rec)
            del a
    calls = pd.concat(per_cell_out)
    calls.to_csv(INTERIM / "cnv" / f"{study}.calls.csv", index=False)
    pd.DataFrame(per_sample).to_csv(INTERIM / "cnv" / f"{study}.samples.csv", index=False)
    log.write(f"DONE {study} | cells {len(calls)} | consensus " + str(calls.consensus_call.value_counts().to_dict()) +
              f" | total wall {(time.time() - t_all) / 60:.1f} min | peak RSS {utils.peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
