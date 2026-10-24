#!/usr/bin/env python
"""Phase 6 / step 14b (audit remediation 2026-10-20, findings SP-03, SP-17) — rebuild the Visium job inputs against
the NEW cell2location reference gene set (modal/inputs/c2l_reference.h5ad, metaprogram-family classes).

Difference from 14_cell2location.py prepare():
  * the gene set comes from the rebuilt reference;
  * genes are restricted PER COHORT to the intersection of the genes present in every section's own feature list,
    so that ad.concat(join='outer', fill_value=0) can no longer turn genes absent from a section's Space Ranger
    reference into observed zeros (SP-17: Medulla_A/B were quantified against a 17,943-gene reference and lacked
    77 of the union genes, affecting 4,670 spots);
  * results/qc/spatial_sections_used.csv gains n_ref_genes_present and n_cohort_genes.
Old inputs are kept as modal/inputs/visium_<cohort>.pre_audit.h5ad.
"""
from __future__ import annotations
import sys, json, glob, gzip, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc  # noqa: E402
from dmg_atlas.ingest import read_10x_mtx_dir  # noqa: E402
sys.path.insert(0, str(ROOT / "scripts"))
import importlib.util
spec = importlib.util.spec_from_file_location("s14", ROOT / "scripts/14_cell2location.py")
s14 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s14)

LOG = StepLogger(ROOT / "logs/fix_spatial.log", context="14b_prepare_visium_fix")
MODAL_IN = ROOT / "modal/inputs"


def main():
    sec = pd.read_csv(ROOT / "results/qc/spatial_sections.csv")
    sec = sec[sec.species == "human"].copy()
    ref_genes = set(ad.read_h5ad(MODAL_IN / "c2l_reference.h5ad", backed="r").var_names)
    LOG.write(f"{len(sec)} human sections; new reference gene set {len(ref_genes)}")
    used = []
    for cohort, rows in sec.groupby("study"):
        parts, present = [], []
        for _, r in rows.iterrows():
            d = ROOT / r.path
            a = s14.read_counts(d)
            pos = s14.read_positions(d); sf = s14.read_scalefactors(d)
            common = a.obs_names.intersection(pos.index)
            a = a[common].copy(); a.obs = a.obs.join(pos)
            n0 = a.n_obs
            a = a[a.obs.in_tissue == 1].copy()
            a.obs["section_id"] = r.section_id; a.obs["cohort"] = cohort
            a.obs["region_or_diagnosis"] = r.region_or_diagnosis
            a.obs["n_counts"] = np.asarray(a.X.sum(1)).ravel(); a.obs["n_genes"] = np.asarray((a.X > 0).sum(1)).ravel()
            a.obsm["spatial"] = a.obs[["pxl_col", "pxl_row"]].values.astype(float)
            g = a.var_names.intersection(sorted(ref_genes))
            present.append(set(map(str, g)))
            a = a[:, g].copy()
            a.X = sp.csr_matrix(a.X).astype(np.int32)
            a.obs_names = [f"{r.section_id}:{b}" for b in a.obs_names]
            a.uns["scalefactors"] = sf
            used.append(dict(study=cohort, section_id=r.section_id, region_or_diagnosis=r.region_or_diagnosis,
                             n_spots_matrix=n0, n_spots_in_tissue=a.n_obs, n_ref_genes=len(ref_genes),
                             n_ref_genes_present=a.n_vars, median_counts=float(np.median(a.obs.n_counts)),
                             median_genes=float(np.median(a.obs.n_genes)), hires_scalef=sf.get("tissue_hires_scalef"),
                             spot_diameter_fullres=sf.get("spot_diameter_fullres")))
            LOG.write(f"  {cohort}/{r.section_id}: {a.n_obs} spots, {a.n_vars} ref genes, "
                      f"median counts {np.median(a.obs.n_counts):.0f}")
            parts.append(a)
        shared = sorted(set.intersection(*present))
        union = sorted(set.union(*present))
        LOG.write(f"{cohort}: union {len(union)} genes, intersection across {len(parts)} sections {len(shared)} "
                  f"(dropped {len(union)-len(shared)} section-specific genes; SP-17)")
        vis = ad.concat([p[:, shared].copy() for p in parts], join="inner", merge="unique", uns_merge=None)
        vis.X = sp.csr_matrix(vis.X).astype(np.int32)
        vis.uns["scalefactors"] = {p.obs.section_id.iloc[0]: p.uns["scalefactors"] for p in parts}
        vis.uns["gene_selection"] = {"n_union": len(union), "n_shared": len(shared), "rule": "intersection across sections (SP-17)"}
        for c in ("section_id", "cohort", "region_or_diagnosis"):
            vis.obs[c] = vis.obs[c].astype(str)
        out = MODAL_IN / f"visium_{cohort}.h5ad"
        pre = MODAL_IN / f"visium_{cohort}.pre_audit.h5ad"
        if out.exists() and not pre.exists():
            shutil.move(str(out), str(pre))
        vis.write_h5ad(out, compression="gzip", compression_opts=4)
        LOG.write(f"wrote {out} ({out.stat().st_size/1e6:.1f} MB, {vis.shape}, nnz {vis.X.nnz})")
    u = pd.DataFrame(used)
    p = ROOT / "results/qc/spatial_sections_used.csv"
    if p.exists() and not (ROOT / "results/qc/spatial_sections_used.pre_audit.csv").exists():
        shutil.copy(p, ROOT / "results/qc/spatial_sections_used.pre_audit.csv")
    u.to_csv(p, index=False)
    LOG.write(f"done | {len(u)} sections | peak RSS {peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
