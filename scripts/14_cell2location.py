#!/usr/bin/env python
"""Phase 6 / step 14 — Visium deconvolution with cell2location (human sections only).

  prepare : read every human section listed in results/qc/spatial_sections.csv (Space Ranger filtered matrix h5 or MTX,
            tissue positions, scale factors), keep in-tissue spots, restrict genes to the cell2location reference gene set
            (modal/inputs/c2l_reference.h5ad var) and write one compressed h5ad per cohort: modal/inputs/visium_<cohort>.h5ad
            (X = int32 counts; obs section_id, cohort, region_or_diagnosis, array_row/col, pxl_row/col, n_counts, n_genes;
            obsm['spatial'] = (pxl_col, pxl_row) full-resolution pixels). Also writes results/qc/spatial_sections_used.csv.
  collect : copy harvested out/<section>_abundance.csv into data/processed/spatial/<cohort>/ and assemble
            results/tables/spatial_abundance_long.csv (one row per spot x class, q05 and mean abundance).
  package : (after 15_spatial_tests.py) one AnnData per cohort, data/processed/spatial/spatial_deconvolution_<cohort>.h5ad =
            visium_<cohort>.h5ad + obsm['q05_cell_abundance_w_sf'] / obsm['means_cell_abundance_w_sf'] + obs['niche']
            (DMG spots) + uns['c2l_run']. Skips existing files unless --force.
Remote part: modal/phase6_c2l.py deconvolve (Cell2location, N_cells_per_location 5, detection_alpha 20, A10G).
"""
from __future__ import annotations

import sys
import json
import gzip
import glob
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas.utils import StepLogger, peak_rss_gb  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import anndata as ad  # noqa: E402
import scanpy as sc  # noqa: E402
from dmg_atlas.ingest import read_10x_mtx_dir  # noqa: E402

LOG = StepLogger(ROOT / "logs/phase6_spatial.log", context="14_cell2location")
SECTIONS = ROOT / "results/qc/spatial_sections.csv"
MODAL_IN = ROOT / "modal/inputs"
TABLES = ROOT / "results/tables"
SPATIAL_OUT = ROOT / "data/processed/spatial"


def read_positions(d: Path) -> pd.DataFrame:
    cands = glob.glob(f"{d}/**/tissue_positions*.csv*", recursive=True) + glob.glob(f"{d}/*tissue_positions*.csv*")
    cands = sorted(set(cands))
    if not cands:
        raise FileNotFoundError(f"no tissue_positions in {d}")
    p = cands[0]
    cols = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    opener = gzip.open if p.endswith(".gz") else open
    with opener(p, "rt") as fh:
        first = fh.readline()
    has_header = first.lower().startswith("barcode")  # some "_list" files (GSE280990) carry a header anyway
    df = pd.read_csv(p, header=0 if has_header else None)
    df.columns = cols
    for c in cols[1:]:
        df[c] = pd.to_numeric(df[c])
    return df.set_index("barcode")


def read_scalefactors(d: Path) -> dict:
    cands = sorted(set(glob.glob(f"{d}/**/scalefactors_json.json*", recursive=True) + glob.glob(f"{d}/*scalefactors_json.json*")))
    if not cands:
        return {}
    p = cands[0]
    opener = gzip.open if p.endswith(".gz") else open
    with opener(p, "rt") as fh:
        return json.load(fh)


def read_counts(d: Path) -> ad.AnnData:
    h5 = sorted(set(glob.glob(f"{d}/**/filtered_feature_bc_matrix.h5", recursive=True) + glob.glob(f"{d}/*filtered_feature_bc_matrix.h5")))
    if h5:
        a = sc.read_10x_h5(h5[0])
        a.var_names_make_unique()
        return a
    return read_10x_mtx_dir(d)  # Huchedé: prefixed MTX triplet in the GSM folder


def prepare():
    sec = pd.read_csv(SECTIONS)
    sec = sec[sec.species == "human"].copy()
    ref_var = ad.read_h5ad(MODAL_IN / "c2l_reference.h5ad", backed="r").var_names
    ref_genes = set(ref_var)
    LOG.write(f"{len(sec)} human sections; reference gene set {len(ref_genes)}")
    used = []
    for cohort, rows in sec.groupby("study"):
        parts = []
        for _, r in rows.iterrows():
            d = ROOT / r.path
            a = read_counts(d)
            pos = read_positions(d)
            sf = read_scalefactors(d)
            common = a.obs_names.intersection(pos.index)
            a = a[common].copy()
            a.obs = a.obs.join(pos)
            n0 = a.n_obs
            a = a[a.obs.in_tissue == 1].copy()
            a.obs["section_id"] = r.section_id; a.obs["cohort"] = cohort; a.obs["region_or_diagnosis"] = r.region_or_diagnosis
            a.obs["n_counts"] = np.asarray(a.X.sum(1)).ravel(); a.obs["n_genes"] = np.asarray((a.X > 0).sum(1)).ravel()
            a.obsm["spatial"] = a.obs[["pxl_col", "pxl_row"]].values.astype(float)
            genes_present = a.var_names.intersection(sorted(ref_genes))
            a = a[:, genes_present].copy()
            a.X = sp.csr_matrix(a.X).astype(np.int32)
            a.obs_names = [f"{r.section_id}:{b}" for b in a.obs_names]
            a.uns["scalefactors"] = sf
            used.append(dict(study=cohort, section_id=r.section_id, region_or_diagnosis=r.region_or_diagnosis, n_spots_matrix=n0,
                             n_spots_in_tissue=a.n_obs, n_genes_total=len(ref_genes), n_ref_genes_present=a.n_vars,
                             median_counts=float(np.median(a.obs.n_counts)), median_genes=float(np.median(a.obs.n_genes)),
                             hires_scalef=sf.get("tissue_hires_scalef"), spot_diameter_fullres=sf.get("spot_diameter_fullres")))
            LOG.write(f"  {cohort}/{r.section_id}: {a.n_obs} spots, {a.n_vars} ref genes, median counts {np.median(a.obs.n_counts):.0f}")
            parts.append(a)
        # union of genes across sections of the cohort (sections of one cohort share a reference; zero-fill missing)
        vis = ad.concat(parts, join="outer", fill_value=0, merge="unique", uns_merge=None)
        vis.X = sp.csr_matrix(vis.X).astype(np.int32)
        vis.uns["scalefactors"] = {p.obs.section_id.iloc[0]: p.uns["scalefactors"] for p in parts}
        for c in ("section_id", "cohort", "region_or_diagnosis"):
            vis.obs[c] = vis.obs[c].astype(str)
        out = MODAL_IN / f"visium_{cohort}.h5ad"
        vis.write_h5ad(out, compression="gzip", compression_opts=4)
        LOG.write(f"wrote {out} ({out.stat().st_size/1e6:.1f} MB, {vis.shape}, nnz {vis.X.nnz})")
    pd.DataFrame(used).to_csv(ROOT / "results/qc/spatial_sections_used.csv", index=False)
    LOG.write(f"done | peak RSS {peak_rss_gb():.2f} GB")


def collect(harvest: dict[str, str]):
    """harvest = {cohort: harvested out/ directory}"""
    sec = pd.read_csv(SECTIONS); sec = sec[sec.species == "human"]
    long = []
    runs = {}
    for cohort, hdir in harvest.items():
        hdir = Path(hdir)
        dest = SPATIAL_OUT / cohort; dest.mkdir(parents=True, exist_ok=True)
        for f in sorted(hdir.glob("*_abundance.csv")):
            df = pd.read_csv(f, index_col=0)
            (dest / f.name).write_bytes(f.read_bytes())
            q = [c for c in df.columns if c.startswith("q05_")]
            m = [c for c in df.columns if c.startswith("mean_")]
            a = df[q].copy(); a.columns = [c[4:] for c in q]
            b = df[m].copy(); b.columns = [c[5:] for c in m]
            a = a.stack().rename("q05").reset_index(); a.columns = ["spot", "class", "q05"]
            b = b.stack().rename("mean").reset_index(); b.columns = ["spot", "class", "mean"]
            l = a.merge(b, on=["spot", "class"])
            l.insert(0, "section_id", df["section_id"].iloc[0]); l.insert(0, "cohort", cohort)
            l = l.merge(df[["array_row", "array_col"]].rename_axis("spot").reset_index(), on="spot")
            long.append(l)
        rj = hdir / f"{cohort}_run.json"
        if rj.exists():
            runs[cohort] = json.load(open(rj))
    long = pd.concat(long, ignore_index=True)
    long = long.merge(sec[["section_id", "region_or_diagnosis"]], on="section_id", how="left")
    long.to_csv(TABLES / "spatial_abundance_long.csv", index=False)
    json.dump(runs, open(TABLES / "c2l_deconvolution_runs.json", "w"), indent=1)
    LOG.write(f"spatial_abundance_long.csv: {long.shape}, sections {long.section_id.nunique()}, spots {long.spot.nunique()}")


def package(force: bool = False):
    """Per cohort: modal/inputs/visium_<cohort>.h5ad + collected abundances (+ niche labels from 15_spatial_tests.py)
    -> data/processed/spatial/spatial_deconvolution_<cohort>.h5ad (cell2location naming: obsm['q05_cell_abundance_w_sf'],
    obsm['means_cell_abundance_w_sf']; uns['c2l_run'], uns['c2l_classes'])."""
    if pd.__version__ >= "2.1":  # pandas 3 infers Arrow-backed strings, which anndata 0.12 cannot write
        pd.set_option("future.infer_string", False)
    runs = json.load(open(TABLES / "c2l_deconvolution_runs.json"))
    niche_p = SPATIAL_OUT / "spatial_niche_labels.csv.gz"
    niche = pd.read_csv(niche_p, index_col=0)["niche"] if niche_p.exists() else None
    for cohort in sorted(runs):
        out = SPATIAL_OUT / f"spatial_deconvolution_{cohort}.h5ad"
        if out.exists() and not force:
            LOG.write(f"skip {out.name} (exists)"); continue
        a = ad.read_h5ad(MODAL_IN / f"visium_{cohort}.h5ad")
        ab = pd.concat([pd.read_csv(f, index_col=0) for f in sorted((SPATIAL_OUT / cohort).glob("*_abundance.csv"))])
        missing = a.obs_names.difference(ab.index)
        assert not len(missing), f"{cohort}: {len(missing)} spots without abundances"
        ab = ab.loc[a.obs_names]
        classes = [str(c[4:]) for c in ab.columns if c.startswith("q05_")]
        a.obs_names = pd.Index([str(x) for x in a.obs_names], dtype=object)  # anndata cannot write Arrow-backed strings
        for k, pre in (("q05_cell_abundance_w_sf", "q05_"), ("means_cell_abundance_w_sf", "mean_")):
            a.obsm[k] = pd.DataFrame(ab[[pre + c for c in classes]].to_numpy(np.float32), index=a.obs_names,
                                     columns=pd.Index(classes, dtype=object))
        if niche is not None:
            v = niche.reindex(list(a.obs_names))
            a.obs["niche"] = pd.Categorical([str(int(x)) if pd.notna(x) else "NA" for x in v])
        for df in (a.obs, a.var):  # Arrow-backed strings (pandas string inference) -> object, which anndata can write
            for c in df.columns:
                col = df[c]
                if isinstance(col.dtype, pd.CategoricalDtype):
                    df[c] = pd.Categorical(np.asarray(col.astype(str), dtype=object), categories=pd.Index([str(x) for x in col.cat.categories], dtype=object))
                elif pd.api.types.is_string_dtype(col.dtype):
                    df[c] = np.asarray(col.astype(str), dtype=object)
        a.var_names = pd.Index([str(x) for x in a.var_names], dtype=object)
        a.uns["c2l_run"] = {k: (v if not isinstance(v, list) else np.array(v, dtype=str)) for k, v in runs[cohort].items()}
        a.uns["c2l_classes"] = np.array(classes, dtype=str)
        a.write_h5ad(out, compression="gzip")
        LOG.write(f"wrote {out.name}: {a.shape}, {len(classes)} classes, niche labels for {int((a.obs.get('niche', pd.Series(dtype=str)) != 'NA').sum())} spots")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["prepare", "collect", "package"])
    ap.add_argument("--harvest", nargs="*", default=[], help="cohort=dir pairs")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if a.mode == "prepare":
        prepare()
    elif a.mode == "collect":
        collect(dict(x.split("=", 1) for x in a.harvest))
    else:
        package(a.force)
