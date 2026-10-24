#!/usr/bin/env python
"""Build results/qc/spatial_sections.csv from unpacked Visium outputs (shape-only checks, no full matrix load)."""
import os, re, gzip, glob, csv, json, sys, time
import h5py, pandas as pd

import pathlib as _pl
ROOT = str(_pl.Path(__file__).resolve().parents[1])  # audit DQ-14: no hard-coded machine path
SP = f"{ROOT}/data/raw/spatial"
LOG = f"{ROOT}/logs/phase1_spatial_refs.log"

def mtx_shape(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("%"): continue
            g, s, nnz = map(int, line.split()); return g, s
def h5_shape(path):
    with h5py.File(path, "r") as f:
        s = f["matrix"]["shape"][:]; return int(s[0]), int(s[1])
def n_positions(path, in_tissue_only=True):
    df = pd.read_csv(path, header=None if "list" in os.path.basename(path) else 0)
    if "list" in os.path.basename(path): df.columns = ["barcode","in_tissue","row","col","pxr","pxc"]
    return int(df.in_tissue.sum()) if in_tissue_only else len(df)

def soft_meta(study, gse):
    meta = {}; cur = None
    for line in open(f"{SP}/{study}/{gse}_family_soft.txt", errors="replace"):
        line = line.rstrip("\r\n")
        if line.startswith("^SAMPLE = "): cur = line.split("= ")[1]; meta[cur] = {}
        elif cur and line.startswith("!Sample_title = "): meta[cur]["title"] = line.split("= ",1)[1]
        elif cur and line.startswith("!Sample_characteristics_ch1 = "):
            k, v = line.split("= ",1)[1].split(": ",1); meta[cur][k.strip()] = v.strip()
    return meta

rows = []
def add(study, section_id, gsm, region, n_spots, n_genes, img, pos, path, extra=""):
    rows.append(dict(study=study, section_id=section_id, gsm=gsm, region_or_diagnosis=region, n_spots=n_spots,
                     n_genes=n_genes, has_image=img, has_positions=pos, path=os.path.relpath(path, ROOT), notes=extra))

# huchede: flat GEO files per GSM
m = soft_meta("huchede", "GSE268577")
for d in sorted(glob.glob(f"{SP}/huchede/GSM*/")):
    gsm = os.path.basename(d.rstrip("/"))
    mtx = glob.glob(f"{d}/*_matrix.mtx.gz")[0]; g, s = mtx_shape(mtx)
    pos = glob.glob(f"{d}/*_tissue_positions.csv.gz")
    img = bool(glob.glob(f"{d}/*_tissue_hires_image.png.gz"))
    mm = m[gsm]; add("huchede", mm["title"].replace(" ","_"), gsm, f"{mm.get('disease','')}; {mm.get('genotype','')}; {mm.get('tissue','')}", s, g, img, bool(pos), d, "MTX")

# ren: spaceranger outs per GSM
m = soft_meta("ren", "GSE194329")
for d in sorted(glob.glob(f"{SP}/ren/GSM*/")):
    gsm = os.path.basename(d.rstrip("/"))
    h5 = glob.glob(f"{d}/*/filtered_feature_bc_matrix.h5")[0]; g, s = h5_shape(h5)
    pos = glob.glob(f"{d}/*/spatial/tissue_positions*.csv"); img = bool(glob.glob(f"{d}/*/spatial/tissue_hires_image.png"))
    mm = m[gsm]
    region = f"{mm.get('integrated diagnosis','')}; {mm.get('primary/ recurrent','')}; H3F3A={mm.get('h3f3a mutant','NA')}; HIST1H3B={mm.get('hist1h3b mutant','NA')}; IDH1={mm.get('idh1 mutant','NA')}"
    add("ren", mm["title"].split("_short")[0], gsm, region, s, g, img, bool(pos), d, "filtered_feature_bc_matrix.h5")

# brainstem: whatever GEO deposited per GSM (filled when unpacked)
if os.path.exists(f"{SP}/brainstem/GSE280990_family_soft.txt"):
    m = soft_meta("brainstem", "GSE280990")
    for d in sorted(glob.glob(f"{SP}/brainstem/GSM*/")):
        gsm = os.path.basename(d.rstrip("/"))
        h5s = glob.glob(f"{d}/**/filtered_feature_bc_matrix.h5", recursive=True) + glob.glob(f"{d}/*filtered_feature_bc_matrix.h5")
        mtxs = glob.glob(f"{d}/**/matrix.mtx*", recursive=True) + glob.glob(f"{d}/*matrix.mtx*")
        if h5s: g, s = h5_shape(h5s[0]); src = "h5"
        elif mtxs: g, s = mtx_shape(mtxs[0]); src = "MTX"
        else: g = s = None; src = "no matrix found"
        pos = glob.glob(f"{d}/**/tissue_positions*", recursive=True) + glob.glob(f"{d}/*tissue_positions*")
        img = bool(glob.glob(f"{d}/**/tissue_hires_image*", recursive=True) + glob.glob(f"{d}/*hires*"))
        mm = m[gsm]; add("brainstem", mm["title"].replace(" ","_"), gsm, f"{mm.get('tissue region','')}; DIPG-infiltrated brainstem (autopsy)", s, g, img, bool(pos), d, src)

# collot_st: Zenodo ST.tar.gz -> ST/ST/{Human,IUE}/<sample>/outs ; species from the matrix genome field
def h5_genome(path):
    with h5py.File(path, "r") as f:
        g = f["matrix"]["features"]["genome"][:]
        return sorted(set(x.decode() for x in g))
for d in sorted(glob.glob(f"{SP}/collot_st/ST/ST/*/*/outs")):
    sample = d.split("/")[-2]; group = d.split("/")[-3]
    h5 = f"{d}/filtered_feature_bc_matrix.h5"; g, s = h5_shape(h5); genome = h5_genome(h5)
    pos = glob.glob(f"{d}/spatial/tissue_positions*.csv"); img = bool(glob.glob(f"{d}/spatial/tissue_hires_image.png"))
    species = "human" if any(x.startswith("GRCh") for x in genome) else ("mouse" if any(("mm10" in x or "GRCm" in x) for x in genome) else "unknown")
    region = ("DMG, primary human; " + ("pons" if "pons" in sample else "tumour core")) if group == "Human" else "mouse IUE model, cortex tumour (EXCLUDE: non-human)"
    add("collot_st", sample, "", region, s, g, img, bool(pos), d, f"h5; genome={'/'.join(genome)}; species={species}; group={group}")

out = pd.DataFrame(rows)
out["species"] = out.notes.str.extract(r"species=(\w+)")[0].fillna("human")
os.makedirs(f"{ROOT}/results/qc", exist_ok=True)
out.to_csv(f"{ROOT}/results/qc/spatial_sections.csv", index=False)
print(out[["study","section_id","gsm","n_spots","n_genes","has_image","has_positions"]].to_string())
print("\nper study spots:", out.groupby("study").n_spots.agg(["count","sum"]).to_dict())
