#!/usr/bin/env python
"""Audit remediation MC-03: scripted reference gene sets with recorded source URL + sha256.

Background. results/tables/genesets/phase6_gene_sets.csv (used by scripts/16_metaprograms.py, 19_metaprogram_scores.py and
22_tme_scvi_clean.py) had no generating script. This script establishes the provenance of its 'Liu2022_DE_author_states' block
and produces published, checksummed alternatives:

  A. Liu2022_TableS2a_<state>   Liu et al. 2022 Nat Genet (doi 10.1038/s41588-022-01236-3) Supplementary Table S2a
                                 "Top 50 metaprogram genes" for S, G2M, AC-like, MES-like, OPC-like-1/2/3, OC-like
                                 (file 41588_2022_1236_MOESM3_ESM.xlsx, sheet table_S2a). Symbols are hg19-era and are mapped to
                                 GRCh38-2020-A symbols with dmg_atlas.genes.GeneMapper (alias rule of Phase 2); unmapped symbols
                                 are listed in the provenance table.
  B. Liu2022_DE_author_states_rederived_<state>
                                 The recipe that reproduces the pre-audit 'Liu2022_DE_author_states' sets (47-50 of 50 genes per
                                 set): data/interim/liu.qc.h5ad, cells whose author_annotation is one of the eight malignant states,
                                 genes with var.measured_in_all_samples, normalize_total(1e4)+log1p, scanpy rank_genes_groups
                                 (wilcoxon, reference='rest'), top 50 by score. These sets are therefore DERIVED FROM THE ATLAS'S OWN
                                 LIU DATA using the authors' per-cell labels - they are not the published tables (A).
  C. Filbin2018_*                Filbin et al. 2018 Science (doi 10.1126/science.aao4750) programme tables. Not retrievable by
                                 script on 2026-10-19: www.science.org answers with a Cloudflare JavaScript challenge (403) and PMC
                                 (PMC5949869, author manuscript NIHMS962385) serves its bin/ files behind a proof-of-work page; Europe
                                 PMC has no supplementary files for this non-OA record. If the user places
                                 data/raw/filbin/supp/NIHMS962385-supplement-S3.xlsx (PMC) or aao4750_tables*.xlsx (Science) in
                                 data/raw/filbin/supp/, parse_filbin() makes a best-effort parse (UNTESTED against the real layout) and the
                                 sets are added; otherwise the gap is recorded in the provenance table.

Outputs
  results/tables/genesets/reference_gene_sets_scripted.csv      gene_set, source, gene, gene_original, map_method
  results/tables/genesets/reference_gene_sets_provenance.csv    one row per source file: url, local path, bytes, sha256, sheet, status
  results/tables/genesets/reference_gene_sets_overlap.csv       per set in use (phase6_gene_sets.csv Liu2022_*, states.STATE_GENE_SETS,
                                                                 curated *_markers) vs the scripted sets: Jaccard, genes added/removed
  data/raw/MANIFEST_tumour.csv                                   rows appended for the downloaded supplementary files (if missing)

Phase 6 is NOT re-run here; the metaprogram track consumes the scripted sets.

Usage: python scripts/01d_reference_gene_sets.py [--fetch] [--no-rederive]
  --fetch        download the Liu MOESM3 xlsx if absent (static-content.springer.com); Filbin is attempted but expected to fail
  --no-rederive  skip block B (saves ~30 s and 3.5 GB RSS)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dmg_atlas import utils  # noqa: E402,F401  (thread caps)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from dmg_atlas import RAW, INTERIM, RESULTS, LOGS  # noqa: E402
from dmg_atlas.genes import GeneMapper  # noqa: E402
from dmg_atlas.states import STATE_GENE_SETS  # noqa: E402

GS_DIR = RESULTS / "tables" / "genesets"
MANIFEST = RAW / "MANIFEST_tumour.csv"

LIU_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01236-3/MediaObjects/41588_2022_1236_MOESM3_ESM.xlsx"
LIU_XLSX = RAW / "liu" / "supp" / "41588_2022_1236_MOESM3_ESM.xlsx"
LIU_SHEET = "table_S2a"
LIU_STATES = ["AC-like", "G2M", "MES-like", "OC-like", "OPC-like-1", "OPC-like-2", "OPC-like-3", "S"]

FILBIN_CANDIDATES = [
    ("https://www.science.org/doi/suppl/10.1126/science.aao4750/suppl_file/aao4750_tables1.xlsx", RAW / "filbin" / "supp" / "aao4750_tables1.xlsx"),
    ("https://pmc.ncbi.nlm.nih.gov/articles/instance/5949869/bin/NIHMS962385-supplement-S3.xlsx", RAW / "filbin" / "supp" / "NIHMS962385-supplement-S3.xlsx"),
]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> str:
    """Download url -> dest; returns a status string. Refuses HTML answers (challenge pages)."""
    import requests
    try:
        r = requests.get(url, timeout=180)
    except Exception as e:  # proxy / allowlist
        return f"blocked: {type(e).__name__}"
    ct = str(r.headers.get("Content-Type", ""))
    if r.status_code != 200:
        return f"HTTP {r.status_code}" + (" (challenge page)" if "html" in ct else "")
    if "html" in ct:
        return "HTTP 200 but HTML challenge page, not the file"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(r.content)
    return "downloaded"


def manifest_add(study: str, accession: str, url: str, path: Path, notes: str):
    man = pd.read_csv(MANIFEST)
    rel = str(path.relative_to(ROOT))
    if rel in set(man.local_path):
        return
    row = {"study": study, "accession": accession, "url": url, "local_path": rel, "bytes": path.stat().st_size, "sha256": sha256(path),
           "download_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(path.stat().st_mtime)), "notes": notes}
    pd.concat([man, pd.DataFrame([row])], ignore_index=True).to_csv(MANIFEST, index=False)


# --------------------------------------------------------------------------------------------------- A. Liu 2022 Table S2a
def parse_liu_s2a(mapper: GeneMapper) -> tuple[pd.DataFrame, dict]:
    raw = pd.read_excel(LIU_XLSX, sheet_name=LIU_SHEET, header=1)
    raw = raw[[c for c in raw.columns if not str(c).startswith("Unnamed")]]
    recs, unmapped = [], {}
    for state in raw.columns:
        genes = raw[state].dropna().astype(str).str.strip().tolist()
        canon, how = mapper.map_symbols(pd.Index(genes))
        seen = set()
        for g0, g, h in zip(genes, canon, how):
            if g is None:
                unmapped.setdefault(state, []).append(g0)
                continue
            if g in seen:
                continue
            seen.add(g)
            recs.append({"gene_set": f"Liu2022_TableS2a_{state}", "source": "Liu2022_NatGenet_TableS2a", "gene": g, "gene_original": g0, "map_method": h})
    df = pd.DataFrame(recs)
    assert set(raw.columns) == set(LIU_STATES), raw.columns.tolist()
    return df, unmapped


# --------------------------------------------------------------------------------------------------- B. re-derived in-house DE sets
def rederive_liu_de(log) -> pd.DataFrame:
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    a = ad.read_h5ad(INTERIM / "liu.qc.h5ad")
    mask = a.obs["author_annotation"].astype(str).isin(LIU_STATES).values
    b = a[mask, a.var["measured_in_all_samples"].values.astype(bool)].copy()
    del a
    b.X = sp.csr_matrix(b.layers["counts"]).astype(np.float32)
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    b.obs["state"] = b.obs["author_annotation"].astype(str)
    sc.tl.rank_genes_groups(b, "state", groups=LIU_STATES, reference="rest", method="wilcoxon")
    recs = []
    for s in LIU_STATES:
        top = sc.get.rank_genes_groups_df(b, group=s).head(50)
        recs += [{"gene_set": f"Liu2022_DE_author_states_rederived_{s}", "source": "Liu2022_DE_author_states_rederived", "gene": g,
                  "gene_original": g, "map_method": "atlas_symbol"} for g in top.names]
    log(f"re-derived Liu DE sets from {b.n_obs} malignant cells x {b.n_vars} genes (wilcoxon, top 50)")
    return pd.DataFrame(recs)


# --------------------------------------------------------------------------------------------------- C. Filbin 2018 (best effort)
def parse_filbin(path: Path, mapper: GeneMapper, log) -> pd.DataFrame:
    """Best-effort parser (UNTESTED: the file could not be downloaded on 2026-10-19). Looks through every sheet for columns whose
    header mentions a programme (OPC, AC, OC, cycle/G1S/G2M, MES) and takes the gene symbols below."""
    import re
    xl = pd.ExcelFile(path)
    recs = []
    for sh in xl.sheet_names:
        for header in (0, 1, 2):
            try:
                d = xl.parse(sh, header=header)
            except Exception:
                continue
            hits = [c for c in d.columns if re.search(r"OPC|AC[- _]?like|OC[- _]?like|astro|oligo|cycl|G1.?S|G2.?M|MES", str(c), re.I)]
            if len(hits) >= 2:
                for c in hits:
                    genes = d[c].dropna().astype(str).str.strip()
                    genes = genes[genes.str.match(r"^[A-Za-z0-9][A-Za-z0-9.\-_]*$")].tolist()
                    canon, how = mapper.map_symbols(pd.Index(genes))
                    name = re.sub(r"[^A-Za-z0-9]+", "_", str(c)).strip("_")
                    recs += [{"gene_set": f"Filbin2018_{name}", "source": f"Filbin2018_Science_{path.name}:{sh}", "gene": g, "gene_original": g0, "map_method": h}
                             for g0, g, h in zip(genes, canon, how) if g is not None]
                log(f"filbin: sheet {sh} header row {header}: {len(hits)} programme columns parsed")
                break
    return pd.DataFrame(recs)


# --------------------------------------------------------------------------------------------------- overlap report
def overlap_report(in_use: dict[str, set], scripted: dict[str, set]) -> pd.DataFrame:
    rows = []
    for k, cur in in_use.items():
        for s, new in scripted.items():
            inter = cur & new
            jac = len(inter) / len(cur | new) if (cur | new) else np.nan
            rows.append({"set_in_use": k, "n_in_use": len(cur), "scripted_set": s, "n_scripted": len(new), "n_shared": len(inter),
                         "jaccard": round(jac, 3), "genes_added_by_scripted": ";".join(sorted(new - cur)), "genes_removed_from_in_use": ";".join(sorted(cur - new))})
    df = pd.DataFrame(rows)
    df["best_match_for_set_in_use"] = df.groupby("set_in_use").jaccard.transform("max") == df.jaccard
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--no-rederive", action="store_true")
    args = ap.parse_args()
    log = utils.StepLogger(LOGS / "fix_data_repro.log", context="01d_reference_gene_sets").write
    GS_DIR.mkdir(parents=True, exist_ok=True)
    prov = []

    # --- Liu 2022 supplementary workbook
    status = "present"
    if not LIU_XLSX.exists():
        status = fetch(LIU_URL, LIU_XLSX) if args.fetch else "missing (run with --fetch)"
    if LIU_XLSX.exists():
        manifest_add("liu", "10.1038/s41588-022-01236-3", LIU_URL, LIU_XLSX, "Liu 2022 Nat Genet Supplementary Tables S1-S7 (xlsx); Table S2a = top-50 metaprogram genes used for reference gene sets (01d)")
    prov.append({"source": "Liu2022_NatGenet_TableS2a", "doi": "10.1038/s41588-022-01236-3", "url": LIU_URL, "local_path": str(LIU_XLSX.relative_to(ROOT)),
                 "bytes": LIU_XLSX.stat().st_size if LIU_XLSX.exists() else None, "sha256": sha256(LIU_XLSX) if LIU_XLSX.exists() else None,
                 "sheet": LIU_SHEET, "status": status})

    # --- Filbin 2018
    filbin_path = None
    for url, dest in FILBIN_CANDIDATES:
        st = "present" if dest.exists() else (fetch(url, dest) if args.fetch else "missing (not fetched; see docstring)")
        if dest.exists() and dest.stat().st_size > 0:
            filbin_path = filbin_path or dest
            manifest_add("filbin", "10.1126/science.aao4750", url, dest, "Filbin 2018 Science supplementary table (programme gene lists)")
        prov.append({"source": "Filbin2018_Science", "doi": "10.1126/science.aao4750", "url": url, "local_path": str(dest.relative_to(ROOT)),
                     "bytes": dest.stat().st_size if dest.exists() else None, "sha256": sha256(dest) if dest.exists() else None, "sheet": "", "status": st})

    mapper = GeneMapper()
    parts = []
    unmapped = {}
    if LIU_XLSX.exists():
        liu_df, unmapped = parse_liu_s2a(mapper)
        parts.append(liu_df)
        log(f"Liu Table S2a: {len(liu_df)} mapped genes in {liu_df.gene_set.nunique()} sets; unmapped per state: { {k: len(v) for k, v in unmapped.items()} }")
    if filbin_path is not None:
        parts.append(parse_filbin(filbin_path, mapper, log))
    if not args.no_rederive:
        parts.append(rederive_liu_de(log))
    scripted = pd.concat(parts, ignore_index=True)
    scripted.to_csv(GS_DIR / "reference_gene_sets_scripted.csv", index=False)
    prov_df = pd.DataFrame(prov)
    prov_df["unmapped_symbols"] = [json.dumps(unmapped) if r["source"].startswith("Liu2022") else "" for _, r in prov_df.iterrows()]
    prov_df.to_csv(GS_DIR / "reference_gene_sets_provenance.csv", index=False)

    # --- overlap with the sets currently in use
    cur = pd.read_csv(GS_DIR / "phase6_gene_sets.csv")
    in_use = {f"phase6:{k}": set(g.gene) for k, g in cur.groupby("gene_set") if k.startswith("Liu2022_") or k.endswith("_markers")}
    in_use.update({f"states.py:{k}": set(v) for k, v in STATE_GENE_SETS.items()})
    scr = {k: set(g.gene) for k, g in scripted.groupby("gene_set")}
    ov = overlap_report(in_use, scr)
    ov.to_csv(GS_DIR / "reference_gene_sets_overlap.csv", index=False)
    best = ov[ov.best_match_for_set_in_use].sort_values("set_in_use")
    log("best scripted match per set in use:\n" + best[["set_in_use", "n_in_use", "scripted_set", "n_shared", "jaccard"]].to_string(index=False))
    log(f"wrote {GS_DIR / 'reference_gene_sets_scripted.csv'} ({len(scripted)} rows, {scripted.gene_set.nunique()} sets), provenance, overlap")


if __name__ == "__main__":
    main()
