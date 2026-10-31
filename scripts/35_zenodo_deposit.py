#!/usr/bin/env python
"""Phase 9 / step 35 — assemble the Zenodo deposit of atlas v1: data and models only (folder + zip, with checksums).

Code, figures, result tables and documentation stay in the GitHub repository, which the deposit links to (README.md and
.zenodo.json related_identifiers). Builds deposit/dmg_atlas_v1/ from existing outputs (nothing upstream is modified) and zips it
to deposit/dmg_atlas_v1_zenodo.zip:
  data/      dmg_atlas_v1.h5ad (release object), query_mapping.h5ad, malignant_calls.csv.gz, sample_manifest.csv (local
             file_paths column dropped)
  models/    scvi_reference/ + scanvi_reference/ (models/scvi_seed1, models/scanvi): model.pt, run.json, var_names.txt
  metadata/  data dictionary: obs_columns.csv (generated from the release h5ad), label_ontology.csv
  README.md  data card (generated here; counts read from the release h5ad; links to the repository for everything else)
  LICENSE.md, .zenodo.json (Zenodo deposition metadata), FILES.tsv (path, bytes, md5, description), MD5SUMS
Large files are hard-linked, not copied. The zip stores already-compressed files (h5ad, gz, pt) without recompression.
Creator: Md. Jubayer Hossain (ORCID 0000-0001-9323-4997); no affiliation recorded.
Run: nice -n 10 python scripts/35_zenodo_deposit.py [--force] [--no-zip]
Log: logs/phase9_deposit.log
"""
import sys, os, json, gzip, shutil, hashlib, argparse, subprocess, zipfile, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import pandas as pd, h5py, torch  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); ap.add_argument("--no-zip", action="store_true")
args = ap.parse_args()
log = U.StepLogger(LOGS / "phase9_deposit.log", context="35_zenodo_deposit")
ROOT = RESULTS.parent; TABLES = RESULTS / "tables"; VERSION = "1.0.0"
DEP = ROOT / "deposit"; OUT = DEP / "dmg_atlas_v1"; ZIP = DEP / "dmg_atlas_v1_zenodo.zip"
RELEASE_MD5 = "4a39f172e4db0a2b470de7b6ac7b81e3"   # results/tables/release_manifest.csv
if OUT.exists() and not args.force:
    log.write(f"{OUT} exists; skip (use --force)"); sys.exit(0)

SOURCES = [  # (study, role, accession / record, URL, publication DOI)
    ("Okada 2026", "reference (snRNA)", "UCSC Cell Browser brain-dmg-multiomics; raw EGAS50000001447 (controlled)",
     "https://ega-archive.org/studies/EGAS50000001447", "10.1093/neuonc/noag119"),
    ("Jessa 2022 scRNA/snRNA", "reference", "GSE210567", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210567", "10.1038/s41588-022-01205-w"),
    ("Jessa 2022 multiome", "reference", "GSE210566", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210566", "10.1038/s41588-022-01205-w"),
    ("Collot 2026", "reference + Visium", "Zenodo 17922773", "https://doi.org/10.5281/zenodo.17922773", "10.1016/j.ccell.2025.12.020"),
    ("Zang 2024", "reference", "GSE245674 (GSM7846769)", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245674", "10.1186/s12883-024-03558-7"),
    ("DeSisto 2024", "reference", "GSE231859; ScPCA SCPCP000001", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231859", "10.1093/neuonc/noad207"),
    ("Liu 2022", "query", "GSE184357", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184357", "10.1038/s41588-022-01236-3"),
    ("LaBelle 2025", "query", "GSE227983, GSE271936", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227983", "10.1016/j.xcrm.2025.102095"),
    ("Filbin 2018", "query", "GSE102130", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130", "10.1126/science.aao4750"),
    ("Braun 2023", "normal query (fetal hindbrain)", "CELLxGENE collection 4d8fed08-2d6d-4692-b5ea-464f1d072077",
     "https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077", "10.1126/science.adf1226"),
    ("Siletti 2023", "normal query (adult pons)", "CELLxGENE collection 283d65eb-dd53-496d-adb7-7570c7caa443",
     "https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443", "10.1126/science.add7046"),
]
GH = "https://github.com/hossainlab/dmg_atlas"
GH_DOC = GH + "/blob/main/"


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def place(src, rel, desc, rows, link=False):
    d = OUT / rel; d.parent.mkdir(parents=True, exist_ok=True)
    os.link(src, d) if link else shutil.copy2(src, d)
    rows.append(dict(path=rel, bytes=d.stat().st_size, description=desc))


with log.step("read release obs"):
    assert md5(PROCESSED / "dmg_atlas_v1.h5ad") == RELEASE_MD5, "dmg_atlas_v1.h5ad md5 differs from release_manifest.csv"
    with h5py.File(PROCESSED / "dmg_atlas_v1.h5ad", "r") as f:
        full = read_elem(f["obs"])
        obs = full[["study", "sample_id", "donor_id", "suspension_type"]].copy()
        n_obs_cols = full.shape[1]
        cols = pd.DataFrame([dict(column=c, dtype=str(full[c].dtype), n_unique=int(full[c].nunique()), n_missing=int(full[c].isna().sum()),
                                  examples="; ".join(map(str, full[c].dropna().unique()[:5]))[:200]) for c in full.columns])
        del full
        n_genes = f["var"]["_index"].shape[0] if "_index" in f["var"] else len(read_elem(f["var"]))
    g = obs.groupby("study", observed=True)
    per_study = pd.DataFrame({"cells": g.size(), "samples": g.sample_id.nunique(), "donors": g.donor_id.nunique(),
                              "suspension": g.suspension_type.agg(lambda s: ", ".join(sorted(map(str, s.unique()))))}).sort_values("cells", ascending=False)
    with h5py.File(PROCESSED / "query_mapping.h5ad", "r") as f:
        q = read_elem(f["obs"]["query"]).value_counts()
    sha = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "src", "scripts", "modal"], capture_output=True, text=True).stdout.strip()
    log.write(f"{len(obs)} cells, {obs.sample_id.nunique()} samples, {obs.donor_id.nunique()} donors, {n_genes} genes; queries {q.to_dict()}; code HEAD {sha}"
              + (f"; WARNING uncommitted code changes:\n{dirty}" if dirty else ""))

with log.step("assemble folder"):
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True); rows = []
    place(PROCESSED / "dmg_atlas_v1.h5ad", "data/dmg_atlas_v1.h5ad", "release object: raw counts in X, CELLxGENE 7.1.0 obs fields, embeddings", rows, link=True)
    place(PROCESSED / "query_mapping.h5ad", "data/query_mapping.h5ad", "query cohorts projected onto the reference (scANVI latent, predicted labels, uncertainty)", rows, link=True)
    with open(PROCESSED / "malignant_calls.csv", "rb") as fi, gzip.open(OUT / "data/malignant_calls.csv.gz", "wb", compresslevel=6) as fo:
        shutil.copyfileobj(fi, fo)
    rows.append(dict(path="data/malignant_calls.csv.gz", bytes=(OUT / "data/malignant_calls.csv.gz").stat().st_size,
                     description="per-cell three-tier malignant calls (author / own CNV / consensus)"))
    sm = pd.read_csv(PROCESSED / "sample_manifest.csv").drop(columns=["file_paths"])
    sm.to_csv(OUT / "data/sample_manifest.csv", index=False)
    rows.append(dict(path="data/sample_manifest.csv", bytes=(OUT / "data/sample_manifest.csv").stat().st_size, description="per-sample metadata and source accessions"))
    for src, dst in [("models/scvi_seed1", "models/scvi_reference"), ("models/scanvi", "models/scanvi_reference")]:
        place(ROOT / src / "model.pt", f"{dst}/model.pt", f"scvi-tools 1.4.2 weights ({src}); load with scvi.model.{'SCANVI' if 'scanvi' in dst else 'SCVI'}.load", rows)
        place(ROOT / src / "run.json", f"{dst}/run.json", "training settings and metrics", rows)
        vn = torch.load(ROOT / src / "model.pt", map_location="cpu", weights_only=False)["var_names"]
        (OUT / dst / "var_names.txt").write_text("\n".join(map(str, vn)) + "\n")
        rows.append(dict(path=f"{dst}/var_names.txt", bytes=(OUT / dst / "var_names.txt").stat().st_size, description=f"the {len(vn)} input gene symbols, in model order"))
    (OUT / "metadata").mkdir()
    cols.to_csv(OUT / "metadata/obs_columns.csv", index=False)
    rows.append(dict(path="metadata/obs_columns.csv", bytes=(OUT / "metadata/obs_columns.csv").stat().st_size,
                     description=f"every obs column of data/dmg_atlas_v1.h5ad ({len(cols)}): dtype, unique values, missing, examples"))
    place(ROOT / "docs/schema/label_ontology.csv", "metadata/label_ontology.csv", "label vocabulary and ontology terms", rows)

with log.step("documentation"):
    if dirty:
        log.write("note: the README links to the repository at its default branch; commit and push code changes first")
    today = time.strftime("%Y-%m-%d")
    src_md = "\n".join(f"| {s} | {r} | [{a}]({u}) | [{d}](https://doi.org/{d}) |" for s, r, a, u, d in SOURCES)
    study_md = "\n".join(f"| {s} | {r.cells:,} | {r.samples} | {r.donors} | {r.suspension} |" for s, r in per_study.iterrows())
    q_md = ", ".join(f"{k} {v:,}" for k, v in q.items())
    readme = f"""# DMG Atlas v{VERSION} — single-cell reference atlas of H3K27-altered diffuse midline glioma

Cross-cohort single-cell / single-nucleus / multiome RNA-seq reference of human diffuse midline glioma (DMG), integrated
with scVI/scANVI, with query projections, malignant metaprograms and spatial deconvolution results.

- **Reference:** {len(obs):,} cells / nuclei · {obs.sample_id.nunique()} samples · {obs.donor_id.nunique()} donors · {obs.study.nunique()} studies · {n_genes:,} genes (GENCODE v48)
- **Queries (projection only):** {q_md}
- **Author:** Md. Jubayer Hossain ([ORCID 0000-0001-9323-4997](https://orcid.org/0000-0001-9323-4997))
- **Code, figures, result tables, full documentation:** {GH} (this deposit was built from commit `{sha}` by `scripts/35_zenodo_deposit.py` on {today})

## Contents

| Folder | What |
|---|---|
| `data/` | `dmg_atlas_v1.h5ad` (release object), `query_mapping.h5ad`, `malignant_calls.csv.gz`, `sample_manifest.csv` |
| `models/` | scVI and scANVI reference weights (scvi-tools 1.4.2) for mapping new data with scArches |
| `metadata/` | data dictionary: every obs column of the release object, and the label vocabulary |

Every file with size, md5 and description: `FILES.tsv`; checksums for `md5sum -c`: `MD5SUMS`.

## Reference cohorts

| Study | Cells | Samples | Donors | Suspension |
|---|---:|---:|---:|---|
{study_md}

## Quick start

```python
import anndata as ad, scanpy as sc
adata = ad.read_h5ad("data/dmg_atlas_v1.h5ad")          # X = raw UMI counts (CSR)
adata.obs[["study", "donor_id", "cell_type_atlas", "consensus_call", "mp_family"]].head()
sc.pl.embedding(adata, "X_umap", color="cell_type_atlas")
```

- `X` holds raw counts; normalise yourself. Gene IDs are Ensembl (GENCODE v48); symbols in `var['gene_symbol']`.
- Embeddings: `X_scvi`, `X_scanvi`, `X_umap`; `X_scvi_malignant` (malignant-only scVI, NaN outside the malignant subset).
- Key labels: `cell_type_atlas` / `cell_type_ontology_term_id` (cell type), `consensus_call` (malignant / normal / ambiguous),
  `scanvi_label` + `scanvi_maxprob` (classifier), `mp_family` (malignant metaprogram family), `tme_subcluster`.
- All {n_obs_cols} obs fields of this file, with type and examples: `metadata/obs_columns.csv`; vocabulary:
  `metadata/label_ontology.csv`. Schema and CELLxGENE field mapping: [atlas_schema.md]({GH_DOC}docs/schema/atlas_schema.md),
  [release_v1.md]({GH_DOC}docs/release/release_v1.md) (written before the deposit clean-up, which dropped the `*_pre_audit` and
  sample-level `author_*` columns it still lists).
- Map new data with scArches: build an AnnData whose `var_names` are the 3,000 gene **symbols** in
  `models/scanvi_reference/var_names.txt` (= `var['gene_symbol']` of the `highly_variable` genes), raw counts in
  `layers['counts']`, and obs columns `sample_id`, `study`, `assay`, `pct_mito` and `lineage_label = "Unknown"`; then
  `scvi.model.SCANVI.prepare_query_anndata(q, "models/scanvi_reference")`,
  `m = scvi.model.SCANVI.load_query_data(q, "models/scanvi_reference")`, `m.train(max_epochs=100, plan_kwargs={{"weight_decay": 0.0}})`,
  `m.predict()`. See the limitation below before trusting normal-lineage calls.

## Key limitations (full list: [limitations.md]({GH_DOC}docs/release/limitations.md))

- On a study it has not seen, the reference classifier transfers malignant and myeloid labels but often mislabels rare normal
  lineages as malignant (balanced accuracy 0.40-0.99 per held-out study; [phase9_quality_notes.md]({GH_DOC}docs/phases/phase9_quality_notes.md)). Confirm normal /
  malignant calls on new data with CNV or markers.
- Assay and cohort are confounded (the whole-cell arm is 6 scRNA samples from 3 publications); only 2 of 15 metaprograms
  (cycling, hypoxia) replicate under the pre-registered rule.
- `is_primary_data` is False for every cell: all data are re-processed from the public deposits below.

## Source data and terms

This atlas re-processes data from the studies below. Cite them when you use the atlas; their original data terms continue to apply.

| Study | Role | Data | Publication |
|---|---|---|---|
{src_md}

Okada 2026 data were obtained from the UCSC Cell Browser (dataset `brain-dmg-multiomics`); the raw reads are under controlled
access at EGA.

## License

Data and models: CC BY 4.0 (see `LICENSE.md`). The pipeline code in the repository is MIT-licensed.

## Citation

Cite this record by the DOI shown on its Zenodo page, and cite the source studies above.
"""
    (OUT / "README.md").write_text(readme)
    (OUT / "LICENSE.md").write_text(
        "# License\n\nThe data, models and metadata in this deposit are released under the Creative Commons "
        "Attribution 4.0 International license (CC BY 4.0): https://creativecommons.org/licenses/by/4.0/\n\n"
        f"The pipeline code ({GH}) is released under the MIT license.\n\n"
        "The atlas is derived from the public data sets listed in README.md. Attribution to those studies is required, and their "
        "original data terms continue to apply to the underlying data.\n")
    zen = {"metadata": {
        "upload_type": "dataset", "title": f"DMG Atlas v{VERSION}: a cross-cohort single-cell reference atlas of H3K27-altered diffuse midline glioma",
        "version": VERSION, "publication_date": today, "access_right": "open", "license": "cc-by-4.0",
        "creators": [{"name": "Hossain, Md. Jubayer", "orcid": "0000-0001-9323-4997"}],
        "description": (f"<p>Cross-cohort reference atlas of human H3K27-altered diffuse midline glioma: {len(obs):,} cells / nuclei from "
                        f"{obs.sample_id.nunique()} samples, {obs.donor_id.nunique()} donors and {obs.study.nunique()} public studies, integrated "
                        "with scVI/scANVI. Includes the release h5ad (raw counts, CELLxGENE 7.1.0 fields), query projections of three tumour "
                        "and two normal-brain cohorts, malignant calls, and scVI/scANVI reference model weights for scArches mapping. "
                        f"Code, figures, result tables and documentation: <a href=\"{GH}\">{GH}</a>. See README.md.</p>"),
        "keywords": ["diffuse midline glioma", "DIPG", "H3K27M", "single-cell RNA-seq", "single-nucleus RNA-seq", "reference atlas", "scVI", "scANVI", "cell2location"],
        "related_identifiers": [{"identifier": "https://github.com/hossainlab/dmg_atlas", "relation": "isSupplementedBy", "resource_type": "software"}]
                               + [{"identifier": d, "relation": "isDerivedFrom", "scheme": "doi", "resource_type": "publication-article"}
                                  for d in dict.fromkeys(s[4] for s in SOURCES)]
                               + [{"identifier": "10.5281/zenodo.17922773", "relation": "isDerivedFrom", "scheme": "doi", "resource_type": "dataset"}]}}
    (OUT / ".zenodo.json").write_text(json.dumps(zen, indent=1) + "\n")
    for f_, desc in [("README.md", "data card"), ("LICENSE.md", "licenses"), (".zenodo.json", "Zenodo deposition metadata")]:
        rows.append(dict(path=f_, bytes=(OUT / f_).stat().st_size, description=desc))

with log.step("checksums"):
    ft = pd.DataFrame(rows)
    ft["md5"] = [md5(OUT / p) for p in ft.path]
    assert ft.set_index("path").loc["data/dmg_atlas_v1.h5ad", "md5"] == RELEASE_MD5
    ft = ft[["path", "bytes", "md5", "description"]].sort_values("path")
    ft.to_csv(OUT / "FILES.tsv", sep="\t", index=False)
    (OUT / "MD5SUMS").write_text("".join(f"{m}  {p}\n" for p, m in zip(ft.path, ft.md5)))
    on_disk = {str(p.relative_to(OUT)) for p in OUT.rglob("*") if p.is_file()} - {"FILES.tsv", "MD5SUMS"}
    assert on_disk == set(ft.path), f"FILES.tsv out of sync: {on_disk ^ set(ft.path)}"
    log.write(f"{len(ft)} files, {ft.bytes.sum() / 1e9:.2f} GB")

if not args.no_zip:
    with log.step("zip"):
        STORE = {".h5ad", ".gz", ".pt"}
        tmp = ZIP.with_suffix(".zip.part")
        with zipfile.ZipFile(tmp, "w", allowZip64=True) as z:
            for p in sorted(OUT.rglob("*")):
                if p.is_file():
                    z.write(p, f"dmg_atlas_v1/{p.relative_to(OUT)}", compress_type=zipfile.ZIP_STORED if p.suffix in STORE else zipfile.ZIP_DEFLATED)
        with zipfile.ZipFile(tmp) as z:
            bad = z.testzip()
        assert bad is None, f"zip CRC check failed at {bad}"
        tmp.replace(ZIP)
        (DEP / "dmg_atlas_v1_zenodo.zip.md5").write_text(f"{md5(ZIP)}  {ZIP.name}\n")
        log.write(f"{ZIP} {ZIP.stat().st_size / 1e9:.2f} GB")
