<div align="center">

# DMG Atlas

**A cross-cohort single-cell reference atlas of H3K27-altered diffuse midline glioma**

[![Python](https://img.shields.io/badge/python-3.11-3776AB?logo=python&logoColor=white)](env/versions.csv)
[![scvi-tools](https://img.shields.io/badge/scvi--tools-1.4.2-5A2D82)](https://scvi-tools.org)
[![scanpy](https://img.shields.io/badge/scanpy-1.11.5-1F6FB2)](https://scanpy.readthedocs.io)
[![cell2location](https://img.shields.io/badge/cell2location-0.1.5-2E8B57)](https://cell2location.readthedocs.io)
[![CELLxGENE schema](https://img.shields.io/badge/CELLxGENE%20schema-7.1.0-0A7BBB)](docs/release/release_validation.md)
[![GPU: Modal](https://img.shields.io/badge/GPU-Modal%20A10G-7FEE64)](modal/)
<br>
[![Cells](https://img.shields.io/badge/cells-311%2C618-informational)](#at-a-glance)
[![Patients](https://img.shields.io/badge/patients-52-informational)](#cohorts)
[![Studies](https://img.shields.io/badge/studies-6%20%2B%205%20queries-informational)](#cohorts)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23623335.svg)](https://doi.org/10.5281/zenodo.23623335)
[![Release](https://img.shields.io/badge/release-v1-success)](docs/release/release_v1.md)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)
[![Last commit](https://img.shields.io/github/last-commit/hossainlab/dmg_atlas)](https://github.com/hossainlab/dmg_atlas/commits/main)

</div>

Public single-cell, single-nucleus and multiome RNA-seq cohorts of human DMG are integrated into a **scVI/scANVI
reference**. Smart-seq2/TPM tumour cohorts and normal brain references are projected onto it as queries, Visium
sections are deconvolved with **cell2location**, and malignant cell states are characterised by per-sample
**consensus NMF**. Human data only.

## Contents

[At a glance](#at-a-glance) · [Cohorts](#cohorts) · [Main results](#main-results) · [Quality checks](#quality-checks) ·
[Quick start](#quick-start) · [Repository layout](#repository-layout) · [Data release](#data-release) ·
[Limitations](#limitations) · [Documentation](#documentation) · [License](#license)

## At a glance

| | |
|---|---|
| **Reference** | 311,618 cells / nuclei · 57 samples · 52 patients · 6 studies |
| **Assays** | 229,043 snRNA · 46,288 multiome RNA · 36,287 scRNA |
| **Queries** (projection only) | Filbin 2018, LaBelle 2025, Liu 2022 (tumour) · Braun 2023 fetal hindbrain, Siletti 2023 adult pons (normal) |
| **Spatial** | 27 Visium sections from 4 cohorts |
| **Release** | `dmg_atlas_v1.h5ad` · 18,865 GENCODE v48 genes · CELLxGENE schema 7.1.0 layout |

## Cohorts

| Study | Role | Samples | Patients | Cells | Assay | Data URL |
|---|---|---:|---:|---:|---|---|
| Okada | reference | 22 | 22 | 156,936 | snRNA | [UCSC Cell Browser](https://cells-test.gi.ucsc.edu/brain-dmg-multiomics/) · raw: [EGAS50000001447](https://ega-archive.org/studies/EGAS50000001447) (controlled) |
| Jessa 2022 (scRNA/snRNA) | reference | 17 | 17 | 83,076 | snRNA, scRNA | [GSE210567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210567) |
| Jessa 2022 (multiome) | reference | 10 | 10 | 46,288 | multiome RNA | [GSE210566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210566) |
| Collot | reference | 5 | 4 | 17,562 | snRNA | [Zenodo 17922773](https://doi.org/10.5281/zenodo.17922773) |
| Zang | reference | 1 | 1 | 6,546 | scRNA | [GSE245674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245674) (GSM7846769) |
| De Sisto | reference | 2 | 2 | 1,210 | scRNA | [GSE231859](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231859) |
| Liu 2022 | query | 42 | 41 | 12,438 | Smart-seq2 | [GSE184357](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184357) |
| LaBelle 2025 | query | 9 | 9 | 8,247 | scRNA, Smart-seq2 | [GSE227983](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227983), [GSE271936](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271936) |
| Filbin 2018 | query | 9 | 9 | 3,037 | Smart-seq2 | [GSE102130](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130) |
| Braun 2023 hindbrain | normal query | 17 | 17 | 100,000 | scRNA | [CELLxGENE collection](https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077) |
| Siletti 2023 pons | normal query | 18 | 3 | 100,000 | snRNA | [CELLxGENE collection](https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443) |

Exact file URLs, sha256 checksums and download dates are written by the fetch scripts to `data/raw/MANIFEST_tumour.csv`
and `data/raw/MANIFEST_spatial_refs.csv` (local, not tracked).
Per-study summary (site, H3 variant, age group, malignant fractions): `results/tables/atlas_v1_summary.csv`.
Per-sample metadata: `data/processed/sample_manifest.csv`.

## Main results

| Finding | Result |
|---|---|
| **Metaprograms** | Per-sample cNMF on 50 samples (46 patients) gives 15 metaprograms; **2 replicate** under the pre-registered rule (≥ 3 studies, cells and nuclei): MP1 cycling and MP4 hypoxia |
| **Spatial** | MES-like × myeloid co-localisation supported (ρ = 0.238, q = 0.034); OPC-like × oligodendrocyte not (ρ = 0.178, q = 0.21) |
| **Malignant calling** | Consensus calls agree with author labels at a per-sample median of 0.898 over 97 labelled samples |
| **Reference classifier** | 0.975 accuracy, 0.917 balanced accuracy on held-out patients; on a held-out *study*, balanced accuracy falls to 0.40–0.99 ([quality checks](#quality-checks)) |

<details>
<summary><b>Figures</b> (generated from tables by <code>scripts/30_figures.py</code> and <code>scripts/31_supp_figures.py</code>)</summary>

| Figure | Content |
|---|---|
| 1 | Integrated reference: composition, scANVI lineage labels, scIB benchmark, hallmark CNV events |
| 2 | Zero-shot projection of query cohorts and label transfer |
| 3 | Malignant metaprograms and the pre-registered replication test |
| 4 | Microenvironment subclusters and sample-level composition models |
| 5 | Spatial organisation of malignant states (cell2location) |
| 6 | ATAC, bulk survival and foundation-model anchoring (*exploratory*) |
| 7 | Robustness: classifier on held-out studies, metaprogram recovery, composition label sensitivity (`scripts/36_fig_robustness.py`; Figure 6 of the preprint) |
| S1–S8 | Malignant calls, integration benchmark, mapping confidence, metaprogram robustness, assay masking, TME markers, spatial niches, survival models |

Output: `results/pub_figs/main`, `results/pub_figs/supp`. Captions with input md5s: `docs/figures/figures.md`,
`docs/figures/supp_figures.md`.

</details>

## Quality checks

Robustness analyses run after release (Phase 9; outputs in `results/tables/quality/`, details in
[`docs/phases/phase9_quality_notes.md`](docs/phases/phase9_quality_notes.md)).

| Check | Script | Result |
|---|---|---|
| Classifier on an unseen **publication** (scVI + scANVI refitted without it, scArches mapping) | `32_loso_classifier` | accuracy 0.88–0.99 per held-out study, but balanced accuracy 0.40–0.99 (pooled 0.51 vs 0.92 on held-out patients); malignant and myeloid transfer, rare normal lineages do not; with Okada unseen, 98 % of its neurons are called malignant |
| Metaprograms without each publication, and with ≤ 5 samples per publication | `33_metaprogram_lopo` | MP1 cycling recovered in every scenario; MP4 hypoxia needs Jessa 2022; MES-like and novel MP9 need both Jessa 2022 and Okada |
| Cell-type composition with consensus-only labels for the 29 % ambiguous cells | `34_composition_label_sensitivity` | all 4 significant effects keep sign and q < 0.05 |

## Quick start

```bash
# environment (local: Python 3.11.16, scanpy 1.11.5, scvi-tools 1.4.2, anndata 0.12.19)
conda create -n dmg_atlas --file env/local_conda_list.txt
conda activate dmg_atlas
pip install -r env/local_pip_freeze.txt

# run any step from the repo root (each script documents its usage in its docstring)
nice -n 10 python scripts/03_ingest_qc.py --study jessa_scrna
```

- Scripts add `src/` to `sys.path` themselves; no package install needed.
- Completed steps are skipped unless `--force` is given.
- Runs on a 24 GB laptop with 8 threads (`DMG_ATLAS_THREADS`). scVI/scANVI training and cell2location run on a Modal
  A10G GPU (`modal/`).
- Pinned versions, local and Modal: `env/versions.csv`, `env/environment.md`.

<details>
<summary><b>Full run order</b></summary>

Stage numbers follow the scripts and logs; they differ from `docs/design/build_plan.md`.

| Stage | Scripts |
|---|---|
| 1. Acquisition | `01_fetch_tumour_geo`, `01c_fetch_nongeo`, `01b_build_label_tables`, `01d_reference_gene_sets`, `02_build_sample_manifest`, `03_subset_braun_hindbrain`, `02_qc_spatial_sections` |
| 2. QC | `03_ingest_qc` (per study; batch: `bash scripts/run_phase2_all.sh`), `04_qc_report` |
| 3. CNV / malignant calls | `05_cnv_calls` (per study), `06_cnv_aggregate`, `06b_tier3_composition_fix` |
| 4. Reference integration | `07_assemble_atlas_input`, `07b_make_training_files`, `07c`, `07d`, Modal `phase4_train.py` + `phase4b_scanvi_retrain.py`, `08_integration_metrics`, `09_assemble_atlas_core`, `09b_refresh_atlas_core` |
| 5. Query mapping | `10_prepare_queries`, `11g_remap_queries_zeroshot`, `11b`–`11f`, `12_package_v1` |
| 6. States and spatial | `13_cnmf_per_sample`, `16_metaprograms`, `19_metaprogram_scores`, `20_assay_masked_metaprograms`, `13a`, `13b`, `14b_prepare_visium_fix` → Modal `phase6_c2l.py` → `14c_collect_fix`, `15b_spatial_tests_fix`, `15c_spatial_figures_fix`, `17`, `22`, `18_composition`, `23_package_phase6` |
| 7. Multimodal and bulk (exploratory) | `24_atac_multiome`, `25_openpedcan_survival`, `25b_survival_rankscore_audit`, `26_foundation_embeddings`, `26b_foundation_comparison_rebuild` |
| 8. Release and figures | `27_package_release`, `28_record_environment`, `29_validate_release`, `30_figures`, `31_supp_figures` |
| 9. Quality checks, deposit, preprint | Modal `phase9_loso_scanvi.py` → `32_loso_classifier`, `33_metaprogram_lopo`, `34_composition_label_sensitivity`, `35_zenodo_deposit`, `36_fig_robustness`, `37_graphical_abstract`; preprint: `bash preprint/scripts/build.sh` |

Superseded routes (`11_map_queries`, `14_cell2location`, `15_spatial_tests`, `21_brainstem_gradient_mp`) are kept for
provenance only; see `docs/release/outputs.md`.

</details>

## Repository layout

```
src/dmg_atlas/   pipeline package: ingest, genes, qc, cnv, states, metaprograms, figstyle, utils
scripts/         numbered CLI entry points 01–37
modal/           Modal GPU image and jobs (scVI/scANVI, cell2location, LOPO test)
configs/         qc.yaml (QC thresholds)
docs/            design/, phases/, schema/, release/, figures/, audit/, datasets/, archive/
env/             recorded software versions
results/         tables/, figures/, pub_figs/, qc/
notebooks/       exploratory only; no result depends on a notebook
preprint/        dmg_preprint.tex/.pdf, references.bib, generated tables/, build scripts/
data/ models/ logs/   not tracked; regenerated by the scripts
```

Git tracks code, configs, docs, environment records, figures and small tables. Data, models, logs, `*.h5ad` and
derived tables over ~20 MB are excluded (see `.gitignore`).

## Data release

| File | Description | md5 |
|---|---|---|
| `dmg_atlas_v1.h5ad` | Release object: raw counts in `X`, CELLxGENE 7.1.0 obs fields, metaprogram families, cleaned TME subclusters | `4a39f172e4db0a2b470de7b6ac7b81e3` |
| `atlas_core.h5ad` | Working reference: counts, QC, author/CNV/scANVI labels, `X_scvi`, `X_scanvi`, `X_umap` | `f14dd67c81ec90395c447c358f87c832` |
| `query_mapping.h5ad` | Query cohorts projected onto the reference, with predicted lineage and uncertainty | `9b585728cc8dcb81c0566362ced87eec` |
| `malignant_calls.csv` | Per-cell three-tier malignant calls (author / own CNV / consensus) | — |

Download from Zenodo: [doi:10.5281/zenodo.23623335](https://doi.org/10.5281/zenodo.23623335) (release h5ad, query mapping, malignant calls, reference models). Locally the files are in `data/processed/`. Field mapping: `docs/release/release_v1.md` · checksums:
`results/tables/release_manifest.csv` · schema: `docs/schema/atlas_schema.md` · labels: `docs/schema/label_ontology.csv`.

Local schema validation passes 119 checks with 0 failures (`docs/release/release_validation.md`); the official
`cellxgene-schema validate` has not been run yet.

## Limitations

- **Zero-shot query mapping.** Normal progenitors without a reference counterpart are predicted malignant
  (87.7 % of Braun fetal hindbrain cells, 29.0 % of Siletti OPCs).
- **Assay and cohort are confounded.** The cell arm is 6 scRNA samples from 3 publications; H3.1 rests on 8 Jessa
  samples (all pontine nuclei); autopsy vs resection is not estimable (all autopsies are Okada).
- **Limited spatial power.** The 10 brainstem Visium sections are from one patient (between-section I² 0.952–0.997).
- **Phase 7 is exploratory.** ATAC, bulk survival and foundation-model analyses were outside the approved plan.
- **Open audit findings.** 11 not closed and 12 partially closed (`docs/audit/CLOSURE_SUMMARY.md`).

Full list with audit references: [`docs/release/limitations.md`](docs/release/limitations.md).

## Documentation

| Folder | Content |
|---|---|
| [`docs/design/`](docs/design) | Frozen design decisions (incl. the pre-registered replication rule) and the phase plan |
| [`docs/phases/`](docs/phases) | What each phase did and found |
| [`docs/schema/`](docs/schema) | obs/obsm schema and label vocabulary |
| [`docs/release/`](docs/release) | Release field mapping, validation, output inventory, limitations |
| [`docs/figures/`](docs/figures) | Figure captions with input checksums |
| [`docs/audit/`](docs/audit) | Audit findings, fix records and the closure re-audit |
| [`docs/datasets/`](docs/datasets) | Dataset survey and selection |
| [`docs/archive/`](docs/archive) | Pre-audit snapshots of revised documents |

## License

Code: [MIT](LICENSE). Data derive from the public deposits of the cited studies and remain under their original terms.
