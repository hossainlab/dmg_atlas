# Phase 1 — spatial cohorts, normal references, bulk cohort: acquisition notes (2026-10-17)

Scripts: `scripts/dl_ranged.py` (parallel byte-range downloader; size + sha256 + optional md5 verification; appends `data/raw/MANIFEST_spatial_refs.csv`),
`scripts/dl_tail.py` (explicit byte ranges, used once for Zenodo), `scripts/02_qc_spatial_sections.py` (shape-only section QC -> `results/qc/spatial_sections.csv`),
`scripts/03_subset_braun_hindbrain.py` (chunked backed subset of the Braun h5ad). Wall times per step: `logs/phase1_spatial_refs.log`.
Manifest: 24 rows, 20,982,985,300 bytes (20.98 GB) in total.

## 1. Visium cohorts (`data/raw/spatial/<study>/`) — all 27 human sections open; 4 mouse sections flagged

| study | source | bytes | sections (human) | spots | notes |
|---|---|---|---|---|---|
| huchede | GSE268577_RAW.tar | 32,071,680 | 3 (GSM8294434–36, pDMG Sample-1/2/3) | 1,650 (577 / 856 / 217) | flat GEO files per GSM (MTX + positions + scalefactors + hires/lowres images), moved into per-GSM folders. All three: midline DMG, H3.3 K27M. Features 17,943 (gene-filtered by authors). Dataset table said 1,618 spots. |
| ren | GSE194329_RAW.tar | 2,611,107,840 | 11 (GSM5833528–38) | 27,233; DMG1–5 = 10,889 | per-sample `*_spaceranger_out.tar.gz` unpacked: filtered h5 + MTX, molecule_info.h5, spatial/ (tissue_positions_list.csv, scalefactors, images). GSM characteristics: DMG1–5 "Diffuse midline glioma, H3K27M-mutant", H3F3A mutant = Yes, HIST1H3B = No for all five (DMG1, DMG4 recurrent; DMG2, 3, 5 primary); GBM1 and GBM4 = secondary GBM (IDH1 mutant); GBM2, GBM3, GBM5_1 = GBM IDH-wt; GBM5_2 = peritumour tissue. 36,601 features. |
| collot_st | Zenodo 17922773 ST.tar.gz | 1,683,064,011 | 3 human (BT01_192_core_HE, VUMC_11_pons_HE, VUMC_17_pons_HE) + 4 mouse IUE (IUE_CORTEX_TUMOR_B_C1/B_D1/C_A1/C_B1) | human 5,296 (1,720 / 1,338 / 2,238); mouse 15,382 | full Space Ranger `outs/` per sample (filtered + raw h5 and MTX, molecule_info, spatial/). Species from the h5 genome field: GRCh38 vs custom_refdata-gex-mm10-2020-A. README confirms only newly generated human sections are included (Ren/Kordowski/Van den Broek not redistributed). Human list: `data/raw/spatial/collot_st/human_sections.csv`. No per-section clinical metadata in the deposit (sample names only; "pons" x2, "core" x1). scRNA.tar.gz NOT downloaded here (other track). |
| brainstem | GSE280990_RAW.tar | 419,450,880 | 10 (GSM8608415–24) | 34,939 | per GSM: filtered_feature_bc_matrix.h5, scalefactors, hires/lowres images, tissue_positions_list.csv, plus author `*_updated_LR_scores.csv.gz`. Regions from GSM titles: Midbrain (3,658), Pons A–D (4,919 / 4,235 / 4,992 / 3,127), Medulla A–C (1,958 / 2,712 / 3,171), Spinal Cord A–B (3,423 / 2,744) — ten sections, two of them spinal cord. Features 18,085 (17,943 for Medulla A/B) — author-filtered h5. Single autopsy patient. |

Images (`tissue_hires_image`) and spot positions are present for every section. Spot counts are rows of the filtered matrix (in-tissue spots).

## 2. Normal references (`data/external/`) — details in `data/external/REFERENCES.md`

* DOIs 10.1126/science.adf1226 (Braun) and 10.1126/science.add7046 (Siletti) verified on api.crossref.org.
* Braun 2023, collection 4d8fed08-2d6d-4692-b5ea-464f1d072077. cellxgene has no Pons/Hindbrain per-region dataset, so the full "All" dataset
  (13149914-ea03-4f01-8bf6-b793b667127b, version 52ff91e6-…, 11,932,827,832 bytes, 1,665,937 cells) was downloaded and subset locally:
  `braun2023/braun2023_hindbrain_pons_medulla.h5ad` = **286,973 cells x 58,226 genes** (Hindbrain 63,238; Pons 96,052 incl. metencephalon 33,488; Medulla 127,683),
  17 donors, 5.0–13.0 pcw. Cerebellum, whole-brain and head specimens excluded. Raw counts preserved (integer float32).
* Siletti 2023, collection 283d65eb-dd53-496d-adb7-7570c7caa443: all six "Dissection: Pons (Pn)" datasets (PnRF, PnEN, PN, DTg, PB, PnAN) =
  **204,168 nuclei x 58,232 genes**, 3 donors, 1,831,054,562 bytes. Each opens in backed mode; obs carries roi, dissection, supercluster_term, cluster_id, cell_type.

## 3. Bulk cohort (`data/external/openpedcan_v14/`)

* Bucket reached via the virtual-hosted hostname (generic s3.amazonaws.com is denylisted in the sandbox; access to
  d3b-openaccess-us-east-1-prd-pbta.s3.amazonaws.com was requested and granted). 133 keys listed (`bucket_listing_v14.xml`).
* Downloaded (md5 checked against md5sum.txt where listed): histologies.tsv (21.4 MB; 47,895 rows x 66 cols), gene-expression-rsem-tpm-collapsed.rds
  (285.4 MB; 60,325 genes x 4,124 biospecimens; gene symbols are current HGNC, e.g. H3-3A, H3C2), snv-consensus-plus-hotspots.maf.tsv.gz (2,077,192,289 bytes
  = 1.93 GiB, taken as "under ~2 GB"; GRCh38 MAF, Tumor_Sample_Barcode = BS_ ids), plus md5sum.txt, release-notes.md, independent-specimens.rnaseq.primary-plus-pre-release.tsv, snv-mutation-tmb-coding.tsv.
  Not downloaded: TPM TSV variants (only RDS collapsed TPM exists for PBTA in v14), isoform TPM (1.2 GB), CNV files, GTEx/TCGA matrices.
* `dmg_histologies.csv`: all rows with molecular_subtype starting "DMG, H3 K28" (509 "DMG, H3 K28, TP53" + 385 "DMG, H3 K28"):
  **238 Kids_First_Participant_IDs / 894 biospecimens**; by assay (counted in dmg_histologies.csv): WGS 281, RNA-Seq 225, Methylation 199, WXS 86, Targeted Sequencing 37, Phospho-Proteomics 29, Whole Cell Proteomics 29, miRNA-Seq 8.
  RNA-Seq: **179 participants / 225 biospecimens** (205 solid tissue, 20 derived cell lines; libraries poly-A stranded 90, stranded 62, exome_capture 44, poly-A 29);
  all 225 are columns of the TPM matrix. OS_status available for 813/894 rows (DECEASED 765, LIVING 48). Columns: participant/biospecimen/sample ids, cohort, sample_type,
  tumor_descriptor, composition, assay, RNA_library, primary_site, CNS_region, sex, age_at_diagnosis_days, age_at_event_days, OS_days/OS_status, EFS, resection,
  pathology/cancer_group/molecular_subtype(+methyl, DKFZ v12 class), tumor_fraction. The earlier estimate (144 / 189) is exceeded in v14.
  A broader filter (cancer_group DMG/DIPG or pathology text) gives 943 rows, including 49 "HGG, H3 wildtype / To be classified / IDH" rows that were NOT kept.

## Surprises / deviations
* Zenodo throttles: 4-way ranged download stalled at 0.07 MB/s (only the offset-0 connection received data); a single stream ran at 1–3 MB/s for ~800 MB then
  dropped to ~70 KB/s. Finished by combining the 1.25 GB single-stream prefix with 3 parallel tail ranges once the burst window reopened; md5 of the assembled
  file matches Zenodo (53ae0ade…). The orphaned single-stream process could not be killed from the sandbox and will eventually rewrite `collot_st/ST.tar.gz`
  with identical bytes (and append a duplicate manifest row); the verified copy used for unpacking is `ST.verified.tar.gz`.
* GEO and cellxgene served 10–45 MB/s with 3–4 range streams (not the ~1 MB/min feared).
* Remote backed reading of the Braun 11.9 GB h5ad over HTTP was impractically slow; the full file was downloaded instead (stated as a deviation from
  "subset remotely if >10 GB"). The subset step itself ran in 25k-cell blocks with <1 GB RSS.
* Stream budget: at most 4 concurrent range streams were active at any time (Collot's 3 idle Zenodo connections during the first attempt excepted).
* OpenPedCan consensus SNV is 1.93 GiB (2.08e9 bytes) — downloaded under the "~2 GB" rule; say so if a stricter cut-off was intended.
* Collot deposit lacks per-section clinical metadata; the dataset-selection table lists 10 patients / 11 samples for scRNA — the ST deposit has only 3 human sections.


## Audit corrections (2026-10-19; DQ-15)
`data/raw/MANIFEST_spatial_refs.csv` has 25 rows totalling 22,666,049,311 bytes (22.67 GB), not 24 rows / 20,982,985,300 bytes as stated in the header: the 25th row is the anticipated re-download of `collot_st/ST.tar.gz` (1,683,064,011 bytes, identical in size to `ST.verified.tar.gz`), appended at 08:22Z. The duplicate row is kept as a record of the re-download.
