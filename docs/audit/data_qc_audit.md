# Independent audit — Phases 1–2 (acquisition, manifest, ingest, QC, gene harmonisation)

Auditor: independent track (data / QC). Date: 2026-10-19. Scope as assigned by the lead agent. Written only under `docs/audit/`; no code, data, table, figure or other doc was modified.
Findings table: `docs/audit/data_qc_findings.csv` (32 rows: 0 critical, 4 major, 16 minor, 12 confirmed-sound notes).

## 1. Overall verdict

**Sound with corrections.** The Phase 2 numbers that matter — 119 atlas samples, 420,381 barcodes, 357,061 `qc_pass`, every one of the 156 per-sample rows, the Okada count inversion, the Liu cell-ID recovery, the gene-harmonisation set sizes and the Jessa sanity checks — were re-derived from the obs tables and the underlying matrices and match the report exactly. No circularity was found in Phase 2 (thresholds are fixed a priori in `configs/qc.yaml`; author labels are used only descriptively). Four items are graded major: (i) the LaBelle layer stored and documented as TPM is in fact log2(TPM/10+1) and the note's "nCount_RNA = 1e6" claim is false (downstream phases rescued this with a heuristic, so conclusions stand); (ii) the label tables and three of the downloads that determine sample inclusion exist in no script; (iii) one Okada sample (1734_pons) enters the atlas with 8,486 nuclei of median 427 genes / 468 UMI and no sample-level flag; (iv) 16 Liu samples (5,375 frozen nuclei) were quantified on an annotation without MT genes, so the mitochondrial filter never acted on them, and the notes misdescribe the deposit (three annotations, not two). None of these invalidates a Phase 1–2 deliverable, but the first and fourth are materially misleading as documented and the third should be acted on before the reference model is considered final.

## 2. Scope covered

Scripts: scripts/01_fetch_tumour_geo.py, 02_build_sample_manifest.py (also re-executed to a temporary path), 02_qc_spatial_sections.py, 03_ingest_qc.py, 03_subset_braun_hindbrain.py, 04_qc_report.py, dl_ranged.py, dl_tail.py, run_phase2_all.sh.
Modules: src/dmg_atlas/__init__.py, ingest.py, genes.py (GeneMapper instantiated from the cached tables for alias spot checks), qc.py, utils.py. Config: configs/qc.yaml.
Data: data/raw/MANIFEST_tumour.csv (90 rows), MANIFEST_spatial_refs.csv (25 rows) with size check of all 115 files and sha256 of 9; data/processed/sample_manifest.csv (216 rows); data/interim/labels/ (jessa_sample_labels, desisto_sample_labels, liu_rawfile_to_sample, labelle_geo_metadata_sample_summary; others read via the manifest script); data/interim/gene_reference.csv, gene_universe.csv; obs + var of all 12 data/interim/*.qc.h5ad files (collot, desisto, desisto.outgroup, filbin, jessa_multiome, jessa_scrna, jessa_scrna.outgroup, labelle, liu, okada, scpca, zang) read with h5py; targeted count blocks of okada.qc.h5ad (586_pons), jessa_multiome (P-6253, P-1709), jessa_scrna (P-6251), scpca (SCPCS000637) and the full data array of liu.qc.h5ad (NaN scan); data/interim/collot_this_study_counts.h5ad (obs); raw deposits: Okada brain-dmg-snRNA-Homo-sapiens.h5ad (obs + raw.X block), Liu postFilter metadata and raw-count file row counts, LaBelle three TPM tables (first 6 columns each), Collot scRNA_metadata.csv.gz; data/external/REFERENCES.md, braun2023_hindbrain_pons_medulla.h5ad (obs, X attrs, first 2e6 values), six siletti2023 h5ads (obs).
Results: results/qc/qc_report.md, qc_per_sample.csv, qc_per_sample_<study>.csv (10), qc_per_study.csv, gene_harmonisation.csv, genes_shared_all_umi_studies.csv, sanity_jessa_malignant_qc.csv, spatial_sections.csv, spatial_sections_used.csv, all 10 qc_violin_*.png (viewed); results/tables/atlas_input_cells_per_sample.csv.
Docs: README.md, CLAUDE.md, docs/decisions.md, build_plan.md, phase1_tumour_notes.md, phase1_spatial_refs_notes.md, phase2_qc_notes.md; logs/phase1_tumour.log, phase1_spatial_refs.log, phase2_qc.log, phase2_driver.log.

## 3. Method of verification

1. Read every script and module in scope line by line (ingest loaders, GeneMapper, qc thresholds/scrublet, report generation, downloaders, manifest builder).
2. Re-executed `02_build_sample_manifest.py` to a temporary path and diffed it cell by cell against `data/processed/sample_manifest.csv`.
3. Loaded `obs` and `var` of all 12 QC'd h5ad files with h5py (no matrices) and re-derived every per-sample statistic in `qc_per_sample.csv`, the per-study totals, the Jessa sanity table and the high-mito-nuclei numbers; verified that `qc_pass` equals the threshold rule for every cell and that the stored 99.5th-percentile caps equal `np.quantile`.
4. Targeted matrix checks: Okada 586_pons recomputed from `raw.X` + `nCount_RNA` and compared with stored counts; harmonised counts for two multiome samples, one scRNA sample and one ScPCA sample to quantify the pre-harmonisation QC effect; full NaN/fractional scan of the Liu matrix; first six columns of each LaBelle table to establish the value scale.
5. Gene harmonisation: recomputed the shared-gene intersections from `var`, instantiated `GeneMapper` from the cached tables and spot-checked alias resolution; identified the Okada genes lost.
6. Data integrity: existence and byte size of all 115 manifest files, sha256 of 9 of them; Braun/Siletti obs and X attributes; Collot metadata; Liu raw-file row counts and mapping table; patient overlaps from the manifest.
7. Viewed all 10 QC figures.

Compute: local env `multiverse`, 8 threads, sparse/h5py only, peak RSS 3.0 GB, no downloads, no GPU.

## 4. Re-derived numbers (claimed vs recomputed)

| quantity | claimed (notes / report) | recomputed by auditor | status |
|---|---|---|---|
| Atlas samples / cells before QC / qc_pass / % | 119 / 420,381 / 357,061 / 84.9 % | 119 / 420,381 / 357,061 / 84.94 % | match |
| Out-group cells before / after | 127,949 / 99,360 | 127,949 / 99,360 (37 samples) | match |
| Per-sample QC rows re-derived from obs (before, after, fail counts, medians, 99.5th-pct cap) | 156 rows | 156/156 identical; qc_pass == threshold rule for every cell | match |
| Okada count recovery (586_pons, 23,731 nuclei) | max |x-round(x)| <= 1e-12; rowsum == nCount | 4.5e-13; max |rowsum-nCount| = 0.0; nnz == nFeature_RNA | match |
| Okada stored counts vs recomputation (586_pons) | exact | identical for 6 test genes and cell order; total 42,200,309 vs 42,204,394 (10 duplicate-suffixed symbols dropped, 0.01 %) | minor deviation (DQ-08) |
| Liu postFilter cells recovered | 8,984 / 9,248; missing 129/62/58/10/5 | 8,984 / 9,248; missing 129/62/58/10/5 = 264 (206 in samples with raw files) | match |
| Liu matrix | RSEM expected counts, 25-31 % fractional | 66,294,689 nnz; 28.3 % fractional; 0 NaN; 0 negative | match |
| Liu annotations in deposit | two (32,617 / 23,686) | three: 23,686 (16 samples), 32,617 (20), 55,765 (6) | MISMATCH (DQ-04) |
| Liu samples without MT genes | not stated (only Filbin/LaBelle) | 16 samples / 5,375 cells, pct_mito ≡ 0 | MISMATCH (DQ-04) |
| LaBelle value scale | TPM; P9 CPM with nCount_RNA = 1e6 | log2(TPM/10+1) in all 3 tables (max 8-9, col sums 3.7k-12.6k; P9 median 3,244) | MISMATCH (DQ-01) |
| Filbin value scale | TPM (n_counts 1e6) | median n_counts 1,000,000 in all 9 samples | match |
| Genes shared by all 7 UMI studies / all 10 | 18,882 / 14,119 | 18,882 (identical set) / 14,119 | match |
| Gene universe | 36,591 unique 2020-A symbols | 36,591 (from 36,601 features) | match |
| Jessa author-malignant total | 143,611 | 71,899 + 71,712 = 143,611; normal 33,301 | match |
| Jessa malignant passing our QC (atlas scRNA / multiome / out-group) | 96.8 % / 96.7 % / 87.3 % | 0.968 / 0.967 / 0.873 | match |
| Jessa barcodes without author metadata | 18,210 scRNA + 39,652 multiome | 18,210 (atlas) + 11,298 (out-group) scRNA; 39,652 multiome | match |
| pct_mito vs author percent.mito | equal to 2 decimals | max diff 0.32 pp (scRNA), 0.21 pp (out-group), 0.00 (multiome) | minor mismatch (DQ-12) |
| High-mito nucleus samples (median > 5 %) | 6 samples; P-6254 47.7 % / 10.6 % / 24 % pass | 6 samples; P-6254 47.69 % / 10.57 % / 24.2 %; P-1779 15.24 / 5.07 / 49.5 %; P-6253 11.15 / 3.07 / 67.0 %; P-1709 5.85 / 3.44 / 100 % | match |
| Scrublet per-study median thresholds | 0.19 ScPCA to 0.63 Zang | 0.19 ScPCA ... 0.38 Zang ... 0.51 Collot (max) | MISMATCH (DQ-06) |
| Samples with 0 auto-called doublets | 'most samples' | 36 / 96 (37.5 %) | mismatch (DQ-06) |
| Fixed-threshold doublet rate per study | 0.1-2.6 % (DeSisto 0.1, Zang 0.5, Okada 2.6) | 0.08 / 0.50 / 2.55 % | match |
| Manifest: integrated included samples / patients / cells | 59 / - / 334,439 | 59 / 54 patients / 334,439 (= author-QC cells for Jessa) | match (semantics DQ-17) |
| Patient overlap GSE210567/GSE210566 | 4 (P-1764, P-6253, P-6337, P-6640) | same 4 | match |
| Filbin ∩ Liu patients | 7 | 7 (BCH1126, BCH836, BCH869, MGH66, MUV1, MUV10, MUV5) | match |
| Manifest regenerated from 02_build_sample_manifest.py | identical | 239 cells differ: k27m_status for 10 G34 rows (WT vs unknown) + float formatting | MISMATCH (DQ-05) |
| Manifest species column | human | 'human' in 7/216 rows, blank in 209 | MISMATCH (DQ-07) |
| Raw-file manifests: files present / sizes / sha256 | 90 + 24 rows | 90 + 25 rows; 115/115 present with matching size; 9/9 re-hashed sha256 match | match (row count DQ-15) |
| Collot This_study cells / samples | 18,767 / 5 (of 146,567) | 18,767 / 5; Jessa2022 copies 127,800 / 24 | match |
| Braun hindbrain subset | 286,973 x 58,226; Hindbrain 63,238; Pons 96,052; Medulla 127,683; 17 donors | identical; integer counts | match |
| Siletti pons | 204,168 nuclei, 3 donors | 204,168; 3 donors; 58,232 genes | match |
| ScPCA | SCPCS000636 0 cells >= 300 genes; SCPCS000637 441 | 0 / 441 (median 362 genes); portal scpca_filter keeps 51 / 1,363 | match |
| Okada 1734_pons | not discussed | 15,228 nuclei, median 321 genes; 8,486 pass with median 427 genes / 468 UMI | NEW (DQ-03) |

## 5. Findings by severity

### Critical
_none_

### Major
**DQ-01 (correctness; src/dmg_atlas/ingest.py (iter_labelle); data/interim/labelle.qc.h5ad; docs/phase2_qc_notes.md §Decisions 7)** — LaBelle values stored as layers['tpm'] and described as TPM; note claims 'P9 is a 10x library distributed as per-cell CPM (nCount_RNA = 1e6)'

*Evidence.* All three GSE227983 tables are log2(TPM/10+1), not linear TPM: first 6 columns of tumor_ss2 have column sums 8,257-12,624, max 8.0-8.65, 100% fractional non-zeros; immune_ss2 sums 4,183-5,391; immune_tenx (P9 source) sums 3,661-6,177, max 7.6-9.0. Inverting (2^x-1)*10 gives 2.0e5-3.0e5 per cell (TPM over a gene-filtered table). qc_per_sample.csv median_counts for LaBelle SS2 samples = 4,823-9,527 and for P9 = 3,244, versus exactly 1e6 for every Filbin sample.

*Impact.* The Phase 2 deliverable mislabels a log-scale matrix as TPM and the note is factually wrong (nCount is ~3e3, not 1e6). Downstream phases detected the scale heuristically (cnv.py / 10_prepare_queries.py: rowmax<30 => log2), so atlas conclusions are not affected, but any consumer of labelle.qc.h5ad['tpm'] without that heuristic would be off by orders of magnitude.

*Recommendation.* Rename the layer (e.g. 'log2_tpm10') or invert to linear TPM at ingest; record the scale in uns and var; correct the note (P9 CPM/1e6 claim) and the qc_report 'median counts' column for TPM cohorts. (re-derived: yes)

**DQ-02 (reproducibility; data/interim/labels/*.csv; data/raw/MANIFEST_tumour.csv (collot, scpca, okada rows); data/raw/okada/snrna-seq-data/per_sample_summary.csv)** — README: 'every results/ file is produced by a script in scripts/'; Phase 1 notes describe label tables built from supplementary xlsx and downloads from Zenodo (Collot scRNA), ScPCA portal (API token) and UCSC Cell Browser (Okada)

*Evidence.* grep over scripts/ and notebooks/: no script writes data/interim/labels/ (02_build_sample_manifest.py only reads them); no script references zenodo.org, scpca-portal or cells-test.gi.ucsc.edu; 01_fetch_tumour_geo.py queues cover only GEO. MANIFEST_tumour.csv nevertheless contains 3 collot, 4 scpca and 13 okada rows and per_sample_summary.csv (okada) feeds the manifest. These steps exist only in the sub-agent's kernel.

*Impact.* Sample inclusion (H3 status, assay, site) for Jessa, Collot, DeSisto, LaBelle rests on label tables with no generating code; a third party cannot regenerate the manifest from the raw deposits.

*Recommendation.* Add scripts (or a documented notebook) that parse the supplementary xlsx/docx into data/interim/labels/*.csv and that fetch Collot/ScPCA/Okada, recording checksums; keep the label_source column. (re-derived: no)

**DQ-03 (validity; src/dmg_atlas/qc.py apply_thresholds; configs/qc.yaml (min_genes 300); results/qc/qc_per_sample.csv row okada__1734_pons; results/tables/atlas_input_cells_per_sample.csv)** — Per-cell thresholds flag 8,486 of 15,228 nuclei of okada__1734_pons as qc_pass; the sample enters the atlas with 8,466 nuclei (2.7% of the 311,618-cell reference)

*Evidence.* obs of okada.qc.h5ad: 1734_pons n_genes median 321 (IQR 249-452), n_counts median 351, 6,665 (43.8%) fail low_genes; passing nuclei have median 427 genes / 468 UMI, versus 971-2,805 genes in the other 21 Okada samples (next lowest 586_pons 971). No sample-level quality flag exists in qc_per_sample.csv or qc_report.md; the author QC (nFeature>200) is equally permissive.

*Impact.* A whole low-complexity sample is admitted with no warning; such nuclei typically form a depth-driven cluster and can distort scVI latent space, marker tests and composition models for the 'pons/ONC201' strata.

*Recommendation.* Add a sample-level QC gate (e.g. median genes of passing cells >= 500 or MAD-based per-sample adaptive min_genes) and flag 1734_pons explicitly; report the sample in qc_report.md. (re-derived: yes)

**DQ-04 (validity; src/dmg_atlas/ingest.py iter_liu; data/interim/liu.qc.h5ad; docs/phase2_qc_notes.md §Decisions 3, 6, 7)** — Notes: Liu has 'two annotations in the deposit' (32,617 / 23,686 genes) and only 'Filbin / LaBelle TPM tables contain no MT- genes'; Liu frozen nuclei are QC'd with the smartseq2 class (mito < 20%)

*Evidence.* gene_harmonisation.csv: Liu raw files come in three annotations: 23,686 genes (16 samples, 4,007 unmapped), 32,617 (20 samples, 1,749 unmapped) and 55,765 (6 samples, 22,216 unmapped). All 16 samples on the 23,686-gene annotation (5,375 cells, all author_prep = frozen_nuclei) have pct_mito identically 0 for every cell (max 0.0), i.e. the annotation lacks MT- genes and the mito filter is inoperative (visible as flat markers in qc_violin_liu.png). Among frozen nuclei that do carry MT genes, 138 of 680 passing nuclei (20.3%) have 5-20% mito and would fail the project's own nuclei rule.

*Impact.* QC is heterogeneous within the Liu cohort (mito-filtered vs unfiltered samples; nuclei judged by the cell threshold). Liu is projection-only, so reference models are unaffected, but Liu-based comparisons (fresh vs frozen, pediatric vs adult) carry a QC confound that the notes do not disclose.

*Recommendation.* State the three annotations and MT-gene absence in the notes; add a per-sample 'mt_genes_present' column to qc_per_sample.csv; apply the nuclei mito rule (or author_prep-aware rule) to frozen nuclei where MT genes exist. (re-derived: yes)

### Minor
**DQ-05 (reproducibility; data/processed/sample_manifest.csv vs scripts/02_build_sample_manifest.py)** — Manifest is the output of 02_build_sample_manifest.py

*Evidence.* Re-ran the script to a temporary path and diffed: 239 cells differ. k27m_status differs for the 10 Jessa H3.3-G34R/V samples (disk 'WT', script output 'unknown' because jessa_excl() sets WT only when h3_variant == 'WT'); 229 further cells are float-vs-int formatting (age_years, n_cells_deposited) indicating the file was re-saved via pandas after generation. File mtime 13:44 vs script mtime 13:33; no doc records the edit.

*Impact.* No effect on the atlas (G34 samples are excluded out-group) and 'WT' is a defensible K27 status for G34 tumours, but the on-disk manifest is not reproducible from the code.

*Recommendation.* Encode the G34 -> k27m_status='WT' rule in jessa_excl()/common and regenerate the manifest from the script. (re-derived: yes)

**DQ-06 (statistics; configs/qc.yaml (scrublet); docs/phase2_qc_notes.md §Decisions 2)** — 'per-study medians 0.19 ScPCA to 0.63 Zang'; 'called 0 doublets in most samples'; predicted_doublet_fixed (score > 0.25) 'gives 0.1-2.6% per study'

*Evidence.* qc_per_sample.csv (95 UMI samples with a threshold): per-study median auto-thresholds are ScPCA 0.19, DeSisto 0.28, multiome 0.38, Zang 0.38, scRNA 0.39, Okada 0.42, Collot 0.51 (no study at 0.63); 36/96 samples (37.5%) have 0 auto-called doublets (median per-sample rate 0.0003). Fixed-threshold pooled rates: DeSisto 0.08%, Zang 0.50%, Collot 1.06%, Okada 2.55% (claim verified). These are 3-10x below the expected 10x multiplet rate (~0.8% per 1,000 loaded cells; e.g. 586_pons has 23.7k nuclei). expected_doublet_rate is fixed at 0.06 regardless of cells loaded and the 0.25 cut-off is arbitrary.

*Impact.* Doublets are likely under-called in the atlas input (which excludes predicted_doublet_fixed); residual heterotypic doublets can appear as hybrid malignant/TME states. Two descriptive numbers in the note are wrong.

*Recommendation.* Correct the note; calibrate doublet calls (scDblFinder or per-sample expected rate from cell number) before state discovery; report per-sample expected vs called rates. (re-derived: yes)

**DQ-07 (correctness; scripts/02_build_sample_manifest.py add(); data/processed/sample_manifest.csv column species)** — Manifest carries a species column ('human' for all tumour rows)

*Evidence.* species is 'human' in only 7 of 216 rows (the placeholder rows that pass species explicitly); 209 rows are blank because add() fills every COLS key with '' before r.setdefault('species','human'), so setdefault never fires.

*Impact.* No analytical effect (all ingested data are human), but the species guarantee of decisions.md #1 is not recorded in the manifest.

*Recommendation.* Set species in the default dict (or use r['species'] = kw.get('species','human')). (re-derived: yes)

**DQ-08 (correctness; src/dmg_atlas/genes.py GeneMapper.map_symbols; data/interim/okada.qc.h5ad)** — Harmonisation rule: 'genes that collapse onto the same canonical symbol are summed'

*Evidence.* Okada var names arrive Seurat-deduplicated (TBCE.1, MATR3.1, HSPA14.1, ...); these 10 symbols are 'unmapped' and dropped (gene_harmonisation.csv okada unmapped_dropped = 10) whereas in Zang the same 10 duplicated 2020-A Ensembl IDs are summed (duplicate_symbols_collapsed = 10). Recomputed on 586_pons: stored total 42,200,309 vs raw 42,204,394 counts (4,085 lost, 0.01%; MATR3.1 alone 2,235).

*Impact.* Negligible count loss; inconsistent treatment of the same loci between two studies.

*Recommendation.* Strip the '.N' suffix before symbol mapping (or map Okada via the 2020-A feature order) so duplicates are summed as elsewhere. (re-derived: yes)

**DQ-09 (correctness; src/dmg_atlas/ingest.py iter_labelle (ad.concat join='inner'); docs/phase2_qc_notes.md §Decisions 7)** — 'LaBelle SS2 cells come from two tables (tumour, immune) with different gene lists (14,642 vs 19,357) -> zero-filled union'

*Evidence.* Within each of samples P1-P8 the tumour and immune SS2 tables are concatenated with join='inner' (gene_harmonisation.csv: n_input_genes 14,621 for P1-P8, 19,357 for P9); the immune cells therefore lose 4,736 measured genes before the across-sample outer join. var: 14,487 genes measured in all samples, 4,625 measured only in P9.

*Impact.* Immune SS2 cells of LaBelle lack genes that were measured; minor for a projection-only cohort but contrary to the documented rule.

*Recommendation.* Use join='outer' with fill_value=0 within samples too, and track measured genes per table. (re-derived: yes)

**DQ-10 (validity; scripts/03_ingest_qc.py lines 52-60 (QC before harmonise))** — QC metrics are comparable across studies

*Evidence.* n_genes/n_counts are computed on each deposit's native feature space before gene harmonisation (jessa_multiome 60,658 features; ScPCA 60,319). Recomputed on harmonised counts: P-6253 multiome median n_genes 508 -> 492, 0 author-QC cells cross the 300-gene threshold; P-1709: 0; jessa_scrna P-6251: 0; ScPCA SCPCS000637: 47 of 441 passing cells would fall below 300 genes (10.7%).

*Impact.* Negligible for the atlas studies; only ScPCA (later dropped) is materially affected.

*Recommendation.* Compute QC metrics after harmonisation, or on the 18,882 shared genes, so thresholds mean the same thing in every study. (re-derived: yes)

**DQ-11 (validity; src/dmg_atlas/qc.py apply_thresholds (max_genes_quantile 0.995))** — Per-sample max genes = 99.5th percentile of n_genes

*Evidence.* The quantile is taken over all deposited barcodes (including empty/low-count ones) and applied on top of author caps: in Okada (authors already enforced nFeature < 5000) 849 nuclei fail only high_genes; in Jessa 0 author-QC cells fail only high_genes (the top 0.5% are all non-author barcodes). The cut has no biological justification and removes a fixed 0.5% per sample irrespective of quality.

*Impact.* Small (<=0.5% per sample) but systematic removal of the highest-complexity cells (often cycling) in studies without an author cap.

*Recommendation.* Replace with a doublet-informed or MAD-based upper bound, or document that 0.5% is intentional. (re-derived: yes)

**DQ-12 (reporting; docs/phase2_qc_notes.md §Sanity (b))** — 'Our pct_mito equals the authors' percent.mito to 2 decimals'; 'The author Technology column confirms the nuclear protocol'

*Evidence.* max |pct_mito - author_percent.mito| = 0.32 percentage points over 86,525 jessa_scrna cells (0.21 pp out-group); exact (0.0) only for multiome. author_Technology is 'sc/snRNAseq' for every GSE210567 sample and 'scMultiome' for GSE210566 — it does not distinguish cells from nuclei; the nuclear label comes from Table S1 column scRNAseq.

*Impact.* Conclusions unchanged (differences are small and the threshold interpretation stands); two statements are imprecise.

*Recommendation.* Reword: 'within 0.32 pp' and cite Table S1 for the nuclear protocol. (re-derived: yes)

**DQ-13 (figure; results/qc/qc_violin_filbin.png, qc_violin_labelle.png, qc_violin_liu.png (scripts/04_qc_report.py lines 110-123))** — Per-sample violins of log10 genes, log10 counts, % mito for all cells before QC

*Evidence.* Viewed all 10 figures. Filbin: the counts panel plots floating-point noise around log10(1e6+1) with a '1e-7+6' offset label above the axes (reads as an in-panel title) and the %mito panel is empty (axis -0.04..0.04). LaBelle: 'log10(counts+1)' displays sums of log2 values (meaningless axis), %mito panel empty. Liu: 2,655 px wide with 6-pt rotated labels, 42 samples; samples with pct_mito == 0 appear as flat markers without explanation. No panel shows the applied thresholds. Otherwise: no in-panel titles, axes labelled, legible.

*Impact.* Figures are internal QC output; two panels per TPM figure carry no information and one carries a misleading axis label.

*Recommendation.* Skip counts/%mito panels for TPM cohorts (or label 'sum of log2(TPM/10+1)'), draw threshold lines, split Liu into two rows. (re-derived: no)

**DQ-14 (reproducibility; scripts/01_fetch_tumour_geo.py, 02_build_sample_manifest.py, 02_qc_spatial_sections.py, 03_subset_braun_hindbrain.py, dl_ranged.py)** — Scripts are CLI entry points runnable from the repo

*Evidence.* Five Phase 1 scripts hard-code ROOT = '/Users/jubayer/Projects/dmg_atlas' (01:16, 02_build:13, 02_qc_spatial:6, 03_subset:8, dl_ranged:10) whereas the package resolves paths from __file__. 02_build_sample_manifest.py and 03_subset_braun_hindbrain.py have no argparse/usage and run on import.

*Impact.* Scripts fail on any other machine without editing.

*Recommendation.* Derive ROOT from Path(__file__).resolve().parents[1] as in src/dmg_atlas/__init__.py. (re-derived: no)

**DQ-15 (reporting; docs/phase1_spatial_refs_notes.md header vs data/raw/MANIFEST_spatial_refs.csv)** — 'Manifest: 24 rows, 20,982,985,300 bytes (20.98 GB)'

*Evidence.* File has 25 rows totalling 22,666,049,311 bytes; the extra row is the anticipated re-download of collot_st/ST.tar.gz (1,683,064,011 bytes, identical size to ST.verified.tar.gz) appended at 08:22Z.

*Impact.* None; the note predicted the duplicate.

*Recommendation.* Update the header numbers or drop the duplicate row. (re-derived: yes)

**DQ-16 (data_integrity; scripts/02_qc_spatial_sections.py line 86; results/qc/spatial_sections.csv)** — Species column of spatial_sections.csv ('all 27 human sections')

*Evidence.* species is parsed from the notes string and .fillna('human'); only Collot ST sections carry a genome field check (GRCh38 vs mm10). Huchede (17,943 features), Ren (36,601) and brainstem (18,085) are labelled human by default without reading the genome/feature names.

*Impact.* Low risk (GEO metadata state Homo sapiens), but the human-only guarantee is asserted, not verified, for 24/27 sections.

*Recommendation.* Read the genome attribute of the h5 / first feature IDs (ENSG vs ENSMUSG) for every section. (re-derived: no)

**DQ-17 (reporting; data/processed/sample_manifest.csv column n_cells_deposited; docs/phase1_tumour_notes.md ('334,439 deposited cells/nuclei'))** — n_cells_deposited = cells in the deposit

*Evidence.* Column semantics differ by study: Jessa = author-QC cells (86,525 / 47,326 for included samples) while the GEO matrices hold 104,735 / 86,978 barcodes; DeSisto/Zang/Collot/Okada/Filbin/LaBelle = matrix barcodes (verified equal); ScPCA = portal estimate; Liu = postFilter metadata (9,190) vs 16,676 raw cells. The headline '334,439 deposited' is therefore a mixed count.

*Impact.* Explains the apparent paradox 'before QC exceeds deposited' noted in Phase 2; no analytical effect.

*Recommendation.* Add separate columns n_barcodes_deposited and n_author_qc_cells. (re-derived: yes)

**DQ-18 (correctness; src/dmg_atlas/qc.py compute_metrics (ribo prefixes))** — pct_ribo = share of counts in RPS/RPL genes

*Evidence.* Prefix matching captures 103 genes in the 2020-A universe of which 16 are not ribosomal proteins (RPS6KA1, RPS6KA2, RPS6KB2, RPL34-AS1, RPS10-NUDT3, RPL36AL, ...).

*Impact.* pct_ribo is descriptive only (not thresholded); bias is small.

*Recommendation.* Use the regex ^RP[SL]\d+[A-Z]?\d*$|^RPLP\d$|^RPSA$. (re-derived: yes)

**DQ-19 (data_integrity; scripts/02_build_sample_manifest.py (DeSisto, Zang, ScPCA tissue_timing='resection'))** — Manifest: 'No value is invented: anything not stated ... is unknown'

*Evidence.* tissue_timing is hard-coded 'resection' for all 19 DeSisto rows, Zang and both ScPCA rows although none of the cited label sources (ScPCA metadata, BMC Neurol table) has a biopsy/resection/autopsy field; Jessa, by contrast, is left 'unknown'.

*Impact.* Autopsy-vs-resection is a declared analysis axis (Okada autopsy vs everything else); the non-Okada 'resection' labels are assumptions, not records.

*Recommendation.* Mark as 'resection (assumed, surgical cohort)' or 'unknown'. (re-derived: no)

**DQ-20 (reporting; src/dmg_atlas/genes.py docstring vs docs/phase2_qc_notes.md §Decisions 3)** — Note: 'canonical symbol per Ensembl ID from the newest source'

*Evidence.* genes.py docstring states the opposite: zang GRCh38-2020-A (Ensembl 98) ranks first 'as a completeness rule, not a newest rule'; ScPCA Ensembl 104 (newer) ranks third. The code follows the docstring (ref priority 1 = zang).

*Impact.* Documentation inconsistency only; mapping spot checks (H3-3A->H3F3A, SEPT7->SEPTIN7, ATP5A1->ATP5F1A, RP11-34P13.7->AL627309.1) behave as the docstring says.

*Recommendation.* Align the note with the docstring. (re-derived: yes)

### Confirmed sound (notes)
**DQ-21 (reporting; docs/phase2_qc_notes.md §Sanity (b) decision vs results/tables/atlas_input_cells_per_sample.csv)** — Phase 2 recommendation: for high-mito nuclei samples 'treat author_qc_pass as the QC gate (or relax mito to 10 %)'

*Evidence.* Recomputed from obs: atlas input equals qc_pass & author_qc_pass & ~predicted_doublet_fixed for 54/57 samples; P-6254_S-8499 (631) corresponds to author & ~doublet & n_genes>=300 & pct_mito<20 (not 10%: that gives 356); P-1779_S-1781 (3,022) and P-6253_S-8498 multiome (2,132) equal author & ~doublet with no mito filter at all. These are Phase 4 choices, recorded here because they depart from the Phase 2 text.

*Impact.* Cross-phase consistency; the Phase 4 audit track should confirm these rules are documented there.

*Recommendation.* State the final per-sample gate in phase2/phase4 notes. (re-derived: yes)

**DQ-22 (correctness; results/qc/qc_report.md; qc_per_sample.csv; qc_per_study.csv)** — Headline: 119 atlas samples, 420,381 -> 357,061 qc_pass (84.9%); out-group 127,949 -> 99,360

*Evidence.* Re-derived from the obs tables of all 12 .qc.h5ad files: 119 samples, 420,381 / 357,061 (84.94%); out-group 37 samples 127,949 / 99,360. Per-study totals in qc_per_study.csv match.

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

**DQ-23 (correctness; results/qc/qc_per_sample.csv (156 rows))** — Per-sample before/after, fail counts, medians, max-genes thresholds

*Evidence.* All 156 rows re-derived from obs: cells_before/after, n_fail_low/high/mito, median genes and mito (passing cells) identical; qc_pass equals (n_genes>=min)&(n_genes<=max)&(pct_mito<th) for every cell; stored max threshold equals np.quantile(n_genes,0.995) in every sample; merged table equals the concatenation of the 10 per-study tables.

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

**DQ-24 (correctness; src/dmg_atlas/ingest.py iter_okada; data/interim/okada.qc.h5ad)** — Okada UMI counts recovered as round(expm1(x)*nCount_RNA/1e4); integer to 1e-12, row sums == nCount

*Evidence.* Re-derived for okada__586_pons (23,731 nuclei) from raw.X: max |x-round(x)| = 4.5e-13, min count 1, nnz per row == nFeature_RNA, max |rowsum-nCount_RNA| = 0.0; stored counts in okada.qc.h5ad identical to the recomputation for GFAP, PDGFRA, OLIG2, PTPRC, MT-CO1, MALAT1 and cell order; obs sample sizes equal manifest n_cells for all 22 samples.

*Impact.* Confirmed (see DQ-08 for the 10 dropped duplicate symbols).

*Recommendation.* None. (re-derived: yes)

**DQ-25 (correctness; src/dmg_atlas/ingest.py iter_liu; data/interim/liu.qc.h5ad)** — Liu: 8,984/9,248 postFilter cells recovered; missing SUDIPG55 129, ICGCGBM96 62, SUDIPG67 58, pSCG1 10, S1D2E3H7 5; 42/43 samples; RSEM expected counts 25-31% fractional

*Evidence.* Recomputed from the two postFilter metadata files and obs: 9,248 metadata cells, 8,984 recovered, missing exactly 129/62/58/10/5 (264; 206 in samples with raw files); author_sample == manifest sample for all 8,984 recovered cells; 42 samples; X has 66,294,689 non-zeros, 0 NaN, 0 negative, 28.3% fractional; multi-file samples (BT836, BT869, MUV10, MUV5) share the 32,617-gene list so no NaN from reindexing; MUV82: 96 duplicate IDs dropped (log).

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

**DQ-26 (correctness; src/dmg_atlas/genes.py; results/qc/genes_shared_all_umi_studies.csv; data/interim/gene_universe.csv)** — 18,882 genes measured in every sample of all 7 UMI studies; 14,119 across all 10; universe 36,591 symbols

*Evidence.* Recomputed from var['measured_in_all_samples'] of the 10 atlas files: 18,882 (set identical to the CSV), 14,119; gene_universe.csv has 36,591 unique symbols from 36,601 2020-A features; 13 MT- genes present in every UMI study and in Liu, 0 in Filbin/LaBelle; alias spot checks correct; 16 ambiguous symbols dropped.

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

**DQ-27 (data_integrity; docs/phase1_tumour_notes.md; results/qc/sanity_jessa_malignant_qc.csv; data/processed/sample_manifest.csv)** — Jessa: 143,611 author-malignant cells; 4 patients in both series; Filbin and Liu share 7 patients; 59 integrated samples / 334,439 author-QC cells

*Evidence.* obs: 71,899 + 71,712 = 143,611 Malignant (176,912 - 143,611 = 33,301 Normal, matching decisions.md); malignant pass fractions 0.968 / 0.873 / 0.967, normal 0.967 / 0.963 / 0.852, non-author barcodes passing 0.377 / 0.394 / 0.270 all re-derived. Manifest: patient overlap GSE210567/GSE210566 = {P-1764, P-6253, P-6337, P-6640}; filbin∩liu = 7 patients (BCH1126, BCH836, BCH869, MGH66, MUV1, MUV10, MUV5); 59 included integrated rows = 54 patients, 334,439 cells; Filbin included 3,057 cells over 9 patients, culture/PDX/Oligo columns excluded.

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

**DQ-28 (data_integrity; data/raw/MANIFEST_tumour.csv (90 rows), MANIFEST_spatial_refs.csv (25 rows))** — Every downloaded file recorded with size and sha256

*Evidence.* All 115 local_path entries exist with byte sizes equal to the manifest; 9 files re-hashed (incl. the 1.51 GB Okada h5ad and a Siletti h5ad) match their sha256; the one logged GEO FAIL (GSM6432686 metadata) was re-fetched 6 min later and manifested.

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

**DQ-29 (data_integrity; data/external/braun2023/braun2023_hindbrain_pons_medulla.h5ad; siletti2023/*.h5ad; REFERENCES.md)** — Braun subset 286,973 cells x 58,226 genes (Hindbrain 63,238; Pons 96,052; Medulla 127,683), 17 donors, integer counts; Siletti 204,168 nuclei, 3 donors

*Evidence.* obs/X attrs read with h5py: 286,973 x 58,226; Region counts exactly as stated; 17 donors; first 2e6 X values all integer (min 1). Six Siletti files total 204,168 nuclei, 3 donors, 58,232 genes, dissections all 'Pons (Pn)'.

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

**DQ-30 (correctness; data/interim/*.qc.h5ad obs)** — obs index = <sample_id>:<original cell id>; author columns preserved

*Evidence.* All 12 files: index unique, 100% of rows satisfy index == sample_id + ':' + original_cell_id, no duplicated (sample_id, original_cell_id); author_* columns present as documented (Collot author_celltype vocabulary OPC_like/Immune/MES_like/AC_like/Macroglia/Vascular/OC_like/Neuronal).

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

**DQ-31 (data_integrity; src/dmg_atlas/ingest.py build_collot_cache; data/interim/collot_this_study_counts.h5ad)** — Collot: 146,567 cells in the deposit, 18,767 This_study cells from 5 samples; 24 Jessa re-processed samples excluded; mouse IUE not ingested

*Evidence.* scRNA_metadata.csv.gz: 146,567 rows, Study This_study 18,767 (5 SampleIDs) / Jessa2022 127,800 (24); cache obs 18,767 rows with the expected per-sample counts (PMC-DMG-7 5,824; -8 5,126; -9_core 2,941; -6 2,455; -9_edge 2,421); all 19,248 gene symbols map exactly to the human 2020-A universe.

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

**DQ-32 (validity; configs/qc.yaml; scripts/03_ingest_qc.py; docs/phase2_qc_notes.md)** — No circularity: thresholds fixed a priori; author labels descriptive only; scrublet seeded

*Evidence.* Thresholds live in qc.yaml and are identical across samples within a class (65 nuclei / 32 cells / 59 smartseq2 samples); author labels are used only in sanity tables, not to set thresholds; scrublet random_state 0; the config snapshot is stored in uns['phase2']; results tables are produced by 03/04 scripts whose log (logs/phase2_qc.log, 617 lines) records wall time and peak RSS per step.

*Impact.* Confirmed.

*Recommendation.* None. (re-derived: yes)

## 6. What could not be verified and why

- The full matrices of jessa_scrna / jessa_multiome / okada / desisto.outgroup beyond the targeted row blocks (RSS budget); QC flags were verified from obs for every cell, matrix content for one sample per loader type.
- The supplementary xlsx/docx tables under data/raw/<study>/supp/ against the label CSVs (no parsing script exists — DQ-02); label provenance was checked only for internal consistency (e.g. GEO metadata vs Table S3 totals: 176,912 vs 176,934).
- The ScPCA portal API JSON and Okada Cell Browser metadata beyond what the manifest script reads.
- Raw Visium matrices (02_qc_spatial_sections.py was read; spatial_sections.csv was checked for consistency, not regenerated), the OpenPedCan files, and the 11.9 GB Braun 'All' file (subset verified instead).
- data/interim/labels/collot_*.csv, jessa_*_geo_sample_summary.csv, labelle_sample_labels.csv, liu_geo_metadata_sample_summary.csv were used only through 02_build_sample_manifest.py, not inspected line by line.
- Scrublet was not re-run (seeded, config snapshot verified; thresholds and rates checked from the tables).
- qc_summary.json (generated by the same script from the same tables as qc_report.md; not separately parsed).

## 7. Deviations from the audit brief

- Matrix-level re-computation was restricted to one sample per loader type and to obs-level checks for all cells (brief: targeted re-computations; RSS cap 10 GB respected).
- sha256 re-hashing covered 9 of 115 files (all 115 checked for existence and byte size).
- `02_build_sample_manifest.py` was re-run to `/tmp` (it hard-codes its output path); nothing under the repository was written except `docs/audit/`.
