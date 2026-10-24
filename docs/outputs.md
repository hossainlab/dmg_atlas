# Output inventory

Moved from `README.md`. Per-file md5s: `results/tables/phase6_output_checksums.csv` (Phase 6, also `docs/phase6_outputs.md`)
and `results/tables/release_manifest.csv` (release). Phase numbers follow the scripts/logs, not `docs/build_plan.md`.

## Atlas v1 (Phase 5, `12_package_v1.py`)

Reference: **311,618 cells / nuclei, 57 samples, 52 patients, 6 studies** (snRNA:229043;multiome_RNA:46288;scRNA:36287); consensus-malignant fraction 0.549. Queries (projection only): braun2023 100,000 cells / 17 samples; filbin 3,037 cells / 9 samples; labelle 8,247 cells / 9 samples; liu 12,438 cells / 42 samples (41 patients); siletti2023 100,000 cells / 18 samples.

| file | content |
|---|---|
| `data/processed/atlas_core.h5ad` | reference counts (`layers['counts']`), obs (manifest, QC, author labels, CNV calls, lineage labels, scANVI labels, Leiden), obsm X_scvi / X_scanvi / X_umap / X_scvi_malignant |
| `data/processed/query_mapping.h5ad` | Filbin 2018, LaBelle 2025, Liu 2022, Braun 2023 hindbrain, Siletti 2023 pons mapped onto the scANVI reference (predicted lineage, soft probabilities, two uncertainty measures, X_scanvi, X_umap_ref) |
| `data/processed/sample_manifest.csv` | per-sample metadata (patient, assay, chemistry, site, age, H3 variant, treatment, tissue timing) |
| `data/processed/malignant_calls.csv` | per-cell three-tier malignant calls (author / own CNV / consensus) with hallmark arm flags |
| `models/` | scvi_seed{0,1,2}, scanvi, scvi_malignant, scanvi_query_<query>; sensitivity and control runs |
| `results/tables/` | atlas_v1_summary.csv, integration_metrics.csv, malignant_call_validation.csv, mapping_summary.csv, mapping_uncertainty.csv, mapping_<query>_author_vs_predicted.csv, mapping_normal_reference_specificity.csv, malignant_nearest_fetal_class.csv, query_preparation.csv |
| `results/figures/` | phase4_umap_overview.png, phase5_query_mapping.png, phase5_mapping_uncertainty.png; the publication set is `results/pub_figs/main` (Figures 1–6) and `results/pub_figs/supp` (Supplementary Figures 1–8) with captions and input md5s in `docs/figures.md` / `docs/supp_figures.md` |
| `docs/` | decisions.md, atlas_schema.md, label_ontology.csv, phase2_qc_notes.md, phase3_cnv_notes.md, phase4_integration_notes.md, phase5_mapping_notes.md |
| `logs/` | phase2_qc.log, phase3_cnv.log, phase4_integration.log, phase4_modal_runs_*.jsonl, phase5_mapping.log (wall time + peak RSS per step) |

Label vocabulary and obs schema: `docs/label_ontology.csv`, `docs/atlas_schema.md`. Per-study summary: `results/tables/atlas_v1_summary.csv`.
Checksums: atlas_core.h5ad md5 `f14dd67c81ec90395c447c358f87c832` (regenerated 2026-10-19, reference-label remediation wave 2; previous md5 `ded17c5f9f9d5cfa2dab16c02bb64b92`, file kept as atlas_core.pre_audit.h5ad); query_mapping.h5ad md5 `9b585728cc8dcb81c0566362ced87eec` (regenerated 2026-10-19 wave 2 against the retrained reference; wave-1 md5 `7b2e0fa74b48a567bc43f717e0bb90b8` kept as query_mapping.pre_audit.h5ad; original md5 `914eb0996f43853c681c1220b2398194`).

## Phase 6 — states, metaprograms, spatial

Malignant metaprograms (per-sample cNMF on 50 samples / 46 patients → 15 metaprograms, 2 of them replicated under the pre-registered rule), per-cell metaprogram scores and families, assay-masked re-analysis, microenvironment subclusters (raw and cleaned scVI), composition models, and cell2location deconvolution of 27 Visium sections (4 cohorts) with spatial tests. Results and caveats: `docs/phase6_notes.md`; every file with md5 and structure: `docs/phase6_outputs.md` (`results/tables/phase6_output_checksums.csv`).

| script | key outputs |
|---|---|
| `13_cnmf_per_sample.py` | `k_selection.csv`, `samples_analysed.csv` |
| `13_reference_signatures.py` | `c2l_reference_classes.csv`, `c2l_reference_classes_mp.csv`, `c2l_reference_state_scores.csv.gz`, `c2l_signatures.csv`, `c2l_signatures_mp.csv`, `c2l_signatures_mp_history.csv` |
| `14_cell2location.py` | `spatial_abundance_long.csv`, `c2l_deconvolution_runs.json`, `spatial_deconvolution_brainstem.h5ad`, `spatial_deconvolution_collot_st.h5ad`, `spatial_deconvolution_huchede.h5ad`, `spatial_deconvolution_ren.h5ad` |
| `15_spatial_tests.py` | `spatial_colocalisation.csv`, `spatial_colocalisation_per_section.csv`, `spatial_niches.csv`, `spatial_niches_by_section.csv`, `brainstem_gradient.csv`, `brainstem_section_means.csv`, `ren_dmg_vs_gbm_composition.csv`, `ren_section_composition.csv` … |
| `16_metaprograms.py` | `metaprograms.csv`, `metaprogram_genes.csv`, `metaprogram_members.csv`, `replication_matrix.csv`, `metaprogram_reference_overlap.csv`, `metaprogram_cut_sensitivity.csv`, `phase6_metaprogram_jaccard.png`, `phase6_replication_matrix.png` |
| `17_tme_subclusters.py` | `tme_subclusters.csv`, `tme_subcluster_markers.csv`, `tme_subcluster_marker_scores.csv`, `tme_subcluster_labels.csv.gz`, `phase6_tme_umap.png`, `phase6_tme_marker_scores.png` |
| `18_composition.py` | `composition_estimability.csv`, `composition_estimability.pre_audit.csv`, `composition_per_sample.csv`, `composition_per_sample.pre_audit.csv`, `composition_results.csv`, `composition_results.pre_audit.csv`, `composition_results_vs_pre_audit.csv`, `malignant_state_assignment.csv.gz` … |
| `19_metaprogram_scores.py` | `metaprogram_scores.parquet`, `metaprogram_families.csv`, `metaprogram_family_composition.csv`, `jessa_cnmf_recovery.csv`, `jessa_cnmf_recovery.pre_audit.csv`, `jessa_cnmf_vs_metaprograms.csv`, `jessa_cnmf_vs_metaprograms.pre_audit.csv`, `phase6_jessa_cnmf_vs_metaprograms.png` |
| `20_assay_masked_metaprograms.py` | `assay_biased_genes.csv`, `metaprograms_assay_masked.csv`, `metaprogram_members_assay_masked.csv`, `metaprogram_masking_summary.csv`, `metaprogram_cell_vs_nucleus.csv`, `phase6_metaprogram_masked_replication.png` |
| `21_brainstem_gradient_mp.py` | `brainstem_gradient_mp.csv`, `brainstem_section_means_mp.csv`, `brainstem_mp_vs_fallback.csv`, `c2l_deconvolution_runs_mp.json`, `phase6_brainstem_gradient_mp.png` |
| `22_tme_scvi_clean.py` | `tme_clean_exclusions.csv`, `tme_subclusters_clean.csv`, `tme_subcluster_markers_clean.csv`, `tme_assay_purity.csv`, `tme_subcluster_labels_clean.csv.gz`, `phase6_tme_umap_clean.png` |

## Phase 7 — multimodal and bulk anchoring

Jessa 2022 multiome ATAC per label / metaprogram family (`24_atac_multiome.py`), OpenPedCan v15 DMG bulk deconvolution and Cox models
(`25_openpedcan_survival.py`, corrected rank-score fits in `25b_survival_rankscore_audit.py`), scGPT / Geneformer vs scVI / PCA embeddings (`26_foundation_embeddings.py`, rebuilt on a common evaluation set by `26b_foundation_comparison_rebuild.py`). All ran on Modal. Results and
caveats: `docs/phase7_multimodal_notes.md`; log: `logs/phase7_multimodal.log`.

| script | key outputs |
|---|---|
| `24_atac_multiome.py` | `atac_summary.csv`, `atac_label_correlation.csv`, `atac_separability_per_sample.csv`, `atac_family_da_peaks.csv`, `atac_pseudobulk_cpm.csv.gz`, `atac_union_peaks.csv.gz`, `phase7_atac_separability.png` |
| `25_openpedcan_survival.py` | `openpedcan_dmg_cohort.csv`, `openpedcan_deconvolution.csv`, `survival_results.csv`, `phase7_survival_forest.png` (all pre-audit: the NNLS layer was not re-fitted, audit P7-07) |
| `25b_survival_rankscore_audit.py` | `survival_rankscore_audit.csv`, `openpedcan_rankscore_audit.csv` (corrected rank-score Cox fits with a Schoenfeld PH test; the readout used in Figure 6c) |
| `26_foundation_embeddings.py` / `26b_foundation_comparison_rebuild.py` | `foundation_model_comparison.csv`, `data/processed/foundation_embeddings.npz`, `phase7_foundation_models.png` |

## Release v1 (`27_package_release.py`)

`data/processed/dmg_atlas_v1.h5ad` — 311,618 cells x 18,865 GENCODE-v48 genes in CELLxGENE schema 7.1.0 layout (raw counts in X; ontology fields, metaprogram families and cleaned TME subclusters added to obs). md5 `4a39f172e4db0a2b470de7b6ac7b81e3`. Field mapping and open items: `docs/release_v1.md`; every release file with md5: `results/tables/release_manifest.csv`. Produced by `scripts/27_package_release.py`.

## Provenance conventions

Conventions: thread cap 8 (torch/numba/BLAS); sparse CSR counts; author labels never overwritten. Every file under
results/ and data/interim/labels/ is produced by a numbered script in scripts/ (label tables: 01b; non-GEO downloads with
sha256: 01c; reference gene sets: 01d, with the Filbin 2018 tables recorded as not retrievable by script); raw downloads are
listed with size and sha256 in data/raw/MANIFEST_*.csv. Exception, declared in docs/audit/fixes_data_repro.csv: the
Liu2022_* sets in results/tables/genesets/phase6_gene_sets.csv predate script 01d and are reproduced (not byte-identically)
by its re-derivation recipe. Audit-era predecessors of regenerated *data* files are kept as *.pre_audit.* (three scripts
read them: 06b, 07d and modal/phase4b_scanvi_retrain.py). Superseded *script* versions are not kept as files — they are in
the git history.

Also in `scripts/`:`dl_tail.py` (byte-range download helpers used by step 1), and three superseded routes kept for provenance —
`14_cell2location.py` and `15_spatial_tests.py` (replaced by `14b`/`14c` and `15b`/`15c` after the spatial corrections)
and `21_brainstem_gradient_mp.py` (its output tables are superseded; `15c` draws the current gradient figure).

