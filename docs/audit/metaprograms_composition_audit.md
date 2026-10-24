# Audit — Phase 6 (non-spatial): per-sample cNMF, metaprograms, replication test, assay-masked re-test, metaprogram scores, TME subclustering, composition models

Auditor: independent adversarial review, 2026-10-19. Nothing outside `docs/audit/` was modified. Companion table: `docs/audit/metaprograms_composition_findings.csv` (37 rows: 1 critical, 8 major, 14 minor, 14 confirmed-sound notes).

## 1. Scope covered

Code read line by line: `scripts/13_cnmf_per_sample.py`, `16_metaprograms.py`, `17_tme_subclusters.py`, `18_composition.py`, `19_metaprogram_scores.py`, `20_assay_masked_metaprograms.py`, `22_tme_scvi_clean.py`, `23_package_phase6.py` (header and file groups), `src/dmg_atlas/metaprograms.py`, `src/dmg_atlas/states.py`; `docs/decisions.md` (#9), `docs/build_plan.md` (Phase 5, checkpoints), `docs/phase6_notes.md` §1, 2, 5, 6, 8, 9, 12, `docs/phase6_outputs.md` (header and 13/16 blocks).

Data and tables read: `results/tables/cnmf/k_selection.csv`, `samples_analysed.csv`, all 54 `<sample>_programs.csv`; `metaprograms.csv`, `metaprogram_genes.csv`, `metaprogram_members.csv`, `replication_matrix.csv`, `metaprogram_reference_overlap.csv`, `metaprogram_cut_sensitivity.csv`, `metaprogram_families.csv`, `metaprograms_assay_masked.csv`, `metaprogram_members_assay_masked.csv`, `metaprogram_masking_summary.csv`, `metaprogram_cell_vs_nucleus.csv`, `assay_biased_genes.csv`, `jessa_cnmf_recovery.csv`, `jessa_cnmf_vs_metaprograms.csv`, `composition_results.csv`, `composition_estimability.csv`, `composition_per_sample.csv`, `malignant_state_assignment.csv.gz`, `tme_subclusters.csv`, `tme_subclusters_clean.csv`, `tme_clean_exclusions.csv`, `tme_assay_purity.csv`, `genesets/phase6_gene_sets.csv`, `phase6_output_checksums.csv`; `data/processed/metaprogram_scores.parquet`, `tme_subcluster_labels.csv.gz`, `tme_subcluster_labels_clean.csv.gz`, 6 of 54 `cnmf_usages_*.parquet`; `atlas_core.h5ad` obs (all 311,618 cells) and `layers/counts` rows of two samples (h5py, backed). Logs: `logs/phase6_metaprograms.log`, `phase6_13_stdout.log`.

Figures viewed: `phase6_metaprogram_jaccard.png`, `phase6_replication_matrix.png`, `phase6_metaprogram_masked_replication.png`, `phase6_jessa_cnmf_vs_metaprograms.png`, `phase6_composition_effects.png`, `phase6_tme_umap.png`, `phase6_tme_marker_scores.png`, `phase6_tme_umap_clean.png`.

## 2. Method of verification

Local env `multiverse`, 8 threads, sparse / h5py block reads only, peak RSS 1.6 GB. Every headline number in the audited sections was located in a table and the key ones re-derived from the underlying files: the full metaprogram clustering was re-executed from the 426 per-sample programme CSVs (unmasked and masked, 0.01-step cut scan, random-mask control); the per-sample consensus NMF was re-run for two samples for k = 4..12; the replication matrix, reference overlaps, Jaccard cell-vs-nucleus values, composition tallies, TME compartment memberships and assay-purity values were recomputed; md5 checksums were recomputed for 10 files; an auditor-designed matched-size nearest-neighbour test was run to separate assay structure from cell-arm power.

## 3. Re-derived numbers (claimed vs recomputed)

| quantity | claimed (notes/table) | recomputed | match |
|---|---|---|---|
| samples analysed / programmes | 54 of 57 / 426 | 54 / 426 (6 cell, 48 nuclei samples) | yes |
| chosen k distribution | — | k=5:1, 6:5, 7:15, 8:11, 9:22 (41 % at ceiling); fallback rule used for 1 sample | — |
| silhouette ≥ 0.90 among (sample,k) pairs | criterion | 71.3 % pass; median silhouette 1.000 at k=4–6 | criterion non-discriminative |
| k-selection desisto__1196 / okada__1963_pons | sil, full, err in k_selection.csv | identical to 4 dp at every k; k*=9 both; k=10–12 fail coverage | yes |
| MPs / assigned at cut 0.85 | 20 / 196 of 426 | 20 / 196; membership identical | yes |
| sensitivity MPs, replicated at 0.75/0.80/0.85/0.90 | 5/14/20/19, 1/1/2/5 | identical; additionally 0.86 → 3 replicated, 0.87–0.89 → 4 | yes (but see MC-02) |
| replication matrix 20 × 6 | replication_matrix.csv | identical from metaprogram_members.csv; replicated = {MP1, MP10} | yes |
| MP1 vs Liu2022_G2M | 29/50 shared | overlap 29, J 0.414, p 1.46e-47 | yes |
| MP10 vs Hypoxia | 19 shared | overlap 19, J 0.317, p 1.05e-31 | yes |
| MP13 / MP14 overlaps | 23 / 21 | 23 / 21 | yes |
| MP1 arms | 5 cell + 29 nuclei, 4 publications | 5 + 29, 4 | yes |
| MP10 arms | 4 cell + 4 nuclei, 3 publications | 4 + 4, 3 | yes |
| assay-biased genes | 3,987 / 1,982 / 108 → 6,080 of 13,185 | identical | yes |
| masked fraction of top-50 | median 56 % (nuclei 58, cells 31) | 0.56 (0.58 / 0.31); remaining genes nuclei median 21 (min 0), cells 34.5 | yes |
| masked MPs / replicated at 0.85 | 17 / 4 | 17 / 4; sensitivity 1/2/4/8 | yes |
| random-mask control (same set sizes, 10 seeds) | — | 2–3 MPs, 1 replicated | masking effect is gene-specific |
| MP13 vs MP7, MP14 vs MP17 Jaccard | 0.12→0.20, 0.03→0.09 | 0.124→0.200, 0.031→0.091 | yes |
| malignant cells scored | 237,886 | 237,886 (all qc_pass) | yes |
| malignant_state counts | 116,630 / 51,117 / 41,952 / 16,026 / 7,219 / 4,942 | identical | yes |
| family assignment counts | 59,369 / 42,905 / 27,926 / 27,601 / 26,471 / 22,254 / 18,846 / 12,514 | identical | yes |
| family assignment by assay | not reported | scRNA: OPC 66.7 %, AC 23.8 %, MES 0.0 %, OC 0.0 %; snRNA: OC 29.2 %, MES 17.5 % | **not disclosed** |
| MP16 contamination | myeloid 0.377 vs malignant 0.029 | identical; MP17 1.84×, MP18 1.87×, MP3 1.81× (not excluded) | yes / see MC-08 |
| Jessa recovery (own cNMF ρ ≥ 0.5) | 59 % scRNA, 45 % multiome | by dataset 59.3 / 44.8 %; by assay scRNA 69.2 %, snRNA 56.1 %, multiome 44.8 % | label error |
| MP16-best author programmes | ρ < 0.25 | median 0.206, max 0.532 (4 of 17 ≥ 0.25) | no |
| Jessa heatmap rows | 153 (figure) | 175 (sample, programme) pairs; 22 collapsed by a pivot-key bug | no |
| composition tests / estimable | 348 / 234 | 348 / 234 | yes |
| assay contrast q < 0.05 | 28 / 44 | 28 / 44 | yes |
| neuron / oligodendrocyte Δlogit (assay) | 3.5 / 2.3 | 3.50 / 2.29 | yes |
| myeloid median cells vs nuclei | 37 % vs 8 % | 0.369 vs 0.078 | yes |
| H3.1 vs H3.3 celltype q | 0.002 / 0.03 / 0.03 / 0.03 / 0.057 | 0.002 / 0.030 / 0.034 / 0.034 / 0.057 | yes |
| H3.1 arm composition | 7 vs 23 samples, Jessa only | 8 H3.1 samples (7 patients), all pons, all nuclei; H3.3 13 pons / 5 thalamus, 4 cells | confound not reported |
| ONC201 vascular_5 | 0.03 vs 0.24, q 0.047, p_mixed 0.004 | 0.030 vs 0.239, q 0.047, p_mixed 0.0044; 9 vs 8 samples; denominators 32–953 | yes |
| TME compartments | 31,508 / 792 / 4,507 / 16,458; 17/8/11/10 | identical, membership recomputed from obs | yes |
| assay purity 17 → clean | 0.95→0.02, 1.00→0.33, 0.80→0.82, 0.99→1.00 | 0.947→0.020, 1.000→0.326, 0.797→0.816, 0.995→1.000 | yes |
| manifest | 122 files, md5 | 122 present; 10/10 md5 match; 108 cNMF per-sample tables + gene-set file absent | partial |
| selected-cell provenance | not reported | 29 % of selected 'malignant' cells are consensus-ambiguous, scANVI-rescued (snRNA 33 %, multiome 22 %, scRNA 12 %) | **not disclosed** |

## 4. Findings in prose (by severity)

### Critical

**MC-01 — The per-cell metaprogram-family assignment is assay-determined.** `19_metaprogram_scores.py` z-scores the family scores across all 237,886 cells pooled (lines 127–129) and takes the argmax. The resulting `metaprogram` column in `metaprogram_scores.parquet` assigns 0.0 % of scRNA cells to MES_like, OC_like or cilia and 66.7 % to OPC_like, whereas nuclei are 17–34 % MES_like and 20–29 % OC_like. This column is the input of the §10 cell2location metaprogram reference (max z ≥ 1) and of `metaprogram_family_composition.csv`. The notes describe the assignment as spreading cells thinly with implausible cilia/OC-like shares but do not disclose that it is a near-perfect assay classifier. As delivered, the per-cell family call is not usable across assays.

### Major

**MC-02 — The '2 replicated' headline sits on a cut step.** At the pre-registered cut 0.85, 2 MPs replicate; at 0.86, OC-like MP4 absorbs two cell programmes and 3 replicate; at 0.87, AC-like MP13 gains two nuclei programmes and 4 replicate. The build-plan checkpoint ('fewer than 3') therefore flips with a 0.01 change of the cut while the number of MPs stays at 20. The 0.05-step sensitivity table in §2 does not show this. The pre-registered rule was applied correctly, and the assay structure is independently supported (MC-N12), but the checkpoint language should not rest on this knife-edge.

**MC-03 — Reference gene sets have no provenance.** `phase6_gene_sets.csv` (Liu 2022 'DE author states', Hypoxia, OXPHOS, etc.) is produced by no script, has no source file in `data/raw/liu`, and is not in the checksum manifest; `states.py` states its Filbin/Liu sets were 'recalled from the published tables'. Every MP label and family rests on these lists.

**MC-04 — 29 % of the 'malignant' input is scANVI-rescued ambiguous cells, at assay-dependent rates** (snRNA 33 %, multiome 22 %, scRNA 12 %). For an analysis whose headline is assay dependence this is an undisclosed, assay-asymmetric selection step.

**MC-05 — The assay-masked re-test is circular in interpretation and leaves nuclei programmes with a median of 21 genes** (33 programmes < 10, 12 < 5, min 0) under an unchanged absolute cut. An auditor control masking the same number of random genes gives 2–3 MPs / 1 replicated, so the gain to 4 replicated is specific to the masked gene identities, not a set-size artefact — but the mask is by definition the genes that separate the assays, 'OC-like replicated' rests on 2 cell samples, and the unmasked tree already merges OC-/AC-like across assays at cut 0.86–0.87.

**MC-06 — H3.1 vs H3.3 is confounded with site and assay.** All 8 H3.1 samples are pons and nuclei (7 patients); H3.3 includes 5 thalamic and 4 cell samples; the model adjusts for assay only. 'H3.1 samples contain more normal cells' cannot be separated from a pons / tumour-content effect.

**MC-07 — MP14 ('OPC-like-3', scRNA-only) is an immediate-early / heat-shock programme** (FOS, JUN, EGR1, ATF3, HSPA1A/B, DNAJB1, DUSP1, …). It is a dissociation signature that cannot exist in nuclei; presenting it as a lineage state that fails to replicate inflates the biological reading of the assay split.

**MC-08 — Three nuclei-only lineage MPs score ~1.8× higher in normal myeloid cells than in tumour cells** (MP17 0.373 vs 0.203, MP18 0.292 vs 0.156, MP3 0.210 vs 0.116), just under the arbitrary 2× exclusion threshold. Their top genes carry generic nuclear signal rather than lineage identity.

**MC-09 — `phase6_tme_marker_scores.png` is illegible**: colourbars overprint the y-tick labels of neighbouring panels, titles inside panels, no colourbar label.

### Minor

MC-10 Jessa heatmap pivot bug (175 → 153 rows, 22 pairs averaged). MC-11 recovery statistics reported by dataset label ('scRNA' = 82 snRNA + 26 scRNA programmes; true scRNA 69 %) and a false 'ρ < 0.25' statement for MP16 matches. MC-12 silhouette ≥ 0.90 criterion non-discriminative (71 % pass), 41 % of samples at the k = 9 ceiling, non-converged NMF fits; two-sample extension to k = 12 did not change k* (provisional). MC-13 no version control; script mtimes post-date the logged runs; clocks mixed between Modal and Mac — mitigated by exact re-execution. MC-14 108 per-sample programme tables and the gene-set file missing from the manifest. MC-15 ribosomal regex catches 17 kinases/paralogues, misses RPSA. MC-16 script 13 omits the qc_pass filter (1,202 cells). MC-17 four Jessa tumours appear twice among 54 samples (54 samples = 50 tumours). MC-18 ONC201 vascular_5: unweighted logits with denominators 32–953, 38 % of the cluster from one sample, 24 BH families; notes' caveat adequate. MC-19 composition heatmap clipped at ±2 (max 3.5). MC-20 in-panel titles in four figures; overlapping labels. MC-21 'credible subtypes' myeloid_13 and T_NK_7 are cell-only clusters with cytoskeletal / ribosomal DE lists. MC-22 TME assay purity 0.95 → 0.02 is guaranteed by conditioning scVI on assay; scVI hit max_epochs. MC-23 malignant_state argmax on raw score_genes scores, assay-shifted (MES_like scRNA 13 % vs multiome 28 %).

### Confirmed sound (notes MC-N01 … MC-N14)

The replication matrix, metaprogram membership, cut-sensitivity table, reference overlaps and hypergeometric tests, assay-biased gene counts, masked re-clustering summary, cell-vs-nucleus Jaccards, composition tallies and estimability, malignant-state counts, family counts, TME compartment sizes, exclusions and assay-purity values, and 10 spot-checked checksums were all reproduced exactly from the underlying files. The per-sample consensus NMF is deterministic (two samples reproduced to 4 dp) and its settings match the notes (the nndsvd → random / 3 → 10 restart deviation is declared). The replication rule in code equals decisions.md #9. No OXPHOS/metabolic claim is made. The auditor's matched-size nearest-neighbour test (cell programmes: best Jaccard 0.264 to 5 other cell samples vs 0.120 to 6 random nuclei samples; nuclei programmes 0.187 vs 0.103) shows that the assay structure is not an artefact of the 6-sample cell arm — the headline is supported as a structural finding, with technical vs biological origin unresolved.

## 5. Overall verdict for this scope

**Sound with corrections.** The pre-registered replication test, its tables and the composition tallies are correct and reproducible, and the assay-dependence headline survives a matched-size test. The following must be corrected before the deliverables are defensible: the per-cell metaprogram-family assignment (MC-01; critical — currently an assay classifier), the provenance of the reference gene sets (MC-03), disclosure of the cut step at 0.86–0.87 and of the scANVI-rescue asymmetry (MC-02, MC-04), the interpretation of the masked re-test and of MP14 (MC-05, MC-07), the H3.1 site confound (MC-06), the contamination threshold (MC-08) and the marker-score figure (MC-09).

## 6. Not verified and why

- Full re-run of the 54-sample cNMF (45 min, declared out of scope): k-selection verified on 2 samples only; the k-ceiling test (k up to 12) is therefore provisional for the other 20 ceiling samples.
- `22_tme_scvi_clean.py` scVI training (GPU, Modal): model outputs taken as given; only the exclusion logic, label files and purity arithmetic were recomputed.
- UCell scores in `metaprogram_scores.parquet` were not recomputed from counts (would require reading 237,886 rows of the counts layer); the family/argmax arithmetic on the stored scores was verified.
- `score_genes` malignant-state scores in `malignant_state_assignment.csv.gz` were not recomputed; counts and assay crosstabs were.
- Spatial sections (§3, 4, 7, 10, 11 beyond the purity table) and `23_package_phase6.py` beyond the manifest check are outside this track.
- Whether the gene sets in `phase6_gene_sets.csv` equal the published Liu 2022 tables could not be checked (no source file; no downloads permitted).
