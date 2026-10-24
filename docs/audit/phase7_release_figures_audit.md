# Audit — Phase 7 (multimodal / bulk anchoring), release layer and figures

Independent adversarial audit of everything produced beyond the approved core plan: build-plan Phase 8 (code Phase 7: ATAC, OpenPedCan survival, foundation-model embeddings), the release packaging / validation / environment record (code Phase 8, scripts 27–29) and the publication figure layer (scripts 30–31). Written 2026-10-19 by the audit track; nothing outside `docs/audit/` was modified. Findings table: `docs/audit/phase7_release_figures_findings.csv` (27 rows: 1 critical, 5 major, 15 minor, 6 confirmed-sound notes).

## 1. Scope — files examined

Scripts (read line by line): `scripts/24_atac_multiome.py`, `25_openpedcan_survival.py`, `26_foundation_embeddings.py`, `27_package_release.py`, `28_record_environment.py`, `29_validate_release.py`, `30_figures.py`, `31_supp_figures.py`; `scripts/19_metaprogram_scores.py` (exclusion rule only, grep).

Data: `data/processed/dmg_atlas_v1.h5ad` (h5py: obs, var, uns, obsm, X attrs; X not read), `data/processed/foundation_embeddings.npz` (all arrays), `data/processed/metaprogram_scores.parquet` (two columns), `data/processed/sample_manifest.csv`, `data/raw/jessa_multiome/GSM6432683/*features.tsv.gz` (feature types and coordinates), `env/versions.csv`, `env/environment.md`, `modal/built_env.json`.

Tables: `atac_summary.csv`, `atac_label_correlation.csv`, `atac_separability_per_sample.csv`, `atac_family_da_peaks.csv`, `atac_union_peaks.csv.gz`, `atac_pseudobulk_cpm.csv.gz` (header only), `openpedcan_dmg_cohort.csv`, `openpedcan_deconvolution.csv`, `survival_results.csv`, `foundation_model_comparison.csv`, `release_manifest.csv`, `release_genes_dropped.csv`, `release_ontology_map.csv`, `release_validation.csv`, `release_validation_terms.csv`, `phase6_output_checksums.csv`, plus the caption inputs `integration_metrics_summary.csv`, `mapping_summary.csv`, `mapping_uncertainty.csv`, `mapping_normal_reference_specificity.csv`, `mapping_liu_author_vs_predicted.csv`, `cnv_event_frequencies.csv`, `metaprogram_families.csv`, `metaprogram_masking_summary.csv`, `assay_biased_genes.csv`.

Figures (every file viewed): `results/pub_figs/main/fig1_atlas.png` … `fig6_multimodal.png` (6), `results/pub_figs/supp/figS1_malignant_calls.png` … `figS8_bulk_survival.png` (8), `results/figures/phase7_atac_separability.png`, `phase7_survival_forest.png`, `phase7_foundation_models.png`. The figures referenced by `docs/figures.md` / `docs/supp_figures.md` live in `results/pub_figs/main/` and `results/pub_figs/supp/` (PNG + PDF each; PDFs not viewed — same matplotlib figure object). There is no `figures/` or `manuscript/` folder.

Docs: `docs/phase7_multimodal_notes.md`, `release_v1.md`, `release_validation.md`, `figures.md`, `supp_figures.md`, `README.md`, `build_plan.md` (Phases 8–9), `decisions.md` (grep), `phase1_tumour_notes.md` (ATAC provenance), `phase5_mapping_notes.md` (grep). Logs: `logs/phase7_multimodal.log`, `logs/phase8_release.log`.

Not examined: `models/foundation/` and `models/c2l*` (out of scope per brief); the PDF versions of the figures; OpenPedCan raw files (not kept locally, by design); the Modal sandbox environments (no access; no Modal used). No GPU, no downloads, no pipeline re-runs; recomputations were targeted (see §2).

## 2. Method of verification

- Code review of each script against its docstring and the notes (undeclared deviations, NaN handling, joins, seeds, hard-coded constants).
- Re-derivation of headline numbers from the tables and from `dmg_atlas_v1.h5ad` obs / `foundation_embeddings.npz` in the local `multiverse` env (threads capped at 8, RSS ≤ 2.5 GB). Cox models were re-fit with `statsmodels` PHReg (Efron ties) because `lifelines` is not installed locally; Schoenfeld residuals were used for a PH check. kNN metrics were recomputed with scikit-learn on the stored embeddings.
- md5 of `dmg_atlas_v1.h5ad` and of every table / figure / doc in `release_manifest.csv` and `phase6_output_checksums.csv` recomputed.
- Every figure viewed for the house rule (no titles / text in panels, legends outside), axis labels, overlaps and agreement with its source table.

## 3. Re-derived numbers (claimed vs recomputed)

| item | claimed (where) | recomputed | status |
|---|---|---|---|
| scVI kNN family accuracy | 0.629 (foundation_model_comparison.csv, Fig 6d, notes) | 0.819 on the 10,719 cells with a finite latent; 0.630 if NaN rows are zeroed | **wrong** (P7-01) |
| scVI nuclei→cells transfer | 0.622 | 0.822 (finite cells) | **wrong** |
| scVI same-assay neighbours | 0.678 | 0.931 (finite cells); 0.938 with NaN→0 on all cells | **wrong / irreproducible** |
| scVI same-sample neighbours | 0.295 | 0.413 (finite cells); 0.302 with NaN→0 | **wrong** |
| cells with NaN scVI latent in the 15,000-cell subset | not reported | 4,281 (28.5%) | — |
| PCA / scGPT / Geneformer kNN accuracy | 0.825 / 0.757 / 0.745 | 0.825 / 0.757 / 0.745 | match |
| scGPT / Geneformer same-assay neighbours | 0.995 / 0.997 | 0.995 / 0.997 | match |
| expected same-assay / same-sample | 0.805 / 0.021 | 0.805 / 0.021 | match |
| ATAC median per-sample z (LSI) / RNA | 7.1 / 52.4 | 7.07 / 52.42 | match |
| 'above chance in every sample' | notes §1 | 2 of 10 samples z < 0 (−0.62, −0.92) | **wrong** (P7-02) |
| union peaks / in all 10 / promoter | 279,115 / 22,182 / 45,226 | 279,115 / 22,182 / 45,226 | match |
| malignant–malignant r range | 0.79–0.95 | 0.795–0.953 | match |
| malignant–normal r range | 0.39–0.66 | 0.389–0.771 | **understated** (P7-03) |
| samples behind AC-like × astrocyte r = 0.39 | not reported | 1 | — |
| family ATAC/RNA enrichment (AC, OPC, cyc, OC, cilia, MES, hyp) | 2.11/2.65, 1.38/2.35, 1.36/2.93, 1.29/1.63, 1.21/5.0, 1.09/1.39, 1.02/1.94 | identical from per-sample table | match |
| DA up/down per family | 85/109, 185/437, 51/264, 30/21, 6/190, 1/386, 0/1,321 | identical | match |
| OpenPedCan rows → independent primary → excluded | 225 → 129 → 1 | 225 (log) → 129 → 1 (C714384) | match |
| Cox n / events / median OS | 99 / 94 / 367 d | 99 / 94 / 367 d | match |
| Cox cohort median age / pons | 7.7 y / 74 | 7.59 y / 64 (7.75 / 74 are the 129-specimen values) | **mismatch** (P7-08) |
| rank:cycling adjusted HR, p | 1.240, 0.056 | 1.240, 0.0563 (PHReg) | match |
| rank:AC_like HR, p | 0.823, 0.053 | 0.823, 0.0534 | match |
| mp:malignant_OPC_like HR, p | 1.380, 0.004 | 1.380, 0.0041 | match |
| mp:share_MES_like HR, p | 0.704, 0.016 | 0.704, 0.0163; PH violated (ρ −0.76, p < 0.001) | match / PH not checked |
| rank-score min BH q | 0.20 | 0.197 | match |
| TP53-adjusted rank:cycling HR | not done | 1.25 (p 0.054); TP53 HR 1.66 (p 0.023) | conclusion robust |
| release cells × genes | 311,618 × 18,865 | 311,618 × 18,865 | match |
| genes dropped / remapped | 17 / 16 | 17 / 16 | match |
| dmg_atlas_v1.h5ad md5 | fa34ed4c… | fa34ed4c9368960930eca7b653a2fbb0 | match |
| manifest / phase6 checksums current | — | 332/332 and 122/122 md5 current | match |
| 'malignant cell' (CL:0001064) cells | 240,603 | 240,603 = 169,492 consensus-malignant + 70,924 ambiguous + 187 normal | **undocumented** (P7-05) |
| Fig 1 caption counts | 311,618; 229,043/46,288/36,287; 77.2% | identical | match |
| Fig 2 caption fractions | 93.3 / 87.3 / 89.0 / 97.8 %; 0 normal→malignant | identical | match |
| Fig 3b masked genes | 6,080 (hard-coded) | 6,077 from assay_biased_genes.csv | **mismatch** (P7-12) |
| supp_figures.md input tables | 'every number computed from tables' | 7 of 19 inputs recorded after a `--figs 1 4` re-run | **incomplete** (P7-06) |
| ontology terms | 43 found / non-obsolete / ancestor / label | 43 / 0 / 0 / 0 failures | match |

## 4. Findings, by severity

### Critical

**P7-01 — The scVI row of the foundation-model comparison is an artefact of NaN handling.** `26_foundation_embeddings.py` draws 15,000 'confident malignant' cells from `metaprogram_scores.parquet` and takes their `X_scvi_malignant` latent. 4,281 of them (28.5%) have no latent — they are `malignant_selected` cells whose `consensus_call` is 'ambiguous' (4,269) or 'normal' (12) and were never embedded by the malignant scVI model (the released h5ad has 142,047 such rows). `np.nan_to_num` turns them into zero vectors, `normalize` leaves them at the origin, and 4,281 cells become one point whose cosine neighbours are tie-broken arbitrarily. Recomputed on the 10,719 cells with a finite latent, scVI's balanced family accuracy is 0.819 (table 0.629), nuclei→cells transfer 0.822 (0.622), same-assay neighbours 0.931 (0.678, expectation 0.765) and same-sample 0.413 (0.295). The other three rows reproduce exactly from the npz. Consequently the statements 'scVI is the only embedding that mixes samples and assays, at the cost of lower within-dataset family accuracy', 'cross-assay transfer is similar for all (0.58–0.65; Geneformer best)' (notes §3, README, Fig 6d caption) and the scVI bars of Fig 6d are wrong; on valid cells scVI has the highest family accuracy and transfer and mixes assays only modestly. The surviving conclusion is the one about scGPT/Geneformer assay purity (0.995 / 0.997), which is sound.

### Major

**P7-02 — 'Separability above chance in every sample' is false.** Per-sample ATAC z-scores: −0.92, −0.62, 2.63, 3.78, 4.46, 9.67, 11.15, 20.42, 23.85, 85.42. Two samples are at or below chance and two are weak. The median (7.07) is correct.

**P7-03 — The ATAC label-correlation matrix averages over 0–9 samples per cell and reports none of this.** Reconstructed from the pseudobulk column names: the highlighted 'AC-like is the least astrocyte-like (0.39)' rests on one sample; every vascular pair and the AC-like × neuron pair rest on ≤ 1 sample; the malignant–normal range is 0.39–0.77, not 0.39–0.66 (OC-like × neuron 0.77, from 2 samples). Fig 6b draws n = 1 and n = 8 cells identically.

**P7-04 — The 'OPC-like' family contains technical programmes (cross-scope with Phase 6).** MP3 (ribosomal), MP17 (mitochondrial) and MP14 (immediate-early genes) are OPC_like members. MP3 scores 7.8× and MP17 3.4× higher in T/NK than in malignant cells, and 1.8× higher in myeloid, escaping the 2× myeloid/vascular exclusion rule. Every Phase 7 OPC-like statement (chromatin enrichment 1.38, promoter-concordance failure, NNLS OPC-like hazard, 2,500 OPC-like cells in the embedding test) inherits this mixture.

**P7-05 — Release cell types over-call 'malignant cell'.** 70,924 consensus-'ambiguous' and 187 consensus-'normal' cells are released as CL:0001064 because `cell_type_atlas` falls back to `scanvi_label`; 2,675 `qc_pass = False` nuclei are also in the release. None of this is stated in `release_v1.md`, whose headline (77.2% malignant) differs from the README's consensus-malignant fraction (0.544) for this reason.

**P7-06 — Supplementary-figure provenance was lost by a partial re-run.** `31_supp_figures.py write_doc()` records md5 only for tables read in the current invocation; the last run (`--figs 1 4`) left `docs/supp_figures.md` with 7 of 19 input tables (log 15:58:24 → 15:58:36). `30_figures.py` has the same defect (a `--figs 1` run wrote 4 inputs) but its final run was complete (27).

### Minor

P7-07 Cox specification: HRs reproduce exactly; TP53 status is available (62% of patients) and prognostic (HR 1.66, p 0.023) but unmodelled; no PH check (one NNLS variable violates PH); BH only within sets of 5–7 tests. Rank-score null result unchanged. P7-08 Cox cohort descriptors (7.7 y, 74 pons) are the 129-specimen values; the 99-patient values are 7.6 y and 64. P7-09 Overlap removal with single-cell cohorts covers only the two 'CBTN-' ScPCA patients; DeSisto (Children's Colorado, a CBTN site) patients cannot be checked. P7-10 'Peaks are called per sample' — they are the authors' cellranger-arc peaks from GEO; ATAC fragments were never downloaded; promoter definition (TSS ± 1 kb) verified correct. P7-11 DA calls: binomial sign test on 5–9 samples, no multiple-testing control over 279k peaks, not depth-normalised (acknowledged; table should be labelled exploratory). P7-12 Hard-coded caption numbers (6,080 vs 6,077 in the table; 279,115) contradict the 'computed from tables' statement. P7-13 release_manifest: 219/332 rows have no producing script; manifest predates validation and the figure layer; all md5s current. P7-14 lifelines/gdown versions not recorded; Phase 9.2 clean-environment re-run not evidenced. P7-15 No seed before scGPT gene subsampling; Geneformer tokenisation hand-rolled (no <eos>). P7-16 README: script range stale (…28 vs 31); 'progenitors predicted malignant' contradicts the table (0 predicted malignant; kNN neighbourhood malignant); Phase 7 never labelled exploratory (only 'optional' in build_plan). P7-17 House rule: Fig 4b has text labels inside the panel and overlapping start markers; legends inside axes in Fig 1a, 3b, 3c, 5e, 6a, S5a, S7c; 'T_NK' vs 'T/NK'; internal embedding/metric ids in Fig 1e and S2. P7-18 S6 dot plot shows only nominating-cluster values (no off-cluster expression) and the size legend overlaps the colour bar. P7-19 S8b 4.5-pt raw labels with a clipped CI; phase7_survival_forest.png axis dominated by HR 0.004; phase7_foundation_models.png y-labels overlap panels. P7-20 Selective captions (Fig 2c Endothelia 54% malignant omitted; Fig 3d assay effect omitted; Fig 5d uniformly positive correlations not discussed). P7-21 Local validator records NaN-embedding rows as 'info' where the official validator warns; otherwise thorough (43 terms verified).

### Confirmed sound (notes)

P7-22 Okada disease term disclosed. P7-23 Cohort and Cox numbers reproduce. P7-24 ATAC union-peak, enrichment and DA counts reproduce. P7-25 Release h5ad internally consistent (shape, gene drop list, md5, obsm, uns provenance, ontology map sums). P7-26 Fig 1 / 2 / 6c caption numbers reproduce. P7-27 PCA / scGPT / Geneformer rows reproduce exactly.

## 5. Answers to the brief's key checks

- **ATAC data**: Jessa GSE210566, 10 multiome samples, 42,004 labelled nuclei. No fragments or bigwigs were ever downloaded; the per-sample peak sets are the authors' cellranger-arc peak-barcode matrices (35,304–174,602 peaks per sample) binarised and merged (bedtools-merge logic) into 279,115 union peaks. Pseudobulk = per sample × label CPM of binarised counts (≥ 50 cells). DA = within-sample log2 ratio of accessible-cell fractions, combined across 5–9 samples by sign test (unit = sample; no depth normalisation, no FDR). Label-correlation claims rest on 0–9 samples per pair (P7-03). Separability unit = sample; not above chance in all samples (P7-02).
- **Survival**: 225 H3 K28 RNA-Seq tumour rows → 129 independent primary biospecimens (one per participant, `independent-specimens.rnaseqpanel.primary.tsv`, all 'Initial CNS Tumor', all PBTA); 1 CBTN participant removed; 100 with OS; Cox n = 99, 94 deaths. Deconvolution = NNLS of sum-normalised TPM on cell2location mean signatures (all shared genes or 50 markers/class) — unreliable, as the notes say (median T/NK 0.13; MES/OC/cycling/hypoxia zero in most patients). Cox: HR per SD, univariate and age + pons adjusted, lifelines; no sex/H3 variant/treatment (unavailable), TP53 available but unused; no PH check; BH within set. The null conclusion is not over-stated; the TP53 sensitivity run confirms it.
- **Foundation models**: scGPT whole-human (zero-shot, 512-d) and Geneformer V2-104M (768-d) on an A10G via Modal; vs PCA (3,000 malignant HVGs, computed on the subset) and the atlas malignant scVI latent; same 15,000 cells for all; metrics kNN (k = 15, cosine) 5-fold balanced accuracy, nuclei→cells transfer, silhouette, same-sample/assay neighbour fraction. The comparison is fair in cell set and metric but the scVI row is invalid (P7-01); labels are expression-derived (acknowledged).
- **Release**: 311,618 × 18,865; 17 genes dropped, 16 remapped; md5 and manifest current; 43 ontology terms verified via OLS4; official validator not run (disclosed). Cell-type provenance issue P7-05.
- **Figures**: no panel titles anywhere; text labels inside Fig 4b; legends inside several panels; plotted quantities match their tables wherever checked (Fig 1, 2, 6a–c, S8b); Fig 6d scVI bars are wrong because the table is wrong.
- **Exploratory marking**: Phase 8 is 'optional' in `build_plan.md`; README/notes call the foundation work a 'methods note' and the survival result 'null rather than negative', but nowhere is Phase 7 labelled exploratory/unplanned (P7-16).

## 6. Verdict for this scope

**Sound with corrections** for the ATAC, survival, release and figure layers, provided P7-02/03/05/06 are fixed; **not yet defensible** for the foundation-model comparison (P7-01), whose scVI row and derived text must be regenerated before any of its comparative statements are used. The cross-scope finding P7-04 (OPC-like family composition) should be resolved by the Phase 6 track before Phase 7 OPC-like statements are kept.

## 7. What could not be verified and why

- lifelines Cox fits themselves (package absent locally): verified instead with statsmodels PHReg to 3 decimals.
- scGPT / Geneformer inference (GPU, Modal, weights): the stored embeddings were taken as given; only downstream metrics were recomputed.
- OpenPedCan raw tables (not kept locally by design): cohort logic verified from the derived table and the log (md5s matched the v15 release list at run time).
- The 'all 10 samples' ATAC log lines and the 42,004-cell count: taken from `logs/phase7_multimodal.log` (reconstructed by the authors from two overlapping runs) — not re-run.
- Official `cellxgene-schema validate`: not installed; the local transcription was reviewed instead.
- PDF versions of the figures (same figure objects as the PNGs) were not opened.
