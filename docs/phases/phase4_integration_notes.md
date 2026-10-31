# Phase 4 — integration: scVI / scANVI reference and malignant-only scVI (2026-10-17)

Code: `scripts/07_assemble_atlas_input.py` (cohort + HVGs), `scripts/07b_make_training_files.py` (slim training files, malignant subset + HVGs),
`modal/phase4_train.py` (all GPU training; modes `scvi`, `controls`, `scanvi`), `scripts/08_integration_metrics.py` (scIB metrics + model choice),
`scripts/09_assemble_atlas_core.py` (atlas_core.h5ad, UMAP, Leiden, figure).
Log: `logs/phase4_integration.log` (wall time + peak RSS per local step), `logs/phase4_modal_runs_*.jsonl` (one JSON record per training run: epochs, ELBO, wall, GPU memory), `logs/phase4_0{7,8,9}*_stdout.log`.
Environment: local conda `multiverse` (scanpy 1.11.5, anndata 0.12, scvi-tools 1.4.2; added `scib-metrics 0.5.10` (jax 0.10.2, CPU) and `harmonypy 2.1.0`); Modal A10G image `im-gl8KES0WsH7PoFjf83hK6v` (torch 2.5.1+cu124, scvi-tools 1.4.2, anndata 0.11.4, harmonypy 2.1.0 installed at job time).
Phase 2/3 inputs were not modified.

## 1. Reference cohort (`data/processed/atlas_input.h5ad`, 9.9 GB; `atlas_input_hvg.h5ad`, 1.9 GB)
- Cell rule: Phase-3 working set (`qc_pass`; Jessa additionally `author_qc_pass`; the four cytoplasm-contaminated nuclei samples P-6254, P-6253 multiome, P-1709, P-1779 under the relaxed 20 % mito rule, via `dmg_atlas.cnv.working_mask`) **minus** `predicted_doublet_fixed` (5,777 cells removed: jessa_scrna 1,118, jessa_multiome 330, okada 4,107, collot 188, zang 33, desisto 1).
- **scpca excluded entirely** (SCPCS000637: 441 unassigned low-quality nuclei, median 362 genes; SCPCS000636: 0 qc_pass cells). Recorded in `uns['phase4']['excluded']`.
- **311,618 cells, 57 samples, 18,882 genes** (symbols measured in every sample of all UMI studies; `results/qc/genes_shared_all_umi_studies.csv`, all present in the six studies).
- Counts: int32 CSR, nnz 601,959,167 (density 0.102). `X` and `layers['counts']` hold the same matrix. obs = Phase-2 obs + `malignant_calls.csv` columns (cnv_score, cnv_call, marker_call, consensus_call, consensus_celltype_coarse, marker scores, hallmark arm flags/means) + `lineage_label`; strings stored as category.

| study | assay | samples | cells |
|---|---|---|---|
| collot | snRNA | 5 | 17,562 |
| desisto | scRNA | 2 | 1,210 |
| jessa_multiome | multiome_RNA | 10 | 46,288 |
| jessa_scrna | scRNA | 4 | 28,531 |
| jessa_scrna | snRNA | 13 | 54,545 |
| okada | snRNA | 22 | 156,936 |
| zang | scRNA | 1 | 6,546 |
| **total** | snRNA 229,043 / multiome_RNA 46,288 / scRNA 36,287 | 57 | **311,618** |

Per-sample table: `results/tables/atlas_input_cells_per_sample.csv` (smallest samples: okada 1732_pons 152, jessa_scrna P-6253 281, desisto 1437 297).

**HVGs** (`results/tables/atlas_input_hvg.csv`): 3,000 genes, seurat_v3 statistics per `sample_id` batch (57 batches, span 0.3, on counts).
*Deviation:* scanpy's `flavor='seurat_v3'` ranks genes by **median per-batch rank first**, which with 57 batches (some of 150–300 cells) selected ~25 % genes that were variable in only 1–2 samples (e.g. FOXA1, OPN4, BICRA; 25th percentile of `highly_variable_nbatches` = 2). `flavor='seurat_v3_paper'` uses the same per-batch statistics but ranks by number of batches first (the Seurat / scvi-tools behaviour); its 3,000 genes are variable in ≥ 21 samples (median 37), overlap with the median-rank list 1,465/3,000, and include PDGFRA, OLIG1, GFAP, AQP4, PLP1, MBP, PTPRC, CD74, CLDN5, MKI67, TOP2A, EGFR, BCAS1, GPR17. The `seurat_v3_paper` list is used; the median-rank-first flag is kept in `var['hvg_seurat_v3_medianrank_first']`.

## 2. scANVI label column (`obs['lineage_label']`)
malignant = `consensus_call == malignant`; consensus-normal cells take `consensus_celltype_coarse`; ambiguous/unassigned → Unknown. No cell fell in `other`.

| label | cells |
|---|---|
| malignant | 169,571 |
| Unknown | 91,309 |
| myeloid | 23,052 |
| oligodendrocyte | 14,141 |
| astrocyte | 5,045 |
| neuron | 4,438 |
| vascular | 3,298 |
| T_NK | 764 |

## 3. Training (Modal, A10G, 8 CPU, 32 GB; data read from Volume `dmg-atlas-data`)
Settings (all scVI runs): `batch_key='sample_id'`, `categorical_covariate_keys=['study','assay']`, `continuous_covariate_keys=['pct_mito']`, n_latent 30, n_layers 2, gene_likelihood nb, batch_size 1024, train_size 0.9, early stopping on `elbo_validation` (patience 5, min_delta 0), `max_epochs` = scvi heuristic `get_max_epochs_heuristic(n)` = **26** for 311,618 cells and **47** for 169,571 cells. Models saved with `save_anndata=False`; latents as float16.

| run | Modal job | cells | epochs (max) | best val ELBO | wall (s) | s/epoch | GPU mem (GB, max alloc) |
|---|---|---|---|---|---|---|---|
| scvi_seed0 | 11a58741 | 311,618 | 26 (26) | 1337.95 | 95.4 | 3.67 | 0.21 |
| scvi_seed1 | 11a58741 | 311,618 | 26 (26) | 1340.04 | 95.6 | 3.68 | 0.22 |
| scvi_seed2 | 11a58741 | 311,618 | 26 (26) | 1340.90 | 95.6 | 3.68 | 0.23 |
| pca (50 PCs, log-normalised HVGs, arpack) | 7d7d267d | 311,618 | – | – | 40.6 | – | CPU |
| harmony (harmonypy 2.1.0 torch/GPU, sample_id, converged 6 it.) | 7d7d267d | 311,618 | – | – | 2.0 | – | <0.01 |
| scvi_seed0_e60 (sensitivity, max 60) | 5fe60f4d | 311,618 | 60 (60) | 1317.78 | 238.6 | 3.98 | 0.21 |
| scvi_malignant_seed0 (provisional, superseded) | 5fe60f4d | 169,571 | 47 (47) | 1499.24 | 95.4 | 2.03 | 0.21 |
| scanvi (n_samples_per_label=100) — **rejected**, see §5 | 44a2abd3 | 311,618 | 25 | – | 230.8 | 9.23 | 0.28 |
| **scvi_malignant** (seed 1) | 44a2abd3 | 169,571 | 47 (47) | 1491.95 | 125.4 | 2.67 | 0.21 |
| **scanvi** (all labelled cells) | dd3eaba7 | 311,618 | 25 | – | 228.3 | 9.13 | 0.28 |

Early stopping never triggered: at epoch 26 the validation ELBO was still falling by ≈1 nat/epoch (seed 0: 1523.6 → 1408.2 → 1370.8 → 1353.1 → 1344.6 → 1337.9 at epochs 1/5/10/15/20/26); the 60-epoch sensitivity run reached 1317.8. Its embedding mixes batches better (iLISI study 0.075 vs 0.053, PCR assay 0.894 vs 0.839) and separates Jessa malignant author types better (fine ARI 0.426 vs 0.185) but scores lower on lineage-level bio (NMI 0.600, label ASW 0.563 vs 0.614 / 0.595), total 0.688 vs 0.739 — the 26-epoch protocol was kept; the model is stored as `models/scvi_seed0_e60/` for later comparison.
GPU time: job wall 311 + 57 + 366 + 399 + 251 = **1,384 s = 23.1 GPU-minutes of job execution**; container lifetime (the billing proxy) ≈ 113 GPU-minutes, dominated by input staging (see deviations) and warm idle.

## 4. Integration metrics (`results/tables/integration_metrics.csv`, `integration_metrics_summary.csv`, `integration_model_choice.json`)
scib-metrics 0.5.10 on a **60,023-cell stratified subsample** of the 220,309 labelled cells (per-sample proportional sampling with a floor of 50; all 57 samples), kNN k = 90 (pynndescent). Batch keys: study, assay, sample_id. Bio label: `lineage_label`. Finer check: 30,014 Jessa malignant cells with `author_Cell_type_consensus_Jessa2022` (OPC, Astrocytes, Oligodendrocytes, Glial progenitors, Ependymal, Neurons, Proliferating OPC, RGC, …; 'Uncertain' excluded).
Metrics: iLISI, kBET acceptance (per label), batch ASW, PCR comparison vs PCA, graph connectivity (batch side); cLISI, label ASW, NMI/ARI of Leiden (resolution optimised) (bio side). Each metric min–max scaled across the six embeddings; batch_score = mean of scaled batch metrics (all keys), bio_score = mean of scaled lineage bio metrics, **total = 0.4·batch + 0.6·bio** (scIB weights). **Rule: the scVI seed (0/1/2, 26-epoch protocol) with the highest total; tie-break lower best validation ELBO.** The e60 row is a sensitivity row and not eligible.

| embedding | ilisi__study | ilisi__assay | ilisi__sample_id | kbet__study | kbet__assay | kbet__sample_id | asw_batch__study | asw_batch__assay | asw_batch__sample_id | pcr__study | pcr__assay | pcr__sample_id | graph_connectivity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| scvi_seed0 | 0.053 | 0.028 | 0.013 | 0.168 | 0.194 | 0.183 | 0.882 | 0.883 | 0.825 | 0.779 | 0.854 | 0.557 | 0.905 |
| scvi_seed1 | 0.053 | 0.029 | 0.013 | 0.199 | 0.203 | 0.186 | 0.876 | 0.874 | 0.821 | 0.777 | 0.839 | 0.559 | 0.969 |
| scvi_seed2 | 0.054 | 0.03 | 0.014 | 0.169 | 0.179 | 0.169 | 0.878 | 0.881 | 0.821 | 0.767 | 0.842 | 0.54 | 0.969 |
| harmony | 0.177 | 0.134 | 0.04 | 0.428 | 0.388 | 0.507 | 0.871 | 0.869 | 0.788 | 0.745 | 0.557 | 0.587 | 0.877 |
| pca | 0.0 | 0.0 | 0.0 | 0.16 | 0.16 | 0.145 | 0.849 | 0.797 | 0.821 | 0.0 | 0.0 | 0.0 | 0.969 |
| scvi_seed0_e60 | 0.075 | 0.05 | 0.018 | 0.187 | 0.195 | 0.192 | 0.898 | 0.901 | 0.839 | 0.822 | 0.894 | 0.645 | 0.968 |

| embedding | clisi | asw_label | nmi | ari | clisi__fine | asw_label__fine | nmi__fine | ari__fine |
|---|---|---|---|---|---|---|---|---|
| scvi_seed0 | 1.0 | 0.592 | 0.612 | 0.393 | 0.997 | 0.524 | 0.34 | 0.197 |
| scvi_seed1 | 1.0 | 0.595 | 0.614 | 0.415 | 0.997 | 0.527 | 0.339 | 0.185 |
| scvi_seed2 | 1.0 | 0.59 | 0.6 | 0.398 | 0.997 | 0.524 | 0.322 | 0.151 |
| harmony | 1.0 | 0.566 | 0.496 | 0.228 | 0.989 | 0.464 | 0.332 | 0.204 |
| pca | 1.0 | 0.543 | 0.375 | 0.066 | 0.998 | 0.466 | 0.288 | 0.057 |
| scvi_seed0_e60 | 1.0 | 0.563 | 0.6 | 0.416 | 0.995 | 0.517 | 0.401 | 0.426 |

| embedding | batch_score | bio_score | fine_bio_score | total_score | rank_total |
|---|---|---|---|---|---|
| scvi_seed0 | 0.493 | 0.841 | 0.654 | 0.702 | 3.0 |
| scvi_seed1 | 0.537 | 0.874 | 0.659 | 0.739 | 1.0 |
| scvi_seed2 | 0.523 | 0.825 | 0.583 | 0.704 | 2.0 |
| harmony | 0.736 | 0.478 | 0.198 | 0.581 | 5.0 |
| pca | 0.127 | 0.125 | 0.257 | 0.126 | 6.0 |
| scvi_seed0_e60 | 0.661 | 0.706 | 0.873 | 0.688 | 4.0 |

**Chosen reference: `scvi_seed1`** (total 0.739; seeds 0/2 0.702/0.704; Harmony 0.581; PCA 0.126). Harmony mixes batches more aggressively (iLISI study 0.18, kBET sample 0.51) at the cost of lineage structure (NMI 0.50, ARI 0.23 vs 0.61/0.42). Absolute iLISI values are low for every method because `sample_id` batches coincide with patients: malignant cells are legitimately patient-specific, so perfect sample mixing is neither expected nor desirable. All scVI seeds agree within 0.04 total score; seed choice is not a strong determinant.
Local cost: 6 embeddings × (kNN + metrics) ≈ 26 min, peak RSS 5.4 GB.

## 5. scANVI (`models/scanvi/`)
`SCANVI.from_scvi_model(scvi_seed1, labels_key='lineage_label', unlabeled_category='Unknown')`, 25 epochs, batch 1024, seed 1.
- **First attempt with `n_samples_per_label=100` (as specified) failed**: the classifier sees only 100 cells per class per epoch, overfits that balanced subset (train accuracy 0.998) and relabels 97k of 169k consensus-malignant cells as astrocyte/neuron/vascular/T_NK (agreement on labelled cells **0.538**). The malignant class is too heterogeneous (OPC-/AC-/OC-/MES-like states across 50 tumours) to be summarised by 100 cells. Kept for the record as `models/scanvi_nspl100_rejected/`.
- **Final model uses all labelled cells (`n_samples_per_label=None`)**: agreement with `lineage_label` on the 220,309 labelled cells **0.987** (this pools 90 % training cells; **held-out validation accuracy 0.980**, see the Audit corrections block below) (per study: collot 0.992, desisto 0.946, jessa_multiome 0.988, jessa_scrna 0.995, okada 0.982, zang 0.994); per-class recall malignant 1.000, myeloid 0.999, neuron 0.970, vascular 0.943, oligodendrocyte 0.928, astrocyte 0.815, **T_NK 0.467** (764 cells; mostly absorbed by myeloid — the smallest class, worth a dedicated immune sub-model later).
- Unknown cells (91,309): **96.2 % receive max probability ≥ 0.8**; predicted malignant 70,924, myeloid 10,218, astrocyte 4,020, oligodendrocyte 3,522, vascular 1,435, neuron 1,092, T_NK 98. Outputs: `predicted_labels.txt`, `maxprob.npy`, `soft_probs_f16.npy` (+ `soft_prob_columns.txt`), `latent_f16.npy`.

## 6. Malignant-only scVI (`models/scvi_malignant/`)
169,571 consensus-malignant cells (56 samples), 3,000 HVGs recomputed on the subset (seurat_v3_paper, batch = sample_id with the 3 samples having < 100 malignant cells pooled into one batch: jessa_scrna P-6253, okada 1934_thalamus, desisto 1437, 177 cells; span 0.5 because 0.3 produced a loess singularity; `results/tables/atlas_malignant_hvg.csv`). Same covariates, seed 1 (the chosen seed), 47 epochs, best val ELBO 1491.95 (1796.9 → 1560.8 → 1520.7 → 1504.4 → 1497.3 → 1492.0 at epochs 1/10/20/30/40/47). Training file: `data/processed/atlas_input_malignant_hvg.h5ad`.

## 7. atlas_core.h5ad (`data/processed/atlas_core.h5ad`, 5.18 GB)
311,618 × 18,882. `layers['counts']` = full-gene int32 counts (single on-disk copy; **X is left empty** because writing X and the layer doubles the file to 9.9 GB > 6 GB — set `adata.X = adata.layers['counts']` before normalising). obsm: `X_scvi` (seed 1), `X_scanvi`, `X_umap` (UMAP on X_scanvi, n_neighbors 30, scanpy defaults, random_state 0), `X_scvi_malignant` (NaN rows for the 142,047 non-malignant cells). obs adds `scanvi_label`, `scanvi_maxprob`, `leiden_scanvi_r1` (Leiden, igraph flavour, resolution 1.0 → 28 clusters). Side tables: `results/tables/atlas_core_obs_labels.csv`, `atlas_core_umap.csv`. Figure: `results/figures/phase4_umap_overview.png` (study / assay / lineage_label / consensus_call).
Local wall: read 60 s, neighbors + UMAP + Leiden 166 s, write 5 s; peak RSS 4.9 GB.

## 8. Models directory
`models/scvi_seed0`, `scvi_seed1` (chosen), `scvi_seed2`, `scanvi`, `scvi_malignant` (deliverables); `scvi_seed0_e60` (sensitivity), `scvi_malignant_seed0` (provisional seed-0 run), `scanvi_nspl100_rejected`, `pca`, `harmony` (controls). Each holds `model.pt`, `latent_f16.npy`, `obs_names.txt`, `history.json`, `run.json`. Modal Volume `dmg-atlas-data` keeps copies of the training files and models.

## 9. Deviations and caveats
1. HVG ranking flavour `seurat_v3_paper` instead of scanpy's `seurat_v3` (reason and overlap in §1).
2. scANVI trained on all labelled cells instead of `n_samples_per_label=100` (the specified setting failed; §5).
3. Early stopping (patience 5) never triggered in any scVI run; the heuristic max_epochs is the binding limit and the ELBO was still decreasing slowly (§3). The 26-epoch protocol was kept after the 60-epoch sensitivity run did not score better.
4. Metrics were computed on a 60,023-cell stratified subsample (and 30,014 cells for the Jessa fine check), not all 311,618 cells, for local CPU/thermal limits; kBET/NMI/ARI use scib-metrics' implementations (kBET with diffusion, Leiden-based NMI/ARI), not the R originals.
5. Modal Volume upload from the local kernel is blocked by the sandbox proxy for files > 4 MB; the training files were streamed once as job inputs (0.29–0.43 MiB/s, 52 min in total) and copied onto the Volume from inside the job. This inflated container time (≈113 GPU-min) relative to job execution (23.1 GPU-min).
6. The malignant-only model was first trained with seed 0 before the seed choice was known, then re-trained with the chosen seed 1; only the seed-1 model is used.
7. `joblib`/loky process pools are forbidden in this sandbox; Leiden resolution optimisation in scib-metrics ran single-threaded.
8. `atlas_core.h5ad` stores counts in `layers['counts']` only (X empty) to stay under 6 GB (§7).
9. Doublets flagged by the fixed 0.25 Scrublet threshold were removed (5,777); cycling cells and the four high-mito nuclei samples are kept as in Phase 3.


## Audit corrections (2026-10-19; `docs/audit/integration_mapping_audit.md`, fixes in `docs/audit/fixes_phase5.csv`)
- **IM-03 — scANVI accuracy.** The 0.987 in §5 is agreement over all labelled cells, 90 % of which were training cells. From the model's own random
  90/10 cell split (`scripts/11d_scanvi_validation.py`, `results/tables/scanvi_validation_split.csv`): train 0.988 (198,279 cells), **validation 0.980**
  (22,030 cells); validation per-class recall malignant 0.999, myeloid 0.994, oligodendrocyte 0.904, neuron 0.958, vascular 0.928, **astrocyte 0.663**
  (reported 0.815), **T_NK 0.388** (reported 0.467). Validation accuracy peaked at epoch 8 (0.986) and validation classification loss rose from 0.048
  (epoch 6) to 0.095 at epoch 24 (`scanvi_training_history.csv`): the fixed 25 epochs over-fit; early stopping on validation classification loss should
  be used in any retraining. Every validation patient is also in training, so patient-level generalisation is unmeasured by scANVI; as a bound, a
  patient-grouped 5-fold estimate for light classifiers on `X_scvi` gives accuracy 0.965 (logistic regression) / 0.967 (128-128 MLP) vs 0.972 / 0.977
  under a random cell split (`scanvi_patient_grouped_cv_summary.csv`; not a scANVI cross-validation — see Phase 5 notes §2b).
- **IM-07 — covariates.** `study` and `assay` are nested in `sample_id` (every sample belongs to one study and one assay), so the scVI categorical
  covariates are linear functions of the batch one-hot and add no information; the model does not correct study and assay separately from sample.
  With `encode_covariates=False` (default) none of them enters the encoder (consequence for query mapping: Phase 5 notes §2).
- **IM-08 — label set (deviation, previously undeclared).** The scANVI labels are the 7 coarse lineages + Unknown (malignant as a single class, no
  OPC-/AC-/MES-/OC-like/cycling substates, vascular not split into endothelial/pericyte), not the plan's label list; malignant-state structure is
  therefore unsupervised in the scANVI latent and comes from Phase 6 cNMF. Added as deviation 10 below.
- **IM-15 — UMAP overview caption.** `phase4_umap_overview.png` is computed on `X_scanvi`, a label-supervised latent, so the clean lineage separation in its lineage panel is partly by construction; the shared tab10 palette across panels was not changed in this revision.
- **IM-23 — typo.** The Jessa fine-check subsample is 30,014 cells (`integration_model_choice.json`, log), not 30,018; corrected in §4 and §9.
- **IM-05 (recorded, not changed).** The scaled scIB scores were min-max normalised over six rows including the ineligible `scvi_seed0_e60`
  sensitivity run; over the five eligible rows the totals are seed1 0.766 / seed0 0.732 / seed2 0.732 / harmony 0.598 (auditor's recomputation);
  the ranking and the chosen model (seed 1) are the same either way, and seed 1 ranks first on raw batch and bio metrics. The 0.4 batch / 0.6 bio
  weighting was set before the e60 run existed (log order) but was not pre-registered in `decisions.md`.

10. (added 2026-10-19) scANVI label set = 7 coarse lineages + Unknown; malignant substates and the endothelial/pericyte split from the plan were not
    used as supervised labels (IM-08).

## Audit corrections — wave 2, reference labels (2026-10-19; `docs/audit/fixes_reflabels.csv`, log `logs/fix_reflabels.log`)
Propagates the CNV-track corrections (regenerated `data/processed/malignant_calls.csv`: per-arm smoothing, NaN handling, `chr4q12_gain` removed,
`cnv_confidence`, `consensus_call_eventtrust`, `low_confidence_sample`) into the scANVI reference, the queries and the malignant-only model.
Every regenerated file keeps its predecessor as `<name>.pre_audit.<ext>` (models: `models/<name>_pre_audit/`); previous obs/obsm columns are kept
as `*_pre_audit` inside `atlas_core.h5ad`.

- **Label rebuild (`scripts/07c_rebuild_lineage_labels.py`).** `lineage_label` re-derived with the §2 rule from the corrected `consensus_call`
  (the rule applied to the pre-audit columns reproduces the pre-audit label for 311,618 / 311,618 cells). **25,491 cells (8.18 %) change label**
  (`results/tables/reflabels/lineage_label_change_crosstab.csv`): malignant → Unknown 9,239; Unknown → malignant 10,608; Unknown → normal lineage
  2,830 (myeloid 969, oligodendrocyte 769, astrocyte 593, neuron 239, vascular 223, T_NK 37); normal lineage → Unknown 2,814 (oligodendrocyte 1,041,
  myeloid 766, astrocyte 647, neuron 177, vascular 159, T_NK 24). No cell moves between two named classes. Totals: malignant 169,571 → **170,940**
  (+1,369), Unknown 91,309 → **89,924**, myeloid 23,052 → 23,255, oligodendrocyte 14,141 → 13,869, astrocyte 5,045 → 4,991, neuron 4,438 → 4,500,
  vascular 3,298 → 3,362, T_NK 764 → 777. Changes concentrate in Okada (1890_pons 3,863 cells, 1689_thalamus 2,346, 586_pons 1,879) and Jessa scRNA
  (P-1775 1,642, P-1764 1,509) (`lineage_label_change_by_sample.csv`).
- **okada__1734_pons (decision).** `sample_qc_pass = False` (DQ-03: 1,163 genes, median 427 genes/nucleus), 22 CNV windows, 64.5 % of its 8,466
  reference cells score-only (`cnv_confidence`), `low_confidence_sample = True`. Decision: **kept in the scVI/scANVI reference** (its 8,466 cells are
  already in the chosen `scvi_seed1` latent; removing them would require retraining and re-selecting the reference model) and flagged
  (`obs['low_confidence_sample']`), but **excluded from the malignant-only model** (3,545 corrected-malignant cells) and **to be excluded from
  Phase 6 malignant analyses** (recorded in `docs/design/decisions.md`). Its lineage labels after correction: Unknown 4,919, malignant 3,545, 2 normal.
- **scANVI retrain (`modal/phase4b_scanvi_retrain.py`, Modal A10G, job cca289c6; `models/scanvi`, previous `models/scanvi_pre_audit`).** Same
  recipe as §5 (from `scvi_seed1`, `labels_key lineage_label`, unlabeled `Unknown`, all labelled cells, batch 1,024, max 25 epochs, seed 1) with two
  deliberate changes answering IM-03: (i) **patient-grouped validation** — whole patients added in a seeded random order until ≥ 10 % of labelled
  cells are held out and all seven classes are present: 5 of 52 patients (okada_1958, okada_2054, okada_1890, okada_1656, PMC-DMG-6; 31,944 cells,
  23,827 labelled = 10.7 %; `models/scanvi/split.json`, per-cell `split_cells.csv.gz`, `obs['scanvi_validation_patient']`); (ii) **early stopping on
  `validation_classification_loss`** (patience 5, best-epoch weights restored): loss 0.127 → 0.077 at **epoch 7**, training stopped after 12 of 25
  epochs (`validation_metrics.json`). Training wall 108 s (9.0 s/epoch), peak GPU memory 0.28 GB.
  Results (`results/tables/reflabels/scanvi_heldout_comparison.csv`, `scanvi_heldout_per_class.csv`): **train accuracy 0.9931** (197,867 labelled
  cells of the 47 training patients); **held-out-patient accuracy 0.9751**, balanced accuracy 0.917 (per patient 0.960–0.999); held-out recall
  malignant 0.999, neuron 0.969, myeloid 0.950, oligodendrocyte 0.917, T_NK **0.896** (77 cells; precision 0.56), vascular 0.883, astrocyte 0.803
  (`validation_confusion.csv`: astrocytes leak to malignant 46 / myeloid 40 / neuron 24; oligodendrocytes to myeloid 158). Agreement over all
  221,694 labelled cells 0.991. The pre-audit model on the **same** held-out cells scores 0.977 accuracy but balanced accuracy 0.815 (T_NK recall
  0.078, astrocyte 0.857) — and it was trained on those patients (random cell split), so its number is not a held-out estimate. Caveat: the seeded
  hold-out contains only Okada and Collot patients (two of six studies, both snRNA); a study-balanced or repeated hold-out was not run (deviation 13).
- **Labels change for 8,020 cells (2.6 %)** between the pre-audit and the new `scanvi_label` (`scanvi_label_change_crosstab.csv`); scANVI-malignant
  240,692 → **240,321** (77.1 %); Unknown cells with max probability ≥ 0.8: 96.2 % → 91.2 % (the early-stopped classifier is less over-confident).
- **atlas_core.h5ad (`scripts/09b_refresh_atlas_core.py`; 5.27 GB; previous `atlas_core.pre_audit.h5ad`).** obs: `lineage_label`, `scanvi_label`,
  `scanvi_maxprob`, `leiden_scanvi_r1`, and all CNV-derived call columns replaced from the corrected `malignant_calls.csv` (categorical call columns
  kept as `*_pre_audit`; `chr4q12_gain` → `chr4q12_gain_pre_audit`; new `cnv_confidence`, `consensus_call_eventtrust`, `low_confidence_sample`,
  `gene_dev_PDGFRA`, `gene_dev_KIT`, `scanvi_validation_patient`). obsm: `X_scanvi` replaced (`X_scanvi_pre_audit`), `X_umap` recomputed on the new
  latent (neighbors k = 30, UMAP seed 0; `X_umap_pre_audit`), Leiden 1.0 recomputed (**34** clusters, previously 28; `leiden_scanvi_r1_pre_audit`).
  `uns['phase4']['audit_fix_2026_10_19']` records provenance. Tables `atlas_core_obs_labels.csv` (now with `*_pre_audit`, confidence and
  validation-patient columns) and `atlas_core_umap.csv` regenerated (`*.pre_audit.csv` kept); `atlas_v1_summary.csv` regenerated (consensus-malignant
  fraction 0.544 → 0.549, scANVI-malignant 0.772 → 0.771).
- **Malignant-only scVI (`scripts/07d_rebuild_malignant_input.py`, `modal/phase4_train.py scvi --seeds 1 --tag scvi_malignant`, job b3e41cc8;
  `models/scvi_malignant`, previous `models/scvi_malignant_pre_audit`).** Input rebuilt from corrected `consensus_call == malignant` minus
  okada__1734_pons: **167,395 cells, 55 samples** (previously 169,571 / 56); 3,000 seurat_v3_paper HVGs, batch = sample_id with 4 samples < 100
  malignant cells pooled (190 cells: jessa_scrna P-6253, okada 1934_thalamus, desisto 1437, okada 1732_pons), span 0.5 used directly; HVG overlap with
  the previous set 2,926 / 3,000 (`atlas_malignant_hvg.csv`, previous `*.pre_audit.csv`). Training: seed 1, same covariates, 48 epochs (heuristic),
  validation ELBO 1852.7 → 1526.7 (best 1526.3; the previous model's 1491.95 is not comparable — different cells and HVGs). `atlas_core.h5ad`
  `obsm['X_scvi_malignant']` replaced (`X_scvi_malignant_pre_audit` kept): values for 167,395 cells, **NaN for the 144,223 cells not in the
  malignant-only model** (non-malignant and okada__1734_pons) — the NaN convention is unchanged from Phase 4 and is documented in
  `uns['phase4']['audit_fix_2026_10_19']['X_scvi_malignant']`.
- **Figure `results/figures/phase4_umap_overview.png`** regenerated (previous `*.pre_audit.png`) in house style (`src/dmg_atlas/figstyle.py`: panel
  letters only, legends outside, project palettes, 300 dpi; IM-15 palette issue resolved). Caption: *UMAP of the 311,618 reference cells computed on
  the retrained scANVI latent (`X_scanvi`, label-supervised — lineage separation is partly by construction). (a) study; (b) assay; (c) corrected
  `lineage_label` used to train scANVI (Unknown = ambiguous consensus, grey); (d) corrected Phase-3 `consensus_call`. Legend counts are cell numbers.*
- **GPU use:** two Modal A10G jobs, job wall 143 s + 122 s = 4.4 min (training 108 s + 104 s); container time including cold start and the
  341 MB input upload for the malignant job ≈ 20 min.
- Downstream (not re-run here): Phase 6 malignant selection (`select_malignant`), metaprogram scoring, TME sub-clustering and composition read
  `scanvi_label` / `consensus_call` / `X_scvi_malignant` from `atlas_core.h5ad` and must be re-run on the refreshed file with okada__1734_pons
  excluded; Phase 7 release (`27_package_release.py`) and the publication figures (`30_figures.py`, `31_supp_figures.py`) read the refreshed
  labels; `README.md` checksums updated below.

13. (added 2026-10-19, wave 2) The scANVI patient hold-out is a single seeded draw (5 patients, Okada + Collot only); no repeated or study-stratified
    hold-out was run, so the 0.975 held-out accuracy carries study-composition uncertainty. The pre-audit classifier comparison on the same
    patients is a training-set number for that model.
