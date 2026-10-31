# Phase 5 — query mapping (zero-shot scANVI projection) and atlas v1 packaging (2026-10-17; audit revision 2026-10-19)

Code: `scripts/10_prepare_queries.py` (query HVG files), `scripts/11_map_queries.py train|uncertainty` (load_query_data + training — see §2 for what this did and did not do — kNN uncertainty,
UMAP projection, `query_mapping.h5ad`), **`scripts/11c_fix_reference_classifier.py`** (audit fix: reference-classifier predictions, §2),
`scripts/11b_mapping_summaries.py` (label-transfer tables, specificity check, cell-of-origin kNN), `scripts/11d_scanvi_validation.py` (held-out /
patient-grouped classifier accuracy, §2b), `scripts/11e_filbin_m_sensitivity.py` (§1 pseudo-count sensitivity),
`scripts/12_package_v1.py` (schema, ontology, summary table, README section, figures). `src/dmg_atlas/figstyle.py` holds the figure rcParams.
Log: `logs/phase5_mapping.log` (original run) and `logs/fix_phase5.log` (audit corrections). Environment: local conda `multiverse` (scvi-tools 1.4.2 CPU, torch 2.14.1 as recorded in `models/scanvi_query_*/run.json`,
scanpy 1.11.5, anndata 0.12, pynndescent 0.6.0, umap-learn 0.5.12; `pyarrow` added because the Modal-trained `model.pt` unpickles a pyarrow-backed
category array). Everything ran locally with 8 threads; Modal was not needed (longest step 22 min). No Phase 2–4 file was modified.

## 1. Query preparation (`data/processed/queries/<q>.hvg.h5ad`, `results/tables/query_preparation.csv`)
All queries are restricted to the 3,000 reference HVGs in the reference var order; missing genes are zero-filled. obs carries the model
fields `sample_id` (batch), `study`, `assay`, `pct_mito`, `lineage_label = 'Unknown'`.

| query | source | cells | samples | assay | counts used | HVGs missing | median HVG counts / cell |
|---|---|---|---|---|---|---|---|
| filbin | `filbin.qc.h5ad`, qc_pass, linear TPM | 3,037 | 9 | smartseq2 | pseudo-counts (below) | 39 | 30,632 |
| labelle | `labelle.qc.h5ad`, qc_pass, log2(TPM/10+1) (407 HVGs not measured in all tables) | 8,247 (was 8,240, DQ-09) | 9 | scRNA (P9 10x, 6,366) + smartseq2 (1,874) | pseudo-counts | 41 | 767 |
| liu | `liu.qc.h5ad`, qc_pass, RSEM expected counts (27.9 % non-integer values) | 12,438 (was 13,100, DQ-04) | 42 | smartseq2 (fresh cells 5,326 / frozen nuclei 3,478 / unannotated 4,296) | round(RSEM) | 5 | 33,270 |
| braun | `braun2023_hindbrain_pons_medulla.h5ad`, UMI counts (X integer-valued) | 100,000 of 286,973 | 17 donors | scRNA (10x v2/v3, 5–14 pcw) | UMI | 3 | 1,015 |
| siletti | 6 `siletti2023_pons_*.h5ad`, UMI counts | 100,000 of 204,168 | 18 donor × dissection (3 donors) | snRNA (10x v3) | UMI | 2 | 2,032 |

**Pseudo-counts for TPM cohorts** (`10_prepare_queries.py::tpm_cohort`). Linear TPM (LaBelle: `10·(2^x − 1)`) is rescaled per cell so that the
cell's total over all harmonised genes equals `n_detected_genes × m`, then rounded; `m` = median (counts / detected genes) of the matching
counts cohort in this project: **33.82** (Liu RSEM Smart-seq2, qc_pass) for Smart-seq2 tables and **2.87** (Jessa 10x scRNA) for the LaBelle
P9 10x CPM table. Median total pseudo-count per cell: Filbin 170,250, LaBelle 2,265. The HVG-level depth of Filbin pseudo-counts (30,632)
matches Liu's real RSEM depth (33,270). Note that TPM is gene-length normalised whereas RSEM expected counts and UMIs are not, so depth matching
makes the per-cell totals Liu-like but cannot make the gene-wise distribution Liu-like; no length correction is possible from the deposited tables.
**Sensitivity to m** (`scripts/11e_filbin_m_sensitivity.py`, `results/tables/phase5_filbin_m_sensitivity.csv`; Filbin re-encoded zero-shot at
m × 0.25 / 0.5 / 1 / 2 / 4): kNN-lineage agreement with the m × 1 deliverable 0.991 / 0.995 / 1 / 0.996 / 0.993 and reference-classifier agreement
0.990 / 0.995 / 1 / 0.992 / 0.985; consensus-malignant → kNN-malignant 0.987 / 0.983 / 0.979 / 0.976 / 0.973. The labels are insensitive to m over a
16-fold range, but the position is not: the median 15-NN distance to the reference rises monotonically with m (2.38 / 2.61 / 2.83 / 3.05 / 3.29;
reference self-distance 1.566), so the kNN-uncertainty percentile (≥ 0.99 at every m) partly measures the pseudo-count scale rather than biology.
These are **not measured counts**; Filbin and LaBelle are
projection-only and are never used to train the reference. `pct_mito = 0` for both (the deposited TPM tables contain no MT- genes).
Gene matching for Braun/Siletti: Ensembl ID of the 2020-A universe first, symbol (`var['Gene']`) second.

**Subsampling** (seed 0): Braun 286,973 → 100,000 stratified by `CellClass × donor_id` (proportional, floor 100 per stratum); Siletti
204,168 → 100,000 stratified by `supercluster_term × donor_id × dissection` (floor 50). `pct_mito` computed from the 37 MT- genes of each
file. Braun `sample_id` = donor; Siletti `sample_id` = donor × dissection (3 donors are the biological replicates).

## 2. Mapping protocol — what was run, and what it actually is (audit revision)
**As run** (`11_map_queries.py train`): `SCANVI.prepare_query_anndata` → `SCANVI.load_query_data(query, models/scanvi, accelerator='cpu')` with the
scArches defaults (`freeze_expression`, `freeze_decoder_first_layer`, `freeze_batchnorm_encoder`, `freeze_classifier` = True), then 100 epochs of
training per query (`plan_kwargs={'weight_decay': 0.0}`, seed 1, all cells `Unknown`, batch 128 / 512), wall 2,461 s in total. Models in
`models/scanvi_query_<q>/`.

**What the training did** (independent audit 2026-10-19, `docs/audit/integration_mapping_audit.md` IM-01/IM-02; verified by
`scripts/11c_fix_reference_classifier.py`, `results/tables/phase5_state_dict_diff.csv`, `results/qc/phase5_reference_classifier_verification.json`):
- The reference scANVI was built with scvi-tools defaults, i.e. `encode_covariates=False`: the encoder receives no batch / study / assay input
  (`z_encoder.encoder.n_cat_list == []`) and the scArches gradient hook zeroes its entire gradient. The encoder state dict is **identical** before and
  after training for all five queries (max |Δ| = 0 for every `z_encoder` tensor). Re-encoding Filbin with the reference model through
  `load_query_data` *without any training* reproduces the saved `X_scanvi` to float16 precision (max |Δ| 0.0039, r = 1.000000).
  **Every query latent in `query_mapping.h5ad` is therefore the zero-shot output of the reference encoder.** No query batch effect is removed, and
  the method is a projection, not scArches surgery. The decoder's new first-layer columns and decoder BatchNorm did move, but nothing downstream
  uses the decoder.
- `freeze_classifier=True` did **not** freeze the classifier. In scvi-tools 1.4.2 `_archesmixin._set_params_online_update` tests `'.0.' in key`,
  which matches the module prefix `classifier.classifier.0.`, so every classifier parameter gets `requires_grad=True` (verified live on Filbin:
  all 8 `classifier.*` weight/bias tensors trainable after `load_query_data`). The gradient hooks zero the `nn.Linear` weights, but the two
  BatchNorm affine layers (2 × 128 γ + 2 × 128 β) have no hook and trained for 100 epochs under the unlabelled scANVI objective. Max |Δ| of the
  classifier BatchNorm parameters: Filbin 0.48, LaBelle 0.57, Liu 0.82, Braun 1.36, Siletti 1.52. The five saved `predicted_labels.txt` /
  `soft_probs_f16.npy` were produced by five different, unintentionally re-fitted classifiers; the classifier rebuilt from each query state dict
  reproduces them 100 %.

**Correction applied.** The reference classifier (rebuilt from `models/scanvi/model.pt`: `Classifier(30 → 128 → 128 → 7)`, BatchNorm, logits →
softmax; it reproduces `models/scanvi/predicted_labels.txt` on the reference latent for 311,616 of 311,618 cells, the two exceptions being float16
rounding) was applied to `obsm['X_scanvi']` of every query cell. `predicted_lineage`, `scanvi_maxprob`, `uncertainty_prob`, `classifier_knn_agree`
and `obsm['scanvi_probs']` in `data/processed/query_mapping.h5ad` now hold these values; the drifted-classifier values are kept as
`*_pre_audit` columns / `obsm['scanvi_probs_pre_audit']`. The untrained `load_query_data` model gives the same labels (Filbin: 100 % agreement,
max |Δprob| 0.001), so the corrected deliverable is exactly what the protocol would have produced with no training. All tables below were regenerated
from the corrected file (previous versions: `results/tables/*.pre_audit.csv`). `obs['liu_prep']` (fresh cell / frozen nucleus / unannotated) was added.

| query | cells | predicted lineage counts, reference classifier (zero-shot) | previous (drifted per-query classifier) |
<!-- wave-1 table: superseded 2026-10-19 by IM-01 (retrained reference); current counts in the wave-2 block below -->
|---|---|---|---|
| filbin | 3,037 | malignant 2,776; oligodendrocyte 113; myeloid 85; astrocyte 46; vascular 7; T_NK 6; neuron 4 | malignant 2,835; oligodendrocyte 101; myeloid 93; vascular 5; neuron 2; astrocyte 1 |
| labelle | 8,240 | myeloid 6,665; malignant 901; T_NK 656; astrocyte 13; neuron 3; vascular 2 | myeloid 7,335; malignant 904; vascular 1 |
| liu | 13,100 | malignant 9,816; neuron 1,182; myeloid 998; oligodendrocyte 642; astrocyte 322; vascular 106; T_NK 34 | malignant 11,442; myeloid 935; oligodendrocyte 512; neuron 116; vascular 70; astrocyte 25 |
| braun | 100,000 | malignant 82,910; neuron 13,137; vascular 2,105; myeloid 1,338; astrocyte 385; oligodendrocyte 124; T_NK 1 | astrocyte 97,822; myeloid 1,136; vascular 1,042 |
| siletti | 100,000 | neuron 51,975; oligodendrocyte 31,318; astrocyte 8,124; myeloid 3,882; malignant 3,541; vascular 1,065; T_NK 95 | neuron 51,934; oligodendrocyte 31,231; astrocyte 11,823; myeloid 3,978; vascular 1,034 |

**Uncertainty** (`11_map_queries.py uncertainty`, unchanged): (i) `uncertainty_prob = 1 − max scANVI probability` (now from the reference
classifier); (ii) `uncertainty_knn` = mean Euclidean distance of a query cell to its 15 nearest reference cells in `X_scanvi` (pynndescent index on all
311,618 reference cells), expressed as the percentile of the reference's own 15-NN mean-distance distribution (50,000-cell reference subsample, self
excluded; `results/tables/phase5_reference_knn_distances.csv`, **median 1.566**, p90 2.145). Also stored: `knn_mean_dist`, **`knn_lineage`** =
majority `scanvi_label` of the same 15 reference neighbours with `knn_vote_frac` (HLCA-style kNN transfer). `knn_lineage` is the primary label for
off-manifold queries; `predicted_lineage` is the classifier read-out of the same latent. **UMAP projection**: umap-learn UMAP (n_neighbors 30,
min_dist 0.5, seed 0) fitted on a 100,000-cell reference subsample of `X_scanvi` (`results/tables/phase5_umap_reference_fit.csv`), queries
transformed into it (`obsm['X_umap_ref']`; not the coordinates of `atlas_core.obsm['X_umap']`).
`data/processed/query_mapping.h5ad`: 224,377 cells × 3,000 HVGs (`X = layers['counts']`), 48 obs columns, obsm `X_scanvi`, `scanvi_probs`,
`scanvi_probs_pre_audit` (7 columns in `uns['phase5']['scanvi_prob_columns']`), `X_umap_ref`; `uns['phase5']['audit_fix_2026_10_19']` records the fix.

## 2b. How good is the reference classifier? (audit IM-03; `scripts/11d_scanvi_validation.py`)
The Phase-4 figure "agreement 0.987 on the 220,309 labelled cells" mixes 90 % training data with 10 % validation data.
`results/tables/scanvi_validation_split.csv` (scvi-tools' own random 90/10 cell split; every validation patient is also in training):

| split | labelled cells | accuracy | macro recall | malignant | myeloid | oligodendrocyte | astrocyte | neuron | vascular | T_NK |
|---|---|---|---|---|---|---|---|---|---|---|
| all labelled (as reported in Phase 4) | 220,309 | 0.987 | 0.875 | 0.9995 | 0.999 | 0.928 | 0.815 | 0.970 | 0.943 | 0.467 |
| train split | 198,279 | 0.988 | 0.879 | 0.9997 | 0.999 | 0.931 | 0.832 | 0.972 | 0.945 | 0.477 |
| **validation split (held-out cells)** | 22,030 | **0.980** | **0.833** | 0.999 | 0.994 | 0.904 | **0.663** | 0.958 | 0.928 | **0.388** |

The training history (`results/tables/scanvi_training_history.csv`) shows validation accuracy peaking at 0.986 (epoch 8) and falling to 0.978 by
epoch 24 while validation classification loss rises from 0.048 (epoch 6) to 0.095 and train accuracy climbs to 0.995: the fixed 25-epoch
classifier was over-fitting when training stopped. **Patient-level generalisation** is not measured by this split. As a bound,
`results/tables/scanvi_patient_grouped_cv_summary.csv` reports a 5-fold GroupKFold on `patient_id` (52 patients) for two light classifiers trained on
the *unsupervised* scVI latent (`X_scvi`, per-fold standardised) of the same 220,309 cells — **this is not a cross-validation of scANVI** (five
scANVI retrainings were outside the local CPU budget), and the scVI embedding itself was fitted on all cells (transductive):

| classifier on X_scvi | split | accuracy (mean ± sd) | macro recall | astrocyte | T_NK | oligodendrocyte |
|---|---|---|---|---|---|---|
| logistic regression | patient-grouped 5-fold | 0.965 ± 0.009 | 0.840 | 0.394 | 0.894 | 0.847 |
| logistic regression | random cell 5-fold | 0.972 ± 0.000 | 0.869 | 0.455 | 0.921 | 0.901 |
| MLP 128-128 (scANVI classifier width) | patient-grouped 5-fold | 0.967 ± 0.008 | 0.851 | 0.460 | 0.879 | 0.848 |
| MLP 128-128 | random cell 5-fold | 0.977 ± 0.000 | 0.890 | 0.591 | 0.898 | 0.910 |

Holding out patients costs ≈ 1 point of accuracy and 3–4 points of macro recall relative to a random cell split, with the loss concentrated in
astrocyte and oligodendrocyte (fold sd up to 0.13 for oligodendrocyte: individual tumours carry most of the normal glia). Astrocyte is the weakest
class under every scheme (0.39–0.66), T_NK is weak in scANVI (0.39 held-out) but not in the light classifiers (0.88–0.92), i.e. the T_NK failure is a
property of the scANVI classifier/latent, not of the information in the data. Accuracy figures to quote for the reference classifier: ~~**0.980 held-out cells; ≈ 0.965 as a patient-held-out bound**~~ *[superseded 2026-10-19 by the wave-2 retraining — the reference scANVI was refitted with a patient-grouped split and early stopping, so the figure to quote is **0.9751 accuracy / 0.9166 balanced accuracy on 23,827 labelled cells from 5 held-out patients** (train 0.9931 on 197,867 cells / 47 patients). Held-out-patient recall: malignant 0.999, neuron 0.969, myeloid 0.950, oligodendrocyte 0.917, T_NK 0.896 (77 cells), vascular 0.883, astrocyte 0.803 — `results/tables/reflabels/scanvi_heldout_comparison.csv` and `scanvi_heldout_per_class.csv`. The tables in this section describe the pre-audit model and are kept for the record.]*

## 3. Mapping quality (`results/tables/mapping_uncertainty.csv`, `mapping_summary.csv`)
Per-query medians over cells (reference classifier). *[The table below is the **wave-1** state — superseded 2026-10-19 by IM-01/IM-21
(retrained reference, queries re-encoded zero-shot). Current values, recomputed from `mapping_uncertainty.csv` /
`mapping_summary.csv` on 2026-10-21: per-sample median max probability Filbin 0.997, LaBelle 1.000, Liu 0.999, Braun 0.997,
Siletti 0.999; cell-weighted classifier = kNN agreement Filbin 0.963, LaBelle 0.989, Liu 0.959, Braun 0.948, Siletti 0.969;
per-sample median kNN uncertainty percentile Filbin 0.999, LaBelle 0.997, Liu 0.993, Braun 0.974, Siletti 0.485; predicted-malignant
fraction Filbin 0.895, LaBelle 0.111, Liu 0.771, Braun 0.877, Siletti 0.018. The wave-1 rows are kept for the record.]*

| query | median max prob | cells with 1 − max prob < 0.01 | cells with max prob < 0.8 | median kNN uncertainty | cells with kNN uncertainty > 0.5 | classifier = kNN label (previous) | predicted malignant | kNN malignant |
|---|---|---|---|---|---|---|---|---|
| filbin | 1.000 | 0.925 | 0.023 | 0.999 | 1.000 | 0.976 (0.978) | 0.914 | 0.916 |
| labelle | 1.000 | 0.820 | 0.079 | 0.906 | 0.988 | 0.899 (0.820) | 0.109 | 0.110 |
| liu | 1.000 | 0.851 | 0.051 | 0.992 | 0.997 | 0.959 (0.869) | 0.749 | 0.755 |
| braun | 1.000 | 0.689 | 0.102 | 0.985 | 0.998 | 0.886 (0.028) | 0.829 | 0.934 |
| siletti | 1.000 | 0.931 | 0.029 | 0.542 | 0.527 | 0.985 (0.966) | 0.035 | 0.030 |

With the reference classifier, classifier and kNN labels agree for 0.95–0.99 of cells in every query against the retrained reference
(0.89–0.99 in wave 1; the original 0.028 for Braun was the drifted classifier, not a property of the latent). The classifier is still near-certain almost everywhere (1 − max prob < 0.01 for 81.4 % of all
query cells, 6.4 % below max prob 0.8; the least certain plates are Liu MUV16 22 %, MUV1 18 %, SUDIPG45 17 % of cells with 1 − max prob > 0.2),
whereas the kNN percentile flags all Smart-seq2 cohorts (median ≥ 0.91) and the fetal reference (0.985) as outside the reference's density;
Siletti (adult 10x snRNA, the protocol of 73 % of the reference) is the only query inside it (median 0.54). Because the latent is a zero-shot
projection (§2), the kNN percentile for the TPM/Smart-seq2 cohorts measures uncorrected protocol shift together with biology (and, for the TPM cohorts,
the pseudo-count scale, §1); the `fraction > 0.5` statistic saturates for them and should be read with the agreement columns.

**Liu by group** (`mapping_liu_by_group.csv`; classifier numbers, previous in brackets): fresh cells (5,326 cells, 18 plates) — author-malignant →
malignant 0.974 (0.995), author-normal → matching lineage 0.965 (0.988), classifier/kNN agreement 0.987; frozen nuclei (3,478, 24 plates) — 0.915
(0.990) / 0.762 (0.691) / 0.952; unannotated cells (4,296) — 0.604 (0.890) predicted malignant, agreement 0.929. Adult (3,507, 12 samples)
0.980 / 0.768; paediatric (9,593, 30) 0.945 / 0.810. Pontine 0.955 / 0.815, thalamic 0.967 / 0.818, spinal (541 cells, 2 samples; 37 author-normal
cells) 0.868 / 0.189. Nuclei plates lose both normal-type accuracy and some malignant recall.

## 4. Label transfer vs author / Phase-3 labels (`mapping_<query>_author_vs_predicted.csv`, `mapping_author_agreement.json`, `mapping_agreement_sample_level.csv`)
Pooled cell-level rates, with the per-sample median [IQR] and a 95 % CI from bootstrapping samples (plates / samples / donors; cells within a sample
are not independent) from `mapping_agreement_sample_level.csv`:
- **Liu** (author states, 35 annotated plates): 7,335 author-malignant cells → predicted malignant **0.956** [per-plate median 0.984, IQR 0.922–0.995;
  CI 0.931–0.974] (kNN 0.960, CI 0.935–0.978); per state OPC-like-1 0.960, OPC-like-2 0.996, OPC-like-3 0.975, AC-like 0.929, MES-like 0.963,
  OC-like 0.841, S 0.977, G2M 0.981, frozen-OPC-like-b 0.921. 1,370 author-normal cells → matching lineage **0.800** [0.875, 0.675–1.000; CI
  0.718–0.867] (kNN 0.808): Oligodendrocyte 0.920, Microglia 0.912, Myeloid 0.833, Endothelia 0.446 (kNN 0.508), T cell **0.333** (72 cells,
  10 plates; CI 0.14–0.53; kNN 0.542, CI 0.41–0.70). Author-normal cells predicted malignant: 0.088 (previously 0.170). The previous 0.993 malignant
  recall / 0.000 T-cell recall were the drifted classifier's.
- **Filbin** (no author per-cell states in the deposit; Phase-3 consensus used): 2,325 consensus-malignant → malignant **0.978** [0.989, 0.975–1.000;
  CI 0.956–0.995] (kNN 0.978); 74 consensus-normal → matching type 0.959; 638 consensus-ambiguous → 0.784 malignant; CNV-malignant (tier 2) → 0.938.
- **LaBelle** (immune author labels only; one sample carries them, so no CI): 6,366 author immune cells → matching lineage **0.923** by the classifier
  (Myeloid 1.000, Tcell **0.559**) and **0.981** by kNN (Myeloid 0.997, Tcell 0.908). 752 consensus-malignant tumour-table cells → malignant 0.999;
  74 consensus-malignant cells from the immune SS2 table → myeloid.
T cells are now partly recovered by the classifier (LaBelle 0.56, Liu 0.33) and mostly by kNN (0.91 / 0.54); the remaining gap is the reference
classifier's weak T_NK decision region (held-out recall 0.39, §2b; 764 reference cells). `knn_lineage` remains the label to use for immune sub-typing.

## 5. Specificity check on normal references (`mapping_normal_reference_specificity.csv`)
Classes with a reference counterpart (expected lineage): Braun Oligo→oligodendrocyte, Vascular/Fibroblast→vascular, Immune→myeloid,
Neuron/Neuroblast/Neuronal IPC→neuron; Siletti Oligodendrocyte/COP→oligodendrocyte, Astrocyte/Bergmann/Ependymal→astrocyte, Microglia→myeloid,
Vascular/Fibroblast→vascular, all neuronal superclusters→neuron. No counterpart: Braun Radial glia, Glioblast, Erythrocyte, Placodes, Neural crest;
Siletti Oligodendrocyte precursor, Choroid plexus, Splatter. CIs are over donors (Braun 17) or donor × dissection (Siletti 18).

**Siletti 2023 adult pons** — specific for differentiated types; the adult OPC lineage lands next to malignant cells. Classifier: predicted malignant
**3.5 %** of all nuclei (3,541; CI 2.6–4.4 %; 2.6 % with max prob ≥ 0.8), of which 3,360 are **Oligodendrocyte precursors** (0.580 of the 5,792 OPCs;
0.404 astrocyte) and 62 committed OPCs (0.400); within the 78,137 nuclei of classes with a counterpart only **0.2 %** are called malignant (CI
0.1–0.3 %) and **0.989** receive the expected lineage (CI 0.987–0.992; Oligodendrocyte 1.000, Astrocyte 0.995, Microglia 1.000, Vascular 0.977,
named neuronal superclusters 0.95–1.00; Ependymal 0.897, Bergmann glia 0.705 (0.295 malignant), Fibroblast 0.581 (0.392 astrocyte), COP 0.523,
'Miscellaneous' 0.037). kNN label: malignant for 3.0 % of nuclei (1.5 % with vote ≥ 0.8), 0.2 % within counterpart classes, OPCs 0.481, COPs 0.510.
The previous "0 of 100,000 predicted malignant" was an artifact of the drifted classifier; the corrected classifier and the kNN agree that the only
adult population drawn into the malignant region is the OPC/COP lineage — the lineage without a normal counterpart in the tumour reference and the
putative cell of origin.

**Braun 2023 fetal hindbrain (5–14 pcw)** — fetal progenitors and immature neurons are called malignant by both read-outs. Classifier: predicted
malignant **82.9 %** of cells (CI over donors 77.6–87.6 %; 77.5 % with max prob ≥ 0.8): Radial glia 0.998, Glioblast 1.000, Neuronal IPC 1.000,
Neural crest 0.983, Neuroblast 0.916, Oligo 0.938, Neuron 0.651 (0.348 neuron), Erythrocyte 0.654; Immune → myeloid 0.992 and Vascular → vascular
0.993 are recovered, Fibroblast 0.315 vascular. Expected-lineage recovery over the 62,521 cells with a counterpart **0.259** (CI 0.19–0.34; previous
0.035). kNN label: malignant for 93.4 % of cells (90.2 % with vote ≥ 0.8) — Radial glia 1.000, Glioblast 1.000, Neuroblast 0.989, Neuron 0.927,
Neuronal IPC 1.000, Oligo 0.960 — with Immune (0.990 myeloid) and Vascular (0.994) correct; expected-lineage recovery by kNN 0.093. Median 15-NN
distance to the reference 2.604 (reference self-distance 1.566), kNN uncertainty 0.985: the fetal cells sit off the reference manifold, and both the
classifier and the neighbourhood place them in the malignant (OPC-/AC-like) region. Interpretation: the reference has **no normal progenitor or
immature-neuron class** — its normal lineages are mature, post-natal types — so fetal progenitors / neuroblasts can only land in the malignant region;
this is the expected behaviour of a tumour reference and the basis of the cell-of-origin view below, not evidence that the classifier is wrong about
tumour cells. The earlier narrative (97.8 % "astrocyte", "frozen classifier extrapolating off-manifold", the 0/10/30/100-epoch flip) described the
BatchNorm drift of IM-01; the 0-epoch run in that sensitivity series was the correct read-out and is what the deliverable now contains. For
out-of-distribution queries, scANVI probabilities remain uninformative as an uncertainty (median max prob 1.0 on Braun); `uncertainty_knn` and
`knn_lineage` are the quantities to use.

## 6. Cell-of-origin view (`malignant_nearest_fetal_class.csv`, `malignant_nearest_fetal_class_per_cell.csv.gz`, `malignant_nearest_fetal_class_per_sample_dispersion.csv`)
kNN transfer (k = 15, pynndescent index on the 100,000 zero-shot-projected Braun cells, queried with the 169,571 consensus-malignant reference cells in
`X_scanvi`). Braun `CellClass` of the majority; `cell_type` and `Region` majorities also stored. Uses latents only, so IM-01 does not affect it; the
Braun latent is the zero-shot projection of §2 (not a "surgery-trained" latent). Braun's "Oligo" class at 5–14 pcw is the oligodendroglial lineage
(pre-OPC / OPC / COP), not mature oligodendrocytes.

Overall (median 15-NN distance 3.279 vs 1.566 within the reference): **Glioblast 0.547, Oligo 0.324, Radial glia 0.094**, Neural crest 0.019,
Neuroblast 0.003, Neuron 0.002. Against the Braun subsample class frequencies (Glioblast 0.115, Oligo 0.011, Radial glia 0.243, Neuroblast 0.234,
Neuron 0.316) this is an enrichment of **28.5× for Oligo, 4.8× for Glioblast, 0.39× (depletion) for Radial glia** and ≈ 0 for Neuroblast / Neuron.
Per sample (55 samples with ≥ 50 malignant cells): Glioblast median 0.46 [IQR 0.37–0.55], Oligo 0.40 [0.16–0.57], Radial glia 0.07 [0.02–0.17] — the
Oligo/Glioblast balance varies widely between tumours (it tracks OPC-like vs AC-like composition), the radial-glia depletion is consistent.

| Jessa author state (consensus-malignant cells) | n | samples | top fetal class | second | enrichment Oligo / Glioblast / Radial glia |
|---|---|---|---|---|---|
| OPC | 31,164 | 27 | Oligo 0.628 | Glioblast 0.327 | 55 / 2.9 / 0.16 |
| Proliferating OPC | 313 | 22 | Oligo 0.760 | Radial glia 0.131 | 67 / 0.9 / 0.54 |
| Oligodendrocytes | 2,689 | 22 | Oligo 0.880 | Glioblast 0.080 | 77 / 0.7 / 0.03 |
| Astrocytes | 7,965 | 25 | Glioblast 0.783 | Radial glia 0.194 | 0.1 / 6.8 / 0.80 |
| Ependymal | 1,627 | 16 | Glioblast 0.701 | Radial glia 0.294 | 0 / 6.1 / 1.2 |
| Glial progenitors | 1,821 | 16 | Glioblast 0.506 | Radial glia 0.302 | 15 / 4.4 / 1.2 |
| RGC | 230 | 15 | Radial glia 0.744 | Glioblast 0.170 | 6.5 / 1.5 / 3.1 |
| Neuronal progenitors | 35 | 11 | Radial glia 0.686 | Glioblast 0.143 | 5.0 / 1.2 / 2.8 |
| Neurons | 457 | 21 | Glioblast 0.466 | Radial glia 0.320 | 2.7 / 4.1 / 1.3 |
| Immune (author) | 74 | 12 | Glioblast 0.311 | Immune 0.297 | 13 / 2.7 / 0.9 |
| Vascular & other | 160 | 13 | Radial glia 0.550 | Glioblast 0.263 | 2.8 / 2.3 / 2.3 |

The ordering is the one expected from the OPC-like / AC-like dichotomy: OPC-like states neighbour the fetal oligodendroglial lineage, astrocyte- /
ependymal-like states neighbour glioblasts, and the rare RGC / neuronal-progenitor states are the only ones enriched for radial glia. Per-Leiden-cluster
rows (28 clusters, `leiden_scanvi_r1`) are in the same table. Caveats: this is a nearest-neighbour statement in an off-manifold region (median
query–reference distance 3.279 vs 1.566 within the reference) made with a reference that contains no developmental cells; the ranking and the
enrichments are informative, the absolute fractions are not calibrated, and it is not cell-of-origin evidence in the sense of a developmental
reference integrated into the model.

## 7. Packaging (`scripts/12_package_v1.py`) and figure captions
`docs/schema/atlas_schema.md` (atlas_core obs fields + obsm/var/layers/uns; query_mapping obs fields), `docs/schema/label_ontology.csv`,
`results/tables/atlas_v1_summary.csv` (6 reference studies + 5 queries + total; query `malignant_fraction_scanvi` now from the reference classifier:
braun 0.829, filbin 0.914, labelle 0.109, liu 0.749, siletti 0.035; previous 0.000 / 0.934 / 0.110 / 0.873 / 0.000 in
`atlas_v1_summary.pre_audit.csv`), README section "Atlas v1 contents" (checksums updated). Figures carry no titles or text inside panels; panel
letters sit outside the axes; previous versions are `*.pre_audit.png`.

**Figure `results/figures/phase5_query_mapping.png`.** Reference UMAP (umap-learn refit on a 100,000-cell reference subsample of `X_scanvi`, grey)
with the projected query cells overlaid, one column per query and at most 30,000 randomly drawn cells per panel: (a, f) Filbin 2018, Smart-seq2,
3,037 cells; (b, g) LaBelle 2025, Smart-seq2 + 10x, 8,247 cells; (c, h) Liu 2022, Smart-seq2, 12,438 cells;
(d, i) Braun 2023 fetal hindbrain, 5–14 pcw, 30,000 of 100,000 cells shown; (e, j) Siletti 2023 adult pons, 30,000 of 100,000 shown. Top row (a–e): cells coloured by `predicted_lineage` (reference scANVI
classifier applied to the zero-shot latent); bottom row (f–j): the same cells coloured by `knn_lineage` (majority label of the 15 nearest reference
cells). Colours: malignant red, myeloid blue, T/NK magenta, oligodendrocyte green, astrocyte orange, neuron purple, vascular teal. Fetal
progenitors and neuroblasts (d, i) occupy the malignant region under both labels; adult nuclei (e, j) fall on the normal lineages except the OPC
cluster (orange/red cells at the malignant border).

**Figure `results/figures/phase5_mapping_uncertainty.png`.** Per query sample, ordered by median kNN uncertainty within each query: (a) Filbin 2018,
9 samples; (b) LaBelle 2025, 9; (c) Liu 2022, 42 plates; (d) Braun 2023, 17 donors; (e) Siletti 2023, 18 donor × dissection samples. Top row: violin
of `uncertainty_knn` (percentile of the cell's mean 15-NN distance in the reference's own 15-NN distance distribution; black dot = median; dotted line
0.5). Bottom row: fraction of cells whose reference-classifier label equals the kNN label. The `1 − max probability` row of the previous version was
dropped because the classifier is near-certain for every sample (median 0.000; per-sample values remain in `mapping_uncertainty.csv`).

## 8. Deviations and caveats
1. Filbin has no per-cell author states in the GEO deposit (Phase 2 notes); its crosstab and agreement use Phase-3 `consensus_call` / `cnv_call`.
   LaBelle has immune labels only; its tumour SS2 cells are also judged against the consensus call.
2. Braun subsampled 286,973 → 100,000 and Siletti 204,168 → 100,000 (stratified; §1), as the plan allowed. Specificity and cell-of-origin numbers
   rest on these subsamples.
3. TPM cohorts mapped with pseudo-counts (§1 rule; `m` calibrated on Liu RSEM and Jessa 10x; m-sensitivity in §1). `pct_mito = 0` for Filbin/LaBelle.
4. **The plan's "scArches fine-tuning with frozen reference weights" was not achieved and could not be with this reference** (`encode_covariates=False`;
   §2). Phase 5 is a zero-shot projection by the reference encoder; the 100-epoch training in `11_map_queries.py train` changed nothing that is used
   (and, through the scvi-tools 1.4.2 key-matching bug, re-fitted the classifier BatchNorm — corrected by `11c_fix_reference_classifier.py`).
   Adaptive mapping would require retraining the reference with `encode_covariates=True` (and, per scArches, layer norm instead of batch norm).
5. `knn_lineage`, `knn_vote_frac`, `classifier_knn_agree` were added beyond the plan; `knn_lineage` is the recommended label for off-manifold queries.
6. The Braun epoch-sensitivity runs (0/10/30 epochs, 20k subsample) quoted in the pre-audit §5 were exploratory, never saved, and are superseded: the
   0-epoch result is now the deliverable. The 15-NN distances quoted there (2.60 / 2.66 / 1.52) were from that unsaved run; the wave-1 saved values
   were 2.25 (Braun → reference) and 1.25 (reference self), and the current values after the wave-2 re-encoding are **2.604 and 1.566**.
7. Batch size (128 / 512) was not specified; the reference used 1,024. Training ran on CPU with 8 threads; Modal was not used.
8. UMAP refit on a 100,000-cell reference subsample (umap-learn) because the stored scanpy UMAP has no transform model; `X_umap_ref` ≠ `X_umap`.
9. kNN uncertainty uses an approximate index (pynndescent) and a 50,000-cell reference distance distribution; the percentile saturates at ≈ 1 for the
   Smart-seq2 cohorts (§3).
10. `nice -n 10` is not permitted in the sandbox; thread cap 8 applied. `pyarrow` installed into `multiverse` to unpickle the reference `model.pt`.
11. Liu `loc` axis mixes author locations (pontine/thalamic/spinal) with manifest sites (pons/thalamus) for the 4,296 unannotated cells. All 12,438 Liu
    cells carry `assay = 'smartseq2'` although 3,478 are frozen nuclei; `obs['liu_prep']` records this (not a model covariate — irrelevant to the latent
    given §2).
12. `run.json` of the query models records a hard-coded `freeze` dictionary that does not describe what scvi-tools did (§2); `trainable_params` counts
    `requires_grad=True` tensors including hooked ones. Versions: scvi-tools 1.4.2, torch 2.14.1 (from `run.json`), not torch 2.10.0 as previously written.
13. The patient-grouped estimate in §2b is a bound from light classifiers on `X_scvi`, not a scANVI cross-validation.

## Audit corrections (2026-10-19)
Independent audit `docs/audit/integration_mapping_audit.md` (findings IM-01 … IM-24 in `docs/audit/ALL_findings.csv`); per-finding closure in
`docs/audit/fixes_phase5.csv`. Pre-audit version of these notes: `docs/archive/phase5_mapping_notes.pre_audit.md`.
- **IM-01 (critical)** — query classifier not frozen (BatchNorm γ/β drift). Verified by state-dict diff and live `load_query_data`; all query
  predictions replaced by the reference classifier applied to the saved latents (`scripts/11c_fix_reference_classifier.py`); tables, figures, README and
  §2–7 regenerated; previous columns/tables kept as `*_pre_audit`. Headline reversals: Braun 0 → 82,910 malignant (97,822 → 385 astrocyte); Siletti
  0 → 3,541 malignant (3,360 OPCs); LaBelle T cells → T_NK 0.000 → 0.559; Liu T cell 0.000 → 0.333, normal → match 0.747 → 0.800, malignant → malignant
  0.993 → 0.956; `classifier = kNN` for Braun 0.028 → 0.886.
- **IM-02 (major)** — latent cannot adapt (`encode_covariates=False`). Confirmed (encoder max |Δ| 0; zero-shot Filbin latent = saved, max |Δ| 0.0039);
  Phase 5 re-described as zero-shot projection (§2, deviation 4, README).
- **IM-03 (major)** — 0.987 is training accuracy. Held-out 0.980 and per-class validation recall reported (§2b) with a patient-grouped bound (superseded in wave 2 by the held-out-PATIENT accuracy 0.9751);
  Phase-4 notes §5 amended (`docs/phases/phase4_integration_notes.md`, Audit corrections block).
- **IM-04 (major)** — narrative/figures rewritten (§2–7); figures regenerated from the corrected labels.
- IM-09 saved distances (1.25 / 2.25 / 2.82) now quoted — **superseded**: the wave-2 zero-shot re-encoding against the
  retrained reference changed them to 1.566 (reference self, p90 2.145) / 2.604 (Braun) / 3.279 (malignant → Braun); the body text above
  quotes the current values; IM-10 base-rate enrichment and per-sample dispersion added (§6); IM-12 m-sensitivity added
  (§1); IM-13/14 figures re-drawn with caption-only identification and LaBelle year harmonised to 2025 (`figstyle.QUERY_NAMES`, README); IM-16 versions
  corrected and the unsaved sensitivity runs declared superseded (deviations 6, 12); IM-21 no change needed; IM-22 per-sample IQR and sample-bootstrap
  CIs added (§4–5, `mapping_agreement_sample_level.csv`); IM-24 `obs['liu_prep']` added.

## Audit corrections — wave 2, retrained reference (2026-10-19; `docs/audit/fixes_reflabels.csv`, log `logs/fix_reflabels.log`)
The reference scANVI was retrained with the corrected consensus labels, a patient-grouped validation split and early stopping
(`docs/phases/phase4_integration_notes.md`, wave-2 block; held-out-patient accuracy 0.975, T_NK recall 0.896 vs 0.078 for the pre-audit classifier on the
same cells). All Phase 5 outputs were regenerated against it (`scripts/11g_remap_queries_zeroshot.py`, then `11b_mapping_summaries.py` and
`12_package_v1.py summary figures`). The wave-1 state of every table/figure is kept in `results/tables/phase5_wave1_snapshot/` and
`results/figures/phase5_wave1_snapshot/`; the original pre-audit files remain `*.pre_audit.*`; `data/processed/query_mapping.h5ad` keeps its
predecessor as `query_mapping.pre_audit.h5ad` (wave-1 state, which itself holds the original columns as `*_pre_audit`).
- **Query files.** `10_prepare_queries.py filbin labelle liu` re-run because the Liu and LaBelle QC files changed in wave 1 (DQ-01/04/09) and the
  Filbin pseudo-count scale *m* is the Liu counts-per-gene median: **Liu 13,100 → 12,438 cells** (42 samples), **LaBelle 8,240 → 8,247**, Filbin
  3,037 (unchanged cells; median HVG pseudo-counts per cell 30,632 → 29,484). Braun/Siletti unchanged. `query_preparation.csv` updated (previous
  `*.pre_audit.csv`). query_mapping.h5ad: **223,722 cells** (previously 224,377; 662 Liu and 1 LaBelle cell dropped, 8 LaBelle cells added).
- **Mapping = zero-shot projection, verified.** `SCANVI.load_query_data()` is used only to align each query; no training (`epochs_run 0`). For
  Filbin every z_encoder / classifier tensor of the query model equals `models/scanvi/model.pt` (34 tensors, max |Δ| **0.0**). Per-query outputs in
  `models/scanvi_query_<q>/` (zero-shot; the surgery-trained wave-1 directories moved to `models/scanvi_query_<q>_pre_audit/`). kNN uncertainties,
  kNN label transfer and the reference-UMAP projection recomputed against the new reference latent (`11_map_queries.py uncertainty`). Reference
  self 15-NN median distance 1.57 (p90 2.15; was 1.25 / 1.75); query medians Siletti 1.64, LaBelle 2.21, Braun 2.60, Liu 2.94, Filbin 3.44.
  Per-cell columns of the previous file carried over as `*_pre_audit` (original Phase 5) and `*_fix1` (wave 1) for common cells; `liu_prep` kept.
- **Headline numbers, wave 1 → wave 2 (`mapping_author_agreement.json`, `mapping_agreement_sample_level.csv`, `mapping_uncertainty.csv`):**
  Liu author-malignant → malignant 0.956 → **0.951** (per-plate median 0.975, IQR 0.927–0.996; sample-bootstrap 95 % CI 0.923–0.971); Liu
  author-normal → matching lineage 0.800 → 0.791 (CI 0.709–0.858), author-normal → malignant 0.088 → 0.071; **Liu T cell → T_NK 0.333 → 0.681**
  (72 cells, 10 plates; CI 0.417–0.857; kNN 0.819); **LaBelle T cell → T_NK 0.559 → 0.989** (1,662 cells; CI 0.986–1.000), LaBelle immune → matching
  lineage 0.923 → 0.985; Filbin consensus-malignant → malignant 0.978 → 0.955 (CI 0.903–0.998; kNN 0.967). Classifier = kNN agreement (cell-weighted):
  Filbin 0.976 → 0.963, LaBelle 0.899 → 0.989, Liu 0.959 → 0.959, Braun 0.886 → 0.948, Siletti 0.985 → 0.969. Max probability ≥ 0.8: Filbin 0.977 →
  0.937, Liu 0.949 → 0.914, Siletti 0.971 → 0.942 (the early-stopped classifier is less over-confident, consistent with the reference).
- **Normal references (`mapping_normal_reference_specificity.csv`).** Siletti adult pons: predicted malignant **3.5 % → 1.8 %** (1,830 cells;
  per-dissection median 1.2 %, CI 1.1–2.7 %), of which 1,681 are oligodendrocyte precursors (29.0 % of the 5,792 OPCs, a class without reference
  counterpart); within counterpart classes 0.2 % (0.2 % → 0.2 %; astrocytes 16 / 5,148, oligodendrocytes 1 / 31,146); expected lineage recovered
  for 99.2 % of counterpart-class cells. Braun fetal hindbrain: predicted malignant 82.9 % → **87.7 %** (CI over donors 0.833–0.916; kNN 91.2 %);
  radial glia 99.9 %, glioblasts 100 %, neuroblasts 95.3 %, neuronal IPC 100 %, neurons 80.3 %, Oligo 94.4 %; counterpart classes recovered 25.9 % →
  17.4 % (immune 0.3 % malignant, vascular 0.1 %). The off-manifold interpretation of §5 stands: first-trimester progenitors have no normal
  counterpart in the tumour reference and are absorbed by the malignant class; Braun `predicted_lineage` is not a cell-type annotation.
- **Cell-of-origin view (§6; `malignant_nearest_fetal_class*.csv`).** On the corrected 170,940 consensus-malignant cells and the new latent:
  nearest Braun class Glioblast 0.509 → **0.547**, Oligo 0.359 → 0.324, Radial glia 0.107 → 0.094, Neural crest 0.016 → 0.019 (enrichments and
  per-sample dispersion in the regenerated tables). Same caveat as §6: proximity in a tumour-trained latent, not cell-of-origin evidence.
- **Figures** `phase5_query_mapping.png` and `phase5_mapping_uncertainty.png` regenerated from the new file (captions in §7 unchanged in
  structure; the reference background is the new 100,000-cell UMAP refit `phase5_umap_reference_fit.csv`). `atlas_v1_summary.csv` query rows
  updated (Liu 12,438 cells; Braun scANVI-malignant 0.877; Siletti 0.018).
- Downstream: `docs/figures/figures.md` / `docs/figures/supp_figures.md` are generated by `scripts/30_figures.py` / `31_supp_figures.py`. **Both were
  re-run over all figures on 2026-10-21** (reporting pass), so their captions and md5 provenance tables now describe the wave-2
  deliverable: Figure 2 reports malignant 89.5 % of Filbin, 77.1 % of Liu and 68.0 % myeloid for LaBelle, all computed from
  `mapping_summary.csv`. (The 93.3 / 87.3 / 89.0 % triple recorded here in wave 2 was never what either script wrote; the figure
  scripts aggregate `n_pred_<lineage>` over all cells of each query, which is the row above.)
