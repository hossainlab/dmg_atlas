# Phase 5 — query mapping (scArches / scANVI) and atlas v1 packaging (2026-10-17)

Code: `scripts/10_prepare_queries.py` (query HVG files), `scripts/11_map_queries.py train|uncertainty` (surgery + training, kNN uncertainty,
UMAP projection, `query_mapping.h5ad`), `scripts/11b_mapping_summaries.py` (label-transfer tables, specificity check, cell-of-origin kNN),
`scripts/12_package_v1.py` (schema, ontology, summary table, README section, figures). `src/dmg_atlas/figstyle.py` holds the figure rcParams.
Log: `logs/phase5_mapping.log` (wall time + peak RSS per step). Environment: local conda `multiverse` (scvi-tools 1.4.2 CPU, torch 2.10.0,
scanpy 1.11.5, anndata 0.12, pynndescent 0.6.0, umap-learn 0.5.12; `pyarrow` added because the Modal-trained `model.pt` unpickles a pyarrow-backed
category array). Everything ran locally with 8 threads; Modal was not needed (longest step 22 min). No Phase 2–4 file was modified.

## 1. Query preparation (`data/processed/queries/<q>.hvg.h5ad`, `results/tables/query_preparation.csv`)
All queries are restricted to the 3,000 reference HVGs in the reference var order; missing genes are zero-filled. obs carries the model
fields `sample_id` (batch), `study`, `assay`, `pct_mito`, `lineage_label = 'Unknown'`.

| query | source | cells | samples | assay | counts used | HVGs missing | median HVG counts / cell |
|---|---|---|---|---|---|---|---|
| filbin | `filbin.qc.h5ad`, qc_pass, linear TPM | 3,037 | 9 | smartseq2 | pseudo-counts (below) | 39 | 30,632 |
| labelle | `labelle.qc.h5ad`, qc_pass, log2(TPM/10+1) (all 8,240 cells detected as log scale; 407 HVGs not measured in all tables) | 8,240 | 9 | scRNA (P9 10x, 6,366) + smartseq2 (1,874) | pseudo-counts | 41 | 767 |
| liu | `liu.qc.h5ad`, qc_pass, RSEM expected counts (27.9 % non-integer values) | 13,100 | 42 | smartseq2 (fresh cells 5,326 / frozen nuclei 3,478 / unannotated 4,296) | round(RSEM) | 5 | 33,270 |
| braun | `braun2023_hindbrain_pons_medulla.h5ad`, UMI counts (X integer-valued) | 100,000 of 286,973 | 17 donors | scRNA (10x v2/v3, 5–14 pcw) | UMI | 3 | 1,015 |
| siletti | 6 `siletti2023_pons_*.h5ad`, UMI counts | 100,000 of 204,168 | 18 donor × dissection (3 donors) | snRNA (10x v3) | UMI | 2 | 2,032 |

**Pseudo-counts for TPM cohorts** (`10_prepare_queries.py::tpm_cohort`). Linear TPM (LaBelle: `10·(2^x − 1)`) is rescaled per cell so that the
cell's total over all harmonised genes equals `n_detected_genes × m`, then rounded; `m` = median (counts / detected genes) of the matching
counts cohort in this project: **33.82** (Liu RSEM Smart-seq2, qc_pass) for Smart-seq2 tables and **2.87** (Jessa 10x scRNA) for the LaBelle
P9 10x CPM table. Median total pseudo-count per cell: Filbin 170,250, LaBelle 2,265. The HVG-level depth of Filbin pseudo-counts (30,632)
matches Liu's real RSEM depth (33,270), which is what the rule intends. These are **not measured counts**; Filbin and LaBelle are
projection-only and are never used to train the reference. `pct_mito = 0` for both (the deposited TPM tables contain no MT- genes).
Gene matching for Braun/Siletti: Ensembl ID of the 2020-A universe first, symbol (`var['Gene']`) second.

**Subsampling** (seed 0): Braun 286,973 → 100,000 stratified by `CellClass × donor_id` (proportional, floor 100 per stratum); Siletti
204,168 → 100,000 stratified by `supercluster_term × donor_id × dissection` (floor 50). `pct_mito` computed from the 37 MT- genes of each
file. Braun `sample_id` = donor; Siletti `sample_id` = donor × dissection (3 donors are the biological replicates).

## 2. Mapping protocol (`11_map_queries.py train`)
`SCANVI.prepare_query_anndata` → `SCANVI.load_query_data(query, models/scanvi, accelerator='cpu')` with the scArches defaults
(`freeze_expression`, `freeze_decoder_first_layer`, `freeze_batchnorm_encoder`, `freeze_classifier` = True; `freeze_batchnorm_decoder`,
`freeze_dropout` = False). `extend_categories=True` adds the query `sample_id` batches and the new `study` / `assay` categories as new one-hot
inputs; only the weights attached to new categories receive non-zero gradients (scvi-tools masks the gradients of the pre-existing columns,
so the "trainable parameter" count in `run.json`, 0.85 M of 2.08 M, over-states what actually moves). Training: 100 epochs,
`plan_kwargs={'weight_decay': 0.0}`, `check_val_every_n_epoch=10`, no early stopping, seed 1, train_size 0.9, all cells `Unknown`
(`n_samples_per_label=None`). Batch size: 128 for queries < 20k cells, 512 for the two 100k queries (not specified in the plan).
Outputs per query in `models/scanvi_query_<q>/` (model.pt, latent_f16.npy, soft_probs_f16.npy, predicted_labels.txt, history.json, run.json).

| query | cells | batch | epochs | wall (s) | s/epoch | peak RSS (GB) | median max prob | predicted lineage counts |
|---|---|---|---|---|---|---|---|---|
| filbin | 3,037 | 128 | 100 | 35 | 0.35 | 0.96 | 1.000 | malignant 2,835; oligodendrocyte 101; myeloid 93; vascular 5; neuron 2; astrocyte 1 |
| labelle | 8,240 | 128 | 100 | 89 | 0.89 | 1.03 | 1.000 | myeloid 7,335; malignant 904; vascular 1 |
| liu | 13,100 | 128 | 100 | 142 | 1.42 | 1.18 | 1.000 | malignant 11,442; myeloid 935; oligodendrocyte 512; neuron 116; vascular 70; astrocyte 25 |
| braun | 100,000 | 512 | 100 | 863 | 8.63 | 2.12 | 1.000 | astrocyte 97,822; myeloid 1,136; vascular 1,042 |
| siletti | 100,000 | 512 | 100 | 1,332 | 13.3 | 2.53 | 1.000 | neuron 51,934; oligodendrocyte 31,231; astrocyte 11,823; myeloid 3,978; vascular 1,034 |

Total training wall 2,461 s (41 min). Siletti's s/epoch is inflated by the Braun sensitivity runs (§5) that shared the CPU.

**Uncertainty** (`11_map_queries.py uncertainty`): (i) `uncertainty_prob = 1 − max scANVI probability`; (ii) `uncertainty_knn` = mean
Euclidean distance of a query cell to its 15 nearest reference cells in `X_scanvi` (pynndescent index on all 311,618 reference cells),
expressed as the percentile of the reference's own 15-NN mean-distance distribution (50,000-cell reference subsample, self excluded;
`results/tables/phase5_reference_knn_distances.csv`). Also stored: `knn_mean_dist`, and **`knn_lineage`** = majority `scanvi_label` of the
same 15 reference neighbours with `knn_vote_frac` and `classifier_knn_agree` (added after the Braun finding in §5; HLCA-style kNN transfer).
**UMAP projection**: the stored scanpy UMAP has no transform model, so a umap-learn UMAP (n_neighbors 30, min_dist 0.5, seed 0) was fitted on a
100,000-cell reference subsample of `X_scanvi` (`results/tables/phase5_umap_reference_fit.csv`) and queries were transformed into it
(`obsm['X_umap_ref']`). These coordinates are not those of `atlas_core.obsm['X_umap']`.
`data/processed/query_mapping.h5ad`: 224,377 cells × 3,000 HVGs (`X = layers['counts']`), 43 obs columns (study/sample/assay, author labels,
Phase-3 calls for the tumour cohorts, predicted_lineage, scanvi_maxprob, uncertainty_prob, knn_mean_dist, uncertainty_knn, knn_lineage,
knn_vote_frac, classifier_knn_agree), obsm `X_scanvi`, `scanvi_probs` (7 columns in `uns['phase5']['scanvi_prob_columns']`), `X_umap_ref`.

## 3. Mapping quality (`results/tables/mapping_uncertainty.csv`, `mapping_summary.csv`)
Per-query medians over cells:

| query | median max prob | median 1 − max prob | median kNN uncertainty | cells with kNN uncertainty > 0.5 | classifier = kNN label | predicted malignant | kNN malignant |
|---|---|---|---|---|---|---|---|
| filbin | 1.000 | 0.000 | 0.999 | 1.000 | 0.978 | 0.933 | 0.916 |
| labelle | 1.000 | 0.000 | 0.906 | 0.988 | 0.820 | 0.110 | 0.110 |
| liu | 1.000 | 0.000 | 0.992 | 0.997 | 0.869 | 0.873 | 0.755 |
| braun | 1.000 | 0.000 | 0.985 | 0.998 | 0.028 | 0.000 | 0.934 |
| siletti | 1.000 | 0.000 | 0.542 | 0.527 | 0.966 | 0.000 | 0.030 |

The two measures disagree by design: the frozen classifier is near-certain for every query (1 − max prob < 0.01 for 98.4 % of all query cells; only a few
Liu plates — pSCG1 18 %, pTG5 9 %, S1D2E3H7 4 % of cells — exceed 0.2), whereas the kNN percentile flags all Smart-seq2 cohorts (median ≥ 0.91) and the
fetal reference (0.985) as sitting outside the reference's density. Siletti (adult 10x snRNA, the same protocol as 73 % of the reference) is
the only query whose cells lie inside the reference distribution (median 0.54; oligodendrocytes 0.06, microglia 0.35, astrocytes 0.55,
neurons 0.57–0.80). The `fraction > 0.5` statistic therefore saturates for the Smart-seq2 cohorts and should be read together with the
agreement columns, not as a per-cell rejection rule.

**Liu by group** (`mapping_liu_by_group.csv`): fresh cells (5,326 cells, 18 plates) — author-malignant → malignant 0.995, author-normal →
matching lineage 0.988, classifier/kNN agreement 0.982; frozen nuclei (3,478, 24 plates) — 0.990 / 0.691 / 0.900; unannotated cells (4,296)
— 0.890 predicted malignant, agreement 0.703. Adult (3,507, 12 samples) 0.998 / 0.706; paediatric (9,593, 30) 0.991 / 0.759. Pontine
0.994 / 0.786, thalamic 0.997 / 0.748, spinal (541 cells, 2 samples; 37 author-normal cells) 0.952 / 0.189. Nuclei plates lose normal-type
accuracy, not malignant recall.

## 4. Label transfer vs author / Phase-3 labels (`mapping_<query>_author_vs_predicted.csv`, `mapping_author_agreement.json`)
- **Liu** (author states): 7,335 author-malignant cells → predicted malignant **0.993** (kNN 0.960); per state OPC-like-1 0.999, OPC-like-2
  1.000, OPC-like-3 1.000, AC-like 0.985, MES-like 0.993, OC-like 0.957, S 0.998, G2M 0.998, frozen-OPC-like-b 0.993. 1,370 author-normal
  cells → matching lineage **0.747** (kNN 0.808): Microglia 0.971, Oligodendrocyte 0.886, Myeloid 0.786, Endothelia 0.377 (kNN 0.508),
  T cell **0.000** (kNN 0.542). Author-normal cells predicted malignant: 0.170 (mostly Endothelia and Myeloid from frozen-nuclei plates).
- **Filbin** (no author per-cell states in the deposit; Phase-3 consensus used): 2,325 consensus-malignant → malignant **0.997** (kNN 0.978);
  74 consensus-normal → matching type 0.959; 638 consensus-ambiguous → 0.806 malignant; CNV-malignant (tier 2) → 0.958.
- **LaBelle** (immune author labels only): 6,366 author immune cells → matching lineage **0.827** by the classifier (Myeloid 1.000, Tcell
  0.000: all 1,104 T cells are called myeloid) but **0.981** by kNN (Myeloid 0.997, Tcell 0.908). 752 consensus-malignant tumour-table cells →
  malignant 1.000; 74 consensus-malignant cells from the immune SS2 table → myeloid.
The T-cell failure repeats the Phase-4 finding (T_NK recall 0.467 in the reference; 764 cells): the frozen classifier's T_NK decision region
is too small, while the latent neighbourhood still places T cells next to reference T_NK cells. `knn_lineage` is the label to use for immune
sub-typing until a dedicated immune model exists.

## 5. Specificity check on normal references (`mapping_normal_reference_specificity.csv`)
Classes with a reference counterpart (expected lineage in the table): Braun Oligo→oligodendrocyte, Vascular/Fibroblast→vascular, Immune→myeloid,
Neuron/Neuroblast/Neuronal IPC→neuron; Siletti Oligodendrocyte/COP→oligodendrocyte, Astrocyte/Bergmann/Ependymal→astrocyte, Microglia→myeloid,
Vascular/Fibroblast→vascular, all neuronal superclusters→neuron. No counterpart: Braun Radial glia, Glioblast, Erythrocyte, Placodes,
Neural crest; Siletti Oligodendrocyte precursor, Choroid plexus, Splatter.

**Siletti 2023 adult pons** — clean. Classifier: false-malignant (max prob ≥ 0.8) **0.000** (0 of 100,000 cells predicted malignant);
expected lineage recovered for **0.990** of the 78,137 cells with a counterpart (Oligodendrocyte 1.000, Astrocyte 1.000, Microglia 1.000,
Vascular 0.974, neuronal superclusters 0.996–1.000; Fibroblast 0.536 — the rest astrocyte 0.463; COP 0.374 oligodendrocyte / 0.568 astrocyte). kNN label: malignant for 3.0 % of
all nuclei (1.5 % with vote ≥ 0.8), 0.2 % within counterpart classes; the kNN-malignant nuclei are the adult **OPCs (48.1 %, 24.3 % with
vote ≥ 0.8)** and committed OPCs (51.0 %), i.e. the lineage without a normal counterpart in the tumour reference and the putative cell of
origin — a biologically expected, not spurious, neighbourhood.

**Braun 2023 fetal hindbrain (5–14 pcw)** — the classifier fails, the latent behaves as expected. Classifier: 0 cells predicted malignant
(false-malignant 0.000) but 97.8 % predicted **astrocyte** with max prob 1.0, including 31,639 neurons, 23,401 neuroblasts, 24,330 radial glia
and 11,479 glioblasts; only Immune (0.947 myeloid) and Vascular (0.615 vascular) are recovered, giving expected-lineage recovery of **0.035**.
kNN label: 93.4 % of fetal cells have a malignant majority among their 15 nearest reference cells (90.2 % with vote ≥ 0.8) — Radial glia 1.000,
Glioblast 1.000, Neuroblast 0.989, Neuron 0.927, Neuronal IPC 1.000, Oligo 0.960 — while Immune (0.990 myeloid) and Vascular (0.994) map
correctly. Sensitivity runs on a 20,000-cell Braun subsample (not part of the protocol; results printed to the session, not saved as models):
with **0 training epochs** the surgery model predicts Radial glia / Glioblast / Neuroblast → malignant (0.99–1.00) and Neuron → malignant 0.65 /
neuron 0.35; after **10 and 30 epochs** all of these are malignant; after 100 epochs on the full 100k they flip to astrocyte. The latent
distance to the reference does not change (median 15-NN distance 2.60 untrained vs 2.66 trained, reference median 1.52), so the flip is the
frozen classifier extrapolating in an off-manifold region, not a change of neighbourhood. Interpretation: (i) the reference has **no normal
progenitor or immature-neuron class** — its normal lineages are mature, post-natal types — so fetal progenitors/neuroblasts can only land in
the malignant (OPC-/AC-like) region, which is the expected behaviour of a tumour reference and the basis of the cell-of-origin analysis below;
(ii) scANVI's classifier probabilities are not a usable uncertainty for out-of-distribution queries; `uncertainty_knn` (0.985 for Braun) and
`knn_lineage` are. For Braun, `predicted_lineage` should not be used.

## 6. Cell-of-origin view (`malignant_nearest_fetal_class.csv`, `malignant_nearest_fetal_class_per_cell.csv.gz`)
kNN transfer (k = 15, pynndescent index on the 100,000 mapped Braun cells, queried with the 169,571 consensus-malignant reference cells in
`X_scanvi`). Braun `CellClass` of the majority; `cell_type` and `Region` majorities also stored. Note that Braun's "Oligo" class at 5–14 pcw is the
oligodendroglial lineage (OPC / pre-OPC / COP stages), not mature oligodendrocytes.

Overall: **Glioblast 0.509, Oligo 0.359, Radial glia 0.107**, Neural crest 0.016 (median 15-NN distance 2.82).

| Jessa author state (consensus-malignant cells) | n | top fetal class | second |
|---|---|---|---|
| OPC | 31,164 | Oligo 0.628 | Glioblast 0.327 |
| Proliferating OPC | 313 | Oligo 0.760 | Radial glia 0.131 |
| Oligodendrocytes | 2,689 | Oligo 0.880 | Glioblast 0.080 |
| Astrocytes | 7,965 | Glioblast 0.783 | Radial glia 0.194 |
| Ependymal | 1,627 | Glioblast 0.701 | Radial glia 0.294 |
| Glial progenitors | 1,821 | Glioblast 0.506 | Radial glia 0.302 |
| RGC | 230 | Radial glia 0.744 | Glioblast 0.170 |
| Neuronal progenitors | 35 | Radial glia 0.686 | Glioblast 0.143 |
| Neurons | 457 | Glioblast 0.466 | Radial glia 0.320 |
| Immune (author) | 74 | Glioblast 0.311 | Immune 0.297 |
| Vascular & other | 160 | Radial glia 0.550 | Glioblast 0.263 |

The ordering is the one expected from the OPC-like / AC-like dichotomy: OPC-like states neighbour the fetal oligodendroglial lineage,
astrocyte-/ependymal-like states neighbour glioblasts, and the rare RGC/neuronal-progenitor states neighbour radial glia. Per-Leiden-cluster
rows (28 clusters, `leiden_scanvi_r1`) are in the same table. Caveat: the transfer uses the surgery-trained Braun latent, which §5 shows to be
off-manifold (median query–reference distance 2.7 vs 1.5 within the reference); the ranking of fetal classes is informative, the absolute
fractions are not calibrated.

## 7. Packaging (`scripts/12_package_v1.py`)
`docs/atlas_schema.md` (145 atlas_core obs fields + obsm/var/layers/uns; 43 query_mapping obs fields; generated from the files),
`docs/label_ontology.csv` (20 rows: lineage labels with marker genes and source tier, consensus tiers, query author vocabularies),
`results/tables/atlas_v1_summary.csv` (6 reference studies + 5 queries + total: reference 311,618 cells, 57 samples, 52 patients;
consensus-malignant fraction 0.544, scANVI-malignant 0.772), README section "Atlas v1 contents" appended (md5 of both h5ad files included),
figures `results/figures/phase5_query_mapping.png` (reference UMAP refit in grey, query cells coloured by predicted lineage; the Braun panel
shows the §5 astrocyte collapse) and `phase5_mapping_uncertainty.png` (per-sample violins of both uncertainty measures; dotted line 0.5).
Figure rcParams reproduce the `figure-style` skill mechanics (`src/dmg_atlas/figstyle.py`); no titles inside panels.

## 8. Deviations and caveats
1. Filbin has no per-cell author states in the GEO deposit (Phase 2 notes); its crosstab and agreement use Phase-3 `consensus_call` / `cnv_call`.
   LaBelle has immune labels only; its tumour SS2 cells are also judged against the consensus call.
2. Braun subsampled 286,973 → 100,000 and Siletti 204,168 → 100,000 (stratified; §1), as the plan allowed. Specificity and cell-of-origin numbers
   rest on these subsamples.
3. TPM cohorts mapped with pseudo-counts (§1 rule; `m` calibrated on Liu RSEM and Jessa 10x). `pct_mito = 0` for Filbin/LaBelle.
4. `knn_lineage`, `knn_vote_frac`, `classifier_knn_agree` were added beyond the plan after the Braun classifier collapse; `predicted_lineage` remains
   the protocol-specified scANVI prediction and is the one summarised in `mapping_summary.csv` (kNN fractions sit beside it).
5. Braun epoch-sensitivity runs (0/10/30 epochs, 20k subsample) were exploratory and are reported in §5 only; no model files were written for them.
6. Batch size (128 / 512) was not specified; the reference used 1,024. Training ran on CPU with 8 threads; Modal was not used.
7. UMAP refit on a 100,000-cell reference subsample (umap-learn) because the stored scanpy UMAP has no transform model; `X_umap_ref` ≠ `X_umap`.
8. kNN uncertainty uses an approximate index (pynndescent) and a 50,000-cell reference distance distribution; the percentile saturates at ≈1 for
   the Smart-seq2 cohorts (§3).
9. `nice -n 10` is not permitted in the sandbox; thread cap 8 applied. Siletti's s/epoch (13.3) was inflated by concurrent sensitivity runs.
10. `pyarrow` installed into `multiverse` to unpickle the reference `model.pt`.
11. Liu `loc` axis mixes author locations (pontine/thalamic/spinal) with manifest sites (pons/thalamus) for the 4,296 unannotated cells.


**Correction (2026-10-19, lead-agent review of `mapping_normal_reference_specificity.csv`):** the 0.996–1.000 expected-lineage recovery holds for the named neuronal superclusters (lower/upper rhombic lip, midbrain-derived inhibitory, cerebellar inhibitory, CGE interneurons; n = 25–23,748) but not for Siletti's 'Miscellaneous' neuron catch-all (378 nuclei, 0.034 expected lineage, 0 malignant neighbours). Elevated malignant kNN-neighbour fractions are seen in both adult oligodendrocyte precursors (5,792 nuclei, 0.48; no reference counterpart) and committed oligodendrocyte precursors (155 nuclei, 0.51; 0.39 expected lineage), not only in OPCs.