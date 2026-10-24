# Audit — Phases 4–5 (integration, scANVI reference, query mapping, v1 packaging)

Auditor: independent adversarial review, 2026-10-19. Scope: `scripts/07_assemble_atlas_input.py`, `07b_make_training_files.py`, `08_integration_metrics.py`,
`09_assemble_atlas_core.py`, `10_prepare_queries.py`, `11_map_queries.py`, `11b_mapping_summaries.py`, `12_package_v1.py`; `modal/phase4_train.py`;
`models/*` (attributes, state dicts, run.json, history); `data/processed/atlas_core.h5ad` (obs/obsm/uns, backed), `atlas_input_hvg.h5ad` (obs/var/uns),
`atlas_input_hvg.train.h5ad` (var), `query_mapping.h5ad` (obs/obsm), `queries/filbin.hvg.h5ad`, `data/interim/filbin.qc.h5ad`, `liu.qc.h5ad` (obs);
tables `integration_metrics*.csv`, `integration_model_choice.json`, `atlas_input_hvg.csv`, `atlas_malignant_hvg.csv`, `atlas_core_obs_labels.csv`,
`atlas_core_umap.csv` (existence/shape only), `mapping_*.csv/json`, `query_preparation.csv`, `phase5_reference_knn_distances.csv`,
`phase5_umap_reference_fit.csv` (existence only), `malignant_nearest_fetal_class*.csv`, `atlas_v1_summary.csv`; docs `phase4_integration_notes.md`,
`phase5_mapping_notes.md`, `atlas_schema.md`, `label_ontology.csv`, `decisions.md`, `build_plan.md`, `README.md`; logs `phase4_integration.log`,
`phase4_modal_runs_*.jsonl`, `phase5_mapping.log`; figures `phase4_umap_overview.png`, `phase5_query_mapping.png`, `phase5_mapping_uncertainty.png`
(all three viewed). Findings table: `docs/audit/integration_mapping_findings.csv` (24 rows: 1 critical, 3 major,
13 minor, 7 notes). Supplementary: `integration_mapping_audit_braun_reference_classifier_crosstab.csv`,
`integration_mapping_audit_siletti_reference_classifier_crosstab.csv` (reference-classifier predictions per author class).

## Method of verification
1. Line-by-line reading of every script and module in scope against the notes; comparison with `build_plan.md`/`decisions.md` for undeclared deviations.
2. Model forensics: `torch.load` of every `model.pt` in scope; comparison of reference `models/scanvi` vs each `models/scanvi_query_<q>` state dict
   (per-tensor max |difference|); inspection of `init_params_`, registry, train/validation indices and training history; reading of the installed
   scvi-tools 1.4.2 source for `ArchesMixin._set_params_online_update`, `FCLayers.set_online_update_hooks` and `SCANVAE.__init__`.
3. Targeted recomputation (local `multiverse`, 8 threads, RSS ≤ 4.1 GB): regeneration of the seeded 60,023-cell scIB subsample and five metrics for three
   embeddings with scib-metrics 0.5.10; scANVI agreement split by the stored train/validation indices; scaled scores with and without the e60 row;
   reference-classifier predictions on the saved query latents (classifier rebuilt from the reference state dict; its query-model counterpart
   reproduces the saved labels 100 %); a live `SCANVI.load_query_data` on Filbin without training (latent and predictions); Filbin pseudo-counts
   recomputed from TPM; all per-query uncertainty/agreement medians; atlas_v1_summary vs atlas_core obs; md5 of both h5ad files.
4. All three figures viewed at full resolution.
No file outside `docs/audit/` was written; no model was retrained; no download or GPU/Modal use.

## Re-derived numbers (claimed vs recomputed)
| quantity | claimed (notes/tables) | recomputed | status |
|---|---|---|---|
| Reference cells / samples / patients / genes | 311,618 / 57 / 52 / 18,882 | 311,618 / 57 / 52 / 18,882 (atlas_input_hvg obs, atlas_core obs, layers['counts'] shape) | match |
| lineage_label counts | malignant 169,571; Unknown 91,309; myeloid 23,052; oligo 14,141; astro 5,045; neuron 4,438; vascular 3,298; T_NK 764 | identical | match |
| HVG nbatches min / median; overlap with median-rank list | 21 / 37; 1,465/3,000 | 21 / 37.0; 1,465/3,000 | match |
| Malignant HVG nbatches min / median; span | 19 / 32; 0.5 | 19 / 32.0; 0.5 (log: span 0.3 loess singularity) | match |
| scIB subsample size (lineage / fine) | 60,023 / 30,018 | 60,023 / 30,014 (json + log) | fine subsample mismatch (typo) |
| asw_label scvi_seed1 / harmony / pca | 0.595 / 0.566 / 0.543 | 0.5953 / 0.5660 / 0.5426 | match |
| iLISI study scvi_seed1 / harmony / pca | 0.053 / 0.177 / 0.0 | 0.0534 / 0.1773 / 0.0 | match |
| graph connectivity scvi_seed1 / harmony / pca | 0.969 / 0.877 / 0.969 | 0.9692 / 0.8765 / 0.9687 | match |
| asw_batch study scvi_seed1 / harmony / pca | 0.876 / 0.871 / 0.849 | 0.8765 / 0.8707 / 0.8495 | match |
| total score seed1 / seed0 / seed2 / harmony (6-row scaling) | 0.739 / 0.702 / 0.704 / 0.581 | 0.739 / 0.702 / 0.704 / 0.581; without e60: 0.766 / 0.732 / 0.732 / 0.598; chosen seed1 either way | match; scaled values depend on e60 inclusion |
| scANVI agreement on labelled cells | 0.987 | 0.9873 all labelled; 0.9881 train split (198,279); 0.9804 validation split (22,030) | match, but 90 % training data |
| scANVI per-class recall astrocyte / T_NK | 0.815 / 0.467 | all-labelled 0.815 / 0.467; validation 0.663 / 0.388 | over-stated on held-out |
| Unknown cells with max prob >= 0.8; predicted malignant | 96.2 %; 70,924 | 96.15 %; 70,924 (myeloid 10,218, astro 4,020, oligo 3,522, vascular 1,435, neuron 1,092, T_NK 98) | match |
| atlas_core: NaN rows in X_scvi_malignant; Leiden clusters; X_scvi source | 142,047; 28; scvi_seed1 | 142,047; 28; allclose to models/scvi_seed1/latent_f16.npy | match |
| md5 atlas_core / query_mapping | ded17c5f...; 914eb099... | ded17c5f9f9d5cfa2dab16c02bb64b92; 914eb0996f43853c681c1220b2398194 | match |
| atlas_v1_summary: consensus / scANVI malignant fraction; per-study rows | 0.544 / 0.772 | 0.5442 / 0.7724; all 6 study rows identical to atlas_core obs | match |
| Filbin pseudo-counts (m, median total, HVG median) | 33.82; 170,250; 30,632 | 33.817; 170,250; 30,632; HVG matrix bit-identical to recomputation | match |
| Query cells per cohort | 3,037 / 8,240 / 13,100 / 100,000 / 100,000 | same (query_mapping.h5ad obs) | match |
| Median kNN uncertainty filbin / labelle / liu / braun / siletti | 0.999 / 0.906 / 0.992 / 0.985 / 0.542 | 0.999 / 0.906 / 0.992 / 0.985 / 0.542 | match |
| classifier = kNN agreement per query | 0.978 / 0.820 / 0.869 / 0.028 / 0.966 | 0.978 / 0.820 / 0.869 / 0.028 / 0.966 | match |
| Cells with 1 - max prob < 0.01 | 98.4 % | 98.42 % | match |
| Reference self 15-NN median distance | 1.52 | 1.252 (phase5_reference_knn_distances.csv) | mismatch |
| Braun -> reference median 15-NN distance | 2.60 untrained / 2.66 trained (20k subsample) | 2.25 (saved knn_mean_dist, 100k) | mismatch; exploratory values unsaved |
| Malignant -> Braun kNN: Glioblast / Oligo / Radial glia; median distance | 0.509 / 0.359 / 0.107; 2.82 | 0.5087 / 0.3586 / 0.1068; 2.82; enrichment vs class frequency 4.4x / 31.6x / 0.44x | match (base rates not reported) |
| Braun predicted lineage (saved vs reference classifier on same latents) | astrocyte 97,822; malignant 0 | reference classifier: malignant 82,910; neuron 13,137; vascular 2,105; myeloid 1,338; astrocyte 385 | artifact (IM-01) |
| Braun expected-lineage recovery | 0.035 | 0.259 with the reference classifier; 0.093 kNN | artifact |
| Siletti predicted malignant | 0 of 100,000 | 3,541 (3.5 %) with the reference classifier, 3,360 of them OPCs; kNN 2,970 | artifact |
| LaBelle author T cells -> T_NK (classifier) | 0.000 | 0.559 with the reference classifier (0.908 kNN); immune agreement 0.923 vs 0.827 | artifact |
| Liu author-normal -> matching lineage; T cell | 0.747; 0.000 | 0.800; 0.333 with the reference classifier (kNN 0.808; 0.542) | artifact |
| Liu author-malignant -> malignant | 0.993 | 0.956 with the reference classifier (kNN 0.960) | artifact (the drifted classifier flattered this number) |
| Filbin query latent: untrained surgery model vs saved after 100 epochs | (implicitly: training adapts the latent) | max |diff| 0.0039 (float16 rounding), r = 1.000000; encoder state dict identical for all 5 queries | latent is zero-shot (IM-02) |

## Findings by severity

### Critical
**IM-01 — The query classifier was not frozen; every query-side prediction table and the §5 narrative rest on per-query, unintentionally re-fitted classifiers.**
`SCANVI.load_query_data(..., freeze_classifier=True)` in scvi-tools 1.4.2 sets `requires_grad=True` on all classifier parameters because the test
`'.0.' in key` matches the module prefix `classifier.classifier.0.`; gradient hooks zero the Linear weights but the BatchNorm affine parameters
(2 × 128 γ + 2 × 128 β) have no hook and train under the unlabelled scANVI objective. State-dict diffs confirm it for all five queries (Braun: max |Δγ| 0.99,
|Δβ| 1.36) while the encoder is byte-identical. Re-applying the *reference* classifier to the identical saved latents reverses the headline results:
Braun 82.9 % malignant and 13.1 % neuron (saved: 97.8 % astrocyte, 0 % malignant; expected-lineage recovery 0.259 vs 0.035); Siletti 3.5 % malignant,
58 % of adult OPCs (saved 0 %); LaBelle T cells 0.559 recovered (saved 0.000); Liu T cells 0.333 (saved 0.000), Liu normal→match 0.800 (0.747),
Liu malignant→malignant 0.956 (0.993). The Braun "astrocyte collapse", the "0 of 100,000 false-malignant" specificity claim and the
"frozen classifier T-cell failure" are artifacts of this bug, and the notes' own 0-epoch sensitivity run (progenitors → malignant) was the correct
behaviour. Affected deliverables: `query_mapping.h5ad` (`predicted_lineage`, `scanvi_maxprob`, `uncertainty_prob`, `scanvi_probs`), `mapping_summary.csv`,
`mapping_<q>_author_vs_predicted.csv`, `mapping_normal_reference_specificity.csv`, `mapping_author_agreement.json`, query rows of `atlas_v1_summary.csv`,
both phase-5 figures, notes §2–5, README limitations. `knn_lineage`, `uncertainty_knn`, `X_scanvi`, the cell-of-origin kNN and all Phase-4 outputs are unaffected.

### Major
**IM-02 — Query "scArches surgery" cannot adapt the latent.** The reference was trained with `encode_covariates=False` (default), so the encoder has no
batch/covariate inputs (`z_encoder.encoder.n_cat_list == []`) and the online-update hook zeroes its entire gradient. The untrained surgery model gives the
same Filbin latent as the saved 100-epoch model (max |diff| 0.0039 = float16 rounding). Query `X_scanvi` is therefore a zero-shot projection by the
reference encoder; no query batch effect is removed, which contributes to kNN-uncertainty ≈ 1 for all non-10x cohorts. The 41 min of training only
re-fitted decoder columns and (IM-01) the classifier. The notes' description of the method and the plan's "fine-tune with frozen reference weights"
are not what the deliverable is.

**IM-03 — The scANVI agreement 0.987 is 90 % training-set accuracy.** Train split 0.9881 (198,279 labelled cells), validation split 0.9804 (22,030);
validation recall astrocyte 0.663 (reported 0.815), T_NK 0.388 (0.467). Validation accuracy peaked at epoch 8 (0.986) and fell to 0.977 at epoch 24 while
validation classification loss rose 0.048 → 0.095: the fixed 25-epoch classifier over-fits. The split is by cell, not by patient, so generalisation to new
tumours is unmeasured.

**IM-04 — Reporting built on IM-01/IM-02.** §5's mechanistic interpretation ("frozen classifier extrapolating off-manifold"), §6's "surgery-trained Braun
latent", the README limitation sentence and the Braun panel of `phase5_query_mapping.png` (astrocyte-coloured cells over the malignant region) describe
the artifact. The 2026-10-19 correction addresses only the Siletti 'Miscellaneous' row.

### Minor
IM-05 scaled scIB scores were normalised over six rows including the ineligible e60 sensitivity run (seed1 0.766 → 0.739; ranking unchanged; first run at
19:08 already chose seed1 before e60 existed; 0.4/0.6 weights not pre-registered in decisions.md but outcome is weight-invariant). IM-07 study/assay
covariates are nested in sample_id (redundant, not harmful). IM-08 scANVI label set (single `malignant`, no substates, no endothelial/pericyte split)
deviates from the plan without being declared. IM-09 §5/§6 distance numbers (1.52, 2.60/2.66) do not reproduce from saved files (1.25, 2.25).
IM-10 cell-of-origin fractions are not base-rate corrected (Oligo 31.6×, Glioblast 4.4× enriched; Radial glia 0.44× depleted); the ranking is supportable
only as an off-manifold nearest-neighbour statement. IM-12 pseudo-count rule justified heuristically; no m-sensitivity; TPM has no gene-length
dependence unlike RSEM/UMI. IM-13/14/15 figure house-rule issues (panel headers, in-figure description, uninformative 1 − max-prob row, palette reuse, LaBelle
year). IM-16 provenance inaccuracies (torch 2.14.1 vs 2.10.0; hard-coded freeze dict; unsaved sensitivity runs). IM-22/23 no uncertainty on pooled
cell-level proportions or on the seed comparison; 30,014 vs 30,018 typo. IM-24 Liu nuclei carry assay 'smartseq2'.

### Confirmed sound (notes)
IM-06 five scIB metrics for three embeddings reproduce to 3 decimals on the regenerated subsample. IM-11 Filbin pseudo-counts reproduce bit-exactly.
IM-17 atlas_core/atlas_v1_summary/README md5 fully consistent. IM-18 HVG tables as described, lineage markers present. IM-19 assembly code correct, no
sample-id collisions, 4 Jessa patients correctly kept as separate scRNA/multiome batches. IM-20 scVI under-training acknowledged and recorded. IM-21
all per-query uncertainty/agreement medians in the notes reproduce from `query_mapping.h5ad`.

## Overall verdict for this scope
**Phase 4 (reference): sound with corrections** — the cohort, HVGs, scVI models, scIB comparison, model choice and atlas_core packaging are reproducible
and internally consistent; the scANVI accuracy must be re-reported as held-out/patient-held-out and the label-set deviation declared.
**Phase 5 (query mapping): not yet defensible** — the classifier-derived outputs (`predicted_lineage` and everything built on them) are artifacts of an
unintended classifier update, and the mapping is a zero-shot projection rather than the described scArches adaptation. The latent-based outputs
(`X_scanvi`, `knn_lineage`, `uncertainty_knn`, cell-of-origin kNN) and the query preparation are sound and can be kept; the tables, figures and notes
§2–6 need regeneration with the reference classifier and rewording.

## Not verified and why
- NMI/ARI (Leiden with resolution optimisation) and kBET/PCR were not recomputed (≈ 26 min per full pass single-threaded; five other metrics were).
- Harmony/PCA control embeddings were taken as saved; the Modal training runs themselves were not re-executed (no GPU/Modal per brief); run.json/history
  were cross-checked against the notes instead.
- Pre-specification of the 0.4/0.6 rule cannot be established (no git repository in the project); the log sequence shows the choice preceded the e60 run.
- Braun epoch-sensitivity runs (0/10/30 epochs) exist only in the original session; IM-01/IM-02 explain their reported behaviour without re-running them.
- `X_umap_ref` / umap-learn transform and `phase5_umap_reference_fit.csv` were not recomputed (stochastic; not conclusion-bearing).
- Siletti/Braun pseudo-count-free preparation (UMI reads via h5py) was reviewed by code only; cell counts and subsample sizes were checked.
