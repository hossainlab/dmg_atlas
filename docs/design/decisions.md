# Phase 0 — design decisions (agreed 2026-10-17)

| # | Decision | Choice | Rationale |
|---|---|---|---|
| 1 | Scope | Human data only | Mouse IUE cells and mouse Visium sections in the Collot 2026 deposit are removed at ingestion; no mouse reference. |
| 2 | Priority | Atlas construction first, biological discovery second | Phases 1–4 and 6 are the v1 deliverable; Phases 5, 7, 8 follow. |
| 3 | Liu 2022 (GSE184357) | Projection-only scANVI query | Smart-seq2 plate data (9,248 cells; 18 fresh scRNA + 25 frozen snRNA samples), non-integer deposited matrix. Keeps the reference 10x-homogeneous. Its adult/pediatric and pontine/thalamic/spinal labels are used in downstream tests. |
| 4 | Normal references | Braun 2023 first-trimester hindbrain + Siletti 2023 adult pons (cellxgene) | Developmental and adult references in standard 10x chemistry. Xie 2026 Paired-Tag dropped (adult-only, non-10x, 20 GB). |
| 5 | Reference models | (a) tumour + microenvironment scVI/scANVI reference; (b) malignant-only scVI | (a) supports query mapping and composition for every cell type; (b) resolves malignant states without TME variance. |
| 6 | Okada 2026 postmortem cohort | Attempt: verify UCSC Cell Browser matrix download in Phase 1; include only if counts + sample metadata complete | Only treated-vs-untreated / biopsy-vs-autopsy axis available openly. |
| 7 | Compute | Mac (M4 Pro, 24 GB) for single runs, thread-capped at 8, `nice -n 10`; Modal A10G (image im-gl8KES0WsH7PoFjf83hK6v) for sweeps, cell2location, scGPT | Measured: scVI 30k×3k = 3.6 s/epoch Mac (12 thr), 1.72 s/epoch A10G. |
| 8 | Label tiers | author label; own CNV call; agreement | Validation statistic: agreement with Jessa 2022 consensus (143,611 malignant / 33,301 normal). |
| 9 | Replication rule | a metaprogram is replicated if recovered in ≥3 studies and in both cells and nuclei | Pre-registered before Phase 5. |
| 10 | Project location | /Users/jubayer/Projects/dmg_atlas | All reproducible code, data, figures, tables and docs live here. |


**Data volume correction (2026-10-19):** the approved plan estimated ~5 GB of core deposits; actual Phase 1 downloads were 8.71 GB (tumour cohorts) + 20.98 GB (spatial, references, bulk) = 29.7 GB, driven by the full 11.93 GB Braun 2023 h5ad and the 2.08 GB OpenPedCan SNV file.

## 2026-10-19 — Phase 5 audit remediation (IM-01 … IM-24)
- **Query labels come from the reference classifier applied to the zero-shot latent.** The per-query `models/scanvi_query_*` models are retained for provenance only; their classifier BatchNorm drifted during 'surgery' (scvi-tools 1.4.2 `'.0.' in key` bug) and their `predicted_labels.txt` must not be used. Canonical labels: `query_mapping.h5ad` obs `predicted_lineage` (classifier) and `knn_lineage` (primary for off-manifold queries); previous values kept as `*_pre_audit`.
- **Phase 5 is described as zero-shot encoder projection**, not scArches surgery, because the reference was trained with `encode_covariates=False`. Adaptive mapping, if ever wanted, requires a retrained reference (`encode_covariates=True`, layer norm) — not done in v1.
- **scANVI accuracy is reported as held-out validation accuracy (0.980)** with per-class validation recall; a patient-grouped bound from light classifiers on `X_scvi` (≈0.965) accompanies it, labelled as such. Any retraining uses early stopping on validation classification loss. *[superseded 2026-10-19 by IM-03 / the wave-2 retraining — the reference was retrained with a patient-grouped split, so the figure to quote is the **held-out-patient accuracy 0.9751** (balanced 0.9166) on 23,827 labelled cells from 5 patients, not a held-out-cell accuracy and not the `X_scvi` bound; see the wave-2 block below and `results/tables/reflabels/scanvi_heldout_comparison.csv`]*
- scIB model-choice rule recorded post hoc: total = 0.4 × batch + 0.6 × bio over min-max-scaled metrics; scaling should run over eligible candidates only (the v1 table included the e60 sensitivity row; ranking unaffected).
- Regenerated deliverables keep the previous version as `<name>.pre_audit.<ext>`; notes carry a dated 'Audit corrections' section.

## Audit corrections — metaprogram track (2026-10-19)
- Malignant selection for Phase 6 (§1, 2, 6, 8, 9) = corrected consensus_call == 'malignant' AND qc_pass, okada__1734_pons excluded; scANVI rescue dropped (MC-04). Sensitivity set: consensus_call_eventtrust. Scripted in `scripts/13a_malignant_selection.py`.
- Per-sample cNMF re-run only where the malignant set changed by > 5 % (46 of 50 samples); the four ≤ 5 % samples reuse the 2026-10-17 programmes (`results/tables/cnmf/sample_rerun_plan.csv`).
- Metaprogram families: within-sample z-scoring, mean over member MPs (MC-01); technical programmes (ribosomal, mitochondrial, IEG/dissociation, immune) classified by gene-set rules before lineage labelling and excluded from lineage counts (MC-07, P7-04); contamination exclusion keeps the pre-registered 2× rule, 1.5× reported as sensitivity (MC-08).
- Reference gene sets tiered and scripted (published Liu Table S2a / curated > in-house DE > Filbin-derived untraceable) (MC-03).
- Replication rule unchanged (decisions #9); `replicated_strict` (≥ 2 samples per assay class) and `n_patients` reported alongside; cut sensitivity 0.75–0.92 in 0.01 steps (MC-02, MC-05, MC-17).
- H3.1 vs H3.3 reported with site adjustment and a pons-only fit (MC-06); malignant states standardised within sample (MC-23).

## 2026-10-19 — Reference-label remediation (wave 2)
- **okada__1734_pons stays in the scVI/scANVI reference but is low-confidence.** It fails the sample-level QC gate (DQ-03; 1,163 genes) and 64.5 % of
  its reference cells have score-only CNV calls (CNV-09). It is kept in the latent (already part of the selected `scvi_seed1` model; removal would
  force re-selection), flagged `obs['low_confidence_sample'] = True`, **excluded from `models/scvi_malignant`** (3,545 cells) and to be excluded
  from Phase 6 malignant analyses when that track re-runs.
- **scANVI is validated on held-out patients with early stopping.** `models/scanvi` = retrain with corrected labels, 5 held-out patients (10.7 % of
  labelled cells), early stopping on validation classification loss (best epoch 7/12); accuracy is reported as held-out-patient accuracy (0.975)
  and balanced accuracy (0.917), never as all-cell agreement. The previous model is `models/scanvi_pre_audit`.
- **Query mapping is zero-shot against the retrained reference**; `models/scanvi_query_<q>` now hold zero-shot outputs (no training), the
  surgery-trained directories are `models/scanvi_query_<q>_pre_audit`.
- Wave-2 regeneration keeps the wave-1 state of Phase 5 tables/figures in `results/tables/phase5_wave1_snapshot/` and
  `results/figures/phase5_wave1_snapshot/` in addition to the `*.pre_audit.*` originals.
