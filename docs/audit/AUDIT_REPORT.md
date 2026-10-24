# Independent scientific audit — DMG single-cell atlas (Phases 1–7)

Date: 2026-10-19. Scope: every script (40), module (9), Modal job script, results table (93), figure (18 + publication/supplementary set), model directory and phase note under `/Users/jubayer/Projects/dmg_atlas`. Method: six independent auditors with no prior involvement, each reading code line by line, re-deriving headline numbers from the underlying files, re-executing targeted steps, and viewing every figure in scope. Findings: 163 (4 critical, 33 major, 71 minor, 55 confirmed-sound notes); 84% were verified by recomputation (100% of critical, 88% of major). Per-track reports: `docs/audit/*_audit.md`; full register: `docs/audit/ALL_findings.csv`.

## 1. Verdict by phase

| Phase | Verdict | One-line reason |
|---|---|---|
| 1–2 Acquisition, manifest, QC | **Sound with corrections** | All headline and per-sample QC numbers reproduce; Okada count recovery exact; gaps are reproducibility (unscripted label tables and three downloads) and three QC edge cases. |
| 3 Malignant calling | **Sound with corrections** | Core CNV arithmetic and the 342,213-cell accounting reproduce; one bug falsifies chr17p frequencies; arm-level and 4q12 events are artefacts of chromosome-wide windows; validation is over-stated (tuning set larger than declared, tier-3 agreement computed on a self-selected half). |
| 4 Integration | **Sound with corrections** | scIB metrics reproduce; model choice stands. scANVI "0.987" is training-set accuracy (held-out 0.980; small-class recall 0.08–0.15 lower). |
| 5 Query mapping | **Not yet defensible** | The query classifier was not frozen (scvi-tools key-match bug): all `predicted_lineage`, probabilities, specificity tables and the Braun/Siletti narrative are artefacts. The latent projections themselves are valid zero-shot encoder outputs; the cell-of-origin kNN survives. |
| 6 Metaprograms + composition | **Sound with corrections** | cNMF, clustering, the replication matrix and composition arithmetic reproduce. The "2 replicated" headline flips at a 0.01 change of the cut; family assignment used downstream is assay-determined; two "lineage" MPs are technical programmes. |
| 6 Spatial | **Sound with corrections** | All numbers reproduce, but co-localisation on raw abundances is density-confounded: MES-like–myeloid survives at about half the effect; OPC-like–oligodendroglial does not; two figures are defective. |
| 7 ATAC / survival / foundation models / release (unplanned work) | **ATAC, survival, release: sound with corrections. Foundation-model comparison: not defensible** | scVI row of the embedding comparison is wrong (28.5% of cells had NaN latents zero-filled); on valid cells scVI is the best embedding, reversing the stated conclusion. Release h5ad over-calls "malignant cell" (CL:0001064) for 70,924 ambiguous cells. |

## 2. Status of the headline claims

| Claim as currently written | Status after audit |
|---|---|
| Reference of 311,618 cells / 57 samples / 52 patients; scVI seed 1 chosen by pre-stated rule; Harmony over-corrects | **Stands** (metrics recomputed) |
| scANVI reproduces lineage labels at 0.987 | **Re-state**: 0.980 held-out; astrocyte recall 0.66, T/NK 0.39 on validation; no patient-held-out estimate |
| Query mapping: Liu malignant recall 0.993; Siletti 0 false-malignant; Braun collapses to astrocyte | **Withdraw and regenerate** (IM-01). Latents are valid; labels are not |
| Cell of origin: malignant cells neighbour fetal glioblasts / oligodendroglial lineage | **Stands as a kNN-in-latent observation**, with the off-manifold caveat; must be re-described as zero-shot projection |
| Three-tier malignant calls; tier-2 agreement with Jessa 0.913/0.923 | **Stands at study level**; report sample-level median 0.90 (IQR 0.84–0.96) and the true tuning set (22 samples incl. Okada/Collot); tier-3 "0.95–0.97" is on 51% of cells and must be labelled as such |
| Hallmark CNV frequencies (1q, 14q, 17p, 4q12) | **1q/2/13q/14q stand; 17p and 4q12 withdraw** (NaN→True bug; centromere-spanning windows; 4q12 is 0–2 windows) |
| Only 2 metaprograms replicate across assays (assay-dependence headline) | **Stands at the pre-registered cut, but must be reported with the 0.01-step curve** (2 at 0.84–0.85, 3 at 0.86, 4 at 0.87–0.89, 5 at 0.90); MP14 ("OPC-like-3", cells-only) is an immediate-early/heat-shock dissociation programme and must not be counted as a lineage state |
| Masked re-test: 4 replicate (AC-like, OC-like shared across assays) | **Exploratory only**; OC-like replication rests on 2 cell samples |
| Spatial: MES-like with myeloid ρ 0.30; OPC-like with oligodendrocytes 0.25 | **MES-like–myeloid stands at ρ ≈ 0.15–0.18 on proportions (q 0.023); OPC-like hypotheses do not survive proportions**; AC-like–astrocyte 0.43 → 0.21 |
| Brainstem: OPC-like enriched in tumour core, MES-like at edge | **MES-like-at-edge replicates across both references; the OPC-like effect swaps to OC-like under the metaprogram reference** (what the fallback calls OPC-like is what the MP reference calls OC-like, spot ρ 0.96) |
| Composition: H3.1 tumours carry more normal cells | **Confounded**: all 8 H3.1 samples are pons and nuclei; Jessa only |
| ATAC: chromatin separability above chance in every sample | **False as written**: 8 of 10 samples z > 2, two at chance |
| Foundation models: scVI mixes assays at the cost of accuracy | **Reversed**: on finite cells scVI accuracy 0.82 vs 0.63 stated |
| Survival (OpenPedCan) Cox results | **Reproduce** (HRs/p match to 3 decimals with an independent implementation) |

## 3. Critical findings, plus one release-blocking major (must fix before any use of the affected outputs)

1. **IM-01 — query classifier not frozen** (`scripts/11_map_queries.py`, `models/scanvi_query_*`). State-dict diff shows encoder weights identical but classifier BatchNorm γ/β changed (up to 1.36 for Braun) during "surgery"; applying the reference classifier to the same latents gives Braun 82,910 malignant (saved: 97,822 astrocyte) and Siletti 3,541 malignant (saved: 0). Everything derived from `predicted_lineage` in Phase 5 is invalid. Fix: re-predict with `models/scanvi` on saved latents; regenerate all mapping tables, both figures and notes §2–6; describe Phase 5 as zero-shot projection (IM-02: `encode_covariates=False`, so the 100-epoch surgery could not and did not change the latent).
2. **CNV-01 — NaN treated as True** (`06_cnv_aggregate.py`: `astype(bool).mean()`). chr17p_loss is NaN for 63,694 cells in 8 samples with no 17p window; those cells count as "loss". chr17p frequencies for Okada, Collot, Jessa scRNA and LaBelle, and the "Collot 17p loss 0.74" sentence, are wrong. Fix: `.eq(True)`, report evaluable-cell counts, regenerate both tables.
3. **MC-01 — assay-determined metaprogram-family assignment** (`19_metaprogram_scores.py` L127–129). Argmax after global z-scoring: scRNA cells 66.7% OPC-like / 0% MES-like / 0% OC-like; multiome nuclei 34% MES-like. This column feeds the §10 cell2location reference and `metaprogram_family_composition.csv`. Fix: z-score within sample, mean rather than max over member MPs, report the by-assay crosstab; do not use the current column as a reference.
4. **P7-01 — scVI row of the embedding comparison** (`26_foundation_embeddings.py`). 4,281 of 15,000 cells have NaN `X_scvi_malignant` and were zero-filled; recomputed on 10,719 finite cells: accuracy 0.819 (table 0.629), transfer 0.822 (0.622), same-assay 0.931 (0.678). Fix: same cell set for all four embeddings; regenerate table, Fig 6d, notes.
5. **P7-05 (major, release-blocking) — release ontology** (`27_package_release.py`). 240,603 cells carry CL:0001064 "malignant cell" = 169,492 consensus-malignant + 70,924 ambiguous + 187 normal; 2,675 QC-fail nuclei released. Fix: assign CL:0001064 only to consensus-malignant cells, or document explicitly.

## 4. Major findings by theme

**Circularity and over-stated validation.** CNV thresholds were tuned on 22 samples (9 Jessa scRNA, 7 Okada, 4 multiome, 2 Collot), not the 7 declared, so Okada/Collot agreement is not independent (CNV-06). Tier-3 agreement is computed only on the 51% of cells where marker and CNV calls already agree (CNV-04). scANVI accuracy is 90% training data (IM-03). The masked metaprogram re-test removes 46% of genes and cannot separate technical from biological assay differences (MC-05).

**Unit of replication and confounding.** Agreement, correlation and composition statistics are pooled over cells where the sample/section/patient is the unit (CNV-08). 10 of the 21 "independent" DMG Visium sections are serial sections of one autopsy (66% of DMG spots) — noted by the auditor as a pseudo-replication risk for the random-effects meta-analysis. Autopsy vs resection is unestimable (Okada only); H3.1 is confounded with pons and nuclei (MC-06).

**Assay (cells vs nuclei) is the dominant structure and is under-controlled.** 29% of the "malignant" cells used in Phase 6 are scANVI rescues of ambiguous cells at assay-dependent rates (snRNA 33%, scRNA 12%; MC-04). Three nuclei-only "lineage" MPs score ~1.8× higher in normal myeloid than in tumour cells (MC-08); the OPC-like family mixes ribosomal, mitochondrial and IEG programmes (P7-04). The spatial fallback state labels come from un-standardised score argmax, with MES-like at 17.5% of malignant cells — high for H3K27M DMG (SP-06).

**Methodological artefacts.** CNV windows are chromosome-wide, so short-arm events are diluted and 4q12 is represented by 0–2 windows (CNV-02/03); ~21% of CNV-malignant calls rest on the absolute score floor with no hallmark event (CNV-09). Co-localisation on raw q05 abundances shares a density/UMI component (SP-01). Deconvolution ELBO histories were not harvested; brainstem and Ren ran at 41%/49% of target epochs (SP-03).

**Reproducibility gaps.** Label tables and the Collot/ScPCA/Okada downloads have no generating script (DQ-02); the Filbin/Liu reference gene sets used for every metaprogram label have no source file or checksum (MC-03); a partial figure re-run overwrote the supplementary-figure provenance table (P7-06). okada__1734_pons (median 427 genes) entered the reference with no sample-level flag (DQ-03); 16 Liu samples have no MT genes so the mito filter was inoperative (DQ-04); LaBelle's layer is log2(TPM/10+1) mislabelled as TPM (DQ-01; downstream code detected it heuristically).

**Figures.** Brainstem gradient: legend overlaps panel, labels truncated, no right y-label (SP-04). Abundance maps: array coordinates stretch sections ~1.6×, no section identification or scale bar (SP-05). TME marker-score heatmap illegible (MC-09). Several figures still carry in-panel titles/labels against the house rule.

## 5. What the audit confirmed as sound
Phase 1–2 cell accounting and QC tables; Okada UMI recovery (exact to 1e-12); gene universe and the 18,882-gene intersection; Liu cell-ID recovery (8,984/9,248); CNV arithmetic, determinism and the 342,213-cell reconciliation; held-out circularity check for tier-2 (tuning vs non-tuning Jessa samples); scIB metric values and model choice; HVG selection; pseudo-count construction; cNMF determinism, programme clustering, the replication matrix as implemented vs pre-registered; reference overlaps and hypergeometric tests; cell2location regression convergence and signature sanity; DMG vs GBM Mann-Whitney values; Cox HRs; ATAC headline counts; release h5ad contents vs atlas_core; caption numbers.

## 6. Remediation plan (ordered; estimated effort)
1. Phase 5 regeneration with the reference classifier; rewrite notes as zero-shot projection; redo both figures (IM-01/02/04) — 2–3 h, CPU.
2. CNV: fix NaN handling; per-arm smoothing or gene-level arm means; drop or re-call 4q12; add event-backed confidence flag; report sample-level agreement with CIs and the true tuning set (CNV-01/02/03/06/08/09) — 1 day; re-run 05/06 (~25 min compute).
3. Metaprograms: within-sample z-scoring and mean-over-MPs for families; relabel MP14 as dissociation stress, exclude MP3/MP17/MP14 from the OPC-like family; report 0.01-step sensitivity; sensitivity re-run on consensus-malignant cells only; script the reference gene sets with checksums (MC-01/02/03/04/07/08) — 1 day.
4. Spatial: proportion-based (or CLR) correlations as primary; state the OPC-like→OC-like label swap; harvest ELBO histories or re-run brainstem/Ren to 20k epochs; redo both figures with pixel coordinates and scale bars (SP-01…06) — 1 day plus ~1 GPU-hour.
5. Phase 7: recompute the embedding comparison on a common finite cell set; correct ATAC "every sample" and n-per-pair in the correlation map; fix release ontology and QC-fail cells; restore supplementary provenance (P7-01…06) — half a day.
6. Reproducibility: scripts for label tables and the three non-GEO downloads; sample-level QC gate; MT-gene presence column; LaBelle layer rename (DQ-01…04) — half a day.
7. Reporting: sample-level statistics with CIs throughout; explicit confounding statements (autopsy=Okada, H3.1=pons+nuclei, serial sections); move all in-panel text to captions.

## 7. Limits of this audit
Targeted re-execution only (two small studies for CNV; two samples for NMF; 5 of 9 scIB metrics on 3 embeddings; one query for the live latent-identity check; no GPU runs, no re-deconvolution, no scGPT/Geneformer inference); checksums re-hashed for 9–10 files per track; supplementary xlsx sources not re-parsed; figure PDFs not opened (PNGs of the same objects were). Centromere coordinates were verified for chr1, chr2 and chr17 against external text and for the rest against gene-free gaps.
