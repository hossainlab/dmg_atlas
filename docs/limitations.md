# Limitations of atlas v1

Moved verbatim from `README.md` (consolidated 2026-10-21). Every headline number is the value in the current
`results/tables/` and `data/processed/` files; the audit trail that produced them is `docs/audit/ALL_findings.csv`, with per-track
closure in `docs/audit/fixes_*.csv` and `docs/audit/CLOSURE_SUMMARY.md`.

**Unit of replication.** Statistics are reported at the unit they were estimated on, and the unit is stated wherever a number is quoted.
- *Malignant calling.* Agreement with author labels is a **sample**-level quantity: median 0.898 (IQR 0.793–0.978) over the 97 labelled
  samples, unweighted mean 0.836 (95 % CI 0.790–0.877 over 2,000 sample bootstraps), against a cell-pooled 0.884
  (`malignant_call_validation.csv`, audit CNV-08). Hallmark event frequencies are given per cell *and* per sample
  (`cnv_event_frequencies.csv`, CNV-23); the per-sample fraction is the one comparable with patient-level genomics.
- *Query mapping.* Label-transfer rates carry per-sample medians, IQRs and donor/plate-bootstrap CIs in
  `mapping_agreement_sample_level.csv` (IM-22); the reference classifier is scored on **held-out patients**, not held-out cells.
- *Metaprograms.* The replication rule counts **samples** and **publications**; `metaprograms.csv` also reports `n_patients` per MP
  (50 cNMF samples = 46 patients, MC-17).
- *Composition.* Models are fitted on **sample** proportions with patient-clustered standard errors and report `n_patients_ref` /
  `n_patients_test` per contrast (MC-06).
- *Spatial.* Co-localisation is pooled over the **12 DMG patients**, not the 21 sections (SP-11).
- *Survival.* Cox models are per **patient** (n = 99, 94 deaths).

**Confounds that cannot be removed from these data.**
- **Autopsy vs resection is not estimable.** All 22 autopsy samples come from Okada, so the contrast is collinear with study; it is
  omitted from every composition model rather than reported (MC-06, DQ-19).
- **H3.1 vs H3.3 rests on 8 samples.** All 8 H3.1 samples are Jessa 2022, all pontine and all nuclei, so H3 variant is confounded with
  site and assay. Of the cell types that reach q < 0.05 under assay adjustment alone (astrocyte, myeloid, neuron), only **neuron**
  survives site adjustment (q 0.017) and the pons-only fit (q 0.050); the others are not claimed (MC-06).
- **The 10 brainstem Visium sections are one patient.** They are serial sections of a single autopsy and carry 66 % of the DMG spots, so
  spatial estimates are pooled over patients and the rostral-to-caudal gradient is descriptive only (SP-11, SP-04).
- **The cell arm of every assay comparison is 6 scRNA samples** from 3 publications, so cell-vs-nucleus differences in metaprograms and
  TME subclusters cannot be separated from cohort (MC-01, MC-05); after the corrected selection no cell-derived lineage metaprogram
  exists at all (MC-07).
- **CNV thresholds were tuned on 22 of the 97 labelled samples** (9 Jessa scRNA, 7 Okada, 4 Jessa multiome, 2 Collot). Tuning and
  held-out samples agree to within 0.003 pooled, so there is no evidence of inflation, but the thresholds are not independent of those
  22 samples (CNV-06, CNV-07).
- **Liu QC is heterogeneous**: 16 samples have no MT genes and are unfiltered for mitochondrial content, 8 frozen-nucleus samples are
  filtered at 5 %, and 18 fresh-cell samples at the cell threshold (DQ-04).
- **ATAC, bulk survival and the foundation-model comparison are single-cohort or external**: all 10 multiome samples are Jessa 2022
  nuclei, and the survival cohort is OpenPedCan bulk, not the atlas.

**Method limitations.**
- Query mapping is a **zero-shot projection** by the reference scANVI encoder: the reference was trained with
  `encode_covariates=False`, so scArches surgery cannot adapt the latent and no query batch effect is removed (IM-02). The per-query
  models in `models/scanvi_query_*_pre_audit/` carry an unintentionally re-fitted classifier BatchNorm (scvi-tools 1.4.2) and must not be
  used for prediction.
- The reference classifier scores **0.9751 accuracy / 0.917 balanced accuracy on 23,827 labelled cells from 5 held-out patients**
  (train 0.9931 on 197,867 cells / 47 patients). Held-out-patient recall: malignant 0.999, neuron 0.969, myeloid 0.950,
  oligodendrocyte 0.917, T_NK 0.896 (77 cells), vascular 0.883, **astrocyte 0.803** — astrocyte is the weakest class under every
  scheme (`reflabels/scanvi_heldout_comparison.csv`, IM-03).
- TPM cohorts (Filbin, LaBelle) are mapped with pseudo-counts and are **projection-only**; TPM is gene-length normalised and no length
  correction is possible from the deposited tables. Labels are insensitive to the pseudo-count scale *m* over a 16-fold range, but the
  kNN-uncertainty percentile partly measures *m* rather than biology (IM-12).
- Normal progenitor classes with **no counterpart in the reference** (fetal radial glia, glioblast, neuroblast; adult OPC) are placed in
  the malignant region by both the classifier and the kNN label: 87.7 % of Braun fetal hindbrain cells and 29.0 % of Siletti adult OPCs
  are predicted malignant. `predicted_lineage` is not a cell-type annotation for those cohorts; within classes that do have a
  counterpart the false-malignant rate is 0.2 % (IM-01, IM-04).
- `chr4q12_gain` is **not a hallmark event**: with 100-gene windows the band is covered by at most one ~40-Mb window measuring
  PDGFRA/OPC-lineage expression against an immune reference. It was removed from the hallmark set and exported only as exploratory
  `gene_dev_PDGFRA` / `gene_dev_KIT` (CNV-03). `chr17p_loss` frequencies are much lower than first reported because NaN arm means were
  previously counted as events and windows straddled the centromere (CNV-01, CNV-02).
- **2 of 15 metaprograms replicate** under the pre-registered rule (MP1 cycling, MP4 hypoxia); the count is 2 across the whole
  0.84–0.90 cut range and reaches 3 only at 0.91 (0.01-step curve in `metaprogram_cut_sensitivity.csv`, MC-02). Three programmes are
  technical and carry no family — MP2 ribosomal, MP12 mitochondrial and **MP13 IEG / dissociation stress** (MC-07, MC-08).
- The **assay-purity drop** after per-compartment scVI (myeloid 0.95 → 0.02) is not independent evidence that the assay effect was
  removed: `assay_class` was an encoded covariate and no held-out-gene or covariate-free control was run (MC-22, **open**).
- **Spatial co-localisation is measured on spot proportions.** Of the two pre-stated hypotheses, MES-like × myeloid is supported
  (ρ = 0.238, q = 0.034) and OPC-like × oligodendrocyte is not (ρ = 0.178, q = 0.210); the strongest pair is AC-like × astrocyte
  (ρ = 0.393, q = 1e-4). Between-section heterogeneity is high for every pair (I² 0.952–0.997 over the 42 state × microenvironment
  pairs), so the pooled estimates are not section-generalisable (SP-01, SP-11).
- **Ren DMG vs GBM is null**: tumour fraction 0.528 vs 0.460 (p = 0.056, q = 0.667) and cycling share 0.112 vs 0.214 (p = 0.095,
  q = 0.667); nothing survives BH across the 14 variables at n = 5 vs 5 sections (SP-12).
- **Bulk NNLS deconvolution is unreliable** (median T/NK proportion 0.13; several malignant families zero in most patients), so the
  survival analysis is null rather than negative. The rank-score Cox fits show no association (min BH q = 0.17), and proportional
  hazards is violated for the cycling score in all three specifications (P7-07, Schoenfeld rank test).
- **ATAC peaks are not called by this pipeline** — they are the authors' cellranger-arc per-sample peaks from GEO, and ATAC fragments
  were never downloaded; the differential-accessibility counts are exploratory (binomial sign test over 5–9 samples, no multiple-testing
  control across 279,115 union peaks, no depth normalisation) (P7-10, P7-11).
- **Foundation-model comparison uses expression-derived labels**, so the PCA family accuracy is circular; scGPT and Geneformer
  neighbourhoods are almost assay-pure (0.996 / 0.997) and only scVI mixes assays (P7-01).

**Known-open items (not closed by the remediation).** `docs/audit/CLOSURE_SUMMARY.md` lists 11 not-closed and 12 partially-closed
findings; the substantive ones are MC-22 (no independent test of the TME assay correction), MC-12 (42 % of cNMF samples sit at the
k = 9 ceiling), CNV-11 (score layer not recomputed after MAD refinement), DQ-08/DQ-18 (gene-symbol and ribosomal-prefix issues that
need a full re-ingest), P7-07 (NNLS survival sets not re-fitted) and P7-14 (no clean-environment re-run). 15 findings were never
assessed by any track and are recorded as `not_assessed`, not as fixed.

## Per-layer notes

**Reference and query mapping.** Summary (query mapping is a zero-shot projection, the reference classifier is scored on held-out patients, TPM cohorts are projection-only, and normal progenitor classes without a reference counterpart are absorbed by the malignant class) — see also `docs/phase5_mapping_notes.md` §2–6.

**Phase 6 (states, metaprograms, spatial).** Cell- and nucleus-derived programmes and TME clusters separate by assay (only 6 cell samples, 3 publications), and after the corrected malignant selection no cell-derived lineage metaprogram exists at all; spatial state classes come from the metaprogram-family reference for **all four** cohorts since the spatial remediation (the score_genes-argmax fallback reference was dropped, audit SP-06, so `brainstem_*_mp.csv` and `brainstem_mp_vs_fallback.csv` from step 21 are superseded and no figure reads them); the brainstem gradient is one autopsy (10 serial sections) and is descriptive only; autopsy vs resection composition is not estimable (all 22 autopsies are Okada).

**Phase 7 (multimodal and bulk).** Bulk NNLS deconvolution is unreliable (proportions disagree with rank scores), so the survival analysis is null rather
than negative, and the rank-score fits violate proportional hazards for the cycling score; ATAC peaks are the authors' cellranger-arc peaks
(not called here) and the DA calls are exploratory and not depth-normalised; the foundation-model comparison uses expression-derived (partly
circular) labels. `results/figures/phase7_survival_forest.png` is a pre-audit figure (2026-10-18) built from `survival_results.csv` and
still shows the degenerate NNLS variable sets — Figure 6c supersedes it.

> **Scope note (audit P7-16).** Phase 7 (code) / build-plan Phase 8 — multimodal ATAC, bulk survival anchoring and foundation-model
> embeddings — **was not part of the approved plan**. It is reported as an exploratory methods note, not as an atlas claim, and its
> figures (Figure 6, Supplementary Figure 8) carry that label. See `docs/release_v1.md`, "Scope of the analytical layers".
