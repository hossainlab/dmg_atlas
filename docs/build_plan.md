# Build plan — cross-cohort single-cell atlas of H3K27M diffuse midline glioma

Inputs: dmg_dataset_selection.csv (v3). Compute: MacBook (14 cores, 24 GiB; scVI measured 3.6 s/epoch at 30k cells × 3k genes), Modal GPU.

## Phase 0 — Freeze design decisions (before any download)

1. Scope of the reference: malignant + microenvironment in one scVI model (recommended), with a second malignant-only model for state discovery.
2. Liu 2022 (Smart-seq2): projection-only query (recommended), or integrated from raw RSEM counts with technology covariate.
3. Normal references: Braun 2023 fetal hindbrain (cell-of-origin) + Siletti 2023 adult pons (standard 10x); keep Xie 2026 Paired-Tag as a secondary adult reference or drop.
4. Okada 2026 autopsy cohort: include if processed matrices download from the UCSC Cell Browser; otherwise watch-list.
5. Target outputs: (a) scANVI reference model + harmonised labels published for query mapping; (b) biology paper.
   Deliverable: decisions.md (one line per decision with rationale).

## Phase 1 — Acquisition and per-sample manifest (≈1 day, ~5 GB)

1. Download: Jessa GSE210567 (1.4 GB) + GSE210566 multiome RNA (RNA part only), Collot Zenodo scRNA.tar.gz (1.0 GB), DeSisto GSE231859 (0.2 GB), Zang GSM7846769 (0.07 GB), ScPCA SCPCP000021 two DMG samples, Filbin GSE102130 TPM, Liu GSE184357 (raw counts tar + metadata), LaBelle GSE227983 + GSE271936 TPM, Visium: Huchedé GSE268577, Ren GSE194329 (DMG sections), Collot ST.tar.gz (human), whole-brainstem GSE280990.
2. Build sample_manifest.csv: one row per biological sample with patient_id, study, assay (cell vs nucleus), chemistry, tissue site (pons/thalamus/spinal/other), age, H3 variant (H3.1/H3.2/H3.3/WT/unknown), treatment status, biopsy vs autopsy, source file. H3 status and sc/sn labels for Jessa, DeSisto and LaBelle come from the papers' supplementary tables — record the table and column used.
3. Deduplicate patients across GSE210567/GSE210566 (shared IDs e.g. P-1764_S-1766) and check CBTN-sourced overlap between DeSisto/ScPCA and OpenPedCan.
4. Exclude: non-DMG pHGG samples (kept separately as out-group), ependymoma sample, mouse IUE cells, cell-line/PDX samples.
   Deliverable: sample_manifest.csv; raw AnnData per study (.h5ad, counts layer).

## Phase 2 — Per-study QC and harmonisation (≈1–2 days)

1. Gene-space harmonisation: map all to one annotation (GENCODE symbols; Ensembl IDs where features.tsv provides them); record genes lost per study.
2. QC per sample with assay-specific thresholds (nuclei: mito <5 %; cells: mito <20 %), min genes 200–500, doublet scoring (scrublet) per 10x sample.
3. Ambient RNA check on 10x cell samples (soupX/cellbender optional; at least report ambient markers).
4. Keep author labels where present (Jessa malignant/normal consensus, Liu annotation, Filbin states) as a separate obs column; never overwrite.
   Deliverable: qc_report.md with cells before/after per sample; per-study QC'd h5ad.

## Phase 3 — Malignant vs normal calling (≈1–2 days)

1. inferCNV-style CNV scoring per sample (infercnvpy) with study-internal normal cells (microglia, oligodendrocytes, T cells) as reference; H3K27M DMG hallmarks: chr1q gain, chr2 gain, chr13q/14q loss, PDGFRA amplification where present.
2. Consensus call = CNV score + marker expression (OPC/AC/MES programs vs microglia/oligodendrocyte/endothelial markers) + author label; store three tiers (author, own-CNV, agreement).
3. Report agreement with author labels for Jessa (143,611/33,301) as the validation statistic.
   Deliverable: malignant_calls.csv; per-sample CNV heatmaps (supplementary).

## Phase 4 — Integration (reference model) (Mac: ~3–4 h training total)

1. Assemble UMI cohort: Jessa scRNA/snRNA + multiome RNA, Collot, DeSisto (K27M only), Zang, ScPCA ×2, (Okada if available). Expected ~300k cells.
2. HVGs: 3,000 by seurat_v3 flavor with batch_key = sample, computed on counts.
3. scVI: batch = sample; categorical covariates = study, assay (cell/nucleus); n_latent 30, 2 layers, NB likelihood; default epoch heuristic (~25 epochs); 3 seeds.
4. Model selection: scIB metrics (batch: iLISI, kBET on study; bio: cLISI, NMI/ARI on author labels, cell-type ASW). Compare scVI vs Harmony vs unintegrated PCA as controls.
5. scANVI from the best scVI: labels = coarse lineage (malignant OPC-like / AC-like / MES-like / OC-like / cycling; myeloid; T; oligodendrocyte; astrocyte; neuron; endothelial/pericyte; unknown).
6. Malignant-only model: subset malignant cells, rerun scVI with the same covariates for state discovery.
   Deliverable: scvi_model/, scanvi_model/, atlas_core.h5ad (latent, UMAP, labels), integration_metrics.csv.

## Phase 5 — Cell-state definition and replication (≈1 week)

1. Malignant programs: cNMF per sample (k selected per sample), cluster program gene-sets across samples into consensus metaprograms; compare to Filbin 2018 (OPC/AC/OC-like) and Liu 2022 (OPC-like 1–3, MES-like) and Jessa cNMF programs.
2. Replication test: a metaprogram counts as replicated if recovered in ≥3 studies and both assays (cell and nucleus). Report per-metaprogram study × assay matrix.
3. Microenvironment: subcluster myeloid, T/NK, vascular, oligodendroglia; annotate against reference markers.
4. Composition models: scCODA or per-sample proportions with mixed models for age group (pediatric vs adult), site (pons/thalamus/spinal), H3 variant, biopsy vs autopsy, treatment — study as random effect.
   Deliverable: metaprograms.csv, replication_matrix.csv, composition_results.csv.

## Phase 6 — Query mapping (Mac: <1 h)

1. scArches-style scANVI query mapping of Filbin 2018 and LaBelle 2025 (TPM; fine-tune with frozen reference weights; interpret with caution — report mapping uncertainty).
2. Map Liu 2022 (if projection) — the only adult+pediatric multi-site cohort; test whether its OPC-like 1–3 substates land on distinct reference states.
3. Map normal references (Braun fetal hindbrain, Siletti adult pons, Xie pons) to locate the closest normal counterparts of each malignant metaprogram; report label-transfer entropy.
   Deliverable: query_mapping.h5ad, mapping_uncertainty.csv.

## Phase 7 — Spatial validation (≈1 week)

1. Deconvolve Visium spots with cell2location (or RCTD) using the atlas signatures — four cohorts: Huchedé (3 sections), Ren (5 DMG sections + 5 GBM out-group), Collot human sections, whole-brainstem DIPG (10 sections).
2. Test co-localisation hypotheses derived from Phase 5 (e.g., MES-like with myeloid; OPC-like with oligodendroglial lineage), with section as unit of replication.
3. Infiltration gradient: in GSE280990, test which metaprograms decline from pons core to medulla/spinal cord.
   Deliverable: spatial_deconvolution.h5ad per cohort, colocalisation_results.csv.

## Phase 8 — Multimodal and bulk anchoring (optional, Modal GPU useful)
> **Not part of the approved plan (audit P7-16).** This phase was listed as optional and was never approved as a deliverable; it was
> nevertheless executed as code Phase 7 (`scripts/24`–`26b`). Everything it produced is an exploratory methods note, not an atlas
> claim — see `README.md`, `docs/release_v1.md` ("Scope of the analytical layers") and `docs/phase7_multimodal_notes.md`. Nothing in
> Phases 1–6 depends on it.

1. Jessa multiome ATAC (12 GB): peak accessibility per scANVI label; check that transcriptional metaprograms have distinct chromatin states (Okada snATAC as second cohort if available).
2. OpenPedCan DMG bulk (144 patients): deconvolve with atlas signatures; Cox models for metaprogram abundance vs survival, after removing patients overlapping single-cell cohorts.
3. Foundation-model embeddings (scGPT/Geneformer, GPU): compare state separability with scVI as a methods note.
   Deliverable: atac_summary.csv, survival_results.csv.

## Phase 9 — Packaging and release

1. Final h5ad with standardised obs schema (cellxgene schema fields), model files, label ontology, manifest and code repository.
2. Reproducibility: run the integration end-to-end from raw deposits once in a clean environment; record versions.
3. Figures per house style (no panel titles; descriptions in captions). Manuscript with every claim tied to a number or citation.
   Deliverable: dmg_atlas_v1.h5ad, models, GitHub repo, manuscript.

## Checkpoints where the plan can fail

- Phase 1: H3 status unobtainable for DeSisto or LaBelle samples → those cohorts shrink to the labelled subset.
- Phase 3: own-CNV calls disagree with Jessa author labels by >10 % → revisit reference cell choice before proceeding.
- Phase 4: study dominates the embedding after scVI (iLISI on study near 1) → add study as a batch (not covariate) or try scPoli.
- Phase 5: fewer than 3 metaprograms replicate across assays → the headline becomes an assay-dependence finding; still publishable.
