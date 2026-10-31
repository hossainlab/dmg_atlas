# Phase 6 — malignant metaprograms (cNMF) and Visium deconvolution / spatial tests (2026-10-17; audit corrections to §1, 2, 6, 8, 9 on 2026-10-19 — see `docs/audit/fixes_metaprograms.csv`; pre-audit text in `docs/archive/phase6_notes.pre_audit.md`)

Code: `scripts/13a_malignant_selection.py` (scripted malignant selection, audit), `scripts/13_cnmf_per_sample.py` (per-sample consensus NMF), `scripts/16a_curated_gene_sets.py`, `scripts/16_metaprograms.py` (metaprogram clustering + replication test),
`scripts/13_reference_signatures.py prepare|collect` and `scripts/14_cell2location.py prepare|collect` (cell2location inputs / harvest),
`modal/phase6_c2l.py signatures|deconvolve` (GPU, A10G), `scripts/15_spatial_tests.py` (co-localisation, niches, brainstem gradient, Ren out-group).
`src/dmg_atlas/metaprograms.py` holds the cNMF / Jaccard / UCell helpers. Logs: `logs/phase6_metaprograms.log`, `logs/phase6_spatial.log`,
`logs/phase6_13_stdout.log`, `logs/phase6_14_prepare_stdout.log`, `logs/phase6_cnmf_stdout.log`.
Build-plan mapping: this is build-plan Phase 5 (cell states, replication) + Phase 7 (spatial validation).

**Where things ran.** cNMF, signature/Visium preparation: local Mac. cell2location `signatures` and the four `deconvolve` jobs: Modal A10G
(Claude Science HPC jobs; outputs harvested from `~/.claude-science/.../hpc/<job>/out/`). `16_metaprograms.py`, `14_cell2location.py collect`
and `15_spatial_tests.py`, and later `17_tme_subclusters.py`, `18_composition.py`, `14_cell2location.py package`: Modal CPU sandbox (8 CPU, 32 GB; same image `im-gl8KES0WsH7PoFjf83hK6v`, scanpy 1.11.5) via
`modal/phase6_cpu_run.py`, which stages the needed slice of the repo, runs the scripts unchanged and copies the outputs back; moved off the
Mac because of thermal load. No earlier-phase file was modified.

## 1. Per-sample cNMF
(`results/tables/cnmf/`, `data/processed/cnmf_usages_<sample>.parquet`; pre-audit run kept in `results/tables/cnmf.pre_audit/` and
`cnmf_usages_<sample>.pre_audit.parquet`)

**Audit corrections (2026-10-19; MC-04, MC-16, MC-17).** Malignant selection is now scripted (`scripts/13a_malignant_selection.py`,
`data/processed/malignant_selection_audit.parquet`) and uses the corrected Phase-3 consensus calls only: **primary = `consensus_call ==
'malignant'` (event-backed or score-only CNV confidence) AND `qc_pass`**, excluding `okada__1734_pons` (`sample_qc_pass` False; the ScPCA
samples are not in `atlas_core`). The pre-audit rule added scANVI-rescued consensus-ambiguous cells (maxprob ≥ 0.9) at assay-dependent
rates — 27.9 % of the pre-audit input (cells 11.9 %, nuclei 29.9 %; `malignant_selection_by_assay.csv`). Selection breakdown per assay class
(qc_pass cells, 56 samples):

| assay class | qc_pass | pre-audit rule | of which scANVI rescue | primary | event-backed | score-only | eventtrust (sensitivity) |
|---|---|---|---|---|---|---|---|
| cells (scRNA) | 36,287 | 24,601 | 2,929 (11.9 %) | 21,603 | 19,486 | 2,117 (9.8 %) | 23,479 |
| nuclei (snRNA + multiome) | 264,190 | 204,831 | 61,180 (29.9 %) | 145,409 | 131,664 | 13,745 (9.5 %) | 171,447 |
| all | 300,477 | 229,432 | 64,109 (27.9 %) | 167,012 | 151,150 | 15,862 (9.5 %) | 194,926 |

The score-only share is now assay-balanced (9.8 % vs 9.5 %). Sensitivity set: `consensus_call_eventtrust == 'malignant'` (marker-normal cells in
event-backed CNV-malignant clusters trusted; CNV-05), +27,914 cells.

Samples with ≥ 200 primary cells: **50 of 57** (6 cell samples: Jessa 4, DeSisto 1, Zang 1; 44 nuclei samples; 46 distinct patients — four Jessa
tumours are profiled twice, so sample counts over-count tumours and `metaprograms.csv` now carries `n_patients`). The per-sample cNMF was
re-run for the 46 samples whose malignant set changed by > 5 % of the pre-audit set (`cnmf/sample_rerun_plan.csv`; median change 25 %,
range 7–73 %); the stored pre-audit programmes were reused for the four samples with ≤ 5 % change (`desisto__1196` 1.8 %,
`jessa_scrna__P-6117_S-8370` 3.2 %, `jessa_scrna__P-6640_S-9581` 3.8 %, `jessa_multiome__P-6640_S-9581` 4.0 %). Four pre-audit samples dropped
out (`okada__1734_pons` excluded; `jessa_scrna__P-6254_S-8499`, `okada__1006_thalamus`, `okada__1934_thalamus` now < 200 malignant cells).
Settings unchanged: k in 4..9, 10 random-init restarts, k* = largest k with silhouette ≥ 0.90 and every consensus cluster in ≥ 80 % of restarts
(chosen k: 4:1, 5:1, 6:3, 7:8, 8:16, 9:21; the silhouette criterion remains non-discriminative, MC-12, and 42 % of samples sit at the k = 9
ceiling — not changed in this remediation). **400 programmes.** The eventtrust set (51 samples, 195k cells) was run in full
(`results/tables/cnmf_eventtrust/`, `cnmf_usages_eventtrust_<sample>.parquet`). Both runs: Modal CPU sandboxes (8 CPU, 32 GB, no GPU; four
shards, 521–692 s each, 39 sandbox-minutes in total; `logs/fix_metaprograms.log` holds the remote logs). The 2026-10-17 deviation
(sklearn NMF, random init, 10 restarts instead of nndsvd / 3) stands.

## 2. Metaprograms and the pre-registered replication test (`metaprograms.csv`, `replication_matrix.csv`, `metaprogram_cut_sensitivity.csv`, …)
Method fixed in the `16_metaprograms.py` docstring before the clustering was run and unchanged: Jaccard on top-50 genes, average linkage on
1 − J, cut 0.85 (mean pairwise J ≥ 0.15), MP = cluster with ≥ 5 programmes from ≥ 3 samples. Replication rule (decisions.md #9): ≥ 3
publications and both cells and nuclei (`jessa_scrna` and `jessa_multiome` are both Jessa 2022 and count once; 5 publications).
**The cell arm is 6 samples** (Jessa 4, DeSisto 1, Zang 1; 3 publications; two of the Jessa cell samples are the paired patients P-6337 /
P-6640 that also contribute nuclei).

**Audit corrections (2026-10-19; MC-02, MC-03, MC-05, MC-07, MC-08, MC-15, MC-17, P7-04).**
(a) Reference gene sets are now all scripted and tiered (`16_metaprograms.py` docstring): *published* Liu 2022 Table S2a sets
(`genesets/reference_gene_sets_scripted.csv`, `scripts/01d_reference_gene_sets.py`, sha256 in `reference_gene_sets_provenance.csv`), *curated*
functional / marker sets (`genesets/phase6_curated_gene_sets.csv`, written by `scripts/16a_curated_gene_sets.py` from
`dmg_atlas.states.CURATED_PHASE6_SETS`, sha256 `cdb05c32…`; hand-curated, not downloaded), *in-house DE* sets (the re-derived Liu
author-state Wilcoxon lists, clearly labelled `inhouseDE:`), and the *Filbin-derived* `states.py` lists (`filbin_untraceable:`; the Filbin 2018
tables could not be fetched by script). The former `phase6_gene_sets.csv` "Liu2022_*" sets were in-house DE lists, not the published tables
(Jaccard 0.27–0.55 to Table S2a); the file is kept as `phase6_gene_sets.pre_audit.csv`. Labels take the highest tier with a qualifying set
(BH q < 0.05, ≥ 5 shared genes); `label_provenance` and the best published-set overlap are reported for every MP.
(b) Technical programmes are classified before labelling and excluded from lineage counts: ≥ 15 ribosomal-protein genes (explicit RP regex,
MC-15), ≥ 5 MT- genes, ≥ 10 genes of the curated IEG / heat-shock dissociation set, ≥ 5 pan-immune genes (contamination).
(c) `replicated_strict` additionally requires ≥ 2 samples per assay class; `n_patients` counts tumours; the cut-sensitivity curve runs
0.75–0.92 in 0.01 steps.

**Result at the primary cut (primary selection): 15 metaprograms (161 / 400 programmes assigned), 12 lineage MPs, 2 replicated** — both
replicate under the strict rule as well:
- MP1 cycling (published G2M, 21 / 50 shared; 34 programmes, 32 patients, 4 publications, 6 cell + 28 nuclei samples).
- MP4 hypoxia (curated Hypoxia set, 17 shared; VEGFA, NDRG1, HILPDA, SLC2A1, CA9; 15 programmes, 13 patients, 4 publications, 4 cell + 11
  nuclei samples).

Lineage MPs that are **nuclei-only**: OPC-like MP3 (curated OPC markers; 18 samples, 3 publications) and MP6 (published OPC-like-1; 12 samples),
OC-like MP5 (published OC-like, 10 shared; 12 samples), AC-like MP7 (in-house DE AC-like; 8 samples), MP10, MP14, MP15 (published AC-like;
6 / 5 / 5 samples), MES-like MP8 (curated MES markers; 7 samples, 2 publications), S-phase MP11 (published S; 5 samples) and the cilia
programme MP9 (ADGB/CFAP43/DNAH9; novel; 6 samples). **No lineage programme is cell-only.** Technical / excluded: MP2 ribosomal (46 RP genes
in its top 50; 17 nuclei + 1 cell sample; its best published match is OPC-like-2, which is itself ribosome-rich), MP12 mitochondrial (11 MT-
genes with APOE/PMP2/GPX3; 5 nuclei samples) and MP13 **IEG / dissociation stress** (IER2, FOS, JUNB, EGR1, JUN, FOSB, …; 31 of 50 genes in the
IEG set; 5 cell samples, 0 nuclei; the pre-audit "OPC-like-3" label came from the published OPC-like-3 list, which is IEG-rich: 20 shared). The
pre-audit immune-contamination MP (MP16, DOCK8/APBB1IP; Jessa only) no longer forms with the corrected selection.

**Cut sensitivity (`metaprogram_cut_sensitivity.csv`, 0.01 steps).** Replicated MPs = 1 at 0.75–0.83, **2 at 0.84–0.90** (cycling + hypoxia; MP count
13 → 20), 3 at 0.91 (an AC-like programme APOE/CLU/CST3 joins, 2 cell samples) and 5 at 0.92 (3 under the strict rule). The pre-audit run sat on a
knife-edge (2 / 3 / 4 at 0.85 / 0.86 / 0.87); with the corrected selection the replicated count is flat over 0.84–0.90, so the "fewer than 3"
checkpoint no longer depends on a 0.01 step of the cut. On the eventtrust sensitivity set (§9) it does.

**Honest headline.** At the pre-registered cut, cross-assay replication is limited to cycling and hypoxia; it is stable over cuts 0.84–0.90 on
the primary selection and becomes 3 (AC-like joins) only at ≥ 0.91. Every lineage programme (OPC-, OC-, AC-, MES-like, S-phase, cilia) is
nuclei-only, and the only cell-only programme is a dissociation-stress signature that cannot exist in nuclei. Whether the lineage split is
technical (nuclear-retained transcripts) or a power limitation of the 6-sample cell arm is not decidable here: the matched-size test of the audit
(cell programmes nearest to other cell programmes, 0.26 vs 0.12 Jaccard) shows the assay structure is real, and the masked re-test (§9) is
exploratory. **Build-plan checkpoint remains triggered** ("fewer than 3 metaprograms replicate across assays → assay-dependence headline"),
now stated with its cut sensitivity.

## 3. cell2location (`results/tables/c2l_signatures.csv`, `data/processed/spatial/<cohort>/`, `spatial_abundance_long.csv`)
*[Section 3 as written below is **superseded 2026-10-20 by SP-03 / SP-06** — all four cohorts were re-deconvolved to the 20,000-epoch
target against a single 13-class metaprogram-family reference and the fallback reference was dropped; see 'Audit corrections — spatial
sections 3, 4, 7, 10' below. The pre-audit text is kept for the record.]*

Reference: 59,940 cells, 11 classes (6 normal lineages + 5 malignant states). Malignant states came from the **fallback** assignment
(scanpy `score_genes` with Filbin 2018 / Liu 2022 gene sets, argmax) because the metaprograms did not exist yet; the spatial tests are
therefore about Filbin/Liu-style states, not about the MPs of §2. RegressionModel 250 epochs (334 s, A10G).

Deconvolution (N_cells_per_location 5, detection_alpha 20; epochs adaptive to a 1,800 s budget):

| cohort | sections | spots | epochs (target 20k) | wall |
|---|---|---|---|---|
| brainstem (GSE280990) | 10 | 34,939 | 8,185 | 2,049 s |
| ren | 11 (5 DMG, 5 GBM, 1 peritumour) | 27,233 | 9,750 | 1,938 s |
| collot_st | 3 | 5,296 | 15,729 | 1,094 s |
| huchede | 3 | 1,650 | 18,806 | 1,009 s |

All four ELBO traces had plateaued: the last 5 % of epochs reduced the ELBO by 0.07–0.16 % of the total reduction, so the
budget-capped runs (brainstem, ren) were not re-run. Long table: 69,118 spots × 11 classes, 27 sections.

## 4. Spatial tests (`15_spatial_tests.py`; section = unit of replication; q05 abundances)
*[Section 4 as written below is **superseded 2026-10-20 by SP-01 / SP-11 / SP-12** — the unit of replication is the **patient** (12 DMG
patients, not 21 sections: 10 of the sections are serial sections of one autopsy), the statistics are computed on spot proportions
rather than q05 abundances, and the BH family is the 42 state × microenvironment pairs of the 7 corrected families. Corrected values
are in the 'Audit corrections — spatial sections 3, 4, 7, 10' block below.]*

**Co-localisation** (21 DMG sections, DerSimonian-Laird on Fisher z, BH over 30 state × TME pairs): AC-like–astrocyte ρ = 0.43
(0.34–0.51, 21/21 sections positive); OC-like–oligodendrocyte 0.36 (18/21); OC-like–neuron 0.32; MES-like–myeloid 0.30 (0.14–0.44, 15/21);
OPC-like–T/NK 0.29; OPC-like–oligodendrocyte 0.25. Null: MES-like–oligodendrocyte, MES-like–neuron, OPC-like–astrocyte, OPC-like/AC-like–vascular.
Both pre-stated hypotheses (MES-like with myeloid; OPC-like with oligodendroglial lineage) are supported.

**Niches** (Leiden on spot proportions, DMG spots): 13 niches; largest is mixed AC/OPC-like (21 % of spots, 20 sections); MES-like-dominant
niches 0 and 3 (11 % together, present in 18 and 21 sections; niche 3 is MES-like + T/NK + myeloid).

**Brainstem gradient** (one autopsy, 10 sections, rostral → caudal): tumour fraction peaks in pons A–C (0.71–0.74) and falls toward the spinal
cord (0.40–0.45). Section-level trends (n = 10): OPC-like absolute abundance falls caudally (ρ = −0.70, p = 0.025), as do OC-like (−0.66)
and cycling (−0.65); oligodendrocyte proportion rises (ρ = 0.69, p = 0.029). Within section (RE over 10 sections), higher tumour fraction goes
with a larger OPC-like share (ρ = 0.35, 10/10 sections) and a smaller MES-like share (ρ = −0.36, 0/10 positive), i.e. MES-like is enriched at the
infiltrative edge. Single patient: descriptive, not generalisable.

**Ren DMG vs GBM** (5 vs 5 sections, exact Mann-Whitney): DMG has higher tumour fraction (median 0.52 vs 0.38, p = 0.032) and a lower cycling
share (0.12 vs 0.20, p = 0.008); other state shares do not differ (p >= 0.31). No multiple-testing correction across the 12 variables.

## 5. Microenvironment subclustering (`17_tme_subclusters.py`; Modal CPU, 2.7 min, peak RSS 6.0 GB)
Marker sets: `dmg_atlas.states.TME_SUBTYPE_MARKERS` (added this phase; not part of `all_marker_genes()`, so the cell2location gene set
is unchanged). Compartments (qc_pass, not consensus-malignant): myeloid 31,508, T/NK 792, vascular 4,507, oligodendroglia 16,458 cells
→ 17 / 8 / 11 / 10 Leiden subclusters on the atlas scVI latent (`tme_subclusters.csv`).
- **Assay still dominates.** Most subclusters are almost purely nuclei or purely cells (frac_nuclei ≈ 0 or ≈ 1), even on the
  batch-corrected latent. Cell-derived clusters: myeloid_3/11/13/15, T_NK_5/6/7.
- **Credible subtypes:** proliferating myeloid (myeloid_12, MKI67/TOP2A, 47 samples, 5 publications); APOE/APOC1/HLA-DR TAM
  (myeloid_16); immediate-early microglia (myeloid_15, FOS/EGR1/JUN, likely dissociation response, cells only); capillary endothelium
  (vascular_3, SLC2A1/FLT1/TFRC); pericytes (vascular_6, PDGFRB/CALD1; vascular_9, RGS5/NDUFA4L2); committed OPC / newly formed
  oligodendrocytes (oligodendrocyte_7, TCF7L2/ITPR2/TNS3/BCAS1 score 1.44; 51 % consensus-ambiguous, so possibly mixed with OPC-like
  tumour).
- **Removed from the credible list (audit MC-21, 2026-10-20):** SPP1/TREM2/C1QB TAM (myeloid_13) and CD8 T (T_NK_7). Both are
  cell-derived clusters by the bullet above (myeloid_3/11/13/15, T_NK_5/6/7), so their identity is confounded with assay and no
  marker-score justification was recorded for them; they are descriptive only until a marker-score comparison across assay classes
  is run.
- **Not subtypes (contamination / doublets / stress), flagged rather than removed:** myeloid_10 and oligodendrocyte_9 (T-cell genes
  SKAP1/PTPRC/CAMK4), vascular_10 and oligodendrocyte_0 (myeloid genes DOCK8/SPP1/APBB1IP), vascular_0, myeloid_14 and oligodendrocyte_6
  (neural/glial genes CADM2/LSAMP/NRXN1, ambient RNA or doublets), heat-shock clusters in every compartment (myeloid_1/2/9, T_NK_2,
  vascular_1, oligodendrocyte_1; HSPA1A/HSPH1). Many nuclei-derived myeloid clusters have only low-confidence marker labels because
  nuclear transcripts (long genes such as DPYD, PLXDC2, DOCK4) dominate their DE lists.
- T/NK is small (792 cells; 3–5 publications per cluster); its subclusters are descriptive only.

**Audit corrections (2026-10-20; MC-18, MC-21, MC-22).** Three findings registered against this track carried no closure row until now.

- **MC-21 — closed.** myeloid_13 and T_NK_7 are out of the credible list (bullet above).
- **MC-18 — partially closed.** The one TME treatment hit, `vascular_5:endothelial_tip_angiogenic` lower with ONC201, survives the
  re-run: coefficient −2.174, p 0.0047, mixed-model p 0.0044, BH q 0.0469, median proportion 0.239 → 0.030. Distinct-patient counts
  are now in `composition_results.csv` and the test rests on **17 samples from 17 patients (9 standard vs 8 ONC201)** — below the
  n ≥ 50 floor the audit proposed for TME subcluster sets, so **it is reported as exploratory and carries no interpretation here**.
  It is one of 37 estimable TME × ONC201 tests. Not done: the precision-weighted or binomial-GLM refit the audit also offered as an
  alternative; the reported coefficients remain unweighted.
- **MC-22 — not closed.** The assay-purity drop (myeloid 0.95 → 0.02) still rests on a per-compartment scVI fitted with assay_class as
  a covariate, with no held-out-gene independence test, no covariate-free comparison run and no ELBO trace recorded. Closing it
  requires refitting the per-compartment models, which was not attempted.

## 6. Composition models (`18_composition.py`; local, 25 s, peak RSS 2.4 GB)
**Audit corrections (2026-10-19; MC-04, MC-06, MC-19, MC-23).** Malignant cells = primary selection (167,012; no scANVI rescue). Malignant
state = argmax of the score_genes scores **z-scored within sample** (scoring per sample): OPC-like 42,685, AC-like 34,495, MES-like 31,861,
OC-like 29,954, cycling 19,558, unassigned 8,459 (`malignant_state_assignment.csv.gz`; the raw-score argmax, kept as `state_raw`, gives OPC-like
84,979 / OC-like 1,105 — it favoured the set with the highest mean expression). State shares by assay class are now close (cells vs nuclei:
AC-like 0.21 / 0.21, MES-like 0.17 / 0.19, OC-like 0.17 / 0.18, OPC-like 0.27 / 0.25, cycling 0.15 / 0.11). New set `mp_family`: metaprogram
family of §8 with within-assay-class standardisation (50 samples). Every contrast reports `n_patients` per arm; the H3 contrast is also fitted
with site as a covariate (`coef_site_adj`, `q_site_adj`) and restricted to pons (`coef_pons`, `q_pons`).
Sets: cell type (55 samples), malignant state (50), MP family (50), TME subclusters (myeloid 49, oligodendroglia 39, vascular 26, T/NK 8).
Primary model: empirical-logit proportion ~ contrast + publication fixed effects (+ assay where it varies within a publication),
patient-clustered SE; secondary: MixedLM with publication random intercept (`p_mixed`). **396 tests, 275 estimable.**
**Estimability** (`composition_estimability.csv`, now with patients and site / assay per arm): autopsy vs resection **not estimable**
(autopsy = Okada only); H3.1 vs H3.3 rests on Jessa 2022 only — **H3.1 = 8 samples / 7 patients, all pons, all nuclei; H3.3 = 23 samples / 19 patients,
13 pons + 5 thalamus (cell-type set), 4 cell samples**; ONC201 vs standard care on Okada only (10 vs 11 autopsies, 10 vs 11 patients);
adult vs paediatric: 7 adult samples (6 patients, all nuclei).
- **Assay (nuclei vs cells, within Jessa 2022):** nuclei have more neurons (Δlogit 3.5, q < 0.001) and oligodendrocytes (2.3, q 0.03);
  cells are myeloid-rich (median 37 % vs 8 %; Δ −2.2, q 0.13). With the corrected selection and within-sample states **no malignant-state
  difference remains** (pre-audit: more "unassigned" malignant cells in nuclei, Δ 2.1, q 0.018 — lost). The MP-family assay contrast is not
  interpretable by construction (family scores standardised within assay class) and is blanked in the figure / flagged in the table.
- **H3.1 vs H3.3** (7 vs 23 samples with ≥ 200 cells; 6 vs 19 patients; Jessa only): neurons higher in H3.1 (Δ 2.6, q 0.004; **site-adjusted
  q 0.017; pons-only (7 vs 16 samples) q 0.050**). The astrocyte (q 0.038), myeloid (0.038), vascular (0.051) and malignant-fraction (0.068)
  effects **do not survive site adjustment (q 0.17–0.27) or the pons restriction (q 0.28–0.29)**: every H3.1 sample is a pons nucleus sample,
  so these effects are confounded with site / tumour content and are not claimed. No malignant-state or MP-family difference (q ≥ 0.32).
- **Age (adult vs paediatric):** no cell-type or malignant-state effect (q ≥ 0.42); MP-family AC-like higher in adults (Δ 2.3, q 0.013; 6
  adult samples / 5 patients, all nuclei) — single estimable publication set, exploratory. TME subcluster effects as before rest on 2 adult vascular samples.
- **Site (thalamus vs pons):** nothing at q < 0.1 (15 thalamic samples, 14 patients).
- **ONC201 vs standard care** (Okada): no cell-type, malignant-state or MP-family shift (q ≥ 0.72); the pre-audit TME hit (tip/angiogenic
  endothelium vascular_5, q 0.047) is unchanged (TME labels were not re-derived; MC-18 caveats apply).
Changes vs pre-audit (`composition_results_vs_pre_audit.csv`, cell-type + malignant-state sets): 7 → 6 tests at q < 0.05; one lost (malignant
"unassigned", assay), none gained. Figure `phase6_composition_effects.png`: colour scale spans the full coefficient range (±3.5); glyphs
· q < 0.1, * q < 0.05, ** q < 0.01 (BH within set × contrast); rows are MP family / cell type / malignant state.

## 7. Per-cohort spatial objects (`14_cell2location.py package`)
`data/processed/spatial/spatial_deconvolution_<cohort>.h5ad` for brainstem (34,939 spots), collot_st (5,296), huchede (1,650) and ren
(27,233): Visium counts restricted to the reference genes + `obsm['q05_cell_abundance_w_sf']`, `obsm['means_cell_abundance_w_sf']`,
`obs['niche']` (DMG spots; 'NA' for Ren GBM / peritumour), `uns['c2l_run']`. The Modal image has pandas 3.0.6 + anndata 0.11.4, which
cannot write Arrow-backed strings, so `package()` switches off `future.infer_string`.

## 8. Per-cell metaprogram scores, families, Jessa comparison (`19_metaprogram_scores.py`; local, 75 s, peak RSS 5.6 GB)
UCell of the 15 MP gene sets on the 167,012 primary malignant cells (+ up to 3,000 normal cells per lineage; streamed from disk) →
`data/processed/metaprogram_scores.parquet`.
**Audit corrections (2026-10-19; MC-01 critical, MC-07, MC-08, MC-10, MC-11, P7-04).**
- **Contamination check for every MP** (`metaprogram_families.csv`: `ratio_myeloid`, `ratio_vascular` = median UCell in 3,000 normal myeloid /
  vascular cells over the malignant median; flags at 1.5× and 2×). Exclusion keeps the pre-registered 2× rule, but technical programmes are
  excluded by class regardless of the ratio: MP2 ribosomal (myeloid ratio 1.62), MP12 mitochondrial (1.86), MP13 IEG/dissociation (1.55). Lineage
  MPs: all ratios ≤ 1.43 except the cycling MP1 vascular ratio of 1.50, which is noise (median UCell 0.02 in every compartment) — the reason the
  threshold was not lowered to 1.5× on medians. Pre-audit MP3 / MP17 / MP18 (1.8× under the 2× threshold) correspond to the now-excluded
  ribosomal / mitochondrial programmes.
- **Families** (from the published / curated label): OPC_like (MP3, MP6), AC_like (MP7, MP10, MP14, MP15), OC_like (MP5), MES_like (MP8), cycling
  (MP1, MP11), hypoxia (MP4), novel_MP9 (cilia). Ribosomal, mitochondrial and IEG programmes are in no family (the pre-audit OPC-like family
  mixed MP3 ribosomal, MP17 mitochondrial and MP14 IEG with lineage programmes).
- **Assignment** = argmax of family scores, family score = **mean** (not max) of the member MPs' UCell scores **z-scored within sample**;
  unassigned if max ≤ 0. The pre-audit column (pooled z-scoring, max over members) was an assay classifier: scRNA cells 66.7 % OPC-like, 0.0 %
  MES-/OC-like vs nuclei 17–34 % MES-like, 20–29 % OC-like. Now (`metaprogram_family_by_assay.csv`), cells vs nuclei: AC-like 0.12 / 0.13,
  MES-like 0.13 / 0.15, OC-like 0.17 / 0.15, OPC-like 0.16 / 0.16, cycling 0.13 / 0.10, hypoxia 0.13 / 0.14, cilia 0.12 / 0.13, unassigned 0.05 / 0.04.
  **Caveat:** within-sample standardisation equalises the families inside every sample by construction (each sample's cells are spread over the
  families by relative rank), so this column measures a cell's relative state within its tumour, not tumour composition. A second column,
  `metaprogram_assayz` (z-scores within assay class), keeps between-sample variation and is what §6 uses for the `mp_family` composition set;
  by assay class it gives AC-like 0.15 / 0.14, MES-like 0.15 / 0.21, OC-like 0.18 / 0.15, OPC-like 0.13 / 0.19, cycling 0.13 / 0.10, hypoxia 0.08 / 0.12,
  cilia 0.14 / 0.05.
- `results/tables/metaprogram_family_reference_classes.csv` (cell_id, family, max_family_z, family_assayz, confident = max z ≥ 1; 160,078
  assigned cells, 106,123 confident) is the input the spatial track builds its cell2location metaprogram reference from; the pre-audit
  `c2l_reference_classes_mp.csv` (§10) was built from the assay-determined column and is superseded.
**Jessa 2022 author cNMF** (within-sample Spearman; 25 samples, 168 author programmes; heatmap keyed by full sample_id — the pre-audit pivot
collapsed 22 pairs of the four twice-profiled tumours, MC-10): our per-sample cNMF recovers an author programme (best usage ρ ≥ 0.5) for
**62 % of scRNA (n 26, median ρ 0.61), 64 % of snRNA (n 75, 0.60) and 43 % of multiome (n 67, 0.45)** programmes — reported by assay, not by
GEO-series label (MC-11; the pre-audit "scRNA 59 %" pooled 82 snRNA programmes under `jessa_scrna`). Best-matching MP ρ medians by family:
AC-like 0.61 (n 30), OC-like 0.61 (12), OPC-like 0.57 (30), MES-like 0.52 (8), hypoxia 0.51 (12); cycling author programmes match the cycling family weakly
(0.27, n 29; per-cell cycling usage is sparse). 43 author programmes have their best match in an excluded technical MP (median ρ 0.23).

## 9. Cell vs nucleus: assay-masked metaprograms (`20_assay_masked_metaprograms.py`; local, 4 s) — **exploratory**
**Audit corrections (2026-10-19; MC-05, MC-15).** Primary selection (74,164 Jessa malignant cells, 4 cell + 22 nuclei samples); random-mask
control added; ≥ 2-samples-per-assay rule reported; 0.01-step cuts; explicit ribosomal regex. The analysis is exploratory for a structural
reason: the mask is by definition the set of genes that separate the assays (Welch test over 4 cell vs 23 nuclei samples, two of the four cell
samples being the paired patients), so a cross-assay re-test after removing them is partly tautological, and it leaves nuclei programmes with a
median of ~21 of 50 genes.
Assay-biased genes from Jessa 2022 pseudobulk (|log2FC| ≥ 1 in both paired patients AND Welch q < 0.05, plus MT-/RP genes): 3,964 nuclei-biased,
1,957 cells-biased, 92 mito/ribo → **6,015 masked genes** (of 13,061 tested; 46 %). Masking removes a median 54 % of each programme's top-50
genes (nuclei 58 %, cells 33 %).
Re-clustering with the unchanged rule (`metaprogram_masking_summary.csv`): at the primary cut **12 masked MPs, 3 replicated** (unmasked 15 / 2):
cycling (mMP1), hypoxia (mMP3) and OPC-like mMP4 (TNR, CA10, SMOC1, MYT1, PDGFRA, OLIG1; 11 nuclei samples + **1 cell sample**, so it fails the
≥ 2-samples-per-assay rule: **strict 2**). The AC-like masked programmes stay assay-split (mMP6 / mMP9 nuclei-only; mMP12 3 cell + 2 nuclei samples
but 2 publications), OC-like / COP (mMP2, 20 nuclei samples) and MES-like (mMP5) stay nuclei-only. **Random-mask control** (same number of
genes removed per programme, random identities, 10 seeds; `metaprogram_masking_random_control.csv`): median 2.5 MPs [1–3] and 1 replicated
[1–2] at 0.85 — the real mask produces more MPs and more replication than a size-matched random mask, i.e. the effect is gene-specific, but
only one extra cross-assay programme is gained and it rests on a single cell sample. Cut sensitivity: masked replicated = 1 / 2 / 3 / 4 / 6 at
0.75 / 0.80 / 0.85 / 0.87 / 0.90 (strict 1 / 2 / 2 / 2 / 4; unmasked 1 / 1 / 2 / 2 / 2).
Cell-only vs nucleus-only MP Jaccard (`metaprogram_cell_vs_nucleus.csv`): the only cell-only MP is the IEG programme MP13, whose highest
masked Jaccard to any nuclei MP is 0.125 (to the mitochondrial MP12) — there is no cell-only lineage programme to compare, unlike the
pre-audit run (MP13 AC-like vs MP7 0.12 → 0.20).
**Eventtrust sensitivity set** (`metaprograms_eventtrust.csv`, `metaprogram_cut_sensitivity_eventtrust.csv`; 51 samples, 194,926 cells, 400
programmes re-derived by cNMF): 15 MPs, **1 replicated at 0.85** (cycling; hypoxia splits into a cell-dominated MP12, 4 cell + 1 nucleus
sample, and a nuclei MP6), 3 at 0.86–0.89 (hypoxia merges and an AC-like APOE/AQP4/CLU programme joins; strict 2), 4 at 0.90, 6 at 0.91–0.92.
A cell-only AC-like MP (MP13, 5 cell samples, 22 published AC-like genes) re-appears in this set. The replication count is therefore 1–2 at the
pre-registered cut depending on the selection and 2–3 within ±0.02 of it; the qualitative statement — cycling and hypoxia are shared across
assays, lineage programmes are not at the pre-registered cut — holds for both selections.

## 10. Brainstem gradient on metaprogram families *[superseded 2026-10-20 by SP-06: the comparison against the fallback reference is
no longer meaningful because that reference was dropped, and the deconvolution was re-run to 20,000 epochs; `brainstem_gradient_mp.csv`,
`brainstem_section_means_mp.csv`, `brainstem_mp_vs_fallback.csv` and `c2l_signatures_mp.csv` still carry the pre-audit `novel_MP15`
family and are read by no current figure. Corrected brainstem numbers: 'Audit corrections — spatial sections 3, 4, 7, 10' below.]* (`13_reference_signatures.py --use-metaprograms --tag mp`, `modal/phase6_c2l.py`, `21_brainstem_gradient_mp.py`; Modal A10G)
Reference: 13 classes (6 normal + 7 malignant families incl. hypoxia and cilia), cells with max family z ≥ 1 (105,902 of 169,492
candidates), cap 5,300 per class, 59,993 cells (`c2l_reference_classes_mp.csv`; the h5ad itself, `modal/inputs/c2l_reference_mp.h5ad`, stayed
in the sandbox and is regenerated by the same command). Signatures 250 epochs (`c2l_signatures_mp.csv`, `models/c2l_regression_mp/`).
Brainstem deconvolution: 10,154 epochs in 2,209 s (budget 2,400 s).
- Section tumour fraction agrees with the fallback-reference run (Spearman 0.98 over 10 sections).
- Rostral → caudal: cycling declines (ρ = −0.65, p = 0.043), OC-like (−0.56, p = 0.09) and OPC-like (−0.49) trend down, oligodendrocytes up
  (proportion ρ = 0.72, p = 0.019); none survives BH across the 14 classes (q ≥ 0.34).
- Within section (RE over 10 sections): higher tumour fraction ↔ larger OC-like (COP-like) share (ρ = 0.36, 10/10) and cilia share (0.15);
  lower tumour fraction (infiltrative edge) ↔ larger MES-like (ρ = −0.31, 1/10 positive) and hypoxia (−0.23, 0/10) shares. MES-like-at-the-edge
  replicates the fallback result; hypoxia at the edge is new. Single autopsy: descriptive.

## 11. Cleaned microenvironment clustering (`22_tme_scvi_clean.py`; Modal A10G, 12 min)
Exclusion (pre-specified; **heat-shock changed to flag-only after a first attempt** removed mainly the scRNA clusters, which would have
inflated assay mixing): other-lineage contamination myeloid_10 (T-cell), myeloid_5 (oligodendrocyte), myeloid_8 (astrocyte), vascular_10
(myeloid), oligodendrocyte_0 (myeloid), oligodendrocyte_9 (T-cell), plus small vascular_2 / vascular_7. Per-compartment scVI (batch
sample_id, assay covariate, encode_covariates, 2,000 HVGs after removing the 6,080 assay-biased genes) — `models/scvi_tme_<compartment>/`.
**Assay purity** (share of cells in clusters with < 10 % or > 90 % nuclei; `tme_assay_purity.csv`): myeloid 0.95 → **0.02**, T/NK 1.00 → 0.33,
vascular 0.80 → 0.82, oligodendroglia 0.99 → 1.00 (both almost entirely nuclei, so purity cannot fall).
Clean subclusters (`tme_subclusters_clean.csv`), named by top genes because the marker-set rule labels every myeloid cluster
'lipid_phagocytic_TAM' (APOE-family genes are broadly expressed): myeloid c0 homeostatic microglia (CX3CR1, P2RY12, C3; 11,003 cells, 56
samples), c1 TGFBI/AOAH macrophages, c2 hypoxic / glycolytic TAM (GPNMB, PLIN2, BNIP3L, HK2, PGK1), c3 interferon (MX1, ISG15, IFIT1/3),
c4 cell-only ribosomal residue (4 % nuclei), c5 proliferating (MKI67, TOP2A); T/NK: CD8 (CCL5, CD3D), CD4 (IL7R, LTB, CD4), NK (NKG7, GNLY,
KLRD1/F1), THEMIS/CAMK4 nuclei cluster; vascular: capillary EC (FLT1, ABCG2; TSHZ2/ABCB1), tip / angiogenic EC (NOX4, ANGPT2), pericytes
(PDGFRB, CALD1), and c0 still carrying neural genes (CADM2, LSAMP); oligodendroglia: four mature-oligodendrocyte clusters and committed OPC
(c4: TCF7L2, ITPR2, TNS3). Near-complete myeloid mixing could partly reflect over-correction by the assay covariate; the cell-only c4 shows
assay structure is not erased entirely.

## 12. Documentation (`23_package_phase6.py`)
`docs/phases/phase6_outputs.md` (generated): all 122 Phase-6 files with producing script, size, md5 and structure;
`results/tables/phase6_output_checksums.csv`; README section 'Phase 6 outputs' appended (existing text untouched). `docs/schema/atlas_schema.md`
is not regenerated because Phase 6 does not modify atlas_core.h5ad / query_mapping.h5ad.
Infrastructure changes made for the Modal runs (results unaffected): `dmg_atlas.utils` disables pandas-3 Arrow string inference (anndata
could not write those strings); `13_reference_signatures.py` reads counts with h5py blocks instead of anndata `sparse_dataset` row indexing
(broken with the image's scipy/anndata) and gained `--use-metaprograms` / `--tag`; `modal/phase6_cpu_run.py` gained `--gpu`, `--attach`,
volume links and extra staging.

## Outputs
Tables (`results/tables/`): `malignant_selection_by_assay.csv`, `cnmf/sample_rerun_plan.csv`, `metaprograms.csv`, `metaprogram_genes.csv`, `metaprogram_members.csv`, `replication_matrix.csv`,
`metaprogram_reference_overlap.csv`, `metaprogram_cut_sensitivity.csv`, the `*_eventtrust.csv` counterparts, `metaprogram_families.csv`, `metaprogram_family_composition.csv`, `metaprogram_family_composition_assayz.csv`, `metaprogram_family_by_assay.csv`, `metaprogram_family_reference_classes.csv`, `metaprogram_masking_summary.csv`, `metaprogram_masking_random_control.csv`, `composition_results_vs_pre_audit.csv`, `genesets/phase6_curated_gene_sets.csv`, `c2l_deconvolution_runs.json`, `spatial_abundance_long.csv`,
`spatial_colocalisation.csv`, `spatial_colocalisation_per_section.csv`, `spatial_niches.csv`, `spatial_niches_by_section.csv`,
`brainstem_gradient.csv`, `brainstem_section_means.csv`, `ren_dmg_vs_gbm_composition.csv`, `ren_section_composition.csv`,
`tme_subclusters.csv`, `tme_subcluster_markers.csv`, `tme_subcluster_marker_scores.csv`, `composition_per_sample.csv`,
`composition_results.csv`, `composition_estimability.csv`, `malignant_state_assignment.csv.gz`.
Figures (`results/figures/`): `phase6_metaprogram_jaccard.png`, `phase6_replication_matrix.png`, `phase6_spatial_abundance_maps.png`,
`phase6_colocalisation_matrix.png`, `phase6_brainstem_gradient.png`, `phase6_tme_umap.png`, `phase6_tme_marker_scores.png`,
`phase6_composition_effects.png`. Per-spot / per-cell data: `data/processed/spatial/` (incl. `spatial_deconvolution_<cohort>.h5ad`),
`data/processed/tme_subcluster_labels.csv.gz`. Logs: `logs/phase6_tme.log`, `logs/phase6_composition.log` (these also hold the
lines of failed attempts; the last START…done block of each is the run that produced the outputs).

**Modal runs for §5–7.** `atlas_core.h5ad` uploaded to volume `dmg-atlas-data` (sha256 01b7b7be…a9809, verified inside the sandbox).
Failed attempts before the final runs: anndata `sparse_dataset` row indexing broken in the image (replaced by h5py block reads);
a pandas `DataFrame.flags` name clash; a missing figures directory; Arrow-string writes. None changed results.

## Open items
- Audit items not closed in this track: k-selection ceiling / non-convergence (MC-12), version control (MC-13), TME subcluster items (MC-18, MC-21, MC-22) and the §5 figure titles (MC-20 for `phase6_tme_umap*.png`, still open — see below); MC-14 was closed by re-running `23_package_phase6.py`; §3, 4, 7, 10 (spatial) were corrected by the spatial track on 2026-10-20.
- **MC-20 (`phase6_tme_umap.png`, `phase6_tme_umap_clean.png`) remains open after the 2026-10-21 reporting pass.** Both figures still
  carry per-panel titles (`myeloid (n=31,508)` …) inside the plotting area and have no axis labels. They cannot be redrawn from saved
  state: `17_tme_subclusters.py` and `22_tme_scvi_clean.py` hold the per-compartment UMAP coordinates only in memory (`umaps[comp]`) and
  never write them, so regeneration requires a full per-compartment scVI refit. Every other figure in `results/figures/` and the whole
  publication set are in house style.
- Cell arm power: 6 scRNA samples / 3 publications; the replication test (§2, §9) cannot be decisive until more cell-based DMG data are added.
- ~~Spatial tests for the other three cohorts (Huchedé, Ren, Collot) still use the fallback-state reference; only the brainstem was
  re-deconvolved with metaprogram families.~~ *[closed 2026-10-20 by SP-06 — all four cohorts now use the single corrected
  metaprogram-family reference and the fallback reference is used nowhere.]*
- TME marker-set annotation needs a more specific myeloid panel (the current one labels every cleaned myeloid cluster lipid_phagocytic_TAM);
  vascular_c0 neural contamination remains.


## Audit corrections — spatial sections 3, 4, 7, 10 (2026-10-20, lead agent)

Findings SP-01 to SP-06, SP-11, SP-12. All four Visium cohorts were re-deconvolved with the corrected
within-sample-z metaprogram-family reference (13 classes: 7 malignant families + 6 normal lineages). Training is
capped by a wall-clock budget (`--budget-s 2400`, `max_epochs = clip(budget_s / probed s-per-epoch, 2000, 20000)`),
so the realised epoch counts are the ones logged in `results/tables/c2l_histories/`: huchede, collot_st and
brainstem each reached the 20,000-epoch target, **Ren stopped at 12,940 epochs (65 % of target) when its budget
expired** — still well above the 8,185/9,750 (41-49 %) of the pre-audit runs, which also used a partly
assay-determined reference. Convergence is therefore argued from the harvested histories rather than from reaching
the target: the last 5 % of logged epochs contribute 0.012 % (brainstem), 0.024 % (Ren), 0.058 % (collot_st) and
0.069 % (huchede) of the total ELBO reduction — i.e. under 0.07 % in every cohort. The score_genes-argmax fallback reference is no longer used anywhere.

**Co-localisation (section 4) is now computed on spot proportions, pooled over PATIENTS (12 DMG patients, 21
sections; the 10 brainstem sections count as one autopsy and the Collot core/edge pair as one patient), BH over the
42 state x TME pairs (7 corrected malignant families x 6 normal lineages; the pre-audit family was 30 pairs over 5 states).** Raw-abundance correlations are kept as a secondary table because every class correlates
with the spot total, which tracks UMI depth.

| pair | raw q05 (old) | proportion (new) | q (BH) | sections positive |
|---|---|---|---|---|
| AC-like x astrocyte | 0.50 | 0.393 (0.234-0.531) | 1e-4 | 13/21 |
| AC-like x oligodendrocyte | - | 0.308 (0.138-0.459) | 0.005 | 11/21 |
| MES-like x vascular | - | 0.253 (0.078-0.412) | 0.034 | 13/21 |
| MES-like x myeloid | 0.33 | 0.238 (0.071-0.392) | 0.034 | 14/21 |
| OC-like x oligodendrocyte | 0.36 | 0.301 (0.036-0.526) | 0.087 | 12/21 |
| OPC-like x oligodendrocyte | 0.26 | 0.178 (-0.028-0.370) | 0.210 | 10/21 |

Of the two pre-stated hypotheses, **MES-like with myeloid survives** (q = 0.034, effect about 70 % of the raw
estimate); **OPC-like with the oligodendroglial lineage does not** (q = 0.21), although the related OC-like/COP
pairing is marginal (q = 0.087). The strongest association in the data is AC-like with astrocytes.

**Brainstem (sections 4 and 10, one autopsy, 10 serial sections - descriptive only).** With the single corrected
reference, within-section tumour fraction rises with the OPC-like share (rho_RE 0.213, q 0.0013, positive in 8/10
sections) and the AC-like share (0.285, q 0.0001), and falls with oligodendrocyte (-0.298), astrocyte (-0.431),
vascular (-0.235) and neuron (-0.239) abundance. The earlier section-4/section-10 discrepancy (OPC-like vs
OC-like) came from using two different references and no longer arises.

**Ren DMG vs GBM (section 4).** Recomputed on the new reference with BH across 14 variables: tumour fraction
0.528 vs 0.460 (p = 0.056, q = 0.667) and cycling share 0.112 vs 0.214 (p = 0.095, q = 0.667). **Neither difference
survives correction**; the earlier uncorrected p-values (0.032 and 0.008) should not be cited.

**Figures.** `phase6_spatial_abundance_maps.png` is redrawn on pixel coordinates with per-panel 1 mm scale bars
(sections pDMG_Sample-2, DMG1, BT01_192_core_HE, Pons_A - named here, not in the panels);
`phase6_brainstem_gradient.png` has a shared legend below both panels and labelled axes;
`phase6_colocalisation_matrix.png` shows the proportion-based estimates.

Code: `scripts/13b_c2l_reference_mpfam.py`, `modal/phase6_c2l.py`, `scripts/14c_collect_fix.py`,
`scripts/15b_spatial_tests_fix.py`, `scripts/15c_spatial_figures_fix.py`. Pre-audit tables are kept as
`*.pre_audit.csv` and the previous abundances under `data/processed/spatial_pre_audit/`.
