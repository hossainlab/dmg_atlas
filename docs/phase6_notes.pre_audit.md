# Phase 6 — malignant metaprograms (cNMF) and Visium deconvolution / spatial tests (2026-10-17)

Code: `scripts/13_cnmf_per_sample.py` (per-sample consensus NMF), `scripts/16_metaprograms.py` (metaprogram clustering + replication test),
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

Deviation from the Phase 6 brief: the per-sample NMF uses sklearn NMF with random initialisation and 10 restarts (consensus over restarts) rather than nndsvd initialisation with 3 restarts; this was chosen for stability of the consensus clustering and was not declared in the sub-agent report. See `scripts/13_cnmf_per_sample.py` for the exact settings.
 (`results/tables/cnmf/`, `data/processed/cnmf_usages_<sample>.parquet`)
Malignant selection: `consensus_call == 'malignant'` OR (`scanvi_label == 'malignant'` AND `scanvi_maxprob >= 0.9`); samples with
>= 200 selected cells analysed: 54 of 57 (`samples_analysed.csv`). k in 4..9, 10 restarts, k* = largest k with silhouette >= 0.90 and every
consensus cluster present in >= 80 % of restarts (`k_selection.csv`). 426 programmes in total. Wall 45.7 min, peak RSS 2.47 GB.

## 2. Metaprograms and the pre-registered replication test (`metaprograms.csv`, `replication_matrix.csv`, …)
Method fixed in the `16_metaprograms.py` docstring before the clustering was run: Jaccard on top-50 genes, average linkage on 1 − J,
cut 0.85 (mean pairwise J >= 0.15), MP = cluster with >= 5 programmes from >= 3 samples. Replication rule (decisions.md #9): >= 3 studies
and both cells and nuclei. **Studies are counted by publication**: `jessa_scrna` and `jessa_multiome` are both Jessa 2022 and share patients,
so they count once (5 publications: Jessa 2022, Okada, Collot, DeSisto, Zang). Cells = scRNA (6 samples: Jessa 4, DeSisto 1, Zang 1);
nuclei = snRNA + multiome RNA (48 samples).

Result at the primary cut: **20 metaprograms (196 / 426 programmes assigned), 2 replicated**:
- MP1 cycling (Liu 2022 G2M, 29 / 50 genes shared; 34 programmes, 4 publications, 5 cell + 29 nuclei samples).
- MP10 hypoxia (19 shared genes; 4 cell + 4 nuclei samples, 3 publications).

The lineage-like programmes split by assay: OPC-like (MP2, MP3, MP5, MP11, MP17), AC-like (MP6, MP7, MP18), OC-like (MP4, MP12) and
MES-like (MP9, MP19, MP20) metaprograms are **nuclei-only**, while the scRNA samples form their own AC-like (MP13) and OPC-like-3 (MP14)
metaprograms. Cell and nucleus versions of the same state therefore cluster apart rather than together. Two unlabelled nuclei-only MPs:
MP15 (ADGB/DNAH9/SPAG17, cilia genes) and MP16 (DOCK8/APBB1IP/RHBDF2, immune genes, Jessa only; likely myeloid contamination of the
malignant selection).

Sensitivity (`metaprogram_cut_sensitivity.csv`): replicated MPs = 1 / 1 / 2 / 5 at cuts 0.75 / 0.80 / 0.85 / 0.90 (5 / 14 / 20 / 19 MPs).

**Build-plan checkpoint triggered** ("fewer than 3 metaprograms replicate across assays → the headline becomes an assay-dependence finding").
Caveat: the cell arm is small (6 samples, 3 publications), so absence of cross-assay replication is partly a power limitation.
Follow-ups not yet done: compare cell vs nucleus MPs directly (e.g. MP13 vs MP6/7, MP14 vs MP2/3/5) after removing assay-specific genes
(mitochondrial, ribosomal, nuclear-retained transcripts); UCell-score the MPs on all malignant cells.

## 3. cell2location (`results/tables/c2l_signatures.csv`, `data/processed/spatial/<cohort>/`, `spatial_abundance_long.csv`)
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
- **Credible subtypes:** proliferating myeloid (myeloid_12, MKI67/TOP2A, 47 samples, 5 publications); SPP1/TREM2/C1QB TAM (myeloid_13);
  APOE/APOC1/HLA-DR TAM (myeloid_16); immediate-early microglia (myeloid_15, FOS/EGR1/JUN, likely dissociation response, cells only);
  capillary endothelium (vascular_3, SLC2A1/FLT1/TFRC); pericytes (vascular_6, PDGFRB/CALD1; vascular_9, RGS5/NDUFA4L2); committed
  OPC / newly formed oligodendrocytes (oligodendrocyte_7, TCF7L2/ITPR2/TNS3/BCAS1 score 1.44; 51 % consensus-ambiguous, so possibly mixed
  with OPC-like tumour); CD8 T (T_NK_7).
- **Not subtypes (contamination / doublets / stress), flagged rather than removed:** myeloid_10 and oligodendrocyte_9 (T-cell genes
  SKAP1/PTPRC/CAMK4), vascular_10 and oligodendrocyte_0 (myeloid genes DOCK8/SPP1/APBB1IP), vascular_0, myeloid_14 and oligodendrocyte_6
  (neural/glial genes CADM2/LSAMP/NRXN1, ambient RNA or doublets), heat-shock clusters in every compartment (myeloid_1/2/9, T_NK_2,
  vascular_1, oligodendrocyte_1; HSPA1A/HSPH1). Many nuclei-derived myeloid clusters have only low-confidence marker labels because
  nuclear transcripts (long genes such as DPYD, PLXDC2, DOCK4) dominate their DE lists.
- T/NK is small (792 cells; 3–5 publications per cluster); its subclusters are descriptive only.

## 6. Composition models (`18_composition.py`; Modal CPU, 2.5 min, peak RSS 12.1 GB)
Sets: cell type (55 samples), malignant state (53; all 237,886 selected qc_pass malignant cells scored: OPC-like 116,630, AC-like 51,117, MES-like 41,952, cycling 16,026,
OC-like 7,219, unassigned 4,942, `malignant_state_assignment.csv.gz`),
TME subclusters per compartment (myeloid 49, oligodendroglia 39, vascular 26, T/NK 8 samples with ≥ 30 compartment cells).
Primary model: empirical-logit proportion ~ contrast + publication fixed effects (+ assay where it varies within a publication),
patient-clustered SE; secondary: MixedLM with publication random intercept (`p_mixed`). 348 tests, 234 estimable.
**Estimability** (`composition_estimability.csv`): autopsy vs resection is **not estimable** (autopsy = Okada only, resection =
DeSisto + Zang only; fully confounded with study). H3.1 vs H3.3 rests on Jessa 2022 only; ONC201 vs standard care on Okada only.
- **Assay (nuclei vs cells, within Jessa 2022):** 28 / 44 tests q < 0.05. Nuclei have more neurons (Δlogit 3.5), oligodendrocytes
  (2.3) and unassigned malignant cells; cells are myeloid-rich (median 37 % vs 8 % of cells). Same assay dominance as §2 and §5.
- **H3.1 vs H3.3** (7 vs 23 samples, Jessa only): H3.1 samples contain more normal cells: neurons (q 0.002), astrocytes (0.03),
  myeloid (0.03), vascular (0.03); malignant fraction lower (median 0.79 vs 0.91, q 0.057). No malignant-state difference.
- **Age (adult vs pediatric):** no cell-type or malignant-state effect (q ≥ 0.22; AC-like higher in adults, Δ 1.1, q 0.22). Seven TME
  subcluster effects reach q < 0.05, but the vascular ones rest on 2 adult samples and are not interpretable.
- **Site (thalamus vs pons):** nothing at q < 0.1.
- **ONC201 vs standard care** (Okada, 10 vs 11 autopsies): no cell-type or malignant-state shift (all q ≥ 0.8); one TME hit
  (tip/angiogenic endothelium vascular_5 lower with ONC201, median 0.03 vs 0.24, q 0.047, p_mixed 0.004) — single test, single study.

## 7. Per-cohort spatial objects (`14_cell2location.py package`)
`data/processed/spatial/spatial_deconvolution_<cohort>.h5ad` for brainstem (34,939 spots), collot_st (5,296), huchede (1,650) and ren
(27,233): Visium counts restricted to the reference genes + `obsm['q05_cell_abundance_w_sf']`, `obsm['means_cell_abundance_w_sf']`,
`obs['niche']` (DMG spots; 'NA' for Ren GBM / peritumour), `uns['c2l_run']`. The Modal image has pandas 3.0.6 + anndata 0.11.4, which
cannot write Arrow-backed strings, so `package()` switches off `future.infer_string`.

## 8. Per-cell metaprogram scores, families, Jessa comparison (`19_metaprogram_scores.py`; Modal CPU, 1.8 min, peak RSS 12.1 GB)
UCell of the 20 MP gene sets on 237,886 selected malignant cells (+ up to 3,000 normal cells per lineage) →
`data/processed/metaprogram_scores.parquet`. **Contamination rule revised after the first attempt** (documented in the script): comparing
against every lineage flagged all 20 MPs (lineage-like tumour programmes score higher in their normal neural counterparts; the 356
T/NK reference cells gave implausible scores), so only myeloid / vascular are used as comparators (non-malignant if median ≥ 2× the
malignant median). Only MP16 (DOCK8/APBB1IP; myeloid median 0.377 vs malignant 0.029) is excluded.
Families (from the MP reference label): OPC_like (MP2, 3, 5, 11, 14, 17), AC_like (MP6, 7, 13, 18), OC_like (MP4, 12), MES_like (MP9, 19, 20),
cycling (MP1), hypoxia (MP8, 10), novel_MP15 (cilia). Argmax of z-scored family scores spreads cells thinly (OC_like 59,369, MES_like
42,905, OPC_like 27,926, hypoxia 27,601, AC_like 26,471, cycling 22,254, cilia 18,846, unassigned 12,514); the cilia and OC-like shares are
implausibly large, so this assignment is used only through the confident subset (max z ≥ 1) for the cell2location reference (§10).
**Jessa 2022 author cNMF** (programme numbering is per sample, so all correlations are within sample; 26 samples, 175 author programmes):
our per-sample cNMF recovers an author programme (best usage Spearman ≥ 0.5) for 59 % of scRNA and 45 % of multiome programmes (median best
ρ 0.61 / 0.45); best-matching MP ρ ≥ 0.5 for 45 % overall, highest for AC-like (median 0.64), OPC-like (0.58), OC-like (0.57), MES-like (0.53);
cycling author programmes match MP1 weakly (0.32, per-cell cycling usage is sparse). Programmes whose best match is the excluded MP16 have
ρ < 0.25 (no real match).

## 9. Cell vs nucleus: assay-masked metaprograms (`20_assay_masked_metaprograms.py`; Modal CPU, < 1 min)
Assay-biased genes from Jessa 2022 pseudobulk of malignant cells: |log2FC| ≥ 1 in the same direction in both patients profiled both ways
(P-6337, P-6640; scRNA vs multiome nuclei) AND Welch q < 0.05 across 4 cell vs 23 nuclei samples, plus MT-/ribosomal genes: 3,987
nuclei-biased, 1,982 cells-biased, 108 mito/ribo → 6,080 masked genes (of 13,185 tested). Masking removes a median 56 % of each programme's
top-50 genes (nuclei programmes 58 %, cell programmes 31 %), i.e. the mask is aggressive and certainly removes some biology.
Re-clustering with the unchanged rule: at the primary cut **17 MPs, 4 replicated** (unmasked 20 / 2): cycling (mMP1), OC-like / committed
OPC (mMP2: BCAS1, GPR17, ENPP6, MYRF; 2 cell + 16 nuclei samples), hypoxia (mMP4: VEGFA, NDRG1, ADM, SLC2A1), AC-like (mMP5: APOE, AQP4, CLU,
GJA1; 4 cell + 5 nuclei samples). Sensitivity: replicated 1 / 2 / 4 / 8 at cuts 0.75 / 0.80 / 0.85 / 0.90 (unmasked 1 / 1 / 2 / 5).
OPC-like (mMP3, mMP8, mMP15) and MES-like (mMP7, mMP9) stay nuclei-only. Cell-only vs nucleus-only MP Jaccard rises with masking for
AC-like (MP13 vs MP7 0.12 → 0.20) but stays low for OPC-like (MP14 vs MP17 0.03 → 0.09).
Interpretation: part of the assay split is technical (nuclear-retained / intronic transcripts); AC-like and OC-like/COP programmes are shared
across assays once these are masked, OPC-like and MES-like are not. The headline remains an assay-dependence finding with this cell arm.

## 10. Brainstem gradient on metaprogram families (`13_reference_signatures.py --use-metaprograms --tag mp`, `modal/phase6_c2l.py`, `21_brainstem_gradient_mp.py`; Modal A10G)
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
`docs/phase6_outputs.md` (generated): all 122 Phase-6 files with producing script, size, md5 and structure;
`results/tables/phase6_output_checksums.csv`; README section 'Phase 6 outputs' appended (existing text untouched). `docs/atlas_schema.md`
is not regenerated because Phase 6 does not modify atlas_core.h5ad / query_mapping.h5ad.
Infrastructure changes made for the Modal runs (results unaffected): `dmg_atlas.utils` disables pandas-3 Arrow string inference (anndata
could not write those strings); `13_reference_signatures.py` reads counts with h5py blocks instead of anndata `sparse_dataset` row indexing
(broken with the image's scipy/anndata) and gained `--use-metaprograms` / `--tag`; `modal/phase6_cpu_run.py` gained `--gpu`, `--attach`,
volume links and extra staging.

## Outputs
Tables (`results/tables/`): `metaprograms.csv`, `metaprogram_genes.csv`, `metaprogram_members.csv`, `replication_matrix.csv`,
`metaprogram_reference_overlap.csv`, `metaprogram_cut_sensitivity.csv`, `c2l_deconvolution_runs.json`, `spatial_abundance_long.csv`,
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
- Cell arm power: 6 scRNA samples / 3 publications; the replication test (§2, §9) cannot be decisive until more cell-based DMG data are added.
- Spatial tests for the other three cohorts (Huchedé, Ren, Collot) still use the fallback-state reference; only the brainstem was
  re-deconvolved with metaprogram families.
- TME marker-set annotation needs a more specific myeloid panel (the current one labels every cleaned myeloid cluster lipid_phagocytic_TAM);
  vascular_c0 neural contamination remains.
