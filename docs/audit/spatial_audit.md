# Independent audit — Phase 6 spatial (cell2location deconvolution and spatial tests)

Auditor role: adversarial, independent of the authors. Date: 2026-10-19. Scope as briefed: `scripts/13_reference_signatures.py`,
`14_cell2location.py`, `15_spatial_tests.py`, `21_brainstem_gradient_mp.py`, `modal/phase6_c2l.py`, `models/c2l_regression*`,
`data/processed/spatial/*`, `data/processed/spatial_mp/*`, `results/tables/c2l_*`, `spatial_*`, `brainstem_*`, `ren_*`,
`results/qc/spatial_sections.csv`, the four Phase 6 spatial figures, and `docs/phases/phase6_notes.md` §3, §4, §7, §10.
Companion table: `docs/audit/spatial_findings.csv` (20 findings: 0 critical, 6 major, 6 minor, 8 notes).

**Overall verdict for this scope: sound with corrections.** The deconvolution pipeline, tables and arithmetic reproduce exactly from
the committed files (every headline number re-derived below matches). Two inferential choices are not defensible as written and change
the conclusions' wording: (i) co-localisation was computed on raw q05 abundances, which share a cell-density/UMI-depth component, so all
30 state × TME correlations are shifted positive; on proportions or density-partialled correlations the second pre-stated hypothesis
(OPC-like with oligodendroglial lineage) is weak-to-null and OPC-like–T/NK disappears, while MES-like–myeloid survives at about half the
reported size; (ii) the within-section "OPC-like share rises with tumour fraction" result does not replicate with the metaprogram
reference (it re-appears under the OC-like label) and §10 does not say so. Two of the four figures have correctness defects (geometric
distortion of the maps; overlapping/truncated legend in the brainstem figure). Deconvolution training histories are missing from the
repository, so the convergence claim for the budget-capped runs is unverifiable.

## 1. Files examined

Scripts/modules read line by line: `scripts/13_reference_signatures.py` (223 lines), `scripts/14_cell2location.py` (209),
`scripts/15_spatial_tests.py` (319), `scripts/21_brainstem_gradient_mp.py` (121), `modal/phase6_c2l.py` (152);
`scripts/02_qc_spatial_sections.py` was listed but not reviewed (upstream of scope; its output table was checked instead).
Tables: `spatial_abundance_long.csv` (760,298 × 9), `spatial_colocalisation.csv` (110 × 18), `spatial_colocalisation_per_section.csv`,
`spatial_niches.csv`, `spatial_niches_by_section.csv`, `brainstem_gradient.csv`, `brainstem_section_means.csv`,
`brainstem_gradient_mp.csv` (via log), `brainstem_section_means_mp.csv`, `brainstem_mp_vs_fallback.csv`,
`ren_dmg_vs_gbm_composition.csv`, `ren_section_composition.csv`, `c2l_reference_classes.csv`, `c2l_reference_classes_mp.csv`,
`c2l_reference_state_scores.csv.gz` (90,000 × 6), `c2l_signatures.csv` (3,592 × 11), `c2l_signatures_mp.csv` (3,592 × 13),
`c2l_signatures_mp_history.csv`, `c2l_signatures_mp_run.json`, `c2l_deconvolution_runs.json`, `c2l_deconvolution_runs_mp.json`,
`results/qc/spatial_sections.csv`, `results/qc/spatial_sections_used.csv`.
Data: `data/processed/spatial/<cohort>/*_abundance.csv` (27 files), `data/processed/spatial_mp/brainstem/*_abundance.csv` (10),
`spatial_niche_labels.csv.gz`, obs/obsm/uns/X-dtype of the four `spatial_deconvolution_<cohort>.h5ad` (h5py, no full load),
obs and uns of `modal/inputs/c2l_reference.h5ad`, `models/c2l_regression/{signatures_history.csv,signatures_run.json}`,
`models/c2l_regression_mp/model.pt` (presence only). Figures (viewed): `phase6_spatial_abundance_maps.png`,
`phase6_colocalisation_matrix.png`, `phase6_brainstem_gradient.png`, `phase6_brainstem_gradient_mp.png`.
Docs/logs: `docs/phases/phase6_notes.md`, `docs/design/build_plan.md` (Phase 7), `docs/design/decisions.md` (grep), `docs/figures/figures.md` (Fig. 5 caption and
checksums), `docs/figures/supp_figures.md` (Supp. Fig. 7), `docs/phases/phase6_outputs.md` (md5), `logs/phase6_spatial.log` (295 lines),
`data/processed/sample_manifest.csv` (Collot patient IDs).

## 2. Method of verification

All recomputation was local (env `multiverse`, 8 threads, peak RSS 0.74 GB), from the committed tables and h5ad metadata only; no GPU,
no Modal, no downloads, no re-deconvolution. Procedure: (a) re-derive every per-spot and per-section quantity from
`spatial_abundance_long.csv` and compare against the published tables; (b) re-implement the DerSimonian–Laird pooling and BH exactly as in
the script and confirm the published numbers, then re-run the same pooling on per-spot proportions, on posterior means, and on
density-partialled ranks; (c) sensitivity of pooled estimates to pseudo-replication (leave-brainstem-out, brainstem collapsed to one
effective section); (d) re-run brainstem trend and within-section tests for both references; (e) compare fallback and MP brainstem runs
spot by spot; (f) inspect reference composition and the fallback state assignment; (g) check h5ad obsm against the long table, X dtype,
spot counts against the QC table, figure md5 against `phase6_outputs.md`; (h) view each figure.

## 3. Re-derived numbers (claimed vs recomputed)

| Quantity (source of claim) | Claimed | Recomputed | Match |
|---|---|---|---|
| Spots × classes, sections (§3) | 69,118 × 11, 27 | 69,118 × 11, 27 (brainstem 34,939; ren 27,233; collot_st 5,296; huchede 1,650) | ✓ |
| Reference size / classes (§3) | 59,940 cells, 11 classes | 59,940 (sum of `c2l_reference_classes.csv`; h5ad obs 59,940 × 8) | ✓ |
| AC-like–astrocyte pooled ρ (§4) | 0.43 (0.34–0.51), 21/21 | 0.4301 (0.3377–0.5143), 21/21, q = 2.9e-15 | ✓ |
| OC-like–oligodendrocyte / OC-like–neuron / MES-like–myeloid / OPC-like–T/NK / OPC-like–oligodendrocyte (§4) | 0.36 (18/21) / 0.32 / 0.30 (0.14–0.44, 15/21) / 0.29 / 0.25 | 0.363 (18) / 0.321 (20) / 0.296 (0.136–0.441, 15) / 0.293 (16) / 0.252 (16) | ✓ |
| Same pairs on per-spot proportions (audit) | — | 0.207 / 0.239 / 0.175 (14/21, q 0.023) / 0.027 (p 0.72) / 0.061 (p 0.36); AC–astro 0.209 (13/21) | n/a |
| Same pairs, partial Spearman given log total q05 (audit) | — | 0.242 / 0.184 / 0.150 (p 0.003) / 0.040 (p 0.55) / 0.127 (p 0.027); AC–astro 0.242 | n/a |
| Per-section ρ file vs recomputation | — | max |Δ| 2.2e-16 over 1,155 DMG rows | ✓ |
| Niches (§4) | 13; largest 21 % in 20 sections; niches 0+3 11 % in 18 and 21 sections | 13; 20.95 % (AC-like dominant), 20; 11.5 %, 18 and 21 | ✓ |
| Brainstem tumour fraction per section (§4, `brainstem_section_means.csv`) | pons A–C 0.71–0.74, spinal cord 0.40–0.45 | Midbrain 0.4585, Pons_A 0.7108, Pons_B 0.7174, Pons_C 0.7376, Pons_D 0.5360, Medulla_A 0.5733, Medulla_B 0.5447, Medulla_C 0.5614, Spinal_Cord_A 0.3995, Spinal_Cord_B 0.4538 (max |Δ| vs table 1e-16) | ✓ |
| Brainstem trends (§4) | OPC-like −0.70 p 0.025; OC-like −0.66; cycling −0.65; oligodendrocyte prop 0.69 p 0.029 | −0.697 p 0.025; −0.661 p 0.038; −0.648 p 0.043; 0.685 p 0.029. BH over 12 abs: min q 0.17; over 12 prop: min q 0.236 | ✓ (uncorrected) |
| Within-section RE (§4) | OPC-like 0.35 (10/10); MES-like −0.36 (0/10) | 0.347 (0.251–0.437), 10/10; −0.364 (−0.474 to −0.242), 0/10; posterior-mean version −0.384 | ✓ |
| MP vs fallback tumour fraction (§10) | Spearman 0.98 over 10 sections | 0.976 (Pearson 0.997); spot-level 0.987 over 34,939 spots | ✓ |
| MP run within-section (§10) | OC-like 0.36 (10/10), MES-like −0.31 (1/10), hypoxia −0.23 (0/10), cilia 0.15 | log: 0.363, −0.308, −0.231, 0.152 — and OPC-like −0.028 (5/10), not reported | ✓ / omission |
| Fallback vs MP per-class spot agreement (audit) | — | OPC-like 0.70, AC-like 0.94, OC-like 0.87, MES-like 0.89, cycling 0.99, TME 0.96–1.00; fallback OPC-like vs MP OC-like 0.96 | n/a |
| Ren DMG vs GBM (§4) | TF 0.52 vs 0.38 p 0.032; cycling 0.12 vs 0.20 p 0.008 | 0.516 vs 0.384, p 0.0317; 0.115 vs 0.195, p 0.0079 (minimum attainable for 5 v 5; BH q 0.095) | ✓ |
| Deconvolution epochs / wall (§3) | 8,185 / 9,750 / 15,729 / 18,806 | run json: 8,185 / 9,750 / 15,729 / 18,806; train_s 1,651 / 1,627 / 1,011 / 962 | ✓ |
| ELBO plateau of the deconvolution runs (§3) | last 5 % of epochs = 0.07–0.16 % of total reduction | histories absent from repo — **not verifiable** | ✗ |
| RegressionModel convergence (audit) | 250 epochs | last 5 % of epochs = 0.05 % (fallback) / 0.06 % (MP) of total ELBO reduction | ✓ |
| Fallback state assignment (log) | OPC 45,293 / AC 19,903 / MES 15,737 / cycling 6,227 / OC 1,420 / NA 1,420 | identical from `c2l_reference_state_scores.csv.gz`; z-scored argmax agrees for 68.7 % of cells | ✓ |
| Figure md5 vs `phase6_outputs.md` | 3e8f363c…, d55b9282…, caad66c4…, c2c1993b… | identical | ✓ |
| Packaged h5ad obsm vs long table | — | max |Δ| 1.1e-7 (huchede); X int32 integer counts | ✓ |

## 4. Findings in prose, by severity

### Major

**SP-01 — Co-localisation on raw abundances is confounded by spot cell density (statistics/validity).** `colocalisation()` correlates
raw q05 abundances across spots. Every class rises with the spot's total q05 (median per-section ρ with total: 0.20–0.55), and total q05
tracks UMI depth. The published matrix is therefore almost uniformly positive (28/30 cells; 25/30 q < 0.05), which is what the figure
shows. Pooling the same 21 sections on per-spot proportions halves the headline effects (AC-like–astrocyte 0.43 → 0.21; OC-like–
oligodendrocyte 0.36 → 0.21; MES-like–myeloid 0.30 → 0.18, still q = 0.023 with 14/21 positive) and removes two of the six "supported"
pairs (OPC-like–T/NK 0.29 → 0.03; OPC-like–oligodendrocyte 0.25 → 0.06, p = 0.36). Partial correlations given log total give a similar
picture (OPC-like–oligodendrocyte 0.13, p = 0.027). Pairs reported as "null" become clear segregations (OPC-like–astrocyte −0.27 on
proportions, 3/21 positive). The sentence "both pre-stated hypotheses are supported" must be weakened to: MES-like–myeloid supported with
a small effect; OPC-like–oligodendroglial lineage weak/inconclusive. Proportions induce some negative dependence of their own, so the
truth lies between the two analyses; neither supports the raw magnitudes.

**SP-02 — The §4 OPC-like within-section result does not replicate and §10 does not say so (reporting).** Fallback: OPC-like share vs
tumour fraction ρ = 0.35, 10/10 sections. MP reference: OPC-like −0.03, 5/10; the 0.36, 10/10 effect appears under OC-like. Spot-level
agreement between runs is 0.70 for OPC-like (lowest of all classes) but 0.96 between fallback OPC-like and MP OC-like — the two references
partition the OPC/OC continuum differently. The "Spearman 0.98 agreement" is for tumour fraction only; it should not be read as
agreement of state maps.

**SP-03 — Deconvolution convergence unverifiable; two runs budget-capped (reproducibility).** `phase6_c2l.py` writes `<cohort>_history.csv`,
but none of the five deconvolution histories is in the repository; only `final_elbo` is. Brainstem and Ren trained 8,185 and 9,750 of the
20,000 target epochs. The plateau statistic quoted in §3 cannot be checked. (The two RegressionModel histories that were kept do show
convergence.)

**SP-04 — `phase6_brainstem_gradient.png` legend overlaps the right panel (figure).** Left-panel legend labels are truncated to two
letters and overprint the right panel's y-tick labels; OC-like is drawn in a near-white colour; the right panel has no y label. The md5
matches the released file.

**SP-05 — `phase6_spatial_abundance_maps.png` is geometrically distorted and unlabelled (figure).** Plotting `array_col` vs `−array_row`
with equal aspect stretches Visium tissue ~1.6× horizontally (DMG1: array extent ratio 1.65 vs pixel ratio 1.06). Pixel coordinates are
available in `obs`/`obsm['spatial']` but unused. No section identifiers, scale bar or caption exist anywhere in `docs/`.

**SP-06 — Fallback malignant-state labels depend on an arbitrary scaling (validity).** Argmax over unstandardised `score_genes` scores
whose means differ >30-fold (OPC-like 0.306 … OC-like 0.009) yields OC-like for 1.6 % of candidates and MES-like for 17.5 %; a z-scored
argmax reassigns 31 % of cells. All state-level spatial results in §4 (for all four cohorts) inherit this choice; gene-centred signature
correlations of 0.70 (AC vs MES) and 0.71 (OPC vs OC) show the affected classes are also the least separable. The MP reference of §10
mitigates this for the brainstem only (already an Open item).

### Minor

SP-07 no spot-level UMI filter (1,633 spots < 500 UMIs; Spinal_Cord_A 10 % of spots < 100 UMIs; minimum 1 UMI); SP-08 brainstem trend
p-values in §4 are uncorrected (BH min q 0.17; the oligodendrocyte "p = 0.019" is uncorrected in §10 as well because BH was applied to
`p_abs` only); SP-09 `N_cells_per_location = 5` is unjustified for hypercellular glioma (median total q05 1.5–5.7 cells/spot) — only
relative statements are defensible; SP-16 mixed-timezone logs (downstream steps logged before their inputs) and harvested Modal outputs
outside the repo; SP-17 string version comparison and zero-filling of 77 genes absent from Medulla_A/B; SP-18 cosmetic figure issues
(palette and label inconsistency between the two brainstem figures, '-0.00', no captions).

### Notes (confirmed sound or immaterial)

SP-10 reference is Okada/snRNA-dominated (51.6 % / 68.6 %) with study and assay as covariates; SP-11 pseudo-replication (10/21 sections
from one autopsy) does not drive the headline pairs — leave-brainstem-out AC–astro 0.55, MES–myeloid 0.29 (p = 0.035); I² ≈ 0.99 so CIs
are between-section; SP-12 Ren statistics reproduce, cycling p = 0.008 is the floor for 5 v 5 and nothing survives BH (declared);
SP-13 MES-like-at-the-edge is robust to q05 vs mean and to the reference; SP-14 all counts, tables, checksums and h5ad contents verified;
SP-15 signature models converged; SP-19 1/(n−3) ignores spatial autocorrelation but is immaterial given I²; SP-20 deviations from the
build plan are declared, and the co-localisation hypotheses are literature-derived (no circularity with the atlas).

## 5. What could not be verified and why

- Deconvolution ELBO histories (all five runs) — files not in the repository (SP-03).
- The MP reference h5ad (`modal/inputs/c2l_reference_mp.h5ad`) — never copied back; class table and log were checked instead.
- Whether the Collot Visium patients (BT01_192, VUMC_11, VUMC_17) overlap the Collot snRNA patients in the reference (P-xxxx / PMC-DMG-x
  identifiers); IDs are not comparable, so reference–test independence for those 3 sections is unknown.
- Absolute calibration of abundances (no histology cell counts); only rank/share statements were assessed.
- The Modal GPU step itself was not re-run (brief: no GPU/Modal); its inputs (`modal/inputs/*.h5ad`, `c2l_signatures*.csv`) and outputs
  were checked for internal consistency instead.
