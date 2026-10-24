# Independent audit — Phase 3 (malignant vs normal calling, inferCNV) — DMG atlas

Auditor: independent adversarial review, 2026-10-19. Scope as briefed: `scripts/05_cnv_calls.py`, `scripts/06_cnv_aggregate.py`, `src/dmg_atlas/cnv.py`, `src/dmg_atlas/states.py`, `data/processed/malignant_calls.csv`, `results/tables/{malignant_calls_per_sample, malignant_call_validation, cnv_event_frequencies, jessa_author_celltype_vs_consensus_coarse, liu_author_vs_consensus_crosstab}.csv`, `results/tables/phase3_summary.json`, `results/qc/cnv/*` (6 heatmaps + 1 scatter), `docs/phase3_cnv_notes.md`. Nothing outside `docs/audit/` was modified. Findings table: `docs/audit/cnv_findings.csv` (23 rows: 1 critical, 7 major, 7 minor, 8 notes).

## 1. Files examined

| file | how |
|---|---|
| docs/phase3_cnv_notes.md | read in full; every number cross-checked against tables |
| scripts/05_cnv_calls.py, scripts/06_cnv_aggregate.py | read line by line; 05 re-executed for two studies |
| src/dmg_atlas/cnv.py | read line by line; compared with infercnvpy 0.6.1 `tl/_infercnv.py` (installed source) |
| src/dmg_atlas/states.py | read; not used by Phase 3 code |
| data/processed/malignant_calls.csv (342,213 × 52) | loaded; all tables re-derived from it |
| results/tables/malignant_calls_per_sample.csv (118 × 59), malignant_call_validation.csv, cnv_event_frequencies.csv, liu_author_vs_consensus_crosstab.csv, jessa_author_celltype_vs_consensus_coarse.csv, phase3_summary.json | recomputed |
| results/qc/cnv/*.png (7) | viewed 4 heatmaps + scatter; 2 heatmaps (jessa_multiome P-1709, okada 1277_pons) only listed, same code path |
| results/qc/qc_per_sample.csv, qc_per_study.csv | count reconciliation |
| results/qc/cnv_cluster_stats_vs_author.csv | threshold-provenance check |
| data/interim/cnv/<study>/<sample>.cnv.npz (117) and .clusters.csv (117), <study>.calls.csv, <study>.samples.csv | window/arm coverage, cluster evidence, sensitivity |
| data/external/gencode/{gene_positions_gencode_v32.csv, arm_boundaries_grch38_validation.csv, README_gencode.txt} | arm/centromere checks |
| logs/phase3_cnv.log (1,935 lines) | tuning history |
| data/interim/{jessa_scrna,filbin}.qc.h5ad (obs + row subsets of three samples) | window reconstruction, PDGFRA check |
| data/processed/sample_manifest.csv; scripts/07_assemble_atlas_input.py (lines 35-100) | leakage and downstream use of `ambiguous` |
| env/environment.md, local_pip_freeze.txt | version record |

## 2. Method of verification

* All tables were recomputed from `malignant_calls.csv` with independent code (agreement, precision, recall, event fractions, crosstabs, per-study counts).
* `scripts/05_cnv_calls.py` was re-executed unchanged (via `runpy`, `INTERIM` redirected to a temporary directory with symlinked inputs) for `desisto` and `filbin` and compared cell-by-cell with the stored outputs.
* The inferCNV arithmetic was checked against the installed infercnvpy 0.6.1 source; the window→arm mapping was reconstructed for three samples from the QC h5ad files.
* The tuning history was reconstructed from `logs/phase3_cnv.log` (which samples were in the "benchmark subset" runs) and tier-2 agreement was split into tuning vs held-out samples.
* Centromere coordinates were compared with UCSC hg38 cytoBand `acen` rows for chr1, chr2, chr17 (found in UCSC mailing-list / GitHub copies of the table via web search; hgdownload itself is not reachable from this sandbox) and with the auditor's own knowledge for the remainder.
* Figures were viewed.

## 3. Re-derived numbers (claimed vs recomputed)

| quantity | claimed (notes/table) | recomputed | status |
|---|---|---|---|
| cells in malignant_calls.csv / samples | 342,213 / 118 | 342,213 / 118 | ✓ |
| per-study rows vs QC `cells_after_qc` (non-Jessa) | equal | equal for all 8 studies | ✓ |
| jessa_scrna working set | 84,194 | 83,739 author-QC ∩ ours + 455 relaxed-mito = 84,194 | ✓ |
| jessa_multiome working set | 46,618 | 44,373 + 708 + 1,537 + 0 = 46,618 | ✓ |
| consensus counts per study (10 studies) | as in notes | identical | ✓ |
| tier-2 agreement jessa_scrna / multiome / collot / okada / liu / labelle | 0.913 / 0.923 / 0.905 / 0.854 / 0.808 / 0.923 | 0.9129 / 0.9234 / 0.9046 / 0.8539 / 0.8079 / 0.9232 | ✓ |
| tier-2 precision / recall jessa_scrna | 0.974 / 0.919 | 0.9744 / 0.9194 | ✓ |
| tier-3 author-independent n, agreement (jessa_scrna) | 43,266, 0.965 | 43,266, 0.9648 | ✓ |
| tier-3 consensus vs author | 1.000 | 0 mismatches / 228,643 (by construction) | ✓ (tautological) |
| samples with tier-2 agreement < 0.90 | 49 / 97 | 49 / 97 (9/17, 3/10, 17/22, 2/5, 18/42, 0/1) | ✓ |
| P-6254 / P-6255 / P-4198 / P-1779 / 1732_pons / 874_pons / 1934_thal / MUV1 / MUV5 / MUV82 | 0.47 / 0.69 / 0.71 / 0.72 / 0.25 / 0.50 / 0.62 / 0.03 / 0.16 / 0.28 | 0.465 / 0.694 / 0.712 / 0.719 / 0.250 / 0.501 / 0.615 / 0.034 / 0.162 / 0.277 | ✓ |
| P-6337 marker tuning: Immune / Oligo / Vascular normal | 2,626/2,963, 52/55, 126/167 | identical | ✓ |
| Liu crosstab | as table | identical | ✓ |
| chr17p loss, okada | 0.30 | 0.030 (0.022 among 83,980 evaluable cells) | ✗ bug |
| chr17p loss, collot | 0.74 | 0.000 | ✗ bug |
| chr17p loss, jessa_scrna | 0.28 | 0.177 | ✗ bug |
| chr17p loss, labelle | — (table 0.803) | 0.707 | ✗ bug |
| all other hallmark fractions (6 events × 6 studies) | as notes | identical to 2 dp | ✓ |
| Filbin 4q12 gain "1 of 9 patients dominates" | — | all 9 samples 0.97–1.00 | ✗ |
| tuning set | 7 Jessa scRNA samples; Okada/Collot not used | 9 jessa_scrna + 7 okada + 4 jessa_multiome + 2 collot samples in benchmark runs; cluster-stats file includes okada/liu/collot | ✗ |
| benchmark recall / precision | 0.986 / 0.973 | 0.913 / 0.975 on the 9 benchmark samples (final calls); best 7-of-9 subset 0.940 / 0.978 | ✗ unverifiable |
| score_floor rationale percentiles | 0.017 / 0.010 | 0.0134 / 0.0114 (all studies), 0.0039 / 0.0110 (Jessa) | ✗ |
| centromeres chr1, chr2, chr17 | typed from memory | match UCSC acen rows | ✓ |
| 05_cnv_calls re-run desisto, filbin | — | bit-identical scores and calls | ✓ |

## 4. Findings (ordered by severity)

### Critical

**CNV-01 — chr17p_loss frequencies are wrong in four studies (NaN→True bug).** `06_cnv_aggregate.py` uses `mal[e].astype(bool).mean()`. The `chr17p_loss` column is NaN for 63,694 cells in 8 samples where no CNV window maps to 17p (okada 586_pons, 1734_pons, 1656_pons; labelle P9; jessa_scrna P-3407, P-4198; collot PMC-DMG-7, PMC-DMG-8), and `NaN.astype(bool)` is `True`. The table/notes report okada 0.30 (true 0.03), collot 0.74 (true 0.00), jessa_scrna 0.28 (true 0.18), labelle 0.80 (true 0.71); the per-sample table shows `frac_malignant_chr17p_loss = 1.0` for all eight affected samples. The other six events are unaffected (NaN only for scpca, which has no malignant cells).

### Major

**CNV-02 — Arm events are computed from centromere-spanning windows.** Smoothing is per chromosome (infercnvpy) with 100-gene windows, and each window is assigned to an arm by its mean position. In a well-covered sample (P-6337, 7,021 genes) chr17 yields 30 windows of which 5 are labelled 17p, containing 82→42 % p-arm genes; chr10 has 2 "10p" windows. 18p is absent in every sample; in low-complexity samples most p arms vanish (okada 1734_pons: 17 arms missing from 22 windows; labelle P9: 13 arms missing). 17p loss — the TP53-locus hallmark — is therefore diluted everywhere and unavailable in 8 samples. Large arms (1q, 13q, 14q, whole chr2/chr7) are little affected.

**CNV-03 — "4q12 gain" is a single 44-Mb window driven by OPC-lineage expression.** 0–2 windows per sample; in filbin BCH1126 the one window spans chr4:37.6–81.7 Mb with 19/100 genes in the band. PDGFRA log-expression is 0.00 in the (immune) reference vs 1.29 in tumour cells; with the reference SD floored at 0.01, one such window gives z ≥ 2. Result: 4q12 gain in 97–100 % of malignant cells in all nine Filbin patients (the notes' "1 of 9 patients dominates" is false), 59 % LaBelle, 34 % Liu — against ~20–30 % PDGFRA amplification/activation in published DMG cohorts. The event also feeds malignant calls.

**CNV-04 — The "author-independent" tier-3 agreement and the Liu crosstab over-state validation.** Tier-3 agreement is computed only on cells where marker and CNV calls agree (51 % of Jessa/Okada cells). The excluded cells are the hard ones: of 42,941 labelled cells with cnv=malignant & marker=normal, 79 % are author-malignant. The Liu consensus crosstab cannot show disagreement because `consensus_call` requires author agreement (0 of 228,643 labelled non-ambiguous cells disagree, by construction); the informative tier-2 crosstab shows e.g. OC-like 380 malignant / 129 normal, Myeloid 114 / 576.

**CNV-05 — The `ambiguous` class is state-biased and propagates downstream.** 100,880 cells (29.5 %) are ambiguous; the dominant cause is cnv=malignant with marker=normal (okada 21,980; jessa_scrna 10,526), mostly oligodendrocyte- or neuron-scoring tumour cells (Jessa: 5,503/5,594 and 3,152/3,249 author-malignant). Jessa author "Oligodendrocytes" are 49 % ambiguous, "Neurons" 45 %, OPC 21 %; Liu OC-like 293/509. `07_assemble_atlas_input.py` maps ambiguous → `Unknown` and admits only consensus-malignant cells to the malignant compartment, so Phases 4–6 under-represent OC-like/neuron-like tumour states.

**CNV-06 — Tuning history is misreported.** The log shows benchmark-subset runs on 9 jessa_scrna samples (not 7), plus 7 Okada, 4 jessa_multiome and 2 Collot samples; `results/qc/cnv_cluster_stats_vs_author.csv`, cited in `cnv.py` as the basis of `score_floor`, contains author-labelled clusters from Okada, Liu and Collot. The claim "Okada and Collot labels were not used for tuning" is therefore not supported. The log records no recall/precision, so the design-table numbers (e.g. "recall 0.986 / precision 0.973") cannot be verified; the delivered pipeline achieves 0.913 / 0.975 on the nine benchmark samples. The marker-threshold tuning sample is P-6337 in the notes but P-6251 in `cnv.py`.

**CNV-08 — Pseudo-replication.** Agreement is reported as pooled cell fractions; the sample-level median is 0.899 (IQR 0.835–0.963), the unweighted sample mean 0.837 (95 % bootstrap CI 0.793–0.876) versus the cell-weighted 0.882; Okada alone is 47 % of labelled cells.

**CNV-09 — One fifth of CNV-malignant calls rest on the absolute score floor alone.** 49,499/237,752 labelled CNV-malignant cells (20.8 %) have no hallmark event; their author-malignant fraction is 0.847 vs 0.954 (collot 0.47, jessa_multiome 0.48, okada 0.83). Raising the floor from 0.008 to 0.016 flips 12.4 % of CNV-malignant cells; lowering it to 0.004 flips 8,716 normal cells. okada 1734_pons (1,160 genes, 22 windows) is almost entirely score-only.

### Minor

CNV-10 (z-scores are effectively absolute deviations because the noise filter collapses reference SD to the 0.01 floor); CNV-11 (best-match baseline and post-hoc MAD refinement shrink signal; conservative); CNV-17 (LaBelle gate is one low-complexity sample with a 7.7 % false-positive rate, shown as a pass); CNV-18/19 (heatmaps lack colour bar and call legend, chromosome labels overlap; scatter legend n ≠ plotted n, overplotting order); CNV-20 (score_floor percentile rationale does not reproduce from the file; file is kernel-only); CNV-23 (hallmark fractions are cell-pooled and not comparable with patient-level literature; per-sample columns exist but are unreported).

### Confirmed sound (notes)

CNV-07: tier-2 agreement on held-out vs tuning samples is similar or higher (jessa_scrna 0.916 vs 0.912; multiome 0.946 vs 0.888; okada 0.850 vs 0.860; collot 0.946 vs 0.852) — no evidence that tuning inflated Jessa agreement. CNV-12: centromere coordinates are correct. CNV-13: the inferCNV arithmetic matches infercnvpy 0.6.1 (ordering, windows, clipping, centring, noise filter; no off-by-one). CNV-14: `05_cnv_calls.py` is bit-reproducible from the CLI (desisto, filbin). CNV-15: all cell counts reconcile with the QC tables. CNV-16: validation tables and crosstabs reproduce exactly. CNV-21: no author-label leakage into references; same-patient pooled reference cells occur for two patients only (benign). CNV-22: `states.py` is unused in Phase 3.

## 5. Overall verdict for Phase 3

**Sound with corrections.** The per-cell calling pipeline is deterministic, internally consistent, correctly reconciled with the QC files, and its tier-2 calls agree with author labels at ~0.85–0.92 per study with no evidence of tuning inflation. However: (i) the hallmark-event table has a wrong column (17p; critical) and a mis-specified event (4q12; major) that should not be cited until fixed; (ii) short-arm events need per-arm smoothing; (iii) the tier-3 validation narrative and the Liu crosstab over-claim and must be rewritten; (iv) the state-dependent `ambiguous` class should be declared as a limitation of every downstream malignant-compartment analysis; (v) the tuning description must be corrected. None of these changes the malignant/normal call of most cells, but the hallmark-frequency deliverable and parts of the notes are currently not defensible as written.

## 6. Not verified and why

* The design-table benchmark numbers in the notes (recall/precision per intermediate design) — not in the log, intermediate outputs overwritten.
* Full re-run of the four large studies (jessa_scrna, jessa_multiome, okada, liu) — not attempted (verification brief: targeted recomputation); desisto and filbin were re-run bit-exactly and all large-study tables were recomputed from the stored per-cell file.
* Centromere coordinates for chromosomes other than 1, 2, 17 against the UCSC file itself (hgdownload blocked); checked against gene-gap validation and auditor knowledge only.
* Heatmaps for jessa_multiome P-1709 and okada 1277_pons were not viewed individually (same code path as the four viewed).
* Whether author labels for Okada/Collot were *looked at* during the benchmark iterations beyond what the cluster-stats file shows — not inferable from the log.
