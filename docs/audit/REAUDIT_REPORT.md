# Independent re-audit of audit closure — H3K27M DMG single-cell atlas

**Date:** 2026-10-21 · **Auditor:** independent re-audit, fresh context, no part in the remediation
**Scope checked:** 90 of the 163 findings — all 37 critical/major findings plus every finding that any
`fixes_*.csv` marks `closed` (89 exact-`closed` rows; CNV-05 is in scope as a major finding whose most
advanced status is `partially_closed`).
**Method:** every numeric claim was recomputed from the file on disk, not read from the closure record's
`new_value` field; every method claim was traced into the executed code path; every wording claim was read
in the document together with a search for surviving contradictions. Log: `logs/reaudit.log`.
Row-level results: `docs/audit/reaudit_closure.csv`.

## 1. Verified-closure rate by severity

| severity | in scope | verified closed | partially closed | regressed | unverifiable | rate |
|---|---|---|---|---|---|---|
| critical | 4 | 4 | 0 | 0 | 0 | 1.000 |
| major | 33 | 28 | 4 | 1 | 0 | 0.848 |
| minor | 41 | 34 | 0 | 1 | 6 | 0.829 |
| note | 12 | 12 | 0 | 0 | 0 | 1.000 |
| **total** | **90** | **78** | **4** | **2** | **6** | **0.867** |

All four critical findings (CNV-01, IM-01, MC-01, P7-01) are independently verified closed.

## 2. Mandated re-derivations

Each of the eight re-derivations named in the brief was computed from disk.

| quantity | recomputed | recorded | verdict |
|---|---|---|---|
| Reference cohort | 311,618 cells, 57 samples, 52 patients, 18,882 genes (`atlas_core.h5ad`) | identical in README.md | match |
| scANVI held-out-patient accuracy | 0.9751 on 23,827 labelled cells of 5 held-out patients (okada_1958/2054/1890/1656, PMC-DMG-6); train 0.9931 on 197,867; balanced 0.9166; 0 patient leakage | 0.9751 / 0.9931 / 0.9166 | match |
| Phase 5 zero-shot projection | `epochs_run` 0 and `trainable_params` 0 in all five `models/scanvi_query_*/run.json`; no query `model.pt` exists; the reference classifier rebuilt from `models/scanvi/model.pt` reproduces the saved labels at 0.99967–1.000 (residual = float16 latent rounding) | zero-shot, encoder frozen | match |
| CNV sample-level agreement | 97 labelled samples: median 0.898, IQR 0.793–0.978, unweighted mean 0.836 (bootstrap 95% CI 0.790–0.877, 2,000 resamples, seed 0), pooled 0.884, study-equal-weight 0.9010, min 0.034 | identical | match |
| CNV hallmark feature set | 7 events (chr1q/chr2/chr7 gain; chr13q/chr14q/chr10q/chr17p loss); `chr4q12_gain` absent from both `malignant_calls.csv` and `cnv.py HALLMARK_EVENTS`; `arm_4q12`, `gene_dev_PDGFRA`, `gene_dev_KIT` retained as exploratory | 4q12 dropped | match |
| Metaprogram count and replication | 15 metaprograms at the 0.85 primary cut; 2 replicated under the pre-registered rule (≥3 publications AND both assays), re-derived independently of the stored `replicated` column and identical to it: MP1 cycling, MP4 hypoxia. Cut curve at 0.01 steps: 1 at 0.75–0.83, 2 at 0.84–0.90, 3 at 0.91, 5 at 0.92 | identical | match |
| Spatial co-localisation (6 headline pairs) | Independent reimplementation (spot proportions, per-section Spearman, Fisher-z averaged within patient, DerSimonian–Laird across 12 patients, BH over 42 DMG state×TME pairs): AC-like×astrocyte 0.393 (0.234–0.531) q 1.0e-4; AC-like×oligodendrocyte 0.308 (0.138–0.459) q 0.005; MES-like×vascular 0.253 (0.078–0.412) q 0.034; MES-like×myeloid 0.238 (0.071–0.392) q 0.034; OC-like×oligodendrocyte 0.301 (0.036–0.526) q 0.087; OPC-like×oligodendrocyte 0.178 (−0.028–0.370) q 0.210 | identical to three decimals including confidence limits and sections-positive | match |
| cell2location realised epochs | huchede 20,000 · collot_st 20,000 · brainstem 20,000 · ren 12,940 (`results/tables/c2l_histories/`) | identical | match |
| Ren DMG vs GBM | Under the table's own section summary (sum of q05 per section): tumour fraction 0.5284 vs 0.4605, p 0.0556, q 0.6667; cycling share 0.1118 vs 0.2143, p 0.0952, q 0.6667; all 14 variables reproduce | identical | match |

Two of these carry an attached defect that does not change the quantity itself; see §3 (SP-03 plateau
statistic) and §4 (composition q-values).

## 3. Items not closed, partially closed or regressed

**MC-06 (major) — REGRESSED.** The remediation added the site-adjusted and pons-restricted H3 models and
they are present in `results/tables/composition_results.csv`; the substantive conclusion survives (neuron
is the only effect that holds: q_site_adj 0.0169, q_pons 0.0495, both < 0.05; astrocyte, myeloid, vascular
and malignant fraction fall to q_site_adj 0.167–0.269 and q_pons 0.279–0.287). But `docs/phase6_notes.md`
L196–198 still states "site-adjusted q 0.011; pons-only q 0.032" and "q 0.16–0.24 … or the pons restriction
(q 0.23–0.26)". None of those four figures matches the table on disk. The section carries no `[superseded]`
marker and `docs/audit/reporting_consistency.csv` has no row for it (its only H3 row, index 60, covers
`docs/figures.md` Figure 4c and records the correct current values). The composition table was regenerated
on 2026-10-19 22:10, after the numbers in the notes were written.

**IM-09 (minor) — REGRESSED.** `docs/phase5_mapping_notes.md` quotes the reference self 15-NN distance as
1.25 (p90 1.75) at L34, L86 and L191–192, and the Braun distance as 2.25 at L191. Recomputed from the
current files: `phase5_reference_knn_distances.csv` median 1.566, p90 2.145; `query_mapping.h5ad`
`knn_mean_dist` for Braun median 2.604. The phase5 track corrected these values against the wave-1 latent;
the reflabels track then re-encoded every query zero-shot against the retrained reference (recording
1.57/2.15 in its own closure row, `fixes_reflabels.csv` IM-21) without updating the notes. The reporting
pass did not cover these three statements.

**CNV-04 (major) — PARTIALLY CLOSED.** Tier-3 coverage reproduces exactly (0.520 overall; jessa_scrna 0.514,
multiome 0.505, collot 0.577, okada 0.505, liu 0.554) and so does the excluded-cell count (42,314). The
author composition does not: recomputed 33,934 author-malignant (80.2%), against 34,596 (82%) stated both in
`fixes_cnv.csv` and in `docs/phase3_cnv_notes.md` L53. L60 of the same document states 33,934 for the
45,305-cell superset, so the document contradicts itself, and the recomputation agrees with L60.

**CNV-05 (major) — PARTIALLY CLOSED (as the record itself declares).** Every subset count reproduces:
ambiguous 99,263 (29.0%), cnv-malignant & marker-normal 45,305 (33,934 author-malignant), of which 30,463
become malignant under `consensus_call_eventtrust`, with jessa_scrna 47,520→55,163 and okada 82,440→94,674.
The record's "eventtrust malignant compartment +34,471 cells" is wrong: the compartment gain is +30,463
(no cell moves in the opposite direction), which is the same 30,463 quoted earlier in the same row. The
declared residual — that Phases 4–6 still use `consensus_call` — stands.

**DQ-01 (major) — PARTIALLY CLOSED.** The LaBelle half is complete: `labelle.qc.h5ad` carries
`layers['log2_tpm10']` and `uns['value_scale']`, and every downstream reader (`qc.py EXPRESSION_LAYERS`,
`genes.py`, `05_cnv_calls.py`, `10_prepare_queries.py`) was updated. The Filbin half is code-only:
`ingest.py iter_filbin` sets `uns['value_scale']` at L305, but `data/interim/filbin.qc.h5ad` has mtime
2026-10-17 15:24 (pre-audit), carries `layers['tpm']` and has no `uns['value_scale']`. The closure record
asserts "uns read back with h5py"; for Filbin that read cannot have been against this file.

**SP-03 (major) — PARTIALLY CLOSED.** The epoch counts are exact and the histories are genuinely harvested.
The convergence statistic is not reproducible: the last 5% of logged epochs contribute 0.0686% (huchede),
0.0576% (collot_st), 0.0124% (brainstem) and 0.0235% (ren) of the total ELBO reduction — range 0.012–0.069%,
not "0.020% (huchede) to 0.022% (Ren)" as stated in `fixes_spatial.csv` and in `docs/phase6_notes.md` L355.
The value is also misattributed: huchede is the largest of the four, not the smallest. The conclusion that
all four traces have plateaued is unaffected and remains conservative (every cohort is below 0.07%).

## 4. New defects introduced by the remediation

The five items above are all of this kind: the remediation ran in sequential waves, and three documents
retain numbers that were correct when written and were invalidated by a later wave regenerating the
underlying table (MC-06, IM-09) or by a transcription error in the closure record itself (CNV-04, CNV-05,
SP-03). No regenerated table was found to disagree with another table, no figure was found to plot values
that disagree with its source (`phase6_brainstem_gradient.png` peak 0.799 at Pons_C and minimum 0.404 at
Spinal_Cord_A match `brainstem_section_means.csv` exactly; Figure 3c, 4c, 5d caption values recompute to
`metaprogram_cell_vs_nucleus.csv` 0.075→0.125, `composition_results.csv` 396 fitted / 278 estimable / 275
with q, and I² 0.952–0.997 over the 42 pairs), and no `.pre_*` backup was found to contradict what the
notes claim was changed. All 30 input md5s in `docs/figures.md` and all 21 in `docs/supp_figures.md` match
the files on disk. The release checksum `885e8dfd070c01f7c492aba8c1b77ca3` and the README checksums for
`atlas_core.h5ad` (`f14dd67c81ec90395c447c358f87c832`) and `query_mapping.h5ad`
(`9b585728cc8dcb81c0566362ced87eec`) all verify.

One closure record uses an undisclosed denominator. CNV-02 states "all 39 arms evaluable in 116 of 117
samples" and "17p evaluable in all 117 samples". Both are true of the 117 samples where inferCNV actually
ran, but `data/interim/cnv/*.samples.csv` holds 118 rows: the 118th is `scpca__SCPCS000637`
(441 cells, `status = 'no_reference'`, `n_genes`/`n_windows` NaN, 0 of 39 arms evaluable, all 441 cells
ambiguous, 0 CNV-malignant) for which no CNV call was ever produced. The per-arm smoothing fix is complete
for every processed sample, so CNV-02 is verified closed, but the 117 denominator silently excludes that
sample and the record does not say so. This is the same 441-cell sample the original audit identified as
NaN for every event.

Three smaller bookkeeping errors in closure records, each immaterial to a conclusion: MC-03 states 9
`published_or_curated` label provenances where `metaprograms.csv` holds 10 (the novel MP carries provenance
`none`); MC-07 gives MP13's IEG count as 31/50 where the programme has 45 genes; P7-05 states 199 obs
columns where the release has 200, and "44 ontology terms" against 17 distinct `cell_type_ontology_term_id`
values and 46 across all seven ontology columns.

One latent defect in code that does not affect any released number: `scripts/26_foundation_embeddings.py`
L76 still filters on the renamed family `novel_MP15`. Step 26 was superseded by
`26b_foundation_comparison_rebuild.py` for the audit fix, so the filter is dead, but re-running 26 as the
documented pipeline step would no longer exclude the family it was written to exclude.

## 5. What could not be verified

- **Six figure-legibility fixes were not visually re-inspected** (CNV-18 six CNV heatmaps, CNV-19 the
  cnv_score scatter, DQ-13 ten QC violins, IM-13 the mapping-uncertainty panel, IM-15 the Phase 4 UMAP
  overview, MC-19 the composition-effects heatmap). Each file exists with a `.pre_audit` copy and a later
  mtime, so regeneration happened; whether the specific legibility complaint is resolved is recorded as
  unverifiable rather than closed. Four figures were inspected
  (`phase6_brainstem_gradient.png`, `phase6_spatial_abundance_maps.png`, `phase6_tme_marker_scores.png`
  and, by value, the gradient source table) and all four are fixed.
- **Claims that are provenance statements rather than computations** (CNV-12 centromere coordinates,
  CNV-13 the infercnvpy attribution, DQ-12 the author pct_mito comparison, P7-10 the ATAC peak provenance)
  were confirmed present in the documents but cannot be recomputed from this repository; DQ-12 was not
  independently re-derived by the remediation either and says so.
- **The pre-audit state** was not reconstructed. Where a closure record quotes an `old_value`, this
  re-audit checked only that the current value is what the record claims, not that the old value was what
  it claims.
- **The 51 findings with no closure row in any `fixes_*.csv`** were outside this scope and are not
  assessed here: 28 are note-severity entries whose recommendation is literally "None.", but 9 minor
  (P7-15, P7-18, P7-19, SP-07, SP-08, SP-09, SP-16, SP-17, SP-18) and 14 note findings carry a real
  recommendation and no declared disposition. Several of these were in fact addressed — the
  `15b_spatial_tests_fix.py` docstring implements SP-07 (a <500-UMI sensitivity run, written to
  `spatial_colocalisation_sensitivity.csv`) and SP-08 (BH on the brainstem metrics, present as
  `q_bh_abs`/`q_bh_prop`/`q_bh_re` in `brainstem_gradient.csv`) — so this is an incomplete closure record
  rather than missing work, but it is not possible to tell from the records which of the 23 are done.

## 6. Residual limitations a reader of this atlas must know

1. **29% of the reference is not consensus-malignant and not normal.** 99,263 cells (29.0%) are
   `ambiguous`, dominated by 45,305 cells that CNV calls malignant and markers call normal, of which
   33,934 carry an author-malignant label. Phases 4–6 use `consensus_call`, so these cells sit in
   `lineage_label` "Unknown"; the release moves 69,322 of them to CL:0000000 "cell". The
   `consensus_call_eventtrust` sensitivity label (which would add 30,463 cells to the malignant
   compartment) exists but no downstream phase uses it. This is CNV-05's declared open item and it is the
   largest single interpretive caveat in the atlas.
2. **18.9% of CNV-malignant cells are called on score elevation with no hallmark event** (47,337 of
   250,221), and their author-malignant agreement is 0.862 against 0.952 for event-backed cells. Two
   samples are flagged `low_confidence_sample` (okada__1734_pons, 1,163 genes; labelle__H3K27M_P9).
3. **Cross-assay replication of malignant metaprograms is limited to cycling and hypoxia** and is
   cut-sensitive: 1 replicated programme at cuts 0.75–0.83, 2 at 0.84–0.90, 3 at 0.91, 5 at 0.92. The
   build-plan assay-dependence checkpoint is triggered and remains so.
4. **Metaprogram family assignment is standardised within sample**, which removes between-sample
   differences in family abundance by construction; per-sample composition from that column is not
   interpretable as abundance.
5. **Query projections are zero-shot** — the reference encoder and classifier are applied unchanged, with
   no query adaptation. Braun fetal hindbrain is 87.7% predicted malignant, which is a statement about the
   classifier's behaviour off its training manifold, not about the fetal tissue.
6. **The spatial layer rests on 12 patients**, 10 of the 21 DMG sections being serial sections of one
   autopsy (66% of DMG spots). Of the two pre-stated hypotheses, MES-like with myeloid survives
   (q 0.034) and OPC-like with the oligodendroglial lineage does not (q 0.210). I² is 0.952–0.997 across
   all 42 pairs: between-patient heterogeneity dominates. The Ren DMG-vs-GBM comparison has n=5 vs 5 and
   nothing survives BH correction.
7. **The Ren deconvolution stopped at 12,940 of 20,000 epochs** when its wall-clock budget expired.
8. **Phase 7 is exploratory and outside the approved build plan** (so labelled in `docs/release_v1.md`).
   Within it: the ATAC differential-accessibility counts are an uncorrected sign test over 5–9 samples;
   14 of 35 label-correlation pairs are masked at n<3; one of ten multiome samples is at chance for ATAC
   separability; the survival layer's NNLS fit was never refitted (P7-07 open) and
   `phase7_survival_forest.png` is a stale pre-audit figure, flagged as superseded in
   `docs/audit/figure_checklist.csv` but still present on disk.
9. **Two figures are declared stale and were not regenerated**: `phase6_brainstem_gradient_mp.png` and
   `phase7_survival_forest.png`. A third, `figS8_bulk_survival.png`, is declared a partial match. Three
   step-21 tables (`brainstem_gradient_mp.csv`, `brainstem_section_means_mp.csv`,
   `brainstem_mp_vs_fallback.csv`) and `c2l_signatures_mp.csv` still carry the pre-audit `novel_MP15`
   family name and are marked superseded; no figure reads them.
10. **`phase6_tme_umap.png` and `phase6_tme_umap_clean.png` still carry in-panel titles** (MC-20, declared
    open, regeneration blocked).
11. **One sample has no CNV call at all.** `scpca__SCPCS000637` (441 cells) has `status = 'no_reference'`
    in the per-sample CNV tables: inferCNV never ran for it, no arm is evaluable, and all 441 cells are
    `consensus_call` ambiguous. Per-sample CNV statistics quoted as "of 117 samples" throughout the notes
    and closure records are over the 117 processed samples, not the 118 rows of the table.
12. **The Filbin 2018 reference gene sets remain untraceable by script** (science.org and PMC both block
    automated retrieval); the sets used for labelling are flagged `filbin_untraceable` and the gap is
    recorded (MC-03 / DQ-17).

## 7. Release decision

The atlas is **releasable once the five text defects in §3 are corrected.** Nothing in the data, the models
or the result tables was found to be wrong: all four critical findings verify, every quantity the brief
required to be re-derived reproduces from disk, the release file's ontology assignment is internally
consistent (CL:0001064 holds for 170,869 cells, every one of them `consensus_call == "malignant"`, with
zero violations), and all 51 recorded figure-input checksums and all three headline file checksums match.

The defects that block release as written are documentation defects, and all five are in documents a reader
will use as the authority for a number:

1. `docs/phase6_notes.md` L196–198 — four H3 composition q-values that contradict
   `composition_results.csv` (MC-06).
2. `docs/phase5_mapping_notes.md` L34, L86, L191–192 — reference self 15-NN distance 1.25/1.75 and Braun
   2.25, superseded by 1.566/2.145 and 2.604 (IM-09).
3. `docs/phase3_cnv_notes.md` L53 — 34,596 (82%) author-malignant, contradicted by L60 of the same
   document and by recomputation (33,934, 80.2%) (CNV-04).
4. `docs/phase6_notes.md` L355 — the ELBO plateau range 0.020–0.022%, recomputed as 0.012–0.069% and
   misattributed between cohorts (SP-03).
5. `data/interim/filbin.qc.h5ad` — either regenerate it so `uns['value_scale']` is present as the closure
   record asserts, or correct the DQ-01 record to say the Filbin change is code-only (DQ-01).

None of the five changes any conclusion. Items 1–4 are arithmetic corrections to sentences; item 5 is
either a 1-minute re-ingest or a one-line record correction. The 23 findings in §5 with no declared
disposition should be given closure rows before release, because at present the closure record cannot
distinguish work that was done from work that was skipped.

---

## Response to this report (2026-10-20, lead agent)

All five release-blocking documentation defects are corrected; per-item disposition, recomputed values and the
evidence string for each are in `docs/audit/reaudit_response.csv`.

1. **MC-06** — the H3.1-vs-H3.3 paragraph in `docs/phase6_notes.md` now quotes `composition_results.csv`: neuron
   q 0.004, site-adjusted 0.017, pons-only 0.050; the non-surviving effects 0.038 / 0.038 / 0.051 / 0.068 with
   site-adjusted 0.17–0.27 and pons-only 0.28–0.29. The conclusion (neuron is the only surviving effect) is unchanged.
2. **IM-09** — all five stale kNN distances in `docs/phase5_mapping_notes.md` replaced by the post-re-encoding values
   (reference self 1.566, p90 2.145; Braun 2.604; malignant → Braun 3.279), and the fetal-class fractions and
   enrichments restated from the current `malignant_nearest_fetal_class.csv`. The wave-1 record is marked superseded.
3. **CNV-04** — this was a table defect, not only a document contradiction. `06_cnv_aggregate.py` grouped tier-3
   exclusions on `author_state` and took `author_call` from each group's first cell, which mislabelled 666
   Jessa-multiome cells. The script is patched, the table regenerated by the new `06b_tier3_composition_fix.py`
   (old file kept as `.pre_reaudit.csv`), and the notes now state 33,934 (80.2 %) — matching the per-cell calls.
4. **CNV-05** — compartment gain corrected to +30,463 (182,924 → 213,387), recomputed per cell; no cell moves the
   other way.
5. **SP-03** — the plateau statistic is now the recomputed per-cohort value (0.012 % brainstem, 0.024 % Ren,
   0.058 % collot_st, 0.069 % huchede), replacing the misattributed 0.020–0.022 % in both records.
6. **DQ-01** — `uns['value_scale']` written into `data/interim/filbin.qc.h5ad` (matrix unchanged,
   `.pre_reaudit.h5ad` kept), so the record now holds for both TPM studies.
7. **CNV-02 denominator** — disclosed in `fixes_cnv.csv`: the per-sample table has 118 rows and
   `scpca__SCPCS000637` (441 cells) has no CNV call, so "116 of 117" counts processed samples only.

Still open and **not** closed by this response: **MC-20** (the two TME UMAP figures cannot be regenerated without a
per-compartment scVI refit) and the **23 findings with no closure row**, which remain undeclared rather than
demonstrably done or skipped.
