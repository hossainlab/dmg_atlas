# Audit closure summary — DMG atlas

Generated 2026-10-20 from `docs/audit/ALL_findings.csv` (163 findings) joined to every `docs/audit/fixes_*.csv`. Machine-readable form: `docs/audit/CLOSURE_SUMMARY.csv` (adds `closure_status`, `closure_source`, `blocked_reason`). Where a finding carries rows in two tracks' tables, the most advanced status wins and every source is listed in `closure_source`.

## Status definitions

| status | meaning |
|---|---|
| closed | a fixes row records the change, with old vs new values and how it was verified |
| partially_closed | part of the recommendation was implemented; the remainder is declared in the row |
| not_closed | a fixes row explicitly declares it open, with a reason |
| not_assessed | **no fixes row anywhere** — neither fixed nor declared |
| verification_note | `note`-severity row: an auditor confirmation whose `recommendation` is 'None.' — nothing to close |

## By severity

| severity | closed | partially_closed | not_closed | not_assessed | verification_note | total |
|---|---|---|---|---|---|---|
| critical | 4 | 0 | 0 | 0 | 0 | 4 |
| major | 26 | 1 | 0 | 6 | 0 | 33 |
| minor | 41 | 11 | 10 | 9 | 0 | 71 |
| note | 10 | 0 | 1 | 0 | 44 | 55 |

Actionable findings (severity other than `note`): **108** — 71 closed, 12 partially_closed, 10 not_closed, 15 not_assessed.

## By track

| track | closed | partially_closed | not_closed | not_assessed | verification_note | total |
|---|---|---|---|---|---|---|
| cnv | 16 | 1 | 2 | 0 | 4 | 23 |
| data_qc | 13 | 1 | 6 | 0 | 12 | 32 |
| integration_mapping | 21 | 3 | 0 | 0 | 0 | 24 |
| metaprograms_composition | 18 | 2 | 3 | 0 | 14 | 37 |
| phase7_release_figures | 13 | 5 | 0 | 3 | 6 | 27 |
| spatial | 0 | 0 | 0 | 12 | 8 | 20 |

All 4 critical findings are closed (CNV-01, IM-01, MC-01, P7-01).

## Not assessed — no closure row anywhere

15 findings, all actionable. These are **not** reported as fixed, and no count or absence derived from them should be read as a result.

| id | track | severity | claim | why not assessed |
|---|---|---|---|---|
| P7-15 | phase7_release_figures | minor | Zero-shot embeddings reproducible | code fix would require re-running the GPU embeddings to take effect |
| P7-18 | phase7_release_figures | minor | Dot plot of top-4 markers per cleaned TME subcluster | requires recomputing panel data / per-panel restyling; not attempted |
| P7-19 | phase7_release_figures | minor | Legibility / overlaps in the survival and foundation-model figures | requires recomputing panel data / per-panel restyling; not attempted |
| SP-01 | spatial | major | Co-localisation rho computed on raw q05 abundances; notes conclude 'bo | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-02 | spatial | major | §4: within section, higher tumour fraction goes with a larger OPC-like | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-03 | spatial | major | 'All four ELBO traces had plateaued: the last 5 % of epochs reduced th | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-04 | spatial | major | Two-panel brainstem gradient figure. | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-05 | spatial | major | Dominant-class maps of four sections (pDMG_Sample-2, DMG1, BT01_192_co | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-06 | spatial | major | Fallback malignant states = argmax of scanpy score_genes scores over 5 | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-07 | spatial | minor | All in-tissue spots are deconvolved; no spot-level UMI filter. | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-08 | spatial | minor | Section-level trends (n=10): OPC-like rho=-0.70 p=0.025, OC-like -0.66 | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-09 | spatial | minor | Hyper-parameters stated but not justified. | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-16 | spatial | minor | Provenance of Modal-sandbox steps. | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-17 | spatial | minor | Small code issues. | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |
| SP-18 | spatial | minor | Remaining figure checks. | spatial track not executed - its sub-agent was refused under content-safety filters; not retried or rerouted |

### The spatial track

Nothing was done for the spatial track: of its 20 registered findings, the **12 actionable ones (6 major, 6 minor) are unaddressed**, and the remaining 8 are auditor verification notes that need no closure. The sub-agent assigned to them was refused under content-safety filters rather than failing technically. That is a boundary: the track was not retried, reworded, split into smaller requests, or rerouted through another tool. Consequences recorded elsewhere in the repo: Figure 5 and Supplementary Figure 7 captions are marked *pending spatial corrections*, `docs/release_v1.md` states the spatial layer is pending, and the survival NNLS models (P7-07) could not be re-fitted because `c2l_signatures_mp.csv` is a spatial-track output (SP-06) still carrying the pre-audit families.

## Declared open (not_closed)

| id | track | severity | recommendation left unimplemented |
|---|---|---|---|
| CNV-11 | cnv | minor | Compute category means leave-one-out for reference cells or exclude a cell's own |
| CNV-21 | cnv | note | Exclude same-patient samples from pools for cleanliness. |
| DQ-08 | data_qc | minor | Strip the '.N' suffix before symbol mapping (or map Okada via the 2020-A feature |
| DQ-10 | data_qc | minor | Compute QC metrics after harmonisation, or on the 18,882 shared genes, so thresh |
| DQ-11 | data_qc | minor | Replace with a doublet-informed or MAD-based upper bound, or document that 0.5%  |
| DQ-16 | data_qc | minor | Read the genome attribute of the h5 / first feature IDs (ENSG vs ENSMUSG) for ev |
| DQ-18 | data_qc | minor | Use the regex ^RP[SL]\d+[A-Z]?\d*$/^RPLP\d$/^RPSA$. |
| DQ-19 | data_qc | minor | Mark as 'resection (assumed, surgical cohort)' or 'unknown'. |
| MC-12 | metaprograms_composition | minor | Document that the silhouette criterion is non-discriminative; extend the range t |
| MC-13 | metaprograms_composition | minor | Initialise version control and tag the state that produced each phase's outputs; |
| MC-22 | metaprograms_composition | minor | Report a cluster x assay independence test on held-out genes, or run without the |

## Partially closed

| id | track | severity | source |
|---|---|---|---|
| CNV-05 | cnv | major | fixes_cnv.csv:partially_closed; fixes_reflabels.csv:partially_closed |
| DQ-17 | data_qc | minor | fixes_data_repro.csv:partially_closed |
| IM-05 | integration_mapping | minor | fixes_phase5.csv:partially_closed |
| IM-16 | integration_mapping | minor | fixes_phase5.csv:partially_closed |
| IM-23 | integration_mapping | minor | fixes_phase5.csv:partially_closed |
| MC-18 | metaprograms_composition | minor | fixes_metaprograms.csv:partially_closed |
| MC-20 | metaprograms_composition | minor | fixes_metaprograms.csv:partially_closed |
| P7-07 | phase7_release_figures | minor | fixes_phase7.csv:partially_closed |
| P7-13 | phase7_release_figures | minor | fixes_phase7.csv:partially_closed |
| P7-14 | phase7_release_figures | minor | fixes_phase7.csv:partially_closed |
| P7-17 | phase7_release_figures | minor | fixes_phase7.csv:partially_closed |
| P7-20 | phase7_release_figures | minor | fixes_phase7.csv:partially_closed |
