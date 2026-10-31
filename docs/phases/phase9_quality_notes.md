# Phase 9 — atlas quality assessment follow-ups (2026-10-31)

Code: `modal/phase9_loso_scanvi.py` (GPU, via `modal/phase6_cpu_run.py --gpu A10G`), `scripts/32_loso_classifier.py`,
`scripts/33_metaprogram_lopo.py`, `scripts/34_composition_label_sensitivity.py`.
Outputs: `results/tables/quality/` (each script writes a `*_inputs.json` with input md5s), `models/loso_scanvi/<publication>/`.
Log: `logs/phase9_quality.log` (local steps), `logs/phase9_loso_modal_stdout.log` + `models/loso_scanvi/runs.jsonl` (GPU job).
No earlier-phase output was modified. Publications: jessa_scrna + jessa_multiome = jessa2022 (4 shared patients) throughout.

Why: a dataset-curation review of `dmg_atlas_v1.h5ad` flagged (1) the reference classifier had only been tested on held-out
patients of studies it was trained on, (2) the cell-type composition set labels the 29 % ambiguous-consensus cells with scANVI,
(3) one or two publications could drive the metaprograms. Recommendations (4) more whole-cell data and (5) thin strata need new
data and are not addressed here (see "Not done").

## 1. Leave-one-publication-out classifier (`32_loso_classifier.py`)

For each publication P, scVI (settings of `scvi_seed1`) and scANVI (settings of `models/scanvi`, patient-grouped inner
validation over the other publications) were refitted without P, and P was mapped by scArches surgery with its labels
hidden (`modal/phase9_loso_scanvi.py`, Modal A10G, 5 folds, ~19 min GPU in total). Scored against the corrected
`lineage_label`; balanced accuracy averages classes with ≥ 20 true cells in P.

| Held out | Labelled cells | Accuracy | Balanced acc. | Malignant recall / precision | Worst class (recall) | Patient accuracy (min / median) |
|---|---:|---:|---:|---:|---|---:|
| Okada | 106,655 | 0.893 | 0.402 | 1.000 / 0.913 | T/NK, neuron (0.00) | 0.10 / 0.92 |
| Jessa 2022 | 94,965 | 0.974 | 0.743 | 0.988 / 0.997 | T/NK (0.00) | 0.77 / 0.97 |
| Collot | 13,945 | 0.987 | 0.704 | 1.000 / 0.996 | T/NK (0.03) | 0.97 / 0.99 |
| Zang | 5,026 | 0.992 | 0.988 | 0.999 / 0.993 | myeloid (0.98) | 0.99 |
| DeSisto | 1,103 | 0.880 | 0.715 | 0.861 / 0.973 | astrocyte (0.00) | 0.71 / 0.82 |
| *pooled, all folds* | 221,694 | 0.936 | 0.506 | | | |
| *baseline: 5 held-out patients* (`models/scanvi`) | 23,827 | 0.975 | 0.917 | | | |

- **Malignant and myeloid transfer to an unseen study** (malignant recall 0.86–1.00, myeloid 0.98–1.00 in every fold).
- **Rare normal lineages do not.** T/NK recall is 0.00 (Okada, Jessa 2022) and 0.03 (Collot), mostly predicted myeloid;
  only DeSisto's 96 T/NK cells are recovered (1.00). Astrocyte recall is 0.00–0.56 across folds.
- **The failure mode is over-calling malignant.** With Okada unseen, 30 % of its labelled normal cells are predicted
  malignant: 3,838 / 3,898 neurons (98 %), 1,943 / 3,065 astrocytes (63 %), 1,328 / 2,592 vascular (51 %). Okada is the
  only autopsy cohort and holds most of the reference's neurons, so without it the model has no nucleus-level normal brain
  context — the same pattern as the zero-shot normal queries (87.7 % of Braun hindbrain, 29.0 % of Siletti OPCs called
  malignant, `docs/release/limitations.md`).
- Confidence does not flag the errors: at scANVI maxprob ≥ 0.9 (coverage 0.95–0.99) accuracy is 0.91–0.99, Okada 0.912.
- Within-reference inner validation accuracy was 0.957–0.984 in every fold, so the drop is cross-study shift, not underfitting.

Verdict: the 0.975 / 0.917 held-out-patient figure overstates performance on a new study. For a new cohort the reference
classifier is reliable for malignant vs myeloid; its normal-lineage calls (and therefore its malignant calls on normal-rich,
autopsy or nucleus samples) need independent support (CNV, markers). This affects query mapping (Phase 5) and any future
cohort, not the reference's own labels, which use consensus calls.

## 2. Composition vs labels of ambiguous cells (`34_composition_label_sensitivity.py`)

The `celltype` set of `18_composition.py` was rebuilt with four labelings and refitted with its primary model unchanged
(reproduction of `composition_per_sample.csv` and `composition_results.csv` under the scANVI labeling is asserted).

- Per-sample malignant fraction is stable to the labeling: Spearman vs scANVI 0.979 (consensus only), 0.992 (scANVI
  maxprob ≥ 0.9); 0.764 when 'ambiguous' is kept as its own category (it then absorbs a sample-dependent share).
- Of 40 estimable category × contrast tests, the 4 with q < 0.05 under scANVI labels (astrocyte, myeloid and neuron higher
  in H3.1; neuron higher in nuclei) keep their sign and q < 0.05 under consensus-only labels. Consensus-only adds 3
  (malignant lower in H3.1, q 0.024; vascular higher in H3.1, q 0.038; oligodendrocyte higher in nuclei, q 0.022), all with the
  same sign under scANVI (q 0.05–0.16). Signs flip between labelings only for tests that are far from significant under
  scANVI labels (q ≥ 0.46: T/NK and malignant vs ONC201, malignant nuclei vs cells with 'ambiguous' as a category, neuron
  thalamus vs pons).
- Caveat unchanged from Phase 6: H3.1 = 7 samples, all Jessa 2022, pontine nuclei; these are the primary-model numbers
  without the site adjustment / pons restriction that `18_composition.py` reports alongside.

Verdict: the classifier's placement of ambiguous cells does not drive any cell-type composition result.

## 3. Publication dependence of the metaprograms (`33_metaprogram_lopo.py`)

Programmes re-clustered with the unchanged `16_metaprograms.py` rule (reproduces 15 MPs / 2 replicated exactly; asserted).
Scenarios: drop each publication (5); at most 5 samples per publication (17 samples, 20 seeds); size-matched random 17-sample
pools ignoring publication (20 seeds). A primary MP is recovered if a scenario MP's gene set has Jaccard ≥ 0.5 with it.

| | MPs | replicated |
|---|---|---|
| all 50 samples | 15 | 2 (MP1, MP4) |
| drop one publication | 6–14 (median 13) | 1–2 (median 2) |
| balanced, 17 samples | 2–5 (median 4) | 1–2 (median 1) |
| random, 17 samples | 3–6 (median 4.5) | 0–2 (median 0) |

- **MP1 (cycling)** survives every drop and is recovered in 20/20 balanced and 20/20 random pools; it replicates in 20/20
  balanced pools but 5/20 random pools (random pools hold 0–3 of the 6 whole-cell samples, balanced pools always include
  DeSisto and Zang). The most robust finding of the metaprogram track.
- **MP4 (hypoxia)** is lost only when Jessa 2022 is dropped; recovered in 50 % (balanced) / 60 % (random) of 17-sample pools.
- **MP3, MP7** (OPC-like, AC-like) need Okada; **MP8 (MES-like), MP9 (novel ADGB/CFAP43/DNAH9), MP10, MP11, MP12, MP15** are
  lost when either Jessa 2022 or Okada is dropped, and are not recovered in any (or nearly any) 17-sample pool.
- Balanced and random pools give similar MP counts, so the loss of MPs in small pools is a pool-size effect, not one
  publication dominating; but most of the 15 MPs are joint Jessa 2022 + Okada structures and should be described that way.

Verdict: consistent with the pre-registered headline (only cycling and hypoxia replicate); MES-like and novel_MP9 rest on two
publications (both mostly nuclei) and need independent data before they are claimed as DMG states.

## Not done

- **More whole-cell data / thin strata (recommendations 4–5).** Needs new ingestion and a rebuilt reference (Phases 2–6),
  which would change `atlas_core`. Candidates: the 4 DeSisto/ScPCA SCPCP000001 midline samples with unrecorded H3 status
  (834, 1107, 968 thalamus; 1431 pons/cerebellum; excluded pending K27M confirmation), Sussman 2026 HTAN (dbGaP, controlled),
  LaBelle raw counts (withheld). Decision for the user.
- **Age banding in the public h5ad** (privacy): would change the release file; not applied.
