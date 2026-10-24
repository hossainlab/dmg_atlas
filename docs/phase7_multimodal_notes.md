# Phase 7 (code) = build-plan Phase 8 — multimodal and bulk anchoring (2026-10-17/18)

> **Scope (audit P7-16).** This phase was **not part of the approved plan**: build-plan Phase 8 is listed as optional, and nothing in
> Phases 1-6 depends on it. Everything below is an **exploratory methods note**, not an atlas claim. The same statement appears in
> `README.md` and in `docs/release_v1.md` ('Scope of the analytical layers'), and the two figures drawn from this phase (Figure 6,
> Supplementary Figure 8) carry it in their captions.

Code: `scripts/24_atac_multiome.py` (Jessa multiome ATAC), `scripts/25_openpedcan_survival.py` (OpenPedCan bulk deconvolution + Cox),
`scripts/26_foundation_embeddings.py` (scGPT / Geneformer vs scVI / PCA). Log: `logs/phase7_multimodal.log` (reconstructed from the
per-run stdout because two parallel Modal runs overwrote each other's copy; every line is preserved, deduplicated, time-ordered).
Everything ran on Modal (image `im-gl8KES0WsH7PoFjf83hK6v`) through `modal/phase6_cpu_run.py`: steps 24 and 25 on an 8-CPU / 64 GB
sandbox, step 26 on an A10G. Inputs on the Modal volume `dmg-atlas-data`: `atlas_core.h5ad`, `raw/jessa_multiome/` (uploaded from
`data/raw/jessa_multiome`, 2.1 GB). OpenPedCan files are downloaded by step 25 inside the sandbox (md5 checked against the release
list; not kept locally). Extra pip packages inside the sandbox: pyreadr, lifelines (step 25); gdown, transformers==4.46.3 (step 26;
newer transformers need torch >= 2.6, the image has 2.5.1). Pretrained weights cached in `models/foundation/` (441 MB).
Deliverables of the build plan: `atac_summary.csv` ✓, `survival_results.csv` ✓.

## 1. Jessa 2022 multiome ATAC (`24_atac_multiome.py`; 3 min, peak RSS 9.1 GB)
42,004 labelled cells from the 10 multiome samples (malignant cells by metaprogram family from step 19; consensus-normal cells by
scANVI label). Peaks are called per sample (35k–175k each): 1,022,263 sample peaks merged into **279,115 union peaks** (22,182 present in
all 10 samples; 45,226 overlap a TSS ± 1 kb).
- **Label correlation** *[superseded 2026-10-20 by P7-03 — see Audit corrections below]* (within-sample pseudobulk, top 20k variable peaks; `atac_label_correlation.csv`): malignant families correlate
  0.79–0.95 with each other and 0.39–0.66 with normal lineages; MES-like ~ hypoxia 0.95, OPC-like ~ OC-like 0.95, AC-like is the most
  distinct malignant family (0.79–0.92) and the least astrocyte-like (0.39).
- **Separability** *[superseded 2026-10-20 by P7-02 — see Audit corrections below]* (kNN family agreement in an LSI space vs 200 label permutations, per sample): above chance in every sample (median
  z = 7.1; RNA scVI latent z = 52.4, inflated because the labels are RNA-derived). Same-family neighbours / expected: AC-like 2.11 (RNA 2.65),
  OPC-like 1.38 (2.35), cycling 1.36 (2.93), OC-like 1.29 (1.63), cilia 1.21 (5.0), **MES-like 1.09 (1.39), hypoxia 1.02 (1.94)**.
  Lineage-like states (AC, OPC, OC) and cycling have a chromatin identity; MES-like and hypoxia barely do, consistent with
  microenvironment-induced rather than epigenetically fixed states.
- **Promoter concordance** (promoter accessibility log2FC of the family's own MP genes vs other genes): higher for AC-like (median 0.17 vs
  −0.03), MES-like (0.16 vs 0.03), OC-like (0.18 vs 0.03), cycling (0.39 vs 0.15), hypoxia (−0.07 vs −0.24) and cilia (all p ≤ 0.001);
  **not for OPC-like** (0.006 vs 0.016, p = 0.41).
- **DA peaks** (`atac_family_da_peaks.csv`): up / down: AC-like 85 / 109, OC-like 185 / 437, cycling 51 / 264, MES-like 30 / 21, OPC-like
  6 / 190, hypoxia 1 / 386, cilia 0 / 1,321. The excess of 'down' calls (fraction of cells with an accessible peak) most likely reflects
  lower fragment depth in some families — the test is not depth-normalised — so DA counts are not interpreted further.

**Audit corrections (2026-10-20; P7-02, P7-03, P7-04, P7-10, P7-11).** `24_atac_multiome.py` was re-run on the corrected seven metaprogram
families (run ended 2026-10-19 23:14:45 +06). The two bullets marked above are pre-audit and are superseded by the numbers here; closure
rows are in `docs/audit/fixes_phase7.csv`.

- **Families (P7-04).** The ATAC layer now consumes the corrected families — AC_like, MES_like, OC_like, OPC_like, cycling, hypoxia,
  novel_MP9 — each tested in 8–10 of the 10 multiome samples (`atac_summary.csv`, column `n_samples_tested`).
- **Separability (P7-02).** Per-sample ATAC z (`atac_separability_sample_z.csv`): **9 of 10 samples above chance at z > 2, median 11.70,
  range 0.47–91.97**; `jessa_multiome__P-1694_S-1694` (580 cells) is at chance with z = 0.47. The RNA scVI latent gives median z = 92.2
  (30.2–256.7), inflated because the labels are RNA-derived. Same-family neighbour enrichment over expectation (ATAC / RNA): AC-like
  1.92 / 3.33, cycling 1.86 / 4.78, OPC-like 1.71 / 3.09, MES-like 1.65 / 2.36, OC-like 1.44 / 2.32, hypoxia 1.12 / 2.65, novel_MP9
  0.98 / 1.54. The direction of the pre-audit reading survives — lineage-like states and cycling carry a chromatin identity, hypoxia and
  novel_MP9 (cilia) essentially none — but MES-like now groups with the lineage states instead of sitting at the bottom.
- **Label correlation (P7-03).** Recomputed with per-pair sample counts (`atac_label_correlation_n.csv`, `atac_label_correlation_long.csv`;
  `interpretable` = n ≥ 3): malignant–malignant **0.89–0.96** (median 0.92; all 21 pairs n ≥ 3), malignant–normal **0.46–0.62**
  (median 0.54) over the **21 of 35 pairs with n ≥ 3** — the 14 pairs with n < 3 are masked and carry no interpretation.
  MES-like ~ hypoxia 0.93 (n 9), OPC-like ~ OC-like 0.96 (n 9). **AC-like is no longer the most distinct malignant family**: AC-like ~
  MES-like is the highest malignant pair at 0.96 (n 8), and AC-like ~ astrocyte is 0.61 (n 4), not 0.39.
- **Peak provenance (P7-10).** Peaks are **not** called by this pipeline — they are the authors' cellranger-arc per-sample peaks from GEO,
  and ATAC fragments were never downloaded. The promoter definition (TSS ± 1 kb) was verified correct.
- **DA peaks (P7-11).** `atac_family_da_peaks.csv` is **exploratory**: a binomial sign test over 5–9 samples, with no multiple-testing
  control across the 279,115 union peaks and no depth normalisation. Re-run counts (up / down): AC-like 96 / 112, MES-like 136 / 239,
  OC-like 113 / 74, OPC-like 41 / 83, cycling 80 / 167, hypoxia 2 / 197, novel_MP9 3 / 702. Not interpreted further.

## 2. OpenPedCan v15 bulk RNA-seq and survival (`25_openpedcan_survival.py`; < 1 min after download)
Cohort (`openpedcan_dmg_cohort.csv`): histologies RNA-Seq tumour rows with molecular subtype 'DMG, H3 K28' (94) or 'DMG, H3 K28, TP53'
(131) → 129 independent primary specimens (all PBTA; the plan expected ~144), 1 excluded as already in the atlas (CBTN participant of
ScPCA SCPCP000021), 100 with OS; Cox cohort **n = 99, 94 deaths, median OS 367 d**, median age 7.7 y, 74 pons *[the age and pons figures are the 129-specimen values — superseded 2026-10-20 by P7-08, see Audit corrections below]*.
Deconvolution: NNLS of TPM on the cell2location signatures (fallback states and metaprogram families), all shared genes (3,548) and a
marker-restricted version (550 / 650 genes; added after the all-gene fit looked implausible). **Both fits are unreliable**: median
proportions put T/NK at 10–22 % and give MES-like, OC-like, cycling and hypoxia zero in most patients; NNLS proportions disagree with an
independent rank-based score of the same MP genes (Spearman OPC-like −0.04, AC-like 0.27–0.44, MES-like 0.26–0.48). Likely causes: highly
collinear malignant signatures, single-cell count signatures applied to bulk TPM, no platform correction.
Survival (`survival_results.csv`; HR per SD, adjusted for age + pons, BH within each variable set):
- Rank scores (most trustworthy readout): **no association** for any family (q ≥ 0.2; cycling HR 1.24, p 0.056; AC-like 0.82, p 0.053).
- NNLS, all genes: OPC-like proportion HR 1.38 (q 0.025), MES-like share HR 0.70 (q 0.041) — not supported by the rank scores.
- NNLS, markers: OPC-like 1.31 (p 0.017, q 0.086), AC-like share 0.79 (q 0.079); MES-like degenerate (mostly zeros).
**Conclusion: no robust survival association of metaprogram abundance in 99 OpenPedCan H3 K27-altered DMG.** The recurring
OPC-like > hazard signal in the NNLS fits is not interpretable while the OPC-like proportion does not track OPC-like gene expression.
A proper bulk deconvolution (BayesPrism / platform-aware signatures from TPM-like pseudobulk) is the follow-up.

**Audit corrections (2026-10-20; P7-07, P7-08, P7-09).** Re-run by `scripts/25b_survival_rankscore_audit.py` →
`results/tables/survival_rankscore_audit.csv` (21 fits) and `openpedcan_rankscore_audit.csv` (129 specimens × 7 families).

*Scope limit — the NNLS/deconvolution models were NOT re-run.* They require `results/tables/c2l_signatures_mp.csv`, a
cell2location reference built by `13_reference_signatures.py --use-metaprograms` (GPU regression) that still carries the pre-audit
families — its malignant column is `malignant_novel_MP15` where the corrected families name `novel_MP9`. Regenerating it was spatial
finding **SP-06**; fitting the stale matrix would label one programme as another. Only the rank-score path, which reads
`metaprogram_families.csv` and `metaprogram_genes.csv` directly, was re-run.

*Status update (2026-10-21, reporting pass).* SP-06 has since been closed: the spatial track re-deconvolved all four cohorts against
the corrected metaprogram-family reference and wrote `c2l_signatures.csv` (2026-10-20 21:25). `c2l_signatures_mp.csv` itself was **not**
regenerated and still dates from 2026-10-17, so the NNLS survival sets remain un-refitted and **P7-07 stays partially closed**: the
signature matrix a re-fit would need now exists in corrected form, but the re-fit was not attempted in any track. Everything in
`survival_results.csv`, `openpedcan_deconvolution.csv` and `results/figures/phase7_survival_forest.png` is therefore pre-audit; Figure 6c
and Supplementary Figure 8's rank-score rows use `survival_rankscore_audit.csv` instead.

- **Cohort descriptors (P7-08).** The Cox cohort is **n = 99, 94 deaths, median OS 367 d, median age 7.6 y, 64 pons**. The 7.7 y / 74
  pons figures above describe the 129 independent primary specimens, not the 99 modelled patients. `molecular_subtype` is present for
  all 99 and encodes TP53 for 61 (61.6 %).
- **Proportional hazards (P7-07).** A Schoenfeld test (`proportional_hazard_test`, rank time-transform) is now fitted for every model.
  **`rank:cycling` violates PH in all three specifications** (p 0.021 univariate, 0.026 adjusted, 0.016 TP53-adjusted); `rank:OPC_like`
  violates it univariately (p 0.041) and `rank:MES_like` in the TP53 model (p 0.039). The remaining 16 of 21 fits show no violation.
  TP53 is now carried as a sensitivity covariate rather than left unmodelled. BH is applied across the 7 families within each model.
- **Result (unchanged in substance).** The rank-score association with survival remains **null after multiple-testing correction**:
  no fit reaches q < 0.05 (minimum q 0.145 univariate, 0.171 adjusted, 0.200 TP53-adjusted). Four of 21 fits are nominally
  significant — AC-like HR 0.82 / 0.79 / 0.81 per SD (p 0.042 / 0.024 / 0.039 across the three models) and MES-like HR 0.78
  univariate (p 0.021) — both protective in direction. Against the pre-audit run the direction is preserved and the ranking shifts:
  MES-like univariate was p 0.011 (q 0.076) and AC-like adjusted p 0.053, so AC-like strengthened and MES-like weakened once the
  corrected gene sets were used. The pre-audit rank set also scored the superseded `novel_MP15` programme in place of `novel_MP9`.
- **Atlas overlap (P7-09).** Removal covers only specimens flagged `in_atlas_scpca` — **1 of the 129** cohort rows. DeSisto patients
  (Children's Colorado, a CBTN site) still cannot be checked against this cohort, because the deposit carries no CBTN participant id.

## 3. Foundation-model embeddings (`26_foundation_embeddings.py`; A10G, 14 min incl. weight downloads; Geneformer inference 9.5 min)
15,000 malignant cells (2,500 per family: AC, MES, OC, OPC-like, cycling, hypoxia; confident members, max z ≥ 1; 1,647 cells / 13,353
nuclei). scGPT whole-human (zero-shot, 18,635 / 18,882 genes in vocabulary, 512-d); Geneformer V2-104M (18,785 genes in vocabulary,
<cls> + top 4,095 rank tokens, penultimate-layer mean, 768-d). `foundation_model_comparison.csv`, embeddings in
`data/processed/foundation_embeddings.npz`.

*Pre-audit table — superseded 2026-10-20 by P7-01; the corrected table is in the Audit corrections block below.*

| embedding | kNN family balanced acc. | nuclei → cells transfer | silhouette | same-sample neighbours (exp. 0.021) | same-assay neighbours (exp. 0.805) |
|---|---|---|---|---|---|
| PCA (50) | 0.825 | 0.615 | 0.104 | 0.846 | 0.993 |
| scVI (atlas malignant latent, 30) | 0.629 | 0.622 | 0.055 | **0.295** | **0.678** |
| scGPT zero-shot (512) | 0.757 | 0.582 | 0.032 | 0.488 | 0.995 |
| Geneformer V2-104M (768) | 0.745 | **0.649** | −0.043 | 0.603 | 0.997 |

PCA's top family accuracy is circular (labels are scores on the same log-normalised expression). Both foundation models separate the
states about as well as each other (0.75–0.76) without any training, but their neighbourhoods are almost purely assay-specific (0.995–0.997)
and largely sample-specific; scVI is the only embedding that mixes samples and assays, at the cost of lower within-dataset family accuracy.
Cross-assay transfer is similar for all (0.58–0.65; Geneformer best). Methods-note conclusion: zero-shot foundation-model embeddings do not
remove the cell/nucleus assay effect that dominates this atlas; batch-aware integration (scVI) remains necessary.

**Audit corrections (2026-10-20; P7-01).** The pre-audit comparison zero-filled NaN scVI latents and so understated scVI. The comparison was
rebuilt by `scripts/26b_foundation_comparison_rebuild.py` on a single evaluation set: cells with **finite latents in all four embeddings and a
confident family assignment — 8,688 of the 15,000 candidates** (54 samples; 7,662 nuclei / 1,026 cells; `foundation_eval_cells.csv`,
`in_eval`). The scGPT and Geneformer embeddings were re-intersected, not recomputed. 2,557 of the 8,688 cells changed family relative to
pre-audit, so the family labels are the corrected seven as well.

| embedding | kNN family balanced acc. | nuclei → cells transfer | silhouette | same-sample neighbours (exp. 0.023) | same-assay neighbours (exp. 0.792) |
|---|---|---|---|---|---|
| PCA (50) | 0.667 | 0.284 | 0.067 | 0.850 | 0.994 |
| scVI (atlas malignant latent, 30) | 0.658 | **0.495** | **0.087** | **0.388** | **0.937** |
| scGPT zero-shot (512) | 0.595 | 0.290 | −0.004 | 0.458 | 0.996 |
| Geneformer V2-104M (768) | 0.561 | 0.329 | −0.090 | 0.558 | 0.997 |

What the corrected numbers support: **scVI leads on cross-assay label transfer by a wide margin (0.495 against 0.284–0.329)**, has the
highest silhouette (0.087) and is the only embedding that mixes samples (0.388 against 0.458–0.850) and assays (0.937 against 0.994–0.997).
PCA's lead on within-set family accuracy is now **0.0085** (0.666985 vs 0.658485; 0.667 vs 0.658 as displayed above) rather than the
pre-audit 0.196, and remains circular — the labels
are scores on the same log-normalised expression; no confidence interval was computed for that difference. Two pre-audit statements are
reversed: Geneformer is **not** the best at cross-assay transfer (0.329, second to scVI), and scVI's sample/assay mixing does **not** come
"at the cost of lower family accuracy". The zero-shot foundation embeddings rank below both baselines on within-set family accuracy (0.595 / 0.561 vs 0.667 / 0.658) and on
silhouette (−0.004 / −0.090 vs 0.067 / 0.087). On the other metrics they beat PCA but not scVI: cross-assay transfer 0.290 / 0.329
against PCA 0.284 and scVI 0.495, and same-sample neighbours 0.458 / 0.558 against PCA 0.850 and scVI 0.388.
The methods-note conclusion is unchanged in substance: zero-shot foundation-model embeddings do not remove the cell/nucleus assay effect
that dominates this atlas, and batch-aware integration remains necessary.

## Outputs
Tables (`results/tables/`): `atac_summary.csv` (md5 e19658e2…), `atac_label_correlation.csv`, `atac_family_da_peaks.csv`,
`atac_separability_per_sample.csv`, `atac_pseudobulk_cpm.csv.gz`, `atac_union_peaks.csv.gz`, `openpedcan_dmg_cohort.csv`,
`openpedcan_deconvolution.csv`, `survival_results.csv` (md5 083bb6ab…), `foundation_model_comparison.csv` (md5 9216bd7d…).
Data: `data/processed/foundation_embeddings.npz` (md5 70135936…). Models: `models/foundation/` (scGPT_human, Geneformer V2-104M cache).
Figures (`results/figures/`): `phase7_atac_separability.png`, `phase7_survival_forest.png`, `phase7_foundation_models.png`.

## Caveats / open items
- ATAC: Jessa only (10 samples); Okada snATAC not fetched; DA calls not depth-normalised; no motif analysis (needs genome FASTA + JASPAR).
- Survival: deconvolution not validated; 99 patients, 94 events, single consortium (PBTA); treatment not modelled.
- Foundation models: zero-shot only (no fine-tuning); 15k-cell subset; labels expression-derived.
