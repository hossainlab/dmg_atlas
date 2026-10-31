# Peer review: "A cross-cohort single-cell reference atlas of H3K27-altered diffuse midline glioma"

**Material reviewed.** `preprint/dmg_preprint.tex` in full (main text, Methods, back matter, Supplementary Information), all six table sources in `preprint/tables/`, main Figs 1–6 (`fig1`–`fig5`, `fig7_robustness`) and Supplementary Figs S1–S7, all viewed as images. I checked numbers against `results/tables/` (including `results/tables/quality/`, `results/qc/` and `docs/phases/phase5_mapping_notes.md`), and Methods parameters against `scripts/`, `src/dmg_atlas/`, `modal/` and `configs/qc.yaml`. Line numbers refer to the .tex source. I did not load any .h5ad file.

---

## Summary

The authors combine six public DMG single-cell and single-nucleus RNA-seq studies into an scVI/scANVI reference: 311,618 cells and nuclei from 57 samples and 52 patients, of which 88% are nuclei and half come from one post-mortem series. They project five further cohorts onto the reference without retraining (two Smart-seq2 TPM tumour cohorts, one immune-focused cohort, fetal hindbrain and adult pons). Malignant cells are called by combining marker scores with inferCNV-style profiles, and these calls agree with the original authors' labels in a median 0.898 of cells per sample. Per-sample cNMF yields 15 metaprograms. Under a replication rule the authors state was fixed in advance (≥3 publications, and both cells and nuclei), only cycling and hypoxia replicate. Every lineage-like programme is found only in nuclei, but only six whole-cell samples were available. Twenty-one DMG Visium sections from 12 patients, deconvolved with cell2location, support a pre-specified MES-like–myeloid association and do not support an OPC-like–oligodendrocyte association. A leave-one-publication-out analysis shows that held-out-patient accuracy (0.975, balanced 0.917) overstates performance on an unseen study (pooled balanced accuracy 0.506). On an unseen study, normal neurons, astrocytes and vascular cells are often called malignant.

The work is careful in many respects. The paper is open about its negative results, the leave-one-study-out evaluation is a real contribution, and most numbers I checked trace exactly to scripted outputs. The main weaknesses are these:

- The central "replication" result depends on a rule that six whole-cell samples cannot satisfy for lineage programmes.
- Several analyses (the per-sample family composition, the malignant-state composition contrasts and the spatial reference) use a within-sample standardisation that, by the authors' own code comments, removes between-sample signal.
- The spatial statistics overstate consistency.
- I found a few numerical errors, including one stale pre-correction value.
- The Methods leave out parameters needed to reproduce the work, notably how TPM data were fed to a count model and how the copy-number reference was actually built.

## Recommendation

**Major revision.** As a resource, an openly released, benchmarked cross-cohort DMG reference is useful and, as far as I know, new. The robustness section is more candid than most atlas papers. However:

- Some headline claims are not supported as stated: "cycling is the only metaprogram recovered whichever study is removed", "small and consistent in sign", and "99.3% of Liu malignant cells".
- Figure 3d and the malignant-state rows of Fig. 4c show quantities that cannot carry the meaning given to them.
- Key Methods are incomplete or differ from the code.

None of these needs new data. All can be fixed by reanalysis of existing outputs and by rewriting.

| Dimension | Score /10 |
|---|---|
| Science (novelty, soundness, support for conclusions) | 6 |
| Methods reproducibility | 6 |
| Presentation (structure, figures, tables) | 6 |
| Language | 7 |

---

## Major concerns

**M1. The replication rule cannot test lineage programmes, and the masking result that partly answers this is not reported (lines 200–216, 331–336; Table 3; Fig. 3; Supplementary Fig. 5).**

- Only six whole-cell samples (DeSisto 1, Jessa 4, Zang 1) contributed to programme discovery, so the "both cells and nuclei" criterion is close to unattainable for any programme that is not universal.
- Under the rule, a single whole-cell sample is enough to count as "cells". Under the stricter variant computed by the code (`replicated_strict`, ≥2 samples per assay class), the masked analysis gives 2 replicated metaprograms, not 3 (`metaprogram_masking_summary.csv`, cut 0.85).
- Supplementary Fig. 5a shows that the third masked replicator is the OPC-like programme (mMP4 [MP3]). This is the most biologically important finding of the masking control, because it suggests that assay-biased genes, not absence of the programme, kept OPC-like out of whole cells. The main text never says so.
- The random-mask null is not comparable. Random masks collapse the structure to a median of 2.5 metaprograms, against 12 with the assay mask and 15 unmasked. "3 replicated versus a median of 1" therefore compares a structured clustering with a degenerate one.

*Remedy:*
- Report which metaprograms replicate after masking, and report the strict-rule counts.
- Replace the random null with masks matched to the assay-biased genes for mean expression and variance, or with a permutation of assay labels in the pseudobulk differential-expression step.
- Rephrase the headline to: "the rule could be satisfied only by programmes present in the few whole-cell samples".
- State plainly that the data cannot distinguish assay effects from limited power. Line 333 already does this; the abstract and the Fig. 3 section heading should match.

**M2. The abstract and Discussion misstate the leave-one-study-out metaprogram result (lines 62–63 and 332).**

- `quality/metaprogram_lopo_summary.csv` shows that MP2 (ribosomal), MP5 (OC-like) and MP6 (OPC-like) are recovered in 5 of 5 drop-one scenarios, exactly like MP1. Fig. 6d (`fig7_robustness`) shows the same.
- "Cycling is the only metaprogram recovered whichever study is removed" is therefore false.
- The correct statement is about *replicated* metaprograms: of the two replicated metaprograms, only cycling is recovered in every drop-one scenario. Hypoxia is lost without Jessa et al.

*Remedy:* correct the abstract and Discussion. Also report whether each metaprogram still *replicates* in each drop-one scenario, not only whether it is recovered, since recovery uses a Jaccard ≥0.5 match without applying the rule.

**M3. Within-sample standardisation removes the signal that Fig. 3d, the malignant-state rows of Fig. 4c and the spatial reference are interpreted to show (lines 225–226, 473–477; Fig. 3d; Fig. 4c).**

- Families are assigned by z-scoring UCell scores *within each sample*. The authors' own script (`scripts/19_metaprogram_scores.py`, docstring) states that this "removes between-sample differences in family composition … [and] is NOT a composition measurement".
- Nevertheless Fig. 3d (from `metaprogram_family_composition.csv`, the within-sample assignment) is captioned "Per-sample metaprogram-family composition of malignant cells". It shows hypoxia and ciliated cells at roughly 10–20% in nearly every sample, which contradicts Table 3: ciliated MP9 was found in 6 samples and hypoxia in 15.
- The same within-sample construction is used for the malignant-state rows of Fig. 4c (`scripts/18_composition.py`, audit note MC-04), whose near-zero coefficients follow from it.
- The cell2location reference takes "cells with a highest family score of at least 1" (line 477). It therefore includes, for example, "hypoxia" cells from samples whose hypoxia programme is weak, which dilutes the signatures.

*Remedy:*
- Show Fig. 3d using the within-assay-class assignment the code already produces (`metaprogram_family_composition_assayz.csv`), or drop the panel.
- Remove the malignant-state rows from Fig. 4c, or recompute them on a between-sample-comparable scale.
- For cell2location, justify the reference construction, or rebuild it from cells whose programme is present in that sample, and show a sensitivity analysis.

**M4. The spatial statistics overstate consistency (lines 257–265; Table 4; Supplementary Fig. 7).**

- "These effects are small and consistent in sign" is not supported. Only 13–14 of 21 sections are positive for the significant pairs (Table 4, `n_sections_positive`), and I² is 0.95–0.997.
- The pooled ρ (0.393 for AC-like–astrocyte) is three times the plain mean of per-section ρ (0.130; `spatial_colocalisation.csv`, `rho_mean`). This is because averaging within patient turns the ten brainstem sections into a single unit. The choice is defensible, but it must be shown.
- With 12 patients and I² near 1, DerSimonian–Laird confidence intervals are known to be too narrow. A Hartung–Knapp–Sidik–Jonkman (or t-based) interval, or a patient-level permutation, should be the primary analysis.
- Spearman correlations between closed proportions are biased by the closure. The CLR variant has already been computed (`metric == 'clr'`) but is not reported.
- The AC-like–astrocyte and OC-like–oligodendrocyte associations, and all four significant associations involving the ciliated programme (which the authors suspect is ependymal contamination), may come from signature collinearity in deconvolution rather than from tissue organisation. They need a check, such as the correlation between reference signatures or a simulation with known proportions.
- The Results describe the ten brainstem sections as "serial sections". They come from different levels of the neuraxis (midbrain to spinal cord), so "sections from ten levels of one brainstem and spinal cord" is accurate.

*Remedy:*
- Report the HKSJ intervals, the CLR results and the per-patient estimates (a forest plot).
- Delete "consistent in sign".
- Add a deconvolution-collinearity control.
- For the two pre-specified hypotheses, report unadjusted and two-test-adjusted p values alongside the 42-pair BH q values.

**M5. Classifier evaluation is restricted to the easy cells, and the leave-one-study-out result mixes dataset shift with loss of training data (lines 169–175, 286–298, 340–345).**

- Both the held-out-patient accuracy (0.975) and the LOSO metrics are computed only on cells with a confident consensus label. That excludes the 28.9% ambiguous cells (32% in the Okada fold), which are the hardest cases.
- The consensus label requires agreement with the author label where one exists (`cnv.consensus`), so the "truth" is partly the authors' annotation.
- Holding out Okada et al. removes half the training data and every autopsy sample. The resulting balanced accuracy of 0.402 therefore cannot be attributed to "a new study" alone.
- The pooled balanced accuracy of 0.506 is cell-weighted and dominated by that fold.
- For Zang et al., balanced accuracy averages only two classes (0.988).
- "These errors were made with high confidence: accuracy among predictions with probability of at least 0.9 was 0.912" is (i) for the Okada fold only and (ii) not a measure of confident *errors*.

*Remedy:*
- Report the fraction of misclassified cells with maximum probability ≥0.9, per fold.
- Add a size-matched control: drop a random subset of patients of the same size from the remaining publications.
- Report fold-wise (unweighted) means.
- State how many classes enter each balanced accuracy.
- Report how ambiguous cells are labelled by the model (Supplementary Table).

**M6. The specificity of the reference classifier on normal brain is poor, and the paper should say so more prominently (lines 182–187; Fig. 2d).**

- 91.3% (per-sample median) of first-trimester hindbrain cells are classified malignant, including 80% of fetal neurons and 94% of fetal oligodendrocyte-lineage cells.
- 29% of adult OPCs and 30% of adult ependymal cells are also classified malignant.
- Combined with M5, the classifier's "malignant" label is close to "neither a mature reference lineage nor immune".
- The Discussion recommends copy-number evidence for normal-lineage calls; that recommendation belongs in the abstract and the data-release documentation.

*Remedy:*
- Consider adding normal progenitor and ependymal classes (for example from the Braun reference) to the reference, or a rejection option.
- At minimum, report the false-malignant rate per normal class in a table and qualify "lineage labels transfer" in the abstract.

**M7. Methods are incomplete or differ from the code in ways that affect reproducibility.**

- **Copy-number reference.** Methods (lines 414–415) say that samples with fewer than 30 reference cells used a pooled reference from the same study. The code (`scripts/05_cnv_calls.py`) works differently:
  - *every* sample's reference combines its own immune marker cells, a pooled immune anchor (≤300 cells per type, stratified by assay and preparation, not by study), and tumour-purged lineage cells pooled from *other* samples (≤50 per sample per type);
  - the purge uses a stricter first-stage calling (`CALL_PARAMS_STAGEA`: z 5/2.5, floor 0.02);
  - samples with fewer than 30 reference cells are not run (`MIN_REF`).
- **TPM cohorts.** Methods do not say how TPM matrices (Filbin, Liu, part of LaBelle) were turned into inputs for a negative-binomial count model. The code scales TPM by a factor *m* (`scripts/11e_filbin_m_sensitivity.py`), and the median neighbour distance rises steadily with *m* (2.38–3.29; mapping notes §1). This choice directly affects the "high mapping uncertainty" result in Fig. 2b.
- **Composition.** The primary composition analysis uses scANVI labels (`18_composition.py`: "celltype: qc_pass cells by scanvi_label"). This is not stated. It is only implied at line 306.
- **Unreported contrasts.** The site, age and ONC201 contrasts shown in Fig. 4c are not described in Methods or Results, and the ONC201 annotation has no source.
- **Undocumented analyses.** Spot niches (Supplementary Fig. 7a), the brainstem-gradient analysis and the DMG versus glioblastoma test (Mann–Whitney, by the output table) have no Methods.
- **Missing details.** Harmony and PCA settings, UMAP parameters, the "Seurat v3" HVG variant actually used (`seurat_v3_paper`), the additional Jessa author-QC filter, the exclusion of `okada__1734_pons`, the fallback rule for cNMF rank selection, the immune-contamination technical class and the provenance hierarchy used to label metaprograms are all absent.
- **cNMF rank range.** 21 of 50 samples selected k = 9, the top of the searched range (`cnmf/k_selection.csv`). The range is binding, so programme resolution was capped. Extend it or justify it.
- **Cell2location training.** The realised number of epochs per cohort is "recorded" but not reported. The current brainstem run used 10,154 epochs against a target of 20,000.

*Remedy:* rewrite the Methods to match the code, add the missing sections, and report the realised training lengths. A parameter table in the Supplement would help.

**M8. Composition inference: units, few clusters and selective reporting (lines 239–243, 491–497; Fig. 4c).**

- The H3.1 arm contains 7 samples from 6 patients. With so few clusters, cluster-robust standard errors are anti-conservative, and the mixed-model p value for the headline neuron effect is missing (`p_mixed` NaN in `composition_results.csv`).
- Because the model includes publication fixed effects and all H3.1 samples are Jessa et al. nuclei, the contrast is estimated within Jessa et al. It is therefore confounded with assay and site rather than "with study" as written.
- The Results report only the neuron effect. Fig. 4c and the robustness text ("all four significant effects") also rely on astrocyte and myeloid H3.1 effects (q = 0.038) that the reader never sees.
- The same output file has further significant results that are not mentioned: AC-like family adult versus paediatric (q = 0.013), H3.1 mature-oligodendrocyte subcluster (site-adjusted q = 0.008), ONC201 tip endothelium (q = 0.047). The last one contradicts "No other H3-variant difference remained significant after adjustment for site".

*Remedy:*
- Report all tested contrasts in a Supplementary Table.
- Define the family over which BH is applied.
- Use a wild-cluster bootstrap or a mixed model as the primary inference for contrasts with fewer than about 10 clusters per arm.
- Correct the sentence on site adjustment.

**M9. Numerical errors (details in the spot-check table).**

- "99.3% of the cells that Liu et al. annotated as malignant were classified as malignant" (line 180). The current output gives 0.951 (`mapping_author_agreement.json`). 0.993 is the pre-correction value that the authors' own audit withdrew (`docs/audit/integration_mapping_audit.md`).
- "Their median neighbour distance was twice that within the reference (3.28 versus 1.57)" (line 186). 3.279 is the distance from *malignant reference cells to fetal cells*. The distance from fetal cells to the reference is 2.604 (mapping notes, deviation block). The sentence attributes the wrong quantity to the fetal cells.
- The abstract's metaprogram claim (M2).

*Remedy:* correct the values, and regenerate every number in the text from the result tables with a script.

---

## Minor concerns

1. **Line 78.** "No systemic therapy has improved it". Please check the current regulatory status of dordaviprone (ONC201) for H3 K27M-mutant diffuse midline glioma (I understand it received FDA accelerated approval in 2025) and phrase accordingly. Line 85, "nor cellular immunotherapy has produced durable remission", should also be hedged, since GD2-CAR-T studies reported individual durable responses.
2. **Eligibility.** The H3 variant is unknown for 25 of 57 samples. State how H3 K27 alteration was confirmed for every included sample (author metadata, histology, an H3K27M expression proxy), particularly for Okada et al., DeSisto et al. ("paediatric HGG") and Zang et al.
3. **Search.** The data search (line 367) gives counts but no search strings, screening flow or dates. Add a flow diagram or Supplementary Note.
4. **Table 1.** The Patients column sums to 56, but the total is 52, because the two Jessa series share four patients. Add a footnote. The overlap is first mentioned only at line 289.
5. **Line 180.** "42 samples from 41 patients" follows both the Filbin and Liu medians but applies to Liu only (Filbin has 9 samples).
6. **Line 184.** "a median of 91.3%" for fetal cells is a per-sample median over 17 donors. The pooled value is 87.7%. Say which.
7. **Line 359.** "18.9% of malignant calls rely on the aggregate copy-number score" holds for CNV-malignant cells across all ten studies, queries included. Among consensus-malignant reference cells the value is 11.2%. Specify the denominator.
8. **Line 155 and Table 2.** "324,676 author-annotated cells in 97 samples" pools reference and query cohorts (Liu, LaBelle). The abstract presents it as a reference property. LaBelle has precision and recall shown as "nan".
9. **Line 158.** The tuning and held-out comparison (0.882 vs 0.885) is pooled. Also give the sample medians (0.897 vs 0.900).
10. **"Calibrated malignant calls" (lines 107, 327).** No calibration (reliability) analysis is shown. Use "validated" or "benchmarked".
11. **Lines 129–135.** Seed selection among three identical configurations is not really model selection. The 60-epoch sensitivity run has a higher batch score (0.661 vs 0.537) and fine-label conservation (0.873 vs 0.659), but lower lineage-level conservation. Discuss whether the 26-epoch heuristic under-trained the model. iLISI for assay is 0.029, so assay mixing is poor. Say this more directly, since it bears on M1.
12. **Line 210.** "the only programme restricted to whole cells reflected dissociation stress". Correct, but say that MP13 appears in all three whole-cell publications, consistent with tissue dissociation.
13. **Line 216.** The ciliated programme has no significant overlap with any reference set (Supplementary Fig. 4b). If ependymal contamination is the leading explanation, test it directly by scoring an ependymal signature (e.g. *FOXJ1*) and checking CNV status of MP9-high cells.
14. **Lines 233–237.** Two of the 21 "populations" are artefact clusters ("ribosomal residue (cells)", "endothelial (neural genes)"; Fig. 4a). Exclude them or call them artefacts. "Several populations were derived almost entirely from nuclei (Fig. 4b)": Fig. 4b shows assay purity, not nuclear fraction (that is Fig. 4a, right). The striking drop in myeloid purity from 0.95 to 0.02 after per-compartment correction needs a sentence of interpretation.
15. **Line 268.** "tumour fraction increased with the AC-like and OPC-like shares" is a within-section relation. Fig. 5e plots section means of different classes (tumour, MES-like, cycling, oligodendrocyte). Make the figure match the claim.
16. **Line 509.** "core and edge sections from one resection, were each treated as one patient". No core/edge Visium pair exists in the spatial cohort (`spatial_patient_map.csv`). The core/edge pair is Collot DMG-PMC-9 in the snRNA reference, where it is counted as two cNMF samples from one patient. State the patient-level unit for cNMF counts too: "34 samples" for MP1 is not 34 patients.
17. **Missing citations.** LaBelle et al. (GSE227983, GSE271936) are named in the text and figures but never cited and have no bibliography entry. The same applies to GSE280990 (brainstem Visium), UCell, DerSimonian–Laird, Harmony's implementation (harmonypy), scib-metrics and Space Ranger.
18. **Supplementary Table 2 (software)** omits scib-metrics (0.5.10 per the code), harmonypy, scrublet, pynndescent and python-igraph. It lists lifelines, which no analysis in the manuscript uses.
19. **Line 391.** Doublet removal uses a fixed score threshold of 0.25 for every sample (`predicted_doublet_fixed`), not Scrublet's automatic threshold. Say so and give the removed fraction.
20. **Line 410.** The smoothing is infercnvpy's pyramidal (weighted) running mean, not a plain running mean.
21. **Line 435.** Min–max scaling was done across six embeddings, including the 60-epoch sensitivity run, which changes the scaled scores of the candidates. State this.
22. **Data availability.** Specify which files are in the Zenodo record (h5ad, models, gene lists), the licence, and a versioned code release (tag or DOI), not only a GitHub URL.
23. **Unused outputs.** The repository contains a multiome ATAC figure (`fig6_multimodal`) and a bulk-survival figure (`figS8`) that the manuscript does not use. If the paired ATAC data are out of scope, say so. If not, the chromatin evidence could address the assay question in M1.

---

## Language and style

The prose is generally clear and economical, with consistent British spelling. The main problems are: compressed sentences that pack several results together; project-internal vocabulary ("field", "label", "released model", "primary cut", "Tier-2"); results given inconsistently as fractions and as percentages; and abbreviations defined late or not at all.

**Abbreviations.** Several are used before their definition or never defined:
- MES-like is used in Table 3 and Methods before its definition at line 261.
- scANVI, scArches, kBET, cLISI, ASW, NMI, ARI, IEG, TAM, HGG, PDX, CSF, PBMC, SS2 and UMAP (used in the Fig. 1 caption; defined only at line 437) are never expanded or expanded too late.

**Units.** Use one convention for proportions. Line 156 has "a median 0.898 of cells", while line 179 has "median 90.3% per sample".

Specific edits:

| Location | Before | After |
|---|---|---|
| Abstract l.56 | "2 of 15 malignant expression metaprograms, cycling and hypoxia, recurred across three studies and in both cells and nuclei" | "Only 2 of 15 malignant metaprograms (cycling and hypoxia) recurred in at least three studies and in both cells and nuclei" |
| Abstract l.57–58 | "every lineage-like programme was detected only in nuclei, but only six whole-cell samples were available" | "every lineage-like programme was detected only in nuclei; however, only six whole-cell samples contributed to programme discovery" |
| Abstract l.60 | "Refitting the reference without each study showed that…" | "Refitting the reference with each study left out in turn showed that…" |
| l.155–156 | "the copy-number call agreed with the authors' labels in a median **0.898** of cells per sample" | "the copy-number call agreed with the authors' labels for a median of 89.8% of cells per sample" |
| l.159 | "so event support is released as a confidence field" | "we therefore provide arm-level event support as a per-cell confidence annotation" |
| l.165 | "is released as a sensitivity label" | "is provided as an alternative annotation for sensitivity analyses" |
| l.200 (heading) | "Two of 15 metaprograms replicate" | "Only two of 15 metaprograms meet the replication criterion" |
| l.203 | "which we clustered into 15 metaprograms, that is, programmes recurring across samples" | "which we clustered into 15 metaprograms (programmes recurring across samples)" |
| l.212–214 | "Masking 6,013 genes that differ between paired cells and nuclei yielded 12 metaprograms, 3 of which replicated, compared with a median of 2.5 metaprograms and 1 replicated over ten size-matched random masks" | "After removal of 6,013 genes that differ between matched cells and nuclei, 12 metaprograms remained and 3 replicated, including the OPC-like programme. Ten random masks of the same size left a median of 2.5 metaprograms, of which 1 replicated." |
| l.264–265 | "so these effects are small and consistent in sign rather than strong" | "so these pooled effects are modest and vary substantially between patients" |
| l.268 | "with a single patient, this observation is descriptive" | "because all ten sections come from one patient, this observation is descriptive only" |
| l.297–298 | "These errors were made with high confidence: accuracy among predictions with probability of at least 0.9 was 0.912." | "Many errors were confident: X% of misclassified cells had a maximum class probability of at least 0.9." |
| l.331 | "Most of our findings define the limits of what public DMG data can establish." | "Many of our findings mark the limits of what current public DMG data can establish." |
| l.341–342 | "Accuracy on held-out patients, a commonly reported measure, overstated performance on a new study" | "Accuracy on held-out patients, the measure most often reported, overstated performance on an unseen study" |
| l.369 | "Data sets" | "Datasets" (as in the Abstract) |
| Table 2 caption | "Tier-2 comparison of the copy-number call" | "Comparison of the copy-number call" |
| l.513 | "using the settings of the released model (seed 1)" | "using the settings of the final reference model (seed 1)" |

---

## Figures and tables

**General**
- Panel letters are uppercase in the figures (A, B, …) but lowercase in the captions. Use one style; lowercase bold is the Nature-family convention.
- Study names are inconsistent across figures and tables: "Jessa 2022", "Jessa 2022 multiome", "Jessa 2022 (sc/snRNA)", "Zang" and "DeSisto" without years, "LaBelle 2025" with a year but no citation.
- The ciliated programme is called "cilia (MP9)" in figures, "ciliated (MP9)" in tables and "novel-MP9" in Fig. 6d.

**Fig. 1**
- Panel e tick labels are software identifiers ("scvi_seed1", "scvi_seed0_e60", "harmony", "pca"). Use "scVI seed 1 (selected)", "scVI 60 epochs", "Harmony", "PCA".
- Panel b is titled "scANVI lineage label", but the caption says only "lineage". Say which label is shown.
- Panel f mixes reference and query cohorts without a separator.

**Fig. 2**
- Panel c shows author labels "score-too-low" and "frozen-OPC-like-b" without explanation.
- Panel b's x-axis runs 0–1 but is labelled "percentile". Use 0–100, or call it a "quantile".

**Fig. 3**
- The panel b title overlaps the "B" and "C" panel labels.
- In panel c, the caption ("matched cell- and nucleus-derived metaprograms") does not match the data. `metaprogram_cell_vs_nucleus.csv` pairs the cell-only IEG programme MP13 with nucleus programmes, and Jaccard values are at most 0.125. Explain what is matched, or drop the panel.
- Panel d: see M3. Its sample axis has no labels and its assay colour strip has no legend.
- The x-axis label "1 − mean Jaccard" in panel b differs from the Methods ("Jaccard distance").

**Fig. 4**
- Panel c: significance levels "*" and "**" are not defined. The ONC201, age and site contrasts are not described anywhere (M7, M8).
- Panel a includes the artefact clusters (minor 14).

**Fig. 5**
- The panel d title and x-axis label overlap panel e. The "Midbrain" tick in panel e is clipped.
- Panels a–c do not say which section is shown.
- Panel b has a band of high MES-like share along the lower edge of the capture area that looks like an edge artefact. Comment on it or choose another section.
- The pre-specified OPC-like–oligodendrocyte pair is absent from panel d, although it is one of the two tested hypotheses.

**Fig. 6 (`fig7_robustness`)**
- Panel d labels are internal ("MP9 novel-MP9", "MP13 IEG-dissociation-stress").
- The panel e legend uses the same circle symbol for "H3.1 vs H3.3" and "black edge"; use separate legend entries for shape and edge.

**Supplementary figures**
- **Fig. S2.** Row labels are raw identifiers ("asw_batch__sample_id", "ilisi__assay"). The colour bar says "scaled 0–1 by scIB", but the cells show raw metric values, and the summary rows show min–max-scaled scores. Separate the two.
- **Fig. S4c.** The "programmes assigned" series (right axis, dashed) is missing from the legend.
- **Fig. S6.** The colour bar says log2 fold change; the caption says "log fold change".
- **Fig. S7a.** Niche clustering is not described in the Methods. Two niches cover 0% of DMG spots.
- **Fig. S7b.** The caption says "each point is a section", but each pair shows about 63 points, which equals 21 DMG sections × 3 metrics (proportion, CLR, raw q05) in `spatial_colocalisation_per_section.csv`. Plot one metric.

**Tables**
- **Table 1.** Add a footnote on patient overlap (minor 4) and give the age split.
- **Table 2.** "Tier-2" is jargon. Replace "nan" with "–" and explain it. Separate reference rows from query rows, or drop query rows from the "All" row.
- **Table 3.** "ciliated (MP9)" in the Family column is redundant. Add the number of patients and publications, and mark MP1/MP4 as replicated also under the strict rule.
- **Table 4.** Format q values in scientific notation (9.4 × 10⁻⁵, not 9.38e-05). Add the plain section mean and the HKSJ interval (M4). Include the pre-specified OPC-like–oligodendrocyte row even though it is not in the top ten.
- **Supplementary Table 1.** "DIPG single-cell cohort (2026) – EGA – not deposited" is not verifiable. Give a reference.

---

## Number and Methods spot-checks

Source paths are relative to `results/tables/` unless stated otherwise.

### Numbers

| # | Claim or parameter | Manuscript value | Source value (file) | Match |
|---|---|---|---|---|
| 1 | Reference cells / samples / patients | 311,618 / 57 / 52 | 311,618 / 57 / 52 (`atlas_v1_summary.csv`) | yes |
| 2 | Cells before / after QC | 420,381 / 356,406 | 420,381 / 356,406, summing included studies (`results/qc/qc_summary.json`) | yes (but includes query cohorts and ScPCA; the step to 311,618 also includes doublet removal, which the text does not say) |
| 3 | Shared genes | 18,882 | 18,882 (`qc_summary.json`) | yes |
| 4 | Nuclei fraction | 88.4% | (229,043 + 46,288) / 311,618 = 88.4% | yes |
| 5 | Pons / thalamus samples | 42 / 15 | 42 / 15 (`atlas_v1_summary.csv`) | yes |
| 6 | H3 variant known; H3.1 / H3.2 / H3.3 | 32; 8 / 1 / 23 | 32; 8 / 1 / 23 | yes |
| 7 | Jessa cells / patients | 129,364 / 23 | 83,076 + 46,288 = 129,364; 27 − 4 shared = 23 | yes |
| 8 | Two publications' share of cells | 92% | 91.9% | yes |
| 9 | scVI seed 1 total / batch / bio | 0.739 / 0.537 / 0.874 | 0.7395 / 0.5373 / 0.8743 (`integration_model_choice.json`) | yes |
| 10 | Harmony / PCA total | 0.581 / 0.126 | 0.5809 / 0.1258 | yes |
| 11 | iLISI assay, scVI / Harmony | 0.029 / 0.134 | 0.0285 / 0.1344 (`integration_metrics_summary.csv`) | yes |
| 12 | CNV agreement: median (IQR), pooled, precision, recall | 0.898 (0.793–0.978), 0.884, 0.934, 0.909 | identical (`malignant_call_validation_summary.csv`) | yes |
| 13 | Tuning vs held-out agreement | 0.882 vs 0.885 | 0.8820 vs 0.8852 (pooled) | yes |
| 14 | Event-backed vs score-only precision | 0.952 vs 0.862 | 0.9519 vs 0.8624 | yes |
| 15 | Consensus malignant / normal / ambiguous | 170,940 / 50,754 / 89,924 | identical (`atlas_core_obs_labels.csv`) | yes |
| 16 | Permissive variant adds | 28,252 | 199,192 − 170,940 = 28,252 | yes |
| 17 | Malignant calls without arm event | 18.9% | 18.9% of CNV-malignant cells over all 10 studies (`phase3_summary.json`); 11.2% of consensus-malignant reference cells | partly (denominator unstated) |
| 18 | Held-out-patient accuracy / balanced / cells | 0.975 / 0.917 / 23,827 | 0.9751 / 0.9166 / 23,827 (`quality/loso_classifier_summary.csv`) | yes |
| 19 | Filbin / Liu median malignant | 90.3% / 83.6% | 0.9034 / 0.8361 (`mapping_summary.csv`) | yes |
| 20 | Liu author-malignant classified malignant | 99.3% | 0.951 (`mapping_author_agreement.json`); 0.993 is the withdrawn pre-correction value | **no** |
| 21 | LaBelle median / pooled malignant | 49.9% / 11.1% | 0.4986 / 0.1111 | yes |
| 22 | Siletti / Braun median malignant | 1.2% / 91.3% | 0.0124 / 0.9127 | yes |
| 23 | Fetal vs reference neighbour distance | 3.28 vs 1.57 | reference 1.566 (`phase5_reference_knn_distances.csv`); Braun→reference 2.604; 3.279 is malignant→Braun (`docs/phases/phase5_mapping_notes.md`) | **no** |
| 24 | cNMF malignant cells / samples | 167,012 / 50 | 167,012 / 50 (`cnmf/samples_analysed.csv`, `k_selection.csv`) | yes |
| 25 | MP1 / MP4 samples | 34 / 15 | 34 / 15 (`replication_matrix.csv`) | yes |
| 26 | Assay-biased genes; masked MPs / replicated; random median | 6,013; 12 / 3; 2.5 / 1 | 6,013 (`assay_biased_genes.csv`); 12 / 3 (strict: 2); 2.5 / 1 (`metaprogram_masking_summary.csv`) | yes (strict count unreported) |
| 27 | Replicated count stable across cuts 0.84–0.90 | 2 | 2 at every cut 0.84–0.90; 3 at 0.91 | yes |
| 28 | Neuron H3.1 vs H3.3 coefficient; q; site-adjusted q; pons q | 2.60; 0.004; 0.017; 0.050 | 2.596; 0.0040; 0.0169; 0.0495 (`composition_results.csv`) | yes |
| 29 | Neuron nuclei vs cells coefficient; q | 3.40; < 10⁻⁹ | 3.402; 8.5 × 10⁻¹¹ | yes |
| 30 | "No other H3 difference significant after site adjustment" | none | oligodendrocyte_2 H3.1 site-adjusted q = 0.008 (TME set) | **no** |
| 31 | Spatial sections / patients | 27 total; 21 DMG / 12 patients | 27; 21 / 12 (`spatial_patient_map.csv`) | yes |
| 32 | AC-like–astrocyte ρ, q; AC-like–oligodendrocyte ρ, q | 0.393, < 0.001; 0.308, 0.005 | 0.3929, 9.4 × 10⁻⁵; 0.3076, 0.0051 (`spatial_colocalisation.csv`) | yes |
| 33 | MES-like–myeloid ρ (CI), q | 0.238 (0.071–0.392), 0.034 | 0.2380 (0.0707–0.3923), 0.0339 | yes |
| 34 | OPC-like–oligodendrocyte ρ (CI), q | 0.178 (−0.028 to 0.370), 0.210 | 0.1781 (−0.0280 to 0.3697), 0.2096 | yes |
| 35 | I² range | 0.952–0.997 | 0.9524–0.9973 | yes |
| 36 | "Consistent in sign" | — | 11–14 of 21 sections positive for the significant pairs | **no** (claim unsupported) |
| 37 | DMG vs GBM tumour fraction; cycling share; q | 0.528 vs 0.460; 0.112 vs 0.214; 0.667 | 0.5284 vs 0.4605; 0.1118 vs 0.2143; 0.667 both (`ren_dmg_vs_gbm_composition.csv`) | yes |
| 38 | LOSO accuracy range / balanced range / pooled | 0.880–0.992 / 0.402–0.988 / 0.506 | 0.8803–0.9924 / 0.4015–0.9878 / 0.5060 | yes |
| 39 | Malignant / myeloid / astrocyte recall ranges | 0.86–1.00 / 0.98–1.00 / 0.00–0.56 | 0.861–1.000 / 0.977–1.000 / 0.000–0.560 (`quality/loso_classifier_per_class.csv`) | yes |
| 40 | Okada neurons, astrocytes, vascular called malignant | 3,838 of 3,898; 63%; 51% | 3,838 / 3,898; 1,943 / 3,065 = 63.4%; 1,328 / 2,592 = 51.2% (`quality/loso_classifier_confusion.csv`) | yes |
| 41 | Accuracy at probability ≥ 0.9 | 0.912 (implied overall) | 0.912 is the Okada fold only; other folds 0.92–0.99 | partly (scope unstated) |
| 42 | Hypoxia replicated in balanced pools | 4 of 20 | 0.2 × 20 = 4 (`quality/metaprogram_lopo_summary.csv`) | yes |
| 43 | "Cycling is the only metaprogram recovered whichever study is removed" (abstract) | only MP1 | MP1, MP2, MP5, MP6 all recovered 5/5 | **no** |
| 44 | "MES-like, ciliated and four others recovered only with both large series" | 6 MPs | MP8, MP9, MP10, MP11, MP12, MP15 lost when Jessa or Okada is dropped | yes |
| 45 | "All four significant effects remain significant" | 4 | astrocyte, myeloid and neuron H3.1, and neuron assay: all q < 0.05, same sign, in all four labellings (`quality/composition_label_sensitivity_results.csv`) | yes (but three of the four are not reported in Results) |
| 46 | Table 1 patient total | 52 | column sums to 56 (4 shared Jessa patients) | inconsistent (needs a footnote) |

### Methods parameters against code

| # | Parameter | Manuscript | Code / config | Match |
|---|---|---|---|---|
| 1 | Nuclei / cells / Smart-seq2 min genes and max mito | 300 / 5%; 300 / 20%; 1,000 / 20% (frozen nuclei 5%) | identical (`configs/qc.yaml`; `qc.py`, author_prep_aware) | yes |
| 2 | Relaxed mito for four nucleus samples | 20% | `RELAXED_MITO = 20.0`, four samples (`cnv.py`) | yes |
| 3 | Max-genes quantile; Scrublet rate, HVG, threshold | 99.5th; 0.06, 3,000, 0.25 | 0.995; 0.06, 3,000, fixed 0.25 used for removal (`07_assemble_atlas_input.py`) | yes (fixed threshold not stated) |
| 4 | Low-confidence sample | median < 500 genes | `min_median_genes_pass: 500` | yes |
| 5 | Marker thresholds | 0.25 / 0.10 / 0.50 / 2× | `MARKER_THR` 0.25 / 0.10 / 0.50 / 2.0 | yes |
| 6 | CNV genes ≥10% of cells, clip ±3, window 100 / step 10, noise 1.5 SD | as stated | `run_infercnv` defaults identical | yes (pyramidal smoothing not stated) |
| 7 | MAD refinement, Leiden resolution | 3 MAD, 1.0 | `n_mad = 3.0`, `resolution = 1.0` | yes |
| 8 | Cluster calling | floor 0.008, z ≥ 2 single / ≥ 1 double, 20% reference guard | `CALL_PARAMS` identical | yes |
| 9 | Pooled CNV reference | only when < 30 reference cells, same study | always pooled immune anchor (assay-stratified) plus purged lineage cells from other samples; < 30 → sample not run | **no** |
| 10 | scVI architecture and training | 30 latent, 2 layers, NB, batch 1,024, ≤ 26 epochs, patience 5, 10% validation | identical (`modal/phase4_train.py`; 26 = scvi-tools heuristic) | yes |
| 11 | scIB subsample, k | 60,023, 90 | 60,023, 90 (`integration_model_choice.json`) | yes |
| 12 | Score weights and scaling | 0.4 / 0.6, min–max across candidates | identical; scaling includes the 60-epoch sensitivity run (`08_integration_metrics.py`) | yes (detail missing) |
| 13 | scANVI training | ≤ 25 epochs, patience 5, patient-grouped ~10% validation | identical (`modal/phase4b_scanvi_retrain.py`) | yes |
| 14 | Mapping uncertainty reference sample | 50,000 cells | 50,000 rows (`phase5_reference_knn_distances.csv`) | yes |
| 15 | cNMF | ≥ 200 cells, 2,000 HVG, random init, cd, 500 iterations, 10 restarts, k 4–9, silhouette 0.90, 80% restarts, top 50 | identical (`metaprograms.py`) | yes (fallback rule unstated; 21/50 samples at k = 9) |
| 16 | Metaprogram rules | cut 0.85, ≥ 5 programmes, ≥ 3 samples, ≥ 25% members; technical 15 RP / 5 MT / 10 IEG | identical (`16_metaprograms.py`, `TECH`) | yes (immune-contamination class unstated) |
| 17 | UCell max rank; family threshold; contamination | 1,500; 1; 2× | 1,500; `CONFIDENT_Z = 1.0`; `CONTAM_FOLD = 2.0` | yes |
| 18 | TME subclustering | Leiden 0.5, 15 neighbours, ≥ 50 cells; 2,000 HVG, 20 latent, ≤ 400 epochs, patience 20 | identical (`17_tme_subclusters.py`, `22_tme_scvi_clean.py`); batch = sample_id not stated | yes |
| 19 | cell2location | 59,993 cells, ≤ 8,000 per class, 250 epochs; 5 cells per spot, α 20, 2,000–20,000 epochs | 59,993; `CLASS_CAP = 8_000`; 250; 5; 20; clip(budget, 2,000, 20,000) (`modal/phase6_c2l.py`, `c2l_signatures_mp_run.json`) | yes |
| 20 | Spatial pooling; UMI sensitivity | Fisher z, DerSimonian–Laird; < 500 UMI | identical (`15b_spatial_tests_fix.py`, `MIN_UMI = 500`) | yes |
| 21 | Core/edge sections treated as one patient (spatial) | stated | no core/edge Visium pair in `spatial_patient_map.csv` | **no** (inapplicable) |
| 22 | LOSO surgery | default freezing, weight decay 0, 100 epochs | identical (`modal/phase9_loso_scanvi.py`) | yes |
| 23 | Primary labels for composition | not stated | `scanvi_label` (`18_composition.py`) | missing |
| 24 | TPM-to-count conversion for queries | not stated | TPM × m scaling (`11e_filbin_m_sensitivity.py`) | missing |

---

## Questions for the authors

1. How was H3 K27 alteration confirmed for the 25 samples with unknown H3 variant? Were any H3-wildtype or EZHIP-overexpressing midline gliomas included?
2. How were the TPM matrices converted for the negative-binomial encoder, and how sensitive are Fig. 2a–b to the scaling factor?
3. Does MP3 (OPC-like) replicate under the strict rule after masking? Which whole-cell samples carry it?
4. With the within-assay-class family assignment, how does the family composition of Fig. 3d change, and does any family differ between H3.1 and H3.3, or between thalamus and pons?
5. Does the pooled spatial ρ for MES-like–myeloid survive a Hartung–Knapp interval, the CLR metric and exclusion of the ciliated and AC-like classes from the deconvolution reference?
6. What fraction of misclassified cells in each LOSO fold had a maximum probability of at least 0.9? How does balanced accuracy change when Okada-sized random patient subsets are removed instead?
7. How do ambiguous cells (28.9%) distribute across scANVI labels? Are they enriched in specific samples or assays?
8. Were any pre-registration documents (for example the frozen design decisions) time-stamped publicly before analysis? If so, cite them. "Specified before analysis" carries weight only if verifiable.
9. Why was the cNMF rank range capped at 9, given that 21 of 50 samples selected the maximum?
10. Can the paired multiome ATAC data (Jessa et al.) independently test whether lineage programmes are absent from whole cells or only under-detected?
