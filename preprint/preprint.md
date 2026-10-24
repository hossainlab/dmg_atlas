# A cross-cohort single-cell and spatial atlas of human H3K27M-altered diffuse midline glioma

**Authors.** [Author list to be completed]

**Affiliations.** [To be completed]

**Correspondence.** [To be completed]

**Keywords.** diffuse midline glioma; H3K27M; single-cell RNA sequencing; single-nucleus RNA sequencing; spatial transcriptomics; data integration; reference atlas

---

## Abstract

H3K27M-altered diffuse midline glioma (DMG) is uniformly lethal, and the single-cell literature describing it is
fragmented across small cohorts that differ in dissociation protocol, assay and anatomical site. We assembled a
cross-cohort reference from six publicly deposited human DMG studies — 311,618 cells and nuclei from 57 samples and 52 patients across 18,882 genes shared by every study — and integrated them with a conditional variational autoencoder
(scVI), selecting the model by a pre-stated score combining batch mixing and biological conservation (0.739, versus
0.581 for Harmony and 0.126 for uncorrected principal components). Lineage labels transferred by scANVI reach 0.975
accuracy and 0.917 balanced accuracy on five entirely held-out patients. Malignant cells were called by combining
marker scoring with an expression-based copy-number inference; across 324,676 author-annotated cells in 97 samples the call agrees with the original authors in a median of 89.8% of cells per sample (IQR 79.3–97.8%), with no evidence of
tuning inflation (tuning 88.2% versus held-out 88.5%). Projecting five further cohorts into the reference without any
fine-tuning recovers expected labels in tumour and adult-brain data, including 1.2% false-malignant calls in adult pons
nuclei. Applying a replication rule stated before the analysis — a programme must recur in at least three studies and in
both cells and nuclei — only 2 of 15 consensus metaprograms replicate (cycling and hypoxia), and every lineage-like
programme is nuclei-restricted. In 21 Visium sections from 12 patients, astrocyte-like and ciliated malignant states
co-localise with normal astrocytes, and a pre-stated mesenchymal–myeloid association holds (rho 0.238, q = 0.034) while
a pre-stated oligodendroglial association does not (rho 0.178, q = 0.210). The atlas, its labels and the complete audit
trail are released as a reproducible resource.

## Introduction

Diffuse midline glioma with an H3 K27M alteration is a distinct WHO grade 4 entity of the pons, thalamus and spinal
cord, and remains incurable: median survival is under a year and no systemic therapy has altered that
(doi:10.1093/neuonc/noab106; doi:10.1016/j.ccell.2017.08.017). Single-cell studies have established that these tumours
are dominated by an oligodendrocyte-precursor-like (OPC-like) compartment with smaller astrocyte-like (AC-like) and
oligodendrocyte-like (OC-like) fractions, a picture first drawn from 3,321 cells across six patients
(doi:10.1126/science.aao4750) and subsequently refined with developmental anchoring in canonical and non-canonical H3
variants (doi:10.1038/s41588-022-01205-w). The conceptual framework — a hierarchy of cycling progenitor-like cells
differentiating into more mature states — parallels what has been described in adult glioblastoma
(doi:10.1126/science.1254257; doi:10.1016/j.cell.2019.06.024).

Three problems limit how far that literature can be pushed. First, every individual DMG cohort is small; the largest
public deposits contain tens of samples, and few contain more than ten patients, so sample-level inference is
underpowered and patient-level inference is rarer still. Second, the deposits differ in ways that are confounded with
the biology of interest: tumours of the pons are sampled at autopsy as frozen nuclei while thalamic tumours are more
often resected and dissociated as whole cells, so assay, site and clinical context cannot be separated within any one
study. Third, published cell states are typically derived per study, so it is unclear which are reproducible features
of the disease and which are properties of a protocol — a concern made concrete by pan-cancer work showing that
recurrent expression programmes can be recovered consistently only when they are tested across many tumours
(doi:10.1038/s41586-023-06130-4).

An integrated reference addresses the first problem directly and makes the second and third testable. Variational
autoencoder models trained with batch covariates (doi:10.1038/s41592-018-0229-2) and their semi-supervised extension
for label transfer (doi:10.15252/msb.20209620) are now standard for atlas assembly, with transfer-learning projection
allowing new datasets to be placed in an existing reference (doi:10.1038/s41587-021-01001-7) and established metrics
for judging whether batch effect was removed at the cost of biology (doi:10.1038/s41592-021-01336-8).

We applied this machinery to every human H3K27M DMG single-cell deposit we could obtain, and treated the limitations as
first-class results rather than caveats. We pre-registered a replication rule for expression programmes, kept patients
rather than cells as the unit of inference wherever a claim is made about tumours, and submitted the finished analysis
to an adversarial audit followed by an independent verification of every correction. The resulting resource is, to our
knowledge, the largest integrated human DMG single-cell reference to date, and the accompanying negative results
delimit what cross-cohort DMG data can currently support.

## Results

### A 311,618-cell integrated reference from six human DMG cohorts

We screened the public single-cell landscape for human H3K27M DMG and retained six studies as a reference cohort, with
five further cohorts reserved for projection-only validation (Methods). Quality control was applied per assay class, and
of 420,381 cells in atlas-eligible samples 356,406 passed (84.8%); the reference was then restricted to samples with
recoverable unique molecular identifier counts, giving **311,618 cells and nuclei from 57 samples and 52 patients**
measured over the 18,882 genes detected in every contributing study (Figure 1a). The cohort is dominated by one autopsy
series — Okada 2026 contributes 156,936 nuclei from 22 patients (50.4%) — followed by Jessa 2022 single-cell and
single-nucleus data (83,076 cells from 17 patients), the same study's multiome RNA (46,288 nuclei, 10 patients), Collot
2026 (17,562 nuclei, 4 patients), Zang (6,546 cells, 1 patient) and DeSisto (1,210 cells, 2 patients). By assay the
reference is 229,043 single-nucleus, 46,288 multiome RNA and 36,287 whole-cell profiles; by site, 42 pons and 15
thalamic samples; by histone variant, 8 H3.1, 1 H3.2, 23 H3.3 and 25 samples whose variant is not reported.

Integration with scVI (three seeds, 3,000 highly variable genes, sample as batch, study/assay/mitochondrial fraction as
covariates) was compared against Harmony and uncorrected principal components using scIB-style metrics. The pre-stated
selection score — 0.4 x batch mixing + 0.6 x biological conservation — chose scVI seed 1 at **0.739** (batch 0.537,
biology 0.874) over Harmony at 0.581 (0.736, 0.478) and principal components at 0.126 (0.127, 0.125) (Figure 1e,
Supplementary Figure 2). Harmony mixes batches more aggressively but at a cost in biological structure that the
weighting rejects. In the chosen embedding, malignant cells form a single large territory while myeloid,
oligodendrocyte, astrocyte, neuronal, vascular and T/NK compartments separate cleanly, and samples from all six studies
and all three assays intermingle within each compartment (Figure 1b–d).

### Malignant-cell calling agrees with published annotations in nine of ten cells

Malignant status was assigned per cell by combining lineage-marker scoring with an expression-based copy-number profile
computed against within-sample immune and oligodendrocyte references, smoothed over 100-gene windows and summarised
both chromosome-wide and per chromosome arm (Methods). We evaluated the calls against the original authors'
annotations wherever they exist: **324,676 annotated cells across 97 samples**, with agreement computed per sample so
that large cohorts cannot dominate. The median per-sample agreement is **0.898 (IQR 0.793–0.978)**; the unweighted mean
is 0.836 (95% CI 0.790–0.877) and the cell-weighted pooled value 0.884 (precision 0.934, recall 0.909), the gap
reflecting poorer performance in small Smart-seq2 plates and low-complexity nuclei samples (Supplementary Figure 1).
Thresholds were tuned on 22 samples; agreement on those samples (0.882) is indistinguishable from the 75 held-out ones
(0.885), so the calls are not an artefact of tuning. Calls supported by at least one recurrent arm-level event are more
reliable than score-only calls (precision 0.952 versus 0.862), and we therefore ship the event support as a confidence
field rather than a filter.

Recurrent events behave as expected for the disease: chromosome 1q gain is the most frequent, followed by 14q and 13q
loss, with chromosome 17p loss confined to a minority of samples (Figure 1f). We explicitly dropped 4q12 (PDGFRA/KIT)
from the hallmark set: at 100-gene resolution that band is covered by at most one window, and what it measures is
OPC-lineage expression against an immune reference rather than copy number — the apparent 4q12 "gain" frequencies in
the TPM-based cohorts were lineage expression, and gene-level deviations are released for exploratory use only.
Combining marker and copy-number evidence leaves 182,924 cells (53.5%) called malignant, 60,026 normal and 99,263
(29.0%) ambiguous; a permissive variant that trusts event-backed copy-number calls over discordant markers recovers a
further 30,463 cells and is released as a sensitivity label.

### Lineage labels transfer to held-out patients and to five external cohorts

A scANVI classifier trained on all labelled reference cells, with the validation split made by **patient** rather than
by cell, reaches 0.993 training and **0.975 held-out-patient accuracy** (0.917 balanced accuracy) over 23,827 labelled
cells from five patients never seen in training. Per-class recall is 0.999 malignant, 0.969 neuron, 0.950 myeloid,
0.917 oligodendrocyte, 0.896 T/NK, 0.883 vascular and 0.803 astrocyte; T/NK precision is the weakest figure at 0.561,
reflecting how few lymphocytes these tumours contain (Supplementary Figure 2).

We then projected five cohorts that were deliberately excluded from training, encoding each with the frozen reference
encoder and transferring labels with no fine-tuning (Figure 2a–b). Two are DMG cohorts deposited as transcripts per
million rather than counts — Filbin 2018 (3,037 cells, 9 samples) and Liu 2022 (12,438 cells, 42 plates) — and both
land predominantly in the malignant territory (median 90.3% and 83.6% of cells per sample); for Liu, cells the authors
called malignant are recovered as malignant in 99.3% of cases. LaBelle 2025, an immune-focused cohort, projects as
roughly half malignant and half immune, as its composition implies. The two normal-brain references provide the
specificity test: Siletti 2023 adult pons nuclei are called malignant in a median of only **1.2%** of cells per
dissection, with adult oligodendrocytes, astrocytes, microglia and neuronal superclusters recovered as their matching
compartments, whereas Braun 2023 first-trimester hindbrain projects as 91.3% malignant (Figure 2c–d).

That last number is a property of the reference, not an error: the reference contains no normal progenitor class, so
fetal radial glia, glioblasts and neuroblasts have nowhere else to go. Their nearest neighbours are informative
nonetheless — malignant cells sit closest to fetal glioblast (54.7%) and oligodendroglial (32.4%) classes, a 28.5-fold
enrichment for the oligodendroglial class over its frequency in the fetal reference — but the placement is made at a
median 15-nearest-neighbour distance of 3.28 in latent space versus 1.57 within the reference, so this is a
nearest-neighbour statement in an off-manifold region and not cell-of-origin evidence. We release the per-cell
neighbour distance and uncertainty percentile so that any downstream use can be filtered accordingly.

### Only two of fifteen expression metaprograms replicate across studies and assays

We derived expression programmes per sample by consensus non-negative matrix factorisation on 167,012 malignant cells
from 50 samples, then clustered the resulting programmes into metaprograms. At the primary dendrogram cut, 15
metaprograms absorb 161 programmes (Figure 3a, Supplementary Figure 4). Twelve carry a recognisable lineage identity —
four AC-like, two OPC-like, one OC-like, one MES-like, two cycling, one hypoxia and one ciliated programme — and three
are technical: ribosomal (MP2), mitochondrial (MP12) and an immediate-early-gene dissociation-stress programme (MP13).

Before deriving them we fixed a replication rule: a metaprogram counts as replicated only if it recurs in at least
three studies **and** in both whole cells and nuclei. Under that rule **2 of 15 replicate** — cycling (MP1, 34 samples)
and hypoxia (MP4, 15 samples) — and the count is stable across dendrogram cuts from 0.84 to 0.90, rising to 3 only at
0.91 and 5 at 0.92 (Supplementary Figure 4c). The explanation is structural rather than statistical: **every
lineage-identity programme in this cohort is nuclei-only**, and the single programme found only in whole cells is the
dissociation-stress programme. Masking the 5,921 genes that differ significantly between whole cells and nuclei in the two
paired multiome samples (plus 92 mitochondrial and ribosomal genes) and repeating the derivation recovers 12
metaprograms with 3 replicated, against a median of 2.5 metaprograms and 1 replicated over ten size-matched random
gene masks, so the assay dependence is specific to genes that differ between cells and nuclei rather than a generic
effect of removing genes (Figure 3b, Supplementary Figure 5).

Metaprogram families assigned per cell by within-sample standardised scores are distributed similarly in cells and
nuclei (within 0.03 per family), which is the correct comparison; the naive between-assay comparison is determined by
assay and is not interpretable (Figure 3c–d). Microenvironment subclustering resolves 21 populations — homeostatic
microglia, three transcriptional states of lipid-phagocytic and hypoxic/glycolytic tumour-associated macrophages,
interferon-responsive macrophages, CD4/CD8/NK lymphocyte states, five endothelial and pericyte populations and four
mature oligodendrocyte states plus a committed OPC population (Figure 4a). Composition modelling across samples finds
one effect that survives adjustment: neurons are more abundant in H3.1 than H3.3 tumours (q = 0.004; site-adjusted
q = 0.017; pons-only q = 0.050), whereas the astrocyte, myeloid, vascular and malignant-fraction differences do not
survive site adjustment (q = 0.17–0.27) or restriction to pons (q = 0.28–0.29) — all eight H3.1 samples are pons
nuclei samples, so these contrasts are confounded by design (Figure 4c). The strongest compositional signal in the
whole analysis is technical: neurons differ between nuclei and whole cells with a coefficient of 3.40 (q < 1e-9).

### Spatial co-localisation of malignant states with normal lineages

We deconvolved 27 Visium sections from four cohorts with cell2location against a single 13-class reference built from
the corrected metaprogram families (seven malignant families, six normal lineages), and tested co-localisation on spot
**proportions** with **patients** as the unit of inference: 21 DMG sections from 12 patients, Fisher-z averaged within
patient and pooled by random effects, with Benjamini–Hochberg correction across 42 state-by-microenvironment pairs
(Figure 5, Supplementary Figure 7).

Astrocyte-like tumour co-localises with normal astrocytes (rho = 0.393, 95% CI 0.234–0.531, q = 1e-4, positive in 13 of
21 sections) and with oligodendrocytes (rho = 0.308, q = 0.005). The ciliated programme behaves similarly, pairing with
astrocytes (rho = 0.369, q = 1e-4) and showing the most consistent positive sign of any pair with neurons (rho = 0.230,
q = 1e-4, positive in 16 of 21 sections). Of our two pre-stated hypotheses, the mesenchymal one holds — MES-like tumour
co-localises with myeloid cells (rho = 0.238, CI 0.071–0.392, q = 0.034, 14 of 21 sections) and with vasculature
(rho = 0.253, q = 0.034) — and the oligodendroglial one fails: OPC-like tumour and oligodendrocytes are not associated
after correction (rho = 0.178, CI -0.028 to 0.370, q = 0.210, 10 of 21 sections). Between-patient heterogeneity is high
for every pair (I-squared 0.952–0.997), so these are small, consistent-in-sign effects rather than strong ones, and
analysing abundances instead of proportions inflates them by roughly 30% because every class correlates with spot depth.

Ten serial sections through one whole brainstem autopsy allow a within-tumour gradient to be described, though with one
patient it remains descriptive: within sections, tumour fraction rises with the astrocyte-like (rho = 0.285, CI
0.149–0.410, positive in 9 of 10 sections) and OPC-like (rho = 0.213, q = 0.001) share and falls with the
oligodendrocyte (-0.298), astrocyte (-0.431), vascular (-0.235), neuronal (-0.239), mesenchymal (-0.385) and hypoxic
(-0.185) shares. Finally, five DMG and five glioblastoma Visium sections deposited in a single study permit a
like-for-like comparison: DMG sections carry a higher median tumour fraction (0.528 versus 0.460) and a lower cycling
share (0.112 versus 0.214), but with five sections per group neither survives correction (q = 0.667 for both), and we
report them as not established.

## Discussion

This resource makes three contributions. The first is the reference itself: 311,618 cells and nuclei from 52 patients
with harmonised labels, calibrated malignant calls, a frozen encoder for projection and a validated held-out-patient
classifier. For a tumour where individual cohorts contain a handful of patients, a patient-level denominator of 52 is
the difference between describing a sample and describing a disease.

The second is a set of calibrated negative results. Only two of fifteen metaprograms replicate under a rule we fixed in
advance, and the reason is not noise but confounding: lineage-identity programmes appear in nuclei and the
dissociation-stress programme appears in cells, so studies that profile different material will recover different
"states" and agree only on cycling and hypoxia. One of our two pre-stated spatial hypotheses survived patient-level
correction and the other did not. A compositional difference between H3.1 and H3.3 tumours is present for neurons and
disappears for every other cell type once site is adjusted, because H3.1 sampling is pons sampling in these deposits.
We regard these as the most immediately useful outputs for anyone planning a DMG single-cell study: they identify which
comparisons are currently decidable from public data and which require matched material to decide.

The third is methodological. Agreement between an expression-based copy-number call and the original authors'
annotations is not uniform — a median of 0.898 per sample with an interquartile range spanning nearly twenty
percentage points — and reporting it per sample rather than pooled is what exposes that. Event-backed and score-only
calls differ by nine points of precision and should not be treated alike. Zero-shot projection into a frozen reference
is accurate for tumour and adult-brain queries, and its failure mode on fetal tissue is interpretable and detectable
from the neighbour distances we release.

Several biological observations deserve follow-up in matched material. The ciliated programme is the clearest candidate:
it is the most consistently positive partner of normal neurons across sections and pairs with astrocytes as strongly as
the astrocyte-like programme does, yet it does not appear in any established DMG state nomenclature. The co-localisation
of mesenchymal tumour with myeloid cells and vasculature is consistent with the hypoxia-associated mesenchymal
compartment described in glioblastoma, and the brainstem gradient — tumour-rich regions enriched for astrocyte-like and
OPC-like programmes and depleted of normal glia — is what an infiltrative front is expected to look like, but a single
autopsy cannot establish it.

## Limitations

One autopsy series contributes half the reference cells, and every autopsy sample in the atlas comes from it, so
autopsy versus resection is not an estimable contrast. All eight H3.1 samples are pons nuclei samples. The whole-cell
arm of any assay comparison rests on six samples. Copy-number thresholds were tuned on 22 samples, and 18.9% of
malignant calls rest on the aggregate score with no supporting arm-level event. The reference contains no normal
progenitor population, so developmental queries project into the malignant region by construction. The spatial
conclusions come from 21 sections and 12 patients with high between-patient heterogeneity, the brainstem gradient from
one patient, and the DMG-versus-glioblastoma comparison from five sections per group. Two of the eleven cohorts are
deposited as transcripts per million with no recoverable counts and are used for projection only. Metaprogram derivation
was limited to k = 9 programmes per sample in 21 of 50 samples by cell number, so programmes finer than that ceiling
could not be resolved.

## Methods

**Cohort assembly.** Public deposits were screened by GEO E-utilities queries and a scholarly-literature search for
human H3K27M DMG single-cell and spatial data (Supplementary Table 1). Six cohorts with recoverable UMI counts form the
reference (GEO GSE210567 and GSE210566; Zenodo 10.5281/zenodo.17922773; GEO GSE231859; GEO GSE245674, sample GSM7846769
only; and the Okada 2026 post-mortem series, raw data under EGA accession EGAS50000001447 with the processed matrix
obtained from its UCSC Cell Browser instance). Five cohorts are projection-only: GEO GSE102130 and GSE184357 (both
transcripts per million), GSE227983, and the normal-brain references of Braun et al. (doi:10.1126/science.adf1226;
hindbrain, pons and medulla dissections, 286,973 cells) and Siletti et al. (doi:10.1126/science.add7046; six pons
dissections, 204,168 nuclei). Spatial cohorts are GSE268577, GSE194329, GSE280990 and the spatial archive of Zenodo
10.5281/zenodo.17922773. One cohort (ScPCA SCPCP000021) was excluded after quality control. All data are human; no
mouse data were used.

**Quality control.** Thresholds were applied per assay class, with distinct mitochondrial-fraction limits for whole
cells and nuclei, a minimum gene count per cell, and doublet scoring by Scrublet (doi:10.1016/j.cels.2018.11.005) at a
fixed threshold of 0.25. Samples failing a per-sample gate were flagged and excluded from the reference but retained in
the released object with their flag. Of 420,381 cells in atlas-eligible samples, 356,406 passed; the 311,618-cell
reference is the UMI-count subset of those.

**Malignant-cell calling.** For each sample, cells were scored against curated lineage marker sets, and an
inferCNV-style profile was computed from log-normalised expression relative to within-sample immune and oligodendrocyte
references over 100-gene windows ordered by GENCODE v32 coordinates, summarised both chromosome-wide and per arm.
Cluster-level calls combined an aggregate score threshold with recurrent arm-level events (1q, 2, 7 gain; 10q, 13q, 14q,
17p loss). A three-tier consensus — marker only, copy number only, and the agreement of both — was evaluated against
author annotations per sample; the released `consensus_call` requires marker and copy-number agreement, and
`consensus_call_eventtrust` is the permissive sensitivity variant.

**Integration and label transfer.** scVI (doi:10.1038/s41592-018-0229-2) was trained on 3,000 highly variable genes
(Seurat v3 selection, per-sample batch) with sample as batch key, study, assay and mitochondrial fraction as continuous
or categorical covariates, 30 latent dimensions and three random seeds, on an A10G GPU. Embeddings were compared with
scIB-style metrics (doi:10.1038/s41592-021-01336-8) and selected by the pre-stated 0.4/0.6 batch/biology score. scANVI
(doi:10.15252/msb.20209620) was initialised from the chosen scVI model and trained on all labelled cells with the
validation split stratified by patient. Query cohorts were encoded by the frozen reference encoder without
architectural surgery or fine-tuning, which the released run metadata records as zero training epochs; labels were
transferred both by the scANVI classifier and by 15-nearest-neighbour vote in latent space
(doi:10.1038/s41587-021-01001-7 for the general approach), and both are released per cell with a neighbour-distance
uncertainty percentile.

**Metaprograms.** Per-sample consensus non-negative matrix factorisation (doi:10.7554/eLife.43803) was run on malignant
cells from samples with at least 200 such cells, with k chosen per sample by stability and capped by cell number.
Programmes were clustered by Jaccard similarity of their top genes; the primary dendrogram cut and the pre-stated
replication rule (at least three studies and both assay classes) were fixed before inspection, and the sensitivity of
both the metaprogram count and the replication count to the cut was recorded. Per-cell family assignment used
within-sample standardised programme scores. The assay-masking control removed the 6,013 genes called assay-biased from the two
paired multiome samples (5,921 by a two-sample Welch test with consistent sign in both pairs, plus 92 mitochondrial
and ribosomal genes) and was compared against ten size-matched random masks.

**Spatial deconvolution and statistics.** cell2location (doi:10.1038/s41587-021-01139-4) reference signatures were
estimated from the 13-class metaprogram-family reference (250 epochs), and each section was deconvolved with a
wall-clock-budgeted schedule up to 20,000 epochs; realised epoch counts and full ELBO histories are released, and
convergence is argued from those histories rather than from reaching the target. Co-localisation was computed on spot
proportions (with a centred-log-ratio variant as sensitivity), Fisher-z averaged within patient, pooled by
DerSimonian–Laird random effects, and corrected across 42 state-by-microenvironment pairs. Serial sections of one
autopsy count as one patient, as do the core and edge sections of one resection.

**Audit.** The completed analysis was audited by six independent adversarial reviews with no access to the reasoning
that produced it, yielding 163 findings (4 critical, 33 major). Findings were remediated in four waves and the closure
of every critical and major finding was then verified by an independent re-audit that recomputed each quantity from the
released files: 78 of 90 findings in scope verified closed, including 4 of 4 critical. The complete record — findings,
closure claims, recomputed values and the remaining open items — is released with the atlas.

**Computing.** All analysis ran on a 14-core laptop with 24 GB of memory except scVI/scANVI training and cell2location
deconvolution, which ran on a single A10G GPU. Software versions are recorded for both environments.

## Data and code availability

The integrated object (`dmg_atlas_v1.h5ad`, CELLxGENE schema 7.x, validated with zero MUST failures) carries raw counts,
harmonised metadata, lineage and state labels, malignant calls with confidence fields, and all model embeddings.
*Release of the merged count matrix is pending confirmation of redistribution terms for the post-mortem cohort, whose
raw data is under controlled access (EGAS50000001447); if those terms preclude redistribution, the release will ship
embeddings, labels and metadata with per-cohort pointers to the primary deposits.* All primary accessions are listed in
Methods. Analysis code, documentation and the complete audit trail are in the project repository; every figure panel
records the md5 checksum of each input table.

## Author contributions

[To be completed]

## Competing interests

[To be completed]

## Figure legends

**Figure 1 | An integrated reference of human H3K27M diffuse midline glioma.** (a) Cells and nuclei per contributing
study, coloured by assay. (b–d) UMAP of the selected scVI latent space (311,618 cells, 57 samples, 52 patients),
coloured by harmonised lineage label (b), assay (c) and study (d). (e) scIB-style integration metrics for three scVI
seeds, a 60-epoch scVI variant, Harmony and uncorrected principal components, separated into batch-mixing, biological
conservation and overall scores; the asterisk marks the model selected by the pre-stated 0.4/0.6 score. (f) Fraction of
copy-number-malignant cells carrying each recurrent arm-level event, by study, including the two projection-only
cohorts.

**Figure 2 | Zero-shot projection of five external cohorts.** (a) Predicted lineage composition of each query cohort.
(b) Per-sample median neighbour-based uncertainty, expressed as a percentile of the reference's own
15-nearest-neighbour distance distribution. (c) Author-annotated state against predicted lineage for Liu 2022
(row-normalised). (d) Fetal (Braun 2023) and adult (Siletti 2023) reference classes against predicted lineage; the
horizontal rule separates the two references.

**Figure 3 | Metaprograms and their assay dependence.** (a) Metaprogram by study-and-assay occurrence matrix; replicated
metaprograms are marked, technical programmes are labelled by the class they capture. (b) Number of replicated
metaprograms as a function of dendrogram cut, with and without masking the 6,013 assay-biased genes; the dashed
line marks the primary cut. (c) Jaccard similarity between matched cell-derived and nucleus-derived metaprograms, by
family, before and after assay masking. (d) Per-sample metaprogram-family composition of malignant cells, ordered by
study.

**Figure 4 | Microenvironment populations and sample-level composition.** (a) Microenvironment subclusters: cell count
(log scale) and the fraction of each cluster derived from nuclei. (b) Assay purity of microenvironment clusters in the
atlas embedding versus a per-compartment model. (c) Composition contrasts as the change in logit proportion per cell
type and malignant state; asterisks mark Benjamini–Hochberg significance, and the nuclei-versus-cells column is the
technical reference contrast.

**Figure 5 | Spatial organisation of malignant states.** (a–c) One representative section (pixel coordinates, 1 mm
scale bar) showing tumour fraction, the mesenchymal share of tumour and the myeloid fraction. (d) Pooled spot-level
Spearman correlations between malignant states and microenvironment lineages over 12 DMG patients and 21 sections;
filled points pass Benjamini–Hochberg correction across 42 pairs. (e) Mean spot proportions along ten serial sections of
one brainstem autopsy. (f) Tumour fraction and cycling share per section in five DMG and five glioblastoma sections
from one study.

**Supplementary Figures 1–7** cover malignant-call validation, integration metrics and classifier performance, mapping
and nearest-fetal-class detail, metaprogram gene content and robustness, the assay-masking control, microenvironment
markers, and spatial niche and sensitivity analyses.

## References

1. Louis DN, et al. The 2021 WHO Classification of Tumors of the Central Nervous System: a summary. *Neuro-Oncology* (2021). doi:10.1093/neuonc/noab106
2. Mackay A, et al. Integrated molecular meta-analysis of 1,000 pediatric high-grade and diffuse intrinsic pontine glioma. *Cancer Cell* (2017). doi:10.1016/j.ccell.2017.08.017
3. Filbin MG, et al. Developmental and oncogenic programs in H3K27M gliomas dissected by single-cell RNA-seq. *Science* (2018). doi:10.1126/science.aao4750
4. Jessa S, et al. K27M in canonical and noncanonical H3 variants occurs in distinct oligodendroglial cell lineages in brain midline gliomas. *Nature Genetics* (2022). doi:10.1038/s41588-022-01205-w
5. Patel AP, et al. Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. *Science* (2014). doi:10.1126/science.1254257
6. Neftel C, et al. An integrative model of cellular states, plasticity, and genetics for glioblastoma. *Cell* (2019). doi:10.1016/j.cell.2019.06.024
7. Gavish A, et al. Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours. *Nature* (2023). doi:10.1038/s41586-023-06130-4
8. Lopez R, et al. Deep generative modeling for single-cell transcriptomics. *Nature Methods* (2018). doi:10.1038/s41592-018-0229-2
9. Xu C, et al. Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. *Molecular Systems Biology* (2021). doi:10.15252/msb.20209620
10. Lotfollahi M, et al. Mapping single-cell data to reference atlases by transfer learning. *Nature Biotechnology* (2021). doi:10.1038/s41587-021-01001-7
11. Luecken MD, et al. Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods* (2021). doi:10.1038/s41592-021-01336-8
12. Kotliar D, et al. Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq. *eLife* (2019). doi:10.7554/eLife.43803
13. Kleshchevnikov V, et al. Cell2location maps fine-grained cell types in spatial transcriptomics. *Nature Biotechnology* (2022). doi:10.1038/s41587-021-01139-4
14. Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Systems* (2019). doi:10.1016/j.cels.2018.11.005
15. Braun E, et al. Comprehensive cell atlas of the first-trimester developing human brain. *Science* (2023). doi:10.1126/science.adf1226
16. Siletti K, et al. Transcriptomic diversity of cell types across the adult human brain. *Science* (2023). doi:10.1126/science.add7046

*Data citations.* GEO: GSE210567, GSE210566, GSE102130, GSE184357, GSE231859, GSE245674, GSE227983, GSE268577,
GSE194329, GSE280990. Zenodo: 10.5281/zenodo.17922773. EGA: EGAS50000001447. Primary publications for the Collot,
Okada, DeSisto, Zang, LaBelle, Huchede, Ren and whole-brainstem cohorts are cited in Supplementary Table 1, which
records for each deposit the accession, access status, assay, sample and cell counts, and the role it plays here.
