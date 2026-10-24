# Phase 2 QC report — DMG atlas

Generated 2026-10-19 17:56 by scripts/04_qc_report.py from data/interim/*.qc.h5ad. Thresholds: configs/qc.yaml (nuclei: pct_mito < 5 %, cells: < 20 %, min genes 300 for 10x / 1000 for Smart-seq2, per-sample max genes = 99.5th percentile; scrublet expected rate 0.06 on 10x samples, flagged only). Cells are flagged (`obs['qc_pass']`), none removed. Doublets are not removed.

## Per-study totals

| study | set | samples | cells before | cells after | % retained | median genes (pooled) | median counts | median % mito | doublet rate (auto thr) | doublet rate (score>0.25) | genes in var | layers | file GB | wall min | peak RSS GB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| collot | atlas | 5 | 18,767 | 17,750 | 94.6 | 1387 | 2332 | 0.04 | 0.001 | 0.011 | 19,248 | counts | 0.45 | 0.9 | 3.17 |
| desisto | atlas | 2 | 1,829 | 1,211 | 66.2 | 3604 | 13825 | 9.75 | 0.001 | 0.001 | 32,734 | counts | 0.09 | 0.5 | 3.61 |
| filbin | atlas | 9 | 3,057 | 3,037 | 99.3 | 5414 | n/a | n/a |  |  | 19,664 | tpm | 0.25 | 0.1 | 2.54 |
| jessa_multiome | atlas | 10 | 86,978 | 55,064 | 63.3 | 2742 | 6053 | 1.13 | 0.021 | 0.008 | 34,256 | counts | 3.5 | 3.2 | 8.95 |
| jessa_scrna | atlas | 17 | 104,735 | 90,596 | 86.5 | 2063 | 4305 | 1.08 | 0.080 | 0.016 | 30,660 | counts | 3.89 | 6.2 | 5.91 |
| labelle | atlas | 9 | 8,347 | 8,247 | 98.8 | 793 | n/a | n/a |  |  | 19,114 | log2_tpm10 | 0.2 | 0.3 | 1.20 |
| liu | atlas | 42 | 16,676 | 12,438 | 74.6 | 4805 | 150834 | 1.13 |  |  | 32,632 | counts | 1.07 | 0.5 | 3.85 |
| okada | atlas | 22 | 168,597 | 161,043 | 95.5 | 1665 | 2949 | 0.19 | 0.071 | 0.026 | 36,591 | counts | 4.99 | 5.8 | 8.45 |
| scpca | atlas | 2 | 3,985 | 441 | 11.1 | 362 | 436 | 2.22 | 0.057 | 0.011 | 36,457 | counts | 0.02 | 0.0 | 0.66 |
| zang | atlas | 1 | 7,410 | 6,579 | 88.8 | 2536 | 6970 | 7.22 | 0.000 | 0.005 | 36,591 | counts | 0.28 | 0.4 | 2.20 |
| desisto | out-group | 17 | 26,264 | 14,320 | 54.5 | 1884 | 4862 | 9.86 | 0.001 | 0.006 | 32,734 | counts | 0.77 | 0.5 | 3.61 |
| jessa_scrna | out-group | 20 | 101,685 | 85,040 | 83.6 | 2056 | 4008 | 1.49 | 0.001 | 0.008 | 30,660 | counts | 3.58 | 6.2 | 5.91 |

Atlas (included) samples: **119 samples, 420,381 cells before QC, 356,406 after QC (84.8 %)**. Out-group files: 127,949 cells before / 99,360 after.

## Sample-level flags (audit corrections 2026-10-19)

`sample_qc_pass` (flag only, nothing dropped): median detected genes of the sample's qc_pass cells >= 500 (configs/qc.yaml `sample_gate`). **3 sample(s) flagged:**

- scpca__SCPCS000637 (scpca, atlas): 441 passing cells, median 362 genes / 436 counts
- okada__1734_pons (okada, atlas): 8,486 passing cells, median 427 genes / 468 counts
- scpca__SCPCS000636 (scpca, atlas): 0 passing cells, median nan genes / nan counts

`mt_genes_present` = the sample's feature annotation carries MT- genes with signal; where False the mitochondrial filter is inoperative. **34 sample(s) without MT- genes** (filbin, labelle, liu): filbin__BCH1126, filbin__BCH836, filbin__BCH869, filbin__MGH101, filbin__MGH104, filbin__MGH66, filbin__MUV1, filbin__MUV10, filbin__MUV5, liu__AAA007653, liu__BT1179, liu__ICGCGBM27, liu__ICGCGBM60, liu__ICGCGBM96, liu__MUV16, liu__MUV82, liu__R1E1E3G3, liu__S1D2E3H7, liu__SUDIPG31, liu__SUDIPG55, liu__SUDIPG65, liu__XT8698, liu__XT9824, liu__pTG4, liu__pTG5, labelle__H3K27M_P1, labelle__H3K27M_P2, labelle__H3K27M_P3, labelle__H3K27M_P4, labelle__H3K27M_P5, labelle__H3K27M_P6, labelle__H3K27M_P7, labelle__H3K27M_P8, labelle__H3K27M_P9

Author-prep-aware mito rule (nuclei threshold 5.0 % applied to Smart-seq2 frozen-nuclei samples with MT- genes): 8 samples: liu__MGH6666410, liu__SUDIPG30, liu__SUDIPG45, liu__SUDIPG50Aut, liu__SUDIPG50Bio, liu__SUDIPG54, liu__SUDIPGIII, liu__pSCG1

Value scales: UMI studies and Liu = counts; Filbin = linear TPM (`layers['tpm']`, per-cell sum 1e6); LaBelle = log2(TPM/10+1) as deposited (`layers['log2_tpm10']`, `uns['value_scale']`). 'median counts' is therefore not meaningful for Filbin/LaBelle and is shown as n/a.

## Sanity checks

**(a) Jessa author-labelled malignant cells passing our QC** (`results/qc/sanity_jessa_malignant_qc.csv`):

| file | author malignant | pass QC | fraction | author normal pass fraction | author-QC cells passing ours | barcodes not in author metadata | of which pass ours | malignant flagged doublet |
|---|---|---|---|---|---|---|---|---|
| jessa_scrna.qc.h5ad | 71,899 | 69,593 | 0.968 | 0.967 | 0.968 | 18,210 | 0.377 | 0.088 |
| jessa_scrna.outgroup.qc.h5ad | 71,712 | 62,611 | 0.873 | 0.963 | 0.892 | 11,298 | 0.394 | 0.000 |
| jessa_multiome.qc.h5ad | 35,122 | 33,975 | 0.967 | 0.852 | 0.938 | 39,652 | 0.270 | 0.015 |

**(b) Assay label vs median pct_mito before QC** (UMI studies; TPM tables of Filbin/LaBelle carry no MT- genes so pct_mito is 0 there and not informative):

```
                count   min   25%    50%    75%    max
qc_assay_class                                        
cells            31.0  2.69  9.67  13.46  19.03  44.87
nuclei           65.0  0.00  0.13   0.48   1.50  47.69
smartseq2        26.0  0.21  1.77   2.26   2.81  15.56
```

Samples whose median pct_mito contradicts the manifest assay label (6):

- jessa_scrna__P-6254_S-8499 labelled snRNA but median pct_mito before QC = 47.69% (>5%)
- jessa_scrna__P-3200_S-3867 labelled snRNA but median pct_mito before QC = 11.48% (>5%)
- jessa_scrna__P-5613_S-7162 labelled snRNA but median pct_mito before QC = 5.43% (>5%)
- jessa_multiome__P-6253_S-8498 labelled multiome_RNA but median pct_mito before QC = 11.15% (>5%)
- jessa_multiome__P-1709_S-1709 labelled multiome_RNA but median pct_mito before QC = 5.85% (>5%)
- jessa_multiome__P-1779_S-1781 labelled multiome_RNA but median pct_mito before QC = 15.24% (>5%)

**(c) Gene overlap:** 18,882 harmonised symbols are measured in every sample of all 7 UMI studies (jessa_scrna, jessa_multiome, collot, desisto, zang, scpca, okada); 18,301 are shared by all 10 studies including the Smart-seq2/TPM cohorts. Gene universe = GRCh38-2020-A symbols (36,591 unique). Per-study losses: `results/qc/gene_harmonisation.csv`.

## Gene harmonisation per study

| study | input genes | mapped by Ensembl | exact symbol | alias | ambiguous dropped | unmapped dropped | duplicates collapsed | output genes |
|---|---|---|---|---|---|---|---|---|
| collot | 19,248 | 0 | 19,248 | 0 | 0 | 0 | 0 | 19,248 |
| desisto | 33,538 | 32,753 | 0 | 0 | 0 | 785 | 19 | 32,734 |
| filbin | 23,686 | 0 | 18,700 | 979 | 0 | 4,007 | 15 | 19,664 |
| jessa_multiome | 60,658 | 35,088 | 0 | 0 | 0 | 25,570 | 832 | 34,256 |
| jessa_scrna | 32,738 | 30,847 | 0 | 0 | 0 | 1,891 | 187 | 30,660 |
| labelle | 19,378 | 0 | 18,224 | 903 | 0 | 251 | 13 | 19,114 |
| liu | 23,686 | 0 | 18,700 | 979 | 0 | 4,007 | 15 | 19,664 |
| okada | 36,601 | 0 | 36,591 | 0 | 0 | 10 | 0 | 36,591 |
| scpca | 60,319 | 36,665 | 0 | 0 | 0 | 23,654 | 208 | 36,457 |
| zang | 36,601 | 36,601 | 0 | 0 | 0 | 0 | 10 | 36,591 |

## Per-sample QC (`results/qc/qc_per_sample.csv`)

| sample | study | assay | set | before | after | median genes | median counts | median % mito | doublet rate (auto thr) | doublet rate (score>0.25) | fail low genes | fail high genes | fail mito | sample gate | MT genes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| collot__DMG-PMC-9_core | collot | snRNA | atlas | 2,941 | 2,672 | 2012 | 4128 | 0.04 | 0.000 | 0.002 | 91 | 15 | 187 |  | yes |
| collot__DMG-PMC-9_edge | collot | snRNA | atlas | 2,421 | 2,198 | 2325 | 5029 | 0.05 | 0.000 | 0.003 | 68 | 13 | 160 |  | yes |
| collot__PMC-DMG-6 | collot | snRNA | atlas | 2,455 | 2,278 | 1689 | 3090 | 0.33 | 0.000 | 0.020 | 84 | 13 | 118 |  | yes |
| collot__PMC-DMG-7 | collot | snRNA | atlas | 5,824 | 5,623 | 1153 | 1797 | 0.00 | 0.001 | 0.010 | 171 | 30 | 0 |  | yes |
| collot__PMC-DMG-8 | collot | snRNA | atlas | 5,126 | 4,979 | 1068 | 1637 | 0.05 | 0.001 | 0.015 | 121 | 26 | 0 |  | yes |
| desisto__1196 | desisto | scRNA | atlas | 1,433 | 914 | 4089 | 16256 | 9.94 | 0.000 | 0.001 | 52 | 8 | 510 |  | yes |
| desisto__1437 | desisto | scRNA | atlas | 396 | 297 | 2082 | 6133 | 8.93 | 0.003 | 0.000 | 37 | 2 | 92 |  | yes |
| filbin__BCH1126 | filbin | smartseq2 | atlas | 299 | 297 | 5125 | n/a | n/a |  |  | 0 | 2 | 0 |  | no |
| filbin__BCH836 | filbin | smartseq2 | atlas | 527 | 524 | 6644 | n/a | n/a |  |  | 0 | 3 | 0 |  | no |
| filbin__BCH869 | filbin | smartseq2 | atlas | 492 | 489 | 4942 | n/a | n/a |  |  | 0 | 3 | 0 |  | no |
| filbin__MGH101 | filbin | smartseq2 | atlas | 92 | 91 | 4343 | n/a | n/a |  |  | 0 | 1 | 0 |  | no |
| filbin__MGH104 | filbin | smartseq2 | atlas | 65 | 64 | 7108 | n/a | n/a |  |  | 0 | 1 | 0 |  | no |
| filbin__MGH66 | filbin | smartseq2 | atlas | 442 | 439 | 6635 | n/a | n/a |  |  | 0 | 3 | 0 |  | no |
| filbin__MUV1 | filbin | smartseq2 | atlas | 146 | 145 | 5103 | n/a | n/a |  |  | 0 | 1 | 0 |  | no |
| filbin__MUV10 | filbin | smartseq2 | atlas | 286 | 284 | 4472 | n/a | n/a |  |  | 0 | 2 | 0 |  | no |
| filbin__MUV5 | filbin | smartseq2 | atlas | 708 | 704 | 4806 | n/a | n/a |  |  | 0 | 4 | 0 |  | no |
| jessa_multiome__P-1694_S-1694 | jessa_multiome | multiome_RNA | atlas | 5,530 | 4,990 | 3260 | 8678 | 1.02 | 0.000 | 0.017 | 151 | 28 | 514 |  | yes |
| jessa_multiome__P-1701_S-1701 | jessa_multiome | multiome_RNA | atlas | 12,386 | 10,297 | 2861 | 7498 | 0.90 | 0.000 | 0.017 | 877 | 62 | 1907 |  | yes |
| jessa_multiome__P-1709_S-1709 | jessa_multiome | multiome_RNA | atlas | 20,000 | 8,037 | 2755 | 5434 | 3.44 | 0.000 | 0.008 | 7 | 100 | 11926 |  | yes |
| jessa_multiome__P-1764_S-1766 | jessa_multiome | multiome_RNA | atlas | 6,542 | 3,742 | 2763 | 6444 | 1.74 | 0.001 | 0.002 | 938 | 33 | 2717 |  | yes |
| jessa_multiome__P-1779_S-1781 | jessa_multiome | multiome_RNA | atlas | 10,181 | 2,091 | 2135 | 3957 | 2.52 | 0.000 | 0.005 | 3445 | 51 | 8072 |  | yes |
| jessa_multiome__P-1780_S-1780 | jessa_multiome | multiome_RNA | atlas | 8,616 | 8,060 | 2767 | 6134 | 1.00 | 0.144 | 0.003 | 75 | 44 | 486 |  | yes |
| jessa_multiome__P-6253_S-8498 | jessa_multiome | multiome_RNA | atlas | 5,778 | 1,929 | 2173 | 4626 | 1.89 | 0.000 | 0.012 | 2506 | 29 | 3768 |  | yes |
| jessa_multiome__P-6337_S-8821 | jessa_multiome | multiome_RNA | atlas | 5,904 | 4,211 | 2372 | 4423 | 0.50 | 0.000 | 0.000 | 1227 | 30 | 1313 |  | yes |
| jessa_multiome__P-6640_S-9581 | jessa_multiome | multiome_RNA | atlas | 6,789 | 6,587 | 2377 | 3690 | 0.92 | 0.000 | 0.003 | 1 | 34 | 168 |  | yes |
| jessa_multiome__P-6774_S-10146 | jessa_multiome | multiome_RNA | atlas | 5,252 | 5,120 | 2616 | 6690 | 0.19 | 0.000 | 0.007 | 78 | 27 | 28 |  | yes |
| jessa_scrna__P-1713_S-1713 | jessa_scrna | snRNA | atlas | 7,021 | 4,852 | 3423 | 9667 | 1.51 | 0.000 | 0.011 | 26 | 36 | 2139 |  | yes |
| jessa_scrna__P-1741_S-2756 | jessa_scrna | snRNA | atlas | 1,626 | 1,328 | 2414 | 4756 | 1.72 | 0.006 | 0.012 | 6 | 9 | 297 |  | yes |
| jessa_scrna__P-1764_S-1766 | jessa_scrna | snRNA | atlas | 5,561 | 5,024 | 1149 | 1573 | 0.93 | 0.000 | 0.009 | 6 | 28 | 509 |  | yes |
| jessa_scrna__P-1775_S-1775 | jessa_scrna | snRNA | atlas | 6,375 | 6,323 | 1807 | 2485 | 0.31 | 0.000 | 0.016 | 0 | 32 | 20 |  | yes |
| jessa_scrna__P-2687_S-2688 | jessa_scrna | snRNA | atlas | 1,116 | 1,060 | 3189 | 8065 | 0.08 | 0.001 | 0.001 | 15 | 6 | 52 |  | yes |
| jessa_scrna__P-3407_S-3447 | jessa_scrna | snRNA | atlas | 6,100 | 6,067 | 1226 | 2049 | 0.09 | 0.000 | 0.028 | 0 | 31 | 2 |  | yes |
| jessa_scrna__P-4198_S-4459 | jessa_scrna | snRNA | atlas | 4,541 | 4,514 | 822 | 1123 | 0.00 | 0.000 | 0.018 | 0 | 23 | 4 |  | yes |
| jessa_scrna__P-4504_S-4916 | jessa_scrna | snRNA | atlas | 5,876 | 5,840 | 1157 | 1613 | 0.09 | 0.000 | 0.027 | 0 | 30 | 7 |  | yes |
| jessa_scrna__P-6117_S-8370 | jessa_scrna | scRNA | atlas | 8,507 | 7,879 | 2480 | 7009 | 10.82 | 0.000 | 0.015 | 0 | 43 | 585 |  | yes |
| jessa_scrna__P-6251_S-8496 | jessa_scrna | snRNA | atlas | 2,731 | 2,659 | 3976 | 12041 | 0.14 | 0.001 | 0.003 | 1 | 14 | 57 |  | yes |
| jessa_scrna__P-6253_S-8498 | jessa_scrna | snRNA | atlas | 570 | 296 | 2490 | 5271 | 2.34 | 0.000 | 0.000 | 4 | 3 | 272 |  | yes |
| jessa_scrna__P-6254_S-8499 | jessa_scrna | snRNA | atlas | 2,944 | 229 | 1182 | 1883 | 3.02 | 0.000 | 0.000 | 1003 | 15 | 2710 |  | yes |
| jessa_scrna__P-6255_S-8500 | jessa_scrna | snRNA | atlas | 1,064 | 896 | 4384 | 12365 | 0.91 | 0.000 | 0.007 | 1 | 6 | 162 |  | yes |
| jessa_scrna__P-6328_S-8672 | jessa_scrna | snRNA | atlas | 18,284 | 18,174 | 1969 | 3433 | 0.77 | 0.365 | 0.005 | 0 | 92 | 19 |  | yes |
| jessa_scrna__P-6337_S-8821 | jessa_scrna | scRNA | atlas | 11,823 | 7,897 | 2474 | 6661 | 10.74 | 0.019 | 0.021 | 348 | 60 | 3852 |  | yes |
| jessa_scrna__P-6519_S-9084 | jessa_scrna | scRNA | atlas | 10,626 | 9,014 | 2385 | 6367 | 7.84 | 0.024 | 0.026 | 195 | 54 | 1545 |  | yes |
| jessa_scrna__P-6640_S-9581 | jessa_scrna | scRNA | atlas | 9,970 | 8,544 | 3170 | 11860 | 7.20 | 0.022 | 0.023 | 190 | 50 | 1358 |  | yes |
| labelle__H3K27M_P1 | labelle | smartseq2 | atlas | 208 | 205 | 4185 | n/a | n/a |  |  | 1 | 2 | 0 |  | no |
| labelle__H3K27M_P2 | labelle | smartseq2 | atlas | 356 | 351 | 3544 | n/a | n/a |  |  | 3 | 2 | 0 |  | no |
| labelle__H3K27M_P3 | labelle | smartseq2 | atlas | 458 | 435 | 2404 | n/a | n/a |  |  | 20 | 3 | 0 |  | no |
| labelle__H3K27M_P4 | labelle | smartseq2 | atlas | 299 | 290 | 4626 | n/a | n/a |  |  | 7 | 2 | 0 |  | no |
| labelle__H3K27M_P5 | labelle | smartseq2 | atlas | 219 | 217 | 3668 | n/a | n/a |  |  | 0 | 2 | 0 |  | no |
| labelle__H3K27M_P6 | labelle | smartseq2 | atlas | 192 | 190 | 3214 | n/a | n/a |  |  | 1 | 1 | 0 |  | no |
| labelle__H3K27M_P7 | labelle | smartseq2 | atlas | 83 | 65 | 3215 | n/a | n/a |  |  | 17 | 1 | 0 |  | no |
| labelle__H3K27M_P8 | labelle | smartseq2 | atlas | 134 | 128 | 2646 | n/a | n/a |  |  | 5 | 1 | 0 |  | no |
| labelle__H3K27M_P9 | labelle | scRNA | atlas | 6,398 | 6,366 | 709 | n/a | n/a |  |  | 0 | 32 | 0 |  | no |
| liu__AAA007653 | liu | smartseq2 | atlas | 288 | 269 | 4563 | 85510 | n/a |  |  | 17 | 2 | 0 |  | no |
| liu__BT1126 | liu | smartseq2 | atlas | 384 | 364 | 6111 | 445669 | 3.04 |  |  | 18 | 2 | 1 |  | yes |
| liu__BT1179 | liu | smartseq2 | atlas | 384 | 351 | 2835 | 65786 | n/a |  |  | 31 | 2 | 0 |  | no |
| liu__BT1462 | liu | smartseq2 | atlas | 383 | 351 | 7137 | 699473 | 2.58 |  |  | 29 | 2 | 3 |  | yes |
| liu__BT836 | liu | smartseq2 | atlas | 864 | 582 | 7748 | 326909 | 2.03 |  |  | 277 | 5 | 17 |  | yes |
| liu__BT869 | liu | smartseq2 | atlas | 960 | 907 | 5001 | 155159 | 1.88 |  |  | 48 | 5 | 0 |  | yes |
| liu__ICGCGBM27 | liu | smartseq2 | atlas | 384 | 300 | 2809 | 75080 | n/a |  |  | 82 | 2 | 0 |  | no |
| liu__ICGCGBM60 | liu | smartseq2 | atlas | 384 | 324 | 2836 | 95101 | n/a |  |  | 58 | 2 | 0 |  | no |
| liu__ICGCGBM96 | liu | smartseq2 | atlas | 288 | 259 | 2245 | 68767 | n/a |  |  | 27 | 2 | 0 |  | no |
| liu__MGH6666410 | liu | smartseq2 | atlas | 384 | 220 | 4909 | 325851 | 1.76 |  |  | 62 | 2 | 133 |  | yes |
| liu__MUV1 | liu | smartseq2 | atlas | 576 | 195 | 6055 | 477791 | 1.78 |  |  | 377 | 3 | 25 |  | yes |
| liu__MUV10 | liu | smartseq2 | atlas | 960 | 355 | 4896 | 438240 | 2.37 |  |  | 600 | 5 | 2 |  | yes |
| liu__MUV16 | liu | smartseq2 | atlas | 384 | 251 | 2895 | 106141 | n/a |  |  | 131 | 2 | 0 |  | no |
| liu__MUV17 | liu | smartseq2 | atlas | 384 | 314 | 5914 | 207873 | 2.79 |  |  | 68 | 2 | 0 |  | yes |
| liu__MUV31 | liu | smartseq2 | atlas | 288 | 273 | 8550 | 539270 | 2.37 |  |  | 13 | 2 | 2 |  | yes |
| liu__MUV32 | liu | smartseq2 | atlas | 284 | 262 | 7168 | 556210 | 2.21 |  |  | 18 | 2 | 5 |  | yes |
| liu__MUV35 | liu | smartseq2 | atlas | 192 | 186 | 7318 | 376229 | 1.79 |  |  | 5 | 1 | 1 |  | yes |
| liu__MUV5 | liu | smartseq2 | atlas | 1,056 | 787 | 5173 | 401373 | 2.49 |  |  | 262 | 6 | 18 |  | yes |
| liu__MUV54 | liu | smartseq2 | atlas | 384 | 371 | 6472 | 118276 | 2.08 |  |  | 11 | 2 | 0 |  | yes |
| liu__MUV77 | liu | smartseq2 | atlas | 288 | 255 | 6388 | 160599 | 2.14 |  |  | 31 | 2 | 7 |  | yes |
| liu__MUV78 | liu | smartseq2 | atlas | 192 | 185 | 8343 | 217724 | 3.35 |  |  | 6 | 1 | 0 |  | yes |
| liu__MUV82 | liu | smartseq2 | atlas | 384 | 359 | 2871 | 91935 | n/a |  |  | 23 | 2 | 0 |  | no |
| liu__MUV83 | liu | smartseq2 | atlas | 192 | 184 | 5964 | 134815 | 1.99 |  |  | 7 | 1 | 0 |  | yes |
| liu__MUV86 | liu | smartseq2 | atlas | 192 | 180 | 7163 | 137842 | 2.71 |  |  | 11 | 1 | 0 |  | yes |
| liu__MUV87 | liu | smartseq2 | atlas | 192 | 173 | 7675 | 170055 | 2.40 |  |  | 18 | 1 | 0 |  | yes |
| liu__MUV91 | liu | smartseq2 | atlas | 384 | 375 | 7196 | 171012 | 1.73 |  |  | 7 | 2 | 4 |  | yes |
| liu__R1E1E3G3 | liu | smartseq2 | atlas | 384 | 356 | 3228 | 82134 | n/a |  |  | 26 | 2 | 0 |  | no |
| liu__S1D2E3H7 | liu | smartseq2 | atlas | 96 | 23 | 1625 | 17410 | n/a |  |  | 72 | 1 | 0 |  | no |
| liu__SUDIPG30 | liu | smartseq2 | atlas | 384 | 120 | 6698 | 296430 | 2.17 |  |  | 55 | 2 | 248 |  | yes |
| liu__SUDIPG31 | liu | smartseq2 | atlas | 384 | 247 | 1560 | 24129 | n/a |  |  | 135 | 2 | 0 |  | no |
| liu__SUDIPG45 | liu | smartseq2 | atlas | 480 | 82 | 7168 | 1029093 | 2.96 |  |  | 37 | 3 | 394 |  | yes |
| liu__SUDIPG50Aut | liu | smartseq2 | atlas | 384 | 330 | 7162 | 501434 | 0.98 |  |  | 5 | 2 | 47 |  | yes |
| liu__SUDIPG50Bio | liu | smartseq2 | atlas | 384 | 371 | 6060 | 350932 | 0.21 |  |  | 1 | 2 | 10 |  | yes |
| liu__SUDIPG54 | liu | smartseq2 | atlas | 384 | 254 | 9376 | 573406 | 0.91 |  |  | 0 | 2 | 128 |  | yes |
| liu__SUDIPG55 | liu | smartseq2 | atlas | 95 | 61 | 2290 | 82538 | n/a |  |  | 33 | 1 | 0 |  | no |
| liu__SUDIPG65 | liu | smartseq2 | atlas | 384 | 235 | 2406 | 59348 | n/a |  |  | 147 | 2 | 0 |  | no |
| liu__SUDIPGIII | liu | smartseq2 | atlas | 384 | 177 | 3729 | 270511 | 2.01 |  |  | 71 | 2 | 166 |  | yes |
| liu__XT8698 | liu | smartseq2 | atlas | 384 | 329 | 2028 | 45255 | n/a |  |  | 53 | 2 | 0 |  | no |
| liu__XT9824 | liu | smartseq2 | atlas | 384 | 365 | 2704 | 75509 | n/a |  |  | 17 | 2 | 0 |  | no |
| liu__pSCG1 | liu | smartseq2 | atlas | 362 | 189 | 7782 | 543975 | 0.32 |  |  | 25 | 2 | 171 |  | yes |
| liu__pTG4 | liu | smartseq2 | atlas | 384 | 311 | 1850 | 23650 | n/a |  |  | 71 | 2 | 0 |  | no |
| liu__pTG5 | liu | smartseq2 | atlas | 384 | 356 | 2952 | 91803 | n/a |  |  | 26 | 2 | 0 |  | no |
| okada__1006_thalamus | okada | snRNA | atlas | 1,620 | 1,611 | 2018 | 4214 | 1.63 | 0.003 | 0.010 | 0 | 9 | 0 |  | yes |
| okada__1277_pons | okada | snRNA | atlas | 5,716 | 5,675 | 1768 | 3135 | 0.22 | 0.000 | 0.016 | 12 | 29 | 0 |  | yes |
| okada__1368_pons | okada | snRNA | atlas | 3,057 | 3,039 | 1241 | 1989 | 0.03 | 0.001 | 0.009 | 2 | 16 | 0 |  | yes |
| okada__1428_pons | okada | snRNA | atlas | 6,748 | 6,715 | 1835 | 3707 | 0.89 | 0.044 | 0.034 | 1 | 32 | 0 |  | yes |
| okada__1549_thalamus | okada | snRNA | atlas | 3,297 | 3,279 | 2805 | 6386 | 0.47 | 0.001 | 0.012 | 1 | 17 | 0 |  | yes |
| okada__1578_pons | okada | snRNA | atlas | 20,294 | 20,192 | 1872 | 3546 | 0.18 | 0.048 | 0.039 | 0 | 102 | 0 |  | yes |
| okada__1656_pons | okada | snRNA | atlas | 4,047 | 4,026 | 1441 | 2535 | 2.31 | 0.000 | 0.017 | 0 | 21 | 0 |  | yes |
| okada__1689_thalamus | okada | snRNA | atlas | 12,762 | 12,696 | 1831 | 3200 | 0.14 | 0.000 | 0.024 | 3 | 63 | 0 |  | yes |
| okada__1712_pons | okada | snRNA | atlas | 3,672 | 3,653 | 1181 | 1776 | 2.55 | 0.000 | 0.019 | 0 | 19 | 0 |  | yes |
| okada__1732_pons | okada | snRNA | atlas | 153 | 152 | 1461 | 2538 | 1.77 | 0.007 | 0.000 | 0 | 1 | 0 |  | yes |
| okada__1734_pons | okada | snRNA | atlas | 15,228 | 8,486 | 427 | 468 | 0.51 | 0.860 | 0.002 | 6665 | 77 | 0 | FLAG | yes |
| okada__1855_pons | okada | snRNA | atlas | 13,501 | 13,433 | 2377 | 4447 | 0.18 | 0.000 | 0.028 | 0 | 68 | 0 |  | yes |
| okada__1890_pons | okada | snRNA | atlas | 15,969 | 15,890 | 2136 | 3974 | 0.12 | 0.018 | 0.025 | 0 | 79 | 0 |  | yes |
| okada__1934_thalamus | okada | snRNA | atlas | 3,877 | 3,848 | 2756 | 7588 | 0.10 | 0.011 | 0.018 | 9 | 20 | 0 |  | yes |
| okada__1944_pons | okada | snRNA | atlas | 3,615 | 3,596 | 2679 | 6356 | 0.16 | 0.016 | 0.019 | 0 | 19 | 0 |  | yes |
| okada__1958_pons | okada | snRNA | atlas | 3,146 | 3,130 | 2086 | 4502 | 0.25 | 0.000 | 0.010 | 0 | 16 | 0 |  | yes |
| okada__1963_pons | okada | snRNA | atlas | 2,971 | 2,956 | 1894 | 3690 | 1.20 | 0.001 | 0.010 | 0 | 15 | 0 |  | yes |
| okada__2054_thalamus | okada | snRNA | atlas | 7,369 | 7,326 | 1928 | 3840 | 0.48 | 0.019 | 0.021 | 6 | 37 | 0 |  | yes |
| okada__586_pons | okada | snRNA | atlas | 23,731 | 23,611 | 971 | 1489 | 0.00 | 0.093 | 0.036 | 1 | 119 | 0 |  | yes |
| okada__874_pons | okada | snRNA | atlas | 4,281 | 4,256 | 1437 | 2705 | 0.68 | 0.000 | 0.015 | 3 | 22 | 0 |  | yes |
| okada__916_thalamus | okada | snRNA | atlas | 7,753 | 7,714 | 2502 | 4652 | 0.13 | 0.000 | 0.027 | 0 | 39 | 0 |  | yes |
| okada__933_pons | okada | snRNA | atlas | 5,790 | 5,759 | 2212 | 4502 | 0.89 | 0.026 | 0.030 | 2 | 29 | 0 |  | yes |
| scpca__SCPCS000636 | scpca | snRNA | atlas | 1,054 | 0 | nan | nan | nan |  |  | 1049 | 6 | 0 | FLAG | yes |
| scpca__SCPCS000637 | scpca | snRNA | atlas | 2,931 | 441 | 362 | 436 | 2.22 | 0.057 | 0.011 | 2475 | 15 | 111 | FLAG | yes |
| zang__GSM7846769 | zang | scRNA | atlas | 7,410 | 6,579 | 2536 | 6970 | 7.22 | 0.000 | 0.005 | 41 | 38 | 770 |  | yes |
| desisto__1107 | desisto | scRNA | out-group | 1,190 | 595 | 1935 | 5352 | 8.32 | 0.002 | 0.002 | 163 | 6 | 588 |  | yes |
| desisto__1144 | desisto | scRNA | out-group | 4,617 | 3,022 | 1040 | 2495 | 11.48 | 0.000 | 0.007 | 321 | 24 | 1570 |  | yes |
| desisto__1187 | desisto | scRNA | out-group | 864 | 153 | 3193 | 8356 | 12.84 | 0.000 | 0.000 | 92 | 5 | 707 |  | yes |
| desisto__1190 | desisto | scRNA | out-group | 2,092 | 1,364 | 1195 | 2859 | 9.94 | 0.000 | 0.001 | 162 | 11 | 709 |  | yes |
| desisto__1214 | desisto | scRNA | out-group | 582 | 331 | 3118 | 13249 | 8.14 | 0.015 | 0.015 | 68 | 3 | 246 |  | yes |
| desisto__1241 | desisto | scRNA | out-group | 1,018 | 706 | 2922 | 9052 | 7.15 | 0.000 | 0.000 | 110 | 6 | 305 |  | yes |
| desisto__1414 | desisto | scRNA | out-group | 948 | 484 | 2886 | 10830 | 10.02 | 0.006 | 0.029 | 58 | 5 | 459 |  | yes |
| desisto__1421 | desisto | scRNA | out-group | 705 | 444 | 3815 | 15918 | 8.52 | 0.007 | 0.007 | 107 | 4 | 255 |  | yes |
| desisto__1431 | desisto | scRNA | out-group | 641 | 303 | 3043 | 12055 | 8.05 | 0.000 | 0.000 | 153 | 4 | 334 |  | yes |
| desisto__1458 | desisto | scRNA | out-group | 2,189 | 1,670 | 3325 | 9442 | 10.66 | 0.001 | 0.018 | 53 | 11 | 501 |  | yes |
| desisto__1463_B | desisto | scRNA | out-group | 234 | 138 | 1806 | 6393 | 10.66 | 0.000 | 0.000 | 15 | 2 | 93 |  | yes |
| desisto__834 | desisto | scRNA | out-group | 2,790 | 788 | 3552 | 10271 | 12.74 | 0.000 | 0.003 | 216 | 14 | 1989 |  | yes |
| desisto__868 | desisto | scRNA | out-group | 1,725 | 720 | 3424 | 12478 | 8.34 | 0.000 | 0.000 | 214 | 9 | 996 |  | yes |
| desisto__905 | desisto | scRNA | out-group | 1,906 | 586 | 1706 | 4902 | 9.55 | 0.003 | 0.005 | 415 | 10 | 1310 |  | yes |
| desisto__910 | desisto | scRNA | out-group | 2,769 | 2,229 | 1592 | 3012 | 9.27 | 0.002 | 0.002 | 24 | 14 | 524 |  | yes |
| desisto__954 | desisto | scRNA | out-group | 486 | 324 | 1204 | 2940 | 10.02 | 0.003 | 0.003 | 3 | 3 | 159 |  | yes |
| desisto__968 | desisto | scRNA | out-group | 1,508 | 463 | 4362 | 12219 | 12.57 | 0.000 | 0.000 | 63 | 8 | 1037 |  | yes |
| jessa_scrna__BT2016062 | jessa_scrna | scRNA | out-group | 3,415 | 3,396 | 1680 | 4182 | 2.69 | 0.000 | 0.004 | 0 | 18 | 1 |  | yes |
| jessa_scrna__BT2018022 | jessa_scrna | scRNA | out-group | 3,435 | 3,417 | 2401 | 6863 | 4.07 | 0.016 | 0.019 | 0 | 18 | 0 |  | yes |
| jessa_scrna__P-1190_S-1197 | jessa_scrna | snRNA | out-group | 4,775 | 4,669 | 2378 | 3929 | 0.60 | 0.000 | 0.012 | 0 | 24 | 83 |  | yes |
| jessa_scrna__P-1569_S-1569 | jessa_scrna | snRNA | out-group | 1,302 | 1,288 | 2551 | 5150 | 0.42 | 0.005 | 0.000 | 1 | 7 | 7 |  | yes |
| jessa_scrna__P-1969_S-3584 | jessa_scrna | snRNA | out-group | 10,965 | 10,352 | 2623 | 4776 | 1.37 | 0.000 | 0.001 | 0 | 55 | 569 |  | yes |
| jessa_scrna__P-2077_S-2077 | jessa_scrna | snRNA | out-group | 13,380 | 9,742 | 2410 | 4728 | 1.78 | 0.000 | 0.018 | 0 | 67 | 3577 |  | yes |
| jessa_scrna__P-3198_S-3199 | jessa_scrna | snRNA | out-group | 3,294 | 3,277 | 1341 | 2272 | 0.03 | 0.000 | 0.016 | 0 | 17 | 0 |  | yes |
| jessa_scrna__P-3200_S-3254 | jessa_scrna | snRNA | out-group | 2,511 | 2,498 | 1370 | 2560 | 0.11 | 0.001 | 0.021 | 0 | 13 | 0 |  | yes |
| jessa_scrna__P-3200_S-3867 | jessa_scrna | snRNA | out-group | 5,298 | 827 | 1107 | 1477 | 3.80 | 0.000 | 0.005 | 0 | 27 | 4460 |  | yes |
| jessa_scrna__P-3387_S-3411 | jessa_scrna | snRNA | out-group | 2,896 | 2,881 | 720 | 977 | 0.00 | 0.001 | 0.034 | 0 | 15 | 0 |  | yes |
| jessa_scrna__P-5099_S-6218 | jessa_scrna | snRNA | out-group | 1,377 | 1,369 | 1019 | 1483 | 0.00 | 0.005 | 0.015 | 1 | 7 | 0 |  | yes |
| jessa_scrna__P-5613_S-7162 | jessa_scrna | snRNA | out-group | 2,642 | 1,236 | 1468 | 2583 | 2.33 | 0.000 | 0.014 | 0 | 14 | 1400 |  | yes |
| jessa_scrna__P-6166_S-8321 | jessa_scrna | scRNA | out-group | 1,205 | 1,197 | 4237 | 17963 | 3.17 | 0.006 | 0.006 | 0 | 7 | 1 |  | yes |
| jessa_scrna__P-6240_S-8628 | jessa_scrna | scRNA | out-group | 3,585 | 2,184 | 1008 | 1774 | 7.90 | 0.000 | 0.007 | 10 | 18 | 1383 |  | yes |
| jessa_scrna__P-6292_S-8579 | jessa_scrna | scRNA | out-group | 5,922 | 5,336 | 1070 | 2436 | 10.35 | 0.000 | 0.009 | 2 | 30 | 558 |  | yes |
| jessa_scrna__P-6308_S-8632 | jessa_scrna | scRNA | out-group | 3,652 | 1,492 | 758 | 1366 | 11.96 | 0.000 | 0.009 | 34 | 19 | 2121 |  | yes |
| jessa_scrna__P-6329_S-8673 | jessa_scrna | snRNA | out-group | 2,486 | 2,258 | 1988 | 3606 | 1.37 | 0.000 | 0.001 | 0 | 13 | 215 |  | yes |
| jessa_scrna__P-6331_S-8904 | jessa_scrna | snRNA | out-group | 13,370 | 12,667 | 1963 | 3291 | 0.38 | 0.000 | 0.003 | 0 | 67 | 636 |  | yes |
| jessa_scrna__P-6388_S-8770 | jessa_scrna | snRNA | out-group | 3,523 | 3,166 | 1607 | 2436 | 0.28 | 0.001 | 0.003 | 1 | 18 | 339 |  | yes |
| jessa_scrna__P-6431_S-8842 | jessa_scrna | scRNA | out-group | 12,652 | 11,788 | 2494 | 6728 | 8.19 | 0.000 | 0.002 | 46 | 64 | 798 |  | yes |

## Figures

- results/qc/qc_violin_jessa_scrna.png — per-sample violins (all cells before QC) of log10(genes detected+1), log10(counts+1) and % mitochondrial; dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; one row
- results/qc/qc_violin_jessa_multiome.png — per-sample violins (all cells before QC) of log10(genes detected+1), log10(counts+1) and % mitochondrial; dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; one row
- results/qc/qc_violin_collot.png — per-sample violins (all cells before QC) of log10(genes detected+1), log10(counts+1) and % mitochondrial; dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; one row
- results/qc/qc_violin_desisto.png — per-sample violins (all cells before QC) of log10(genes detected+1), log10(counts+1) and % mitochondrial; dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; one row
- results/qc/qc_violin_zang.png — per-sample violins (all cells before QC) of log10(genes detected+1), log10(counts+1) and % mitochondrial; dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; one row
- results/qc/qc_violin_scpca.png — per-sample violins (all cells before QC) of log10(genes detected+1), log10(counts+1) and % mitochondrial; dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; one row
- results/qc/qc_violin_okada.png — per-sample violins (all cells before QC) of log10(genes detected+1), log10(counts+1) and % mitochondrial; dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; one row
- results/qc/qc_violin_filbin.png — per-sample violins (all cells before QC) of log10(genes detected+1); dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; one row; TPM/log-scale cohort: counts and % mito panels omitted (no MT- genes, values not counts)
- results/qc/qc_violin_liu.png — per-sample violins (all cells before QC) of log10(genes detected+1), log10(counts+1) and % mitochondrial; dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; two rows of samples
- results/qc/qc_violin_labelle.png — per-sample violins (all cells before QC) of log10(genes detected+1); dashed segments = per-sample min-genes (and % mito) thresholds, dotted = 99.5th-percentile max-genes cap; one row; TPM/log-scale cohort: counts and % mito panels omitted (no MT- genes, values not counts)
