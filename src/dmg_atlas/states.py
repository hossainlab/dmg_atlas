"""Malignant-state gene sets used as the FALLBACK state assignment for the cell2location reference (Phase 6).

Used only when data/processed/metaprogram_scores.parquet (Phase-6 metaprogram track) is not available.
Sets are curated, abbreviated versions of published signatures (top genes recalled from the published tables, not the
full lists): OPC-like / AC-like / OC-like / cycling from Filbin et al. 2018 Science (H3K27M glioma single-cell programs),
MES-like from the Liu et al. 2022 Nat Genet / Neftel et al. 2019 Cell mesenchymal programs.
PROVENANCE (audit MC-03, 2026-10-19): these lists are NOT traceable to a downloaded table. Checksummed alternatives are
produced by scripts/01d_reference_gene_sets.py (results/tables/genesets/reference_gene_sets_scripted.csv: Liu 2022 Table S2a
top-50 metaprogram genes; the Filbin 2018 tables could not be fetched by script). Overlap of these curated sets with the
scripted ones: results/tables/genesets/reference_gene_sets_overlap.csv (Jaccard 0.18-0.37). Gene symbols are 2020-A
(GRCh38) symbols. Scoring: scanpy.tl.score_genes on log-normalised counts, assignment = argmax; cells whose maximum
score is <= 0 are left unassigned and excluded from the reference.
"""

STATE_GENE_SETS = {
    "OPC_like": ["PDGFRA", "OLIG1", "OLIG2", "SOX10", "CSPG4", "GPR17", "PTPRZ1", "BCAN", "SOX6", "ASCL1", "NKX2-2",
                 "CSPG5", "PCDH15", "LHFPL3", "VCAN", "SOX8", "APOD", "CCND1", "EPN2", "TNR", "FERMT1", "SMOC1",
                 "DLL3", "HES6", "SOX4", "SOX11", "MEG3", "NEU4", "LIMA1", "FIBIN"],
    "AC_like": ["GFAP", "AQP4", "APOE", "CLU", "SPARCL1", "ALDOC", "SLC1A3", "SLC1A2", "ATP1A2", "GJA1", "MLC1", "AGT",
                "FABP7", "S100B", "CST3", "EDNRB", "F3", "ID3", "ID4", "NDRG2", "SOX9", "CRYAB", "HOPX", "GLUL",
                "ALDH1L1", "SPARC", "PTGDS", "MT3", "ITM2C", "ZFP36L1"],
    "OC_like": ["MBP", "PLP1", "MOG", "MAG", "CNP", "CLDN11", "MOBP", "TF", "OPALIN", "ERMN", "UGT8", "ENPP2", "MAL",
                "SIRT2", "BCAS1", "CNTN2", "KLK6", "SEPTIN4", "GPR37", "TMEM144", "CARNS1", "ANLN", "PLLP", "QDPR",
                "APLP1", "TPPP", "NKX6-2", "MYRF", "ST18", "SOX10"],
    "MES_like": ["CD44", "VIM", "CHI3L1", "ANXA1", "ANXA2", "LGALS1", "LGALS3", "S100A10", "S100A11", "S100A6", "TIMP1",
                 "SERPINE1", "TNC", "NAMPT", "MT2A", "HILPDA", "VEGFA", "ADM", "NDRG1", "SOD2", "EMP1", "CLIC1", "FN1",
                 "TAGLN2", "TGFBI", "HMOX1", "CAV1", "IGFBP3", "C1S", "SERPING1"],
    "cycling": ["MKI67", "TOP2A", "CENPF", "PCNA", "MCM2", "CCNB1", "CDK1", "UBE2C", "BIRC5", "NUSAP1", "TPX2", "ASPM", "HMGB2", "AURKB", "CCNA2", "TYMS", "RRM2", "PTTG1", "PCLAF", "CENPE", "SMC4", "H4C3", "CDC20", "PBK", "KIF2C", "DLGAP5", "CENPU", "HMMR", "ZWINT"],
}

# normal-lineage markers appended to the cell2location gene set (same sets as Phase 3, cnv.py, plus a few vascular/neuron genes)
NORMAL_MARKERS = {
    "myeloid": ["PTPRC", "CSF1R", "CX3CR1", "P2RY12", "C1QA", "C1QB", "CD74", "AIF1", "CD163", "CD14", "LYZ", "HLA-DRA"],
    "T_NK": ["CD3D", "CD3E", "CD2", "NKG7", "CD8A", "GZMB", "IL7R", "CCL5", "TRAC", "GNLY"],
    "oligodendrocyte": ["MBP", "PLP1", "MOG", "MAG", "ST18", "MOBP", "CLDN11"],
    "astrocyte": ["AQP4", "GFAP", "ALDH1L1", "SLC1A2", "GJA1"],
    "neuron": ["RBFOX3", "SYT1", "SNAP25", "GAD1", "GAD2", "SLC17A7", "SLC17A6", "STMN2", "NEFL", "NEFM"],
    "vascular": ["CLDN5", "VWF", "PECAM1", "FLT1", "PDGFRB", "RGS5", "ACTA2", "COL1A1", "COL1A2", "DCN", "KDR"],
}


# TME sub-type markers for the Phase-6 subclustering annotation (scripts/17_tme_subclusters.py); canonical markers from
# brain-tumour myeloid (microglia vs BAM vs monocyte-derived; Pombo Antunes 2021, Abdelfattah 2022), T/NK (Mathewson 2021),
# brain vasculature (Garcia 2022, Yang 2022 vascular atlases) and oligodendroglial lineage (Marques 2016) literature.
# Not used in all_marker_genes() so the cell2location gene set is unchanged.
TME_SUBTYPE_MARKERS = {
    "myeloid": {
        "homeostatic_microglia": ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "CSF1R", "GPR34", "SELPLG", "FCRLS", "OLFML3", "SIGLEC8"],
        "activated_microglia": ["CCL3", "CCL4", "IL1B", "TNF", "NFKBIA", "EGR1", "FOS", "JUN", "CD83", "IER3"],
        "lipid_phagocytic_TAM": ["APOE", "APOC1", "SPP1", "GPNMB", "LPL", "FABP5", "TREM2", "CD9", "LGALS3", "CTSB"],
        "border_associated_macrophage": ["MRC1", "LYVE1", "F13A1", "CD163", "SIGLEC1", "STAB1", "CD38", "COLEC12", "PF4", "MS4A7"],
        "monocyte": ["S100A8", "S100A9", "VCAN", "FCN1", "LYZ", "CD14", "SELL", "CCR2", "S100A12", "LILRB2"],
        "interferon_myeloid": ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "RSAD2", "OAS1", "IFI6", "CXCL10", "STAT1"],
        "dendritic_cell": ["CLEC9A", "XCR1", "CD1C", "FCER1A", "CLEC10A", "LAMP3", "CCR7", "IDO1", "FLT3", "HLA-DQA1"],
        "proliferating_myeloid": ["MKI67", "TOP2A", "CENPF", "NUSAP1", "UBE2C", "TPX2", "ASPM", "BIRC5", "CDK1", "PCLAF"],
    },
    "T_NK": {
        "CD8_T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "CD3D", "CD3E", "CRTAM", "CST7", "NKG7"],
        "CD4_T": ["CD4", "IL7R", "CD40LG", "TCF7", "LTB", "CCR7", "MAL", "CD3D", "CD3E", "LEF1"],
        "regulatory_T": ["FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "TNFRSF4", "IKZF2", "BATF", "CD4", "TIGIT", "LAYN"],
        "exhausted_T": ["PDCD1", "HAVCR2", "LAG3", "TOX", "CXCL13", "ENTPD1", "TIGIT", "CTLA4", "LAYN", "ITGAE"],
        "NK": ["NCAM1", "KLRD1", "KLRF1", "GNLY", "NKG7", "FCGR3A", "PRF1", "GZMB", "KLRC1", "TYROBP"],
        "proliferating_T_NK": ["MKI67", "TOP2A", "CENPF", "NUSAP1", "UBE2C", "TPX2", "ASPM", "BIRC5", "CDK1", "PCLAF"],
    },
    "vascular": {
        "endothelial_capillary": ["CLDN5", "SLC2A1", "MFSD2A", "SLC7A5", "CA4", "SPOCK2", "FLT1", "ABCB1", "ITM2A", "PODXL"],
        "endothelial_arterial": ["GJA5", "BMX", "SEMA3G", "EFNB2", "HEY1", "VEGFC", "FBLN5", "ELN", "SOX17", "DLL4"],
        "endothelial_venous": ["ACKR1", "SELP", "SELE", "VCAM1", "NR2F2", "EMCN", "TSHZ2", "PLVAP", "IL1R1", "POSTN"],
        "endothelial_tip_angiogenic": ["ESM1", "ANGPT2", "PGF", "APLN", "KCNE3", "CXCR4", "INSR", "COL4A1", "COL4A2", "PLVAP"],
        "pericyte": ["PDGFRB", "RGS5", "KCNJ8", "ABCC9", "HIGD1B", "NDUFA4L2", "CSPG4", "VTN", "ATP13A5", "GRM8"],
        "smooth_muscle": ["ACTA2", "MYH11", "TAGLN", "DES", "CNN1", "MYL9", "LMOD1", "SMTN", "ACTG2", "PLN"],
        "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "FBLN1", "CEMIP", "SLC6A13", "COL6A1", "MGP"],
    },
    "oligodendrocyte": {
        "OPC": ["PDGFRA", "CSPG4", "VCAN", "PTPRZ1", "OLIG1", "OLIG2", "GPR17", "SOX6", "PCDH15", "LHFPL3"],
        "committed_OPC_newly_formed": ["BCAS1", "ENPP6", "GPR17", "TCF7L2", "ITPR2", "FYN", "TNS3", "SIRT2", "BMP4", "NKX2-2"],
        "mature_oligodendrocyte": ["MOG", "MOBP", "MBP", "PLP1", "MAG", "CLDN11", "OPALIN", "ERMN", "KLK6", "TF"],
        "stressed_oligodendrocyte": ["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "BAG3", "HSPB1", "FOS", "JUN", "ATF3", "CRYAB"],
    },
}


def all_marker_genes() -> list[str]:
    out = []
    for d in (STATE_GENE_SETS, NORMAL_MARKERS):
        for genes in d.values():
            out.extend(genes)
    return sorted(set(out))


# Curated Phase-6 reference / technical gene sets (audit MC-03, 2026-10-19). Previously these lived only in
# results/tables/genesets/phase6_gene_sets.csv, which no script produced; scripts/16a_curated_gene_sets.py now writes
# results/tables/genesets/phase6_curated_gene_sets.csv (with sha256) from this dictionary. The former 'Liu2022_*' sets of that
# file were in-house Wilcoxon DE sets, not the published tables, and are superseded by reference_gene_sets_scripted.csv
# (scripts/01d_reference_gene_sets.py: Liu 2022 Table S2a published sets + the in-house re-derived DE sets, clearly named).
# Hand-curated lists (hallmark-style; canonical marker genes), NOT downloaded tables: hypoxia, glycolysis, cholesterol
# biosynthesis, OXPHOS (nuclear-encoded), interferon response, heat-shock stress, the four 20-gene lineage marker sets,
# cycling markers; added in the audit: IEG_dissociation (immediate-early / heat-shock dissociation response, van den Brink
# et al. 2017 Nat Methods core genes + canonical IEGs) and Immune_pan (pan-leukocyte / myeloid / T-cell genes used to flag
# immune contamination of the malignant selection).
CURATED_PHASE6_SETS = {
    "AC_like_markers": [
        "GFAP", "AQP4", "APOE", "CLU", "SPARCL1", "ALDOC", "S100B", "SLC1A3", "AGT", "MLC1", "GJA1", "ID3", "SPARC",
        "CST3", "ATP1A2", "NDRG2", "SLC1A2", "F3", "FABP7", "HOPX"],
    "Cholesterol_biosynthesis": [
        "HMGCS1", "HMGCR", "MVD", "MVK", "PMVK", "IDI1", "FDPS", "FDFT1", "SQLE", "LSS", "CYP51A1", "MSMO1", "NSDHL",
        "HSD17B7", "DHCR7", "DHCR24", "SC5D", "EBP", "TM7SF2", "LDLR", "INSIG1", "SREBF2", "ACAT2", "STARD4"],
    "Cycling_markers": [
        "MKI67", "TOP2A", "CENPF", "NUSAP1", "TPX2", "UBE2C", "BIRC5", "CCNB1", "CDK1", "PCNA", "MCM2", "MCM3",
        "MCM4", "MCM5", "MCM6", "MCM7", "TYMS", "RRM2", "HMGB2", "CDC20", "PTTG1", "AURKA", "AURKB", "KIF2C", "KIF11",
        "SMC4", "ASPM", "CENPE", "CKS2", "ZWINT"],
    "Glycolysis": [
        "HK1", "HK2", "GPI", "PFKL", "PFKM", "PFKP", "ALDOA", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1",
        "ENO2", "PKM", "LDHA", "SLC2A1", "SLC2A3", "PFKFB3", "PFKFB4", "PGM1", "GYS1"],
    "Hypoxia": [
        "VEGFA", "ADM", "CA9", "BNIP3", "BNIP3L", "NDRG1", "SLC2A1", "PGK1", "LDHA", "ENO1", "P4HA1", "P4HA2",
        "EGLN3", "EGLN1", "ANKRD37", "DDIT4", "HILPDA", "PDK1", "ALDOA", "MXI1", "ERO1A", "PLOD2", "ANGPTL4",
        "IGFBP3", "KDM3A", "PFKFB3", "LOX", "INSIG2", "HK2", "NDUFA4L2", "FAM162A"],
    "IEG_dissociation": [
        "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "IER5", "ATF3", "NR4A1",
        "NR4A2", "NR4A3", "DUSP1", "DUSP6", "KLF2", "KLF4", "KLF6", "KLF9", "ZFP36", "ZFP36L1", "ZFP36L2", "GADD45A",
        "GADD45B", "GADD45G", "BTG1", "BTG2", "SOCS3", "CEBPB", "CEBPD", "NFKBIA", "NFKBIZ", "PPP1R15A", "MCL1",
        "CDKN1A", "CCNL1", "RHOB", "ARC", "MAFF", "TXNIP", "HSPA1A", "HSPA1B", "HSPB1", "DNAJB1", "HSPH1", "HSP90AA1",
        "HSPA8", "HSPA6", "DNAJA1", "BAG3"],
    "Immune_pan": [
        "PTPRC", "DOCK8", "APBB1IP", "INPP5D", "CSF1R", "CD74", "LAPTM5", "SRGN", "LCP2", "SYK", "FYB1", "ITGAX",
        "MSR1", "CD163", "C3", "DOCK2", "ARHGAP15", "SAMSN1", "PIK3R5", "TBXAS1", "SLCO2B1", "MS4A6A", "MS4A7",
        "CD14", "AIF1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
        "CD53", "CD86", "CX3CR1", "P2RY12", "LYN", "HCK", "FGR", "PLCG2", "LCP1", "CORO1A", "CD3E", "CD3D", "CD2",
        "SKAP1", "THEMIS", "CD247", "IKZF1", "RUNX3", "ITGB2", "SPI1", "IRF8", "MERTK", "SLC11A1", "CYBB", "RGS1",
        "LY86", "CD37", "ARHGDIB", "CTSS", "FCGR2A", "FCGR3A", "MPEG1", "SPP1", "APOC1", "LPCAT2", "ADAM28", "STAB1"],
    "Interferon": [
        "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3", "MX1", "MX2",
        "OAS1", "OAS2", "OAS3", "OASL", "STAT1", "STAT2", "IRF7", "IRF9", "BST2", "RSAD2", "XAF1", "HERC5", "HERC6",
        "EPSTI1", "PARP9", "PARP14", "DDX58", "IFIH1", "SAMD9", "SAMD9L", "USP18", "CMPK2", "LY6E", "B2M", "HLA-A",
        "HLA-B", "HLA-C", "PSMB8", "PSMB9", "TAP1"],
    "MES_like_markers": [
        "CD44", "VIM", "CHI3L1", "ANXA1", "ANXA2", "LGALS1", "LGALS3", "S100A10", "S100A11", "TIMP1", "SERPINE1",
        "TNC", "MT2A", "NAMPT", "SOD2", "TAGLN2", "EMP1", "MYL9", "IGFBP5", "S100A6"],
    "OC_like_markers": [
        "MBP", "PLP1", "MOG", "MAG", "CNP", "TF", "APOD", "MOBP", "CLDN11", "ERMN", "BCAS1", "MYRF", "SIRT2",
        "TUBB4A", "NFASC", "ENPP6", "GPR37", "OPALIN", "UGT8", "KLK6"],
    "OPC_like_markers": [
        "PDGFRA", "OLIG1", "OLIG2", "SOX10", "CSPG4", "ASCL1", "DLL3", "SOX6", "BCAN", "PTPRZ1", "CCND1", "APOD",
        "EPN2", "GPR17", "NKX2-2", "SOX8", "FGFR3", "TNR", "LHFPL3", "PCDH15"],
    "OXPHOS_nuclear": [
        "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF",
        "ATP5MG", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C",
        "COX7A2", "COX7B", "COX7C", "COX8A", "CYC1", "NDUFA1", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13", "NDUFA2",
        "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFB1", "NDUFB10",
        "NDUFB11", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFC1", "NDUFC2",
        "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "NDUFV3",
        "SDHA", "SDHB", "SDHC", "SDHD", "UQCR10", "UQCR11", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ"],
    "Stress_HSP": [
        "HSPA1A", "HSPA1B", "HSPA6", "HSPA8", "HSPH1", "HSP90AA1", "HSP90AB1", "HSPB1", "HSPD1", "HSPE1", "DNAJB1",
        "DNAJA1", "DNAJB4", "HSPA4", "HSPA5", "FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3", "IER2", "IER3", "DDIT3",
        "PPP1R15A", "BAG3", "CRYAB", "UBC", "UBB", "GADD45B", "NR4A1"],
}
