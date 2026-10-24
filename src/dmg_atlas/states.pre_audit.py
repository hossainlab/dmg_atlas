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
