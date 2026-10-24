#!/usr/bin/env python
"""Phase 1 (tumour track): build data/processed/sample_manifest.csv.

One row per biological sample across the tumour cohorts (integrated + projection), plus
placeholder rows (file_paths blank) for the spatial / normal-reference / bulk datasets handled
by the other track. Inputs are the raw deposits under data/raw/ and the per-study label tables in
data/interim/labels/ produced in Phase 1. No value is invented: anything not stated in a deposit,
a paper table or a portal record is 'unknown'.
"""
import glob, gzip, json, os, re
import pandas as pd

import pathlib as _pl
ROOT = str(_pl.Path(__file__).resolve().parents[1])  # audit DQ-14: no hard-coded machine path
RAW = f"{ROOT}/data/raw"
LAB = f"{ROOT}/data/interim/labels"
import sys as _sys
OUT = _sys.argv[1] if len(_sys.argv) > 1 else f"{ROOT}/data/processed/sample_manifest.csv"  # optional output path (audit DQ-05/DQ-14)

COLS = ["sample_id", "patient_id", "study", "role", "assay", "chemistry", "site", "age_years", "age_group",
        "h3_variant", "k27m_status", "sex", "treatment", "tissue_timing", "species", "include_in_atlas",
        "exclusion_reason", "counts_type", "n_cells_deposited", "accession", "file_paths", "label_source", "notes"]


def rel(p):
    return os.path.relpath(p, ROOT)


def age_group(a):
    try:
        a = float(a)
    except (TypeError, ValueError):
        return "unknown"
    return "pediatric" if a < 18 else "adult"


def mtx_dims(path):
    with gzip.open(path, "rt") as f:
        l = f.readline()
        while l.startswith("%"):
            l = f.readline()
    return tuple(int(x) for x in l.split())


rows = []


def add(**kw):
    r = {c: kw.get(c, "unknown" if c in ("chemistry", "site", "age_group", "h3_variant", "k27m_status", "sex",
                                               "treatment", "tissue_timing") else "") for c in COLS}
    r.update(kw)
    r["species"] = kw.get("species", "human")  # audit DQ-07: setdefault never fired because every COLS key was pre-filled
    rows.append(r)


# ---------------------------------------------------------------- Jessa 2022 (GSE210567 scRNA/snRNA, GSE210566 multiome)
jl = pd.read_csv(f"{LAB}/jessa_sample_labels.csv")
js = pd.read_csv(f"{LAB}/jessa_scrna_geo_sample_summary.csv")
jm = pd.read_csv(f"{LAB}/jessa_multiome_geo_sample_summary.csv")
JSRC = "Jessa2022 NatGenet (10.1038/s41588-022-01205-w) MOESM3 Table S1 cols Age/Sex/Location/Group/Histone_status/scRNAseq; Table S3 Kit"


def jessa_excl(r):
    if r.is_dmg:
        return True, ""
    g = r.Group
    if g == "HGG-H3.3K27M-Cortex":
        return False, "H3.3K27M but hemispheric (cortex), not midline"
    if g == "PFA-H3.1K27M-PF":
        return False, "posterior fossa ependymoma (PFA) carrying H3.1K27M, not a glioma"
    if "G34" in g:
        return False, "H3.3 G34R/V hemispheric glioma (out-group)"
    if "H3WT" in g:
        return False, "H3-wild-type hemispheric glioma (out-group)"
    if "PFA" in g:
        return False, "PFA ependymoma (EZHIP), not a glioma"
    return False, "non-DMG"


for _, r in jl.iterrows():
    inc, why = jessa_excl(r)
    common = dict(patient_id=r.patient_id, site=r.site, age_years=r.Age if pd.notna(r.Age) else "",
                  age_group=age_group(r.Age), h3_variant=r.h3_variant,
                  # audit DQ-05: H3.3 G34R/V tumours are K27-wild-type (rule previously applied by hand to the on-disk manifest)
                  k27m_status="K27M" if "K27M" in str(r.Histone_status) else ("WT" if r.h3_variant in ("WT", "H3.3_G34") else "unknown"),
                  sex=r.Sex if pd.notna(r.Sex) else "unknown", tissue_timing="unknown",
                  treatment="unknown", include_in_atlas=inc, exclusion_reason=why, counts_type="umi",
                  label_source=JSRC, notes=f"Group={r.Group}; Histone_status={r.Histone_status}; Material={r.Material}")
    if r.in_GSE210567:
        s = js[js["sample"].str.endswith(r.ID)].iloc[0]
        kit = {2.0: "10x 3' v2", 3.0: "10x 3' v3"}.get(r.kit_version, "10x 3' (version not stated)")
        add(sample_id=f"jessa_scrna__{r.ID}", study="jessa_scrna", role="integrated",
            assay=r.assay_scrna if isinstance(r.assay_scrna, str) and r.assay_scrna else "unknown", chemistry=kit,
            n_cells_deposited=int(s.qc_cells), accession=f"GSE210567/{s['sample'].split('_')[0]}",
            file_paths=rel(f"{RAW}/jessa_scrna/samples/{s['sample']}"), **common)
    if r.in_GSE210566:
        s = jm[jm["sample"] == r.ID].iloc[0]
        add(sample_id=f"jessa_multiome__{r.ID}", study="jessa_multiome", role="integrated", assay="multiome_RNA",
            chemistry="10x Multiome (ARC) NEXT-GEM v1", n_cells_deposited=int(s.qc_cells), accession=f"GSE210566/{s.gsm}",
            file_paths=rel(f"{RAW}/jessa_multiome/{s.gsm}"), **{**common, "notes": common["notes"] +
                 "; features.tsv contains Gene Expression + Peaks rows - subset to Gene Expression at ingest" +
                 ("; patient also in GSE210567" if r.in_GSE210567 else "")})

# ---------------------------------------------------------------- Collot 2026 (Zenodo 17922773)
cs = pd.read_csv(f"{LAB}/collot_sample_summary.csv")
cl = pd.read_csv(f"{LAB}/collot_sample_labels.csv")
cl["pid"] = cl["Paper ID"].str.replace("PMC-DMG-0", "PMC-DMG-", regex=False)
CSRC = "Collot2026 CancerCell (10.1016/j.ccell.2025.12.020) mmc2.xlsx Table S1A cols Tumor_subtype/Location/Gender/Age_years/Primary-recurrent; cells from scRNA_metadata.csv.gz"
CH = f"{RAW}/collot/scRNA/scRNA/Human"
cfiles = ";".join(rel(f"{CH}/{x}") for x in ["scRNA_counts.mtx.gz", "scRNA_genes.tsv.gz", "scRNA_cells.tsv.gz", "scRNA_metadata.csv.gz"])
for _, r in cs.iterrows():
    sid = r.SampleID
    if r.Study == "Jessa2022":
        base = re.sub(r"_(RNA_only|Multiome)$", "", sid)
        add(sample_id=f"collot__{sid}", patient_id=base.split("_S-")[0], study="collot", role="integrated", assay="snRNA" if "Multiome" in sid else "unknown",
            n_cells_deposited=int(r.n_cells), accession="zenodo.17922773", file_paths=cfiles, counts_type="umi",
            include_in_atlas=False, exclusion_reason=f"reprocessed copy of Jessa 2022 sample {base} ({'GSE210566' if 'Multiome' in sid else 'GSE210567'}); GEO original is used",
            label_source=CSRC, notes="Collot Table S1A annotations for Jessa samples disagree with Jessa Table S1 for several IDs; Jessa table is authoritative")
        continue
    pid = {"DMG-PMC-9_core": "PMC-DMG-9", "DMG-PMC-9_edge": "PMC-DMG-9"}.get(sid, sid)
    L = cl[cl.pid == pid].iloc[0]
    site = {"Pons": "pons", "Thalamus": "thalamus"}.get(L.Location, "unknown")
    add(sample_id=f"collot__{sid}", patient_id=pid, study="collot", role="integrated", assay="snRNA", chemistry="10x (version not stated)",
        site=site, age_years=L.Age_years if pd.notna(L.Age_years) else "", age_group=age_group(L.Age_years),
        h3_variant="H3.3" if "H3.3" in str(L.Tumor_subtype) else "unknown", k27m_status="K27M", sex=L.Gender if pd.notna(L.Gender) else "unknown",
        treatment="recurrent (post-therapy)" if str(L["Primary/ recurrent"]) == "Recurrent" else ("treatment-naive (primary)" if str(L["Primary/ recurrent"]) == "Primary" else "unknown"),
        tissue_timing="unknown", include_in_atlas=True, exclusion_reason="", counts_type="umi", n_cells_deposited=int(r.n_cells),
        accession="zenodo.17922773", file_paths=cfiles, label_source=CSRC,
        notes=f"Institute={L.Institute}; {'core/edge regions of one tumour' if 'PMC-9' in sid else ''}; assay given as snRNA per Table S1A sheet name; deposit is 19,248-gene filtered integrated matrix")

# ---------------------------------------------------------------- DeSisto 2024 (GSE231859)
dl = pd.read_csv(f"{LAB}/desisto_sample_labels.csv", dtype=str)
DSRC = "ScPCA portal API SCPCP000001 sample metadata (same Colorado cohort; submitter_id == GEO sample name): diagnosis/tissue_location/age/sex/molecular_characteristics/disease_timing. Paper supplement (10.1093/neuonc/noad207) not retrievable (PMC proof-of-work, OUP Cloudflare)"
for d in sorted(glob.glob(f"{RAW}/desisto/samples/*")):
    gsm, name = os.path.basename(d).split("_", 1)
    L = dl[dl.geo_sample == name].iloc[0]
    g, b = mtx_dims(glob.glob(d + "/*matrix.mtx.gz")[0])[:2]
    midline = L.site in ("pons", "thalamus", "spinal")
    if L.k27m_status == "K27M":
        inc, why = True, ""
    elif midline and L.k27m_status == "unknown":
        inc, why = False, "midline location but H3K27M status not recorded - candidate pending K27M verification (e.g. K27M reads in matched bulk RNA, ScPCA SCPCP000001)"
    elif L.k27m_status == "WT":
        inc, why = False, "H3K27M-negative per molecular characteristics"
    else:
        inc, why = False, f"non-midline ({L.tissue_location}); {L.diagnosis}"
    add(sample_id=f"desisto__{name}", patient_id=f"CO-{L.participant_id}", study="desisto", role="integrated", assay="scRNA", chemistry="10x 3' v3",
        site=L.site, age_years=L.age, age_group=age_group(L.age), h3_variant="unknown", k27m_status=L.k27m_status, sex=L.sex,
        treatment="unknown (initial diagnosis)" if str(L.disease_timing).startswith("Initial") else f"recurrent: {L.disease_timing}",
        tissue_timing="resection", include_in_atlas=inc, exclusion_reason=why, counts_type="umi", n_cells_deposited=b,
        accession=f"GSE231859/{gsm}", file_paths=rel(d), label_source=DSRC,
        notes=f"diagnosis={L.diagnosis}; tissue_location={L.tissue_location}; molecular={L.molecular_characteristics}; ScPCA={L.scpca_id}; GEO sample 905 == ScPCA 905_3")

# ---------------------------------------------------------------- Zang 2024 (GSM7846769)
g, b = mtx_dims(f"{RAW}/zang/GSM7846769_matrix1.mtx.gz")[:2]
add(sample_id="zang__GSM7846769", patient_id="zang_DMG_BG", study="zang", role="integrated", assay="scRNA", chemistry="10x (version not stated)",
    site="pons", age_years=10, age_group="pediatric", h3_variant="unknown", k27m_status="K27M", sex="F", treatment="unknown",
    tissue_timing="resection", include_in_atlas=True, exclusion_reason="", counts_type="umi", n_cells_deposited=b, accession="GSE245674/GSM7846769",
    file_paths=rel(f"{RAW}/zang"), label_source="Zang2024 BMC Neurol (10.1186/s12883-024-03558-7) Supplementary Table 1 (12883_2024_3558_MOESM1_ESM.docx): Sex, Age, Sampling site, Pathological diagnosis",
    notes="GSM7846770 (anaplastic ependymoma) not downloaded; sampling site 'brain stem, pontine cerebellum mass'; H3K27M-mutant WHO IV by pathology")

# ---------------------------------------------------------------- ScPCA SCPCP000021
sm = pd.read_csv(f"{RAW}/scpca/unpacked/SCPCP000021_single-cell/single-cell_metadata.tsv", sep="\t", dtype=str)
for _, r in sm.iterrows():
    api = json.load(open(f"{RAW}/scpca/{r.scpca_sample_id}_api_sample.json"))
    am = api["additional_metadata"]
    d = f"{RAW}/scpca/unpacked/SCPCP000021_single-cell/{r.scpca_sample_id}"
    add(sample_id=f"scpca__{r.scpca_sample_id}", patient_id=f"CBTN-{r.participant_id}", study="scpca", role="integrated", assay="snRNA",
        chemistry="10x 3' v3.1", site="thalamus", age_years=float(r.age), age_group=age_group(r.age), h3_variant="unknown", k27m_status="K27M",
        sex=r.sex, treatment=f"{r.disease_timing}; prior chemo/RT ages (days): {am.get('age_at_chemo_start')}/{am.get('age_at_radiation_start')}",
        tissue_timing="resection", include_in_atlas=True, exclusion_reason="", counts_type="umi",
        n_cells_deposited=int(api["sample_cell_count_estimate"]), accession=f"SCPCP000021/{r.scpca_sample_id}/{r.scpca_library_id}",
        file_paths=rel(d), label_source="ScPCA portal single-cell_metadata.tsv + API sample record (subdiagnosis 'H3 K28-mutant;WHO grade IV', tissue_location, age, sex, disease_timing, mutation_status)",
        notes=f"mutation_status={am.get('mutation_status')}; extent_of_tumor_resection={am.get('extent_of_tumor_resection')}; cancer_predispositions={am.get('cancer_predispositions')}; CBTN submitter_id {r.submitter_id} - check overlap with OpenPedCan/PBTA; processed_rna.h5ad keeps only {'51' if r.scpca_sample_id=='SCPCS000636' else '1363'} cells after miQC - use filtered_rna.h5ad with own QC")

# ---------------------------------------------------------------- Okada 2026 (UCSC Cell Browser brain-dmg-multiomics)
ok = pd.read_csv(f"{RAW}/okada/snrna-seq-data/per_sample_summary.csv")
for _, r in ok.iterrows():
    add(sample_id=f"okada__{r['sample']}", patient_id=f"okada_{r['sample'].split('_')[0]}", study="okada", role="integrated", assay="snRNA",
        chemistry="10x (version not stated)", site=r.location.lower(), age_years=r.age, age_group=age_group(r.age), h3_variant="unknown",
        k27m_status="K27-altered (series scope)", sex=r.sex,
        treatment={"ONC201": "imipridone (ONC201) + standard care", "Standard": "standard care, no imipridone", "Untreated": "untreated"}.get(r.treatment, r.treatment),
        tissue_timing="autopsy", include_in_atlas=True, exclusion_reason="", counts_type="umi",
        n_cells_deposited=int(r.n_cells), accession="cellbrowser:brain-dmg-multiomics/snrna-seq-data (raw: EGAS50000001447)",
        file_paths=rel(f"{RAW}/okada/snrna-seq-data/brain-dmg-snRNA-Homo-sapiens.h5ad"),
        label_source="Cell Browser meta.tsv columns sample/location/treatement/age/sex (Okada 2026 Neuro-Oncology 10.1093/neuonc/noag119)",
        notes="h5ad raw.X = Seurat LogNormalize(1e4) over 36,601 genes; UMI counts recovered exactly as round(expm1(x)*nCount_RNA/1e4) (verified on 5,000 nuclei); author QC 200<nFeature<5000, mito<5% already applied")

# ---------------------------------------------------------------- Filbin 2018 (GSE102130, projection, TPM)
fc = pd.read_csv(f"{RAW}/filbin/GSE102130_sample_characteristics.csv")
FILB_ALIAS = {"BCH1126": "BCH1126", "BCH836": "BCH836", "BCH869": "BCH869", "MUV1": "MUV1", "MUV5": "MUV5", "MUV10": "MUV10", "MGH66": "MGH66", "MGH101": "MGH101", "MGH104": "MGH104"}
for pref, grp in fc.groupby("sample_prefix"):
    st = ";".join(sorted(set(grp["sample type"].astype(str))))
    is_tumor = st == "human tumor"
    add(sample_id=f"filbin__{pref}", patient_id=FILB_ALIAS.get(pref, pref), study="filbin", role="projection", assay="smartseq2", chemistry="Smart-seq2 plates",
        site="unknown", h3_variant="unknown", k27m_status="K27M (series scope)" if is_tumor else "n/a",
        treatment="unknown", tissue_timing="unknown", include_in_atlas=is_tumor,
        exclusion_reason="" if is_tumor else f"not patient tissue ({st})", counts_type="tpm", n_cells_deposited=len(grp),
        accession="GSE102130", file_paths=rel(f"{RAW}/filbin/GSE102130_K27Mproject.RSEM.vh20170621.txt.gz"),
        label_source="GEO SOFT sample characteristics ('sample type'); per-patient age/site not in GEO - paper Table S1 not fetched",
        notes="column prefix in the RSEM TPM table; patients BCH1126/BCH836/BCH869/MUV1/MUV5/MUV10/MGH66 re-appear in Liu 2022 GSE184357 (fresh set) - same patient_id used" if is_tumor else "")

# ---------------------------------------------------------------- Liu 2022 (GSE184357, projection)
ls = pd.read_csv(f"{LAB}/liu_geo_metadata_sample_summary.csv")
LIU_ALIAS = {"BT1126": "BCH1126", "BT836": "BCH836", "BT869": "BCH869", "MGH6666410": "MGH66", "SUDIPG50Aut": "SUDIPG50", "SUDIPG50Bio": "SUDIPG50"}
liumap = pd.read_csv(f"{LAB}/liu_rawfile_to_sample.csv")
for _, r in ls.iterrows():
    s = r["sample"]
    cand = [rel(f"{RAW}/liu/raw_counts/{fn}") for fn in liumap[liumap.metadata_sample == s].file]
    add(sample_id=f"liu__{s}", patient_id=LIU_ALIAS.get(s, s), study="liu", role="projection", assay="smartseq2",
        chemistry="Smart-seq2 plates (fresh cells)" if r.prep == "fresh" else "Smart-seq2 plates (frozen nuclei)",
        site={"pontine": "pons", "thalamic": "thalamus", "spinal": "spinal"}.get(r.location, "unknown"), age_group=r.age,
        h3_variant="unknown", k27m_status="K27M (series scope)", treatment="unknown",
        tissue_timing="autopsy" if s.endswith("Aut") else ("biopsy" if s.endswith("Bio") else "unknown"),
        include_in_atlas=True, exclusion_reason="", counts_type="other",
        n_cells_deposited=int(r.n_cells), accession="GSE184357",
        file_paths=";".join([rel(f"{RAW}/liu/GSE184357_metadata_{r.prep}_postFilter.txt.gz")] + cand),
        label_source="GEO GSE184357_metadata_{fresh,frozen}_postFilter.txt.gz columns sample/age/location",
        notes=f"{r.prep} prep ({'cells' if r.prep=='fresh' else 'nuclei'}); raw RSEM expected counts (mostly integer, some fractional) in GSE184357_DMG_raw_counts.tar.gz, mapped by file name - cell-ID formats differ from postFilter metadata; {int(r.n_malig)} malignant-annotated cells" + ("; NO raw count file for this sample" if not cand else "") + ("; same patient as Filbin 2018 " + LIU_ALIAS[s] if s in LIU_ALIAS and not s.startswith("SUDIPG") else ""))

# ---------------------------------------------------------------- LaBelle 2025 (GSE227983 + GSE271936)
gl = pd.read_csv(f"{LAB}/labelle_geo_metadata_sample_summary.csv")
meta_cols = {}
for fn, k in [("GSE227983_tumor_ss2_meta.csv.gz", "tumor_ss2"), ("GSE227983_immune_ss2_meta.csv.gz", "immune_ss2"), ("GSE227983_immune_tenx_meta.csv.gz", "immune_tenx")]:
    m = pd.read_csv(f"{RAW}/labelle/{fn}")
    meta_cols[k] = m.groupby("sampleid").agg(radio=("RadioChemoTreatBeforeSurgery", "first"), ster=("SteroidTreatment", "first"), imm=("ImmunotherapyTreatment", "first"))
for sid, grp in gl.groupby("sampleid"):
    files = sorted(set(grp.file))
    sub = grp.iloc[0]
    age = sub.Age
    age_y = 8 / 365 if str(age) == "8dys" else float(age)
    k27 = "K27M" if sub.Subtype == "H3K27M" else ("WT" if "H3WT" in str(sub.Subtype) or "G34" in str(sub.Subtype) else "n/a")
    inc = sub.Subtype == "H3K27M"
    why = "" if inc else ("H3-wild-type midline glioma (out-group)" if sub.Subtype == "H3WT_midline" else f"non-DMG ({sub.Subtype}, {sub.Location})")
    tr = next((meta_cols[f].loc[sid] for f in files if sid in meta_cols[f].index), None)
    add(sample_id=f"labelle__{sid}", patient_id=re.sub(r"_R\d$", "", sid), study="labelle", role="projection",
        assay="scRNA" if files == ["immune_tenx"] else ("smartseq2" if "immune_tenx" not in files else "smartseq2+scRNA"),
        chemistry="; ".join(("10x 5' (immune, sorted)" if f == "immune_tenx" else "Smart-seq2 (sorted " + ("tumour" if f == "tumor_ss2" else "CD45/CD3") + " cells)") for f in files),
        site="unknown (midline)" if sub.Location == "Midline" else "other", age_years=round(age_y, 3), age_group=age_group(age_y), h3_variant="unknown", k27m_status=k27, sex=sub.Sex,
        treatment=f"before surgery: {tr.radio}; steroid: {tr.ster}; immunotherapy: {tr.imm}" if tr is not None else "unknown",
        tissue_timing="resection", include_in_atlas=inc, exclusion_reason=why, counts_type="tpm", n_cells_deposited=int(grp.n.sum()), accession="GSE227983",
        file_paths=";".join(rel(f"{RAW}/labelle/GSE227983_{f}_tpm.csv.gz") + ";" + rel(f"{RAW}/labelle/GSE227983_{f}_meta.csv.gz") for f in files),
        label_source="GEO GSE227983 per-cell metadata (sampleid/Subtype/Age/Sex/Location/treatment columns); paper Table S1 (mmc2.xlsx) uses different patient numbering and is kept in data/interim/labels/labelle_sample_labels.csv",
        notes=f"files: {','.join(files)}; GEO Subtype={sub.Subtype}; H3.1 vs H3.3 not reported; 'sort' column distinguishes CD45/CD3 sorts in immune_ss2")
add(sample_id="labelle__GSM8389375_H3K27M_P6_Region_1", patient_id="labelle_xenium_H3K27M_P6", study="labelle", role="spatial", assay="xenium",
    chemistry="10x Xenium in situ", site="unknown (midline)", h3_variant="unknown", k27m_status="K27M", treatment="unknown", tissue_timing="unknown",
    include_in_atlas=False, exclusion_reason="imaging-based spatial (Xenium) - not part of the snRNA atlas; available for spatial validation", counts_type="umi",
    n_cells_deposited=mtx_dims(f"{RAW}/labelle/GSE271936/GSM8389375_H3K27M_P6_Region_1-matrix.mtx.gz")[1], accession="GSE271936/GSM8389375",
    file_paths=rel(f"{RAW}/labelle/GSE271936"), label_source="GEO sample title", notes="GSE271936 = 4 Xenium regions (H3K27M, PFA, ST-ZFTA/RELA, H3WT hemispheric); only the H3K27M region's cell x gene matrix fetched; morphology/transcripts skipped")

# ---------------------------------------------------------------- placeholders for the other track
for sid, study, role, assay, acc, note in [
    ("huchede__GSE268577", "huchede", "spatial", "visium", "GSE268577", "3 pontine H3.3K27M sections"),
    ("ren__GSE194329", "ren", "spatial", "visium", "GSE194329", "5 DMG + 5 GBM + 1 peritumour sections"),
    ("collot_st__zenodo17922773", "collot_st", "spatial", "visium", "zenodo.17922773 ST.tar.gz", "human + mouse sections; keep human"),
    ("brainstem__GSE280990", "brainstem_dipg", "spatial", "visium", "GSE280990", "1 DIPG autopsy, 10 sections along the neuraxis"),
    ("braun2023__hindbrain", "braun_2023", "reference", "scRNA", "cellxgene (10.1126/science.adf1226)", "first-trimester hindbrain subset"),
    ("siletti2023__pons", "siletti_2023", "reference", "snRNA", "cellxgene (10.1126/science.add7046)", "adult pons subset"),
    ("openpedcan_v14__dmg", "openpedcan", "bulk", "bulk", "OpenPedCan v14 / PBTA", "DMG subset, TPM + SNV; check overlap with ScPCA CBTN IDs and DeSisto"),
]:
    add(sample_id=sid, patient_id="", study=study, role=role, assay=assay, accession=acc, include_in_atlas=role == "reference",
        exclusion_reason="" if role == "reference" else "placeholder - per-sample rows to be filled by the spatial/reference track",
        counts_type="umi" if assay != "bulk" else "tpm", label_source="placeholder", notes=note, species="human")

man = pd.DataFrame(rows, columns=COLS)
assert man.sample_id.is_unique, man.sample_id[man.sample_id.duplicated()].tolist()
man.to_csv(OUT, index=False)
print(man.groupby(["study", "role"]).agg(n=("sample_id", "size"), included=("include_in_atlas", lambda s: int(s.astype(bool).sum()))).to_string())
print("rows", len(man), "->", OUT)
