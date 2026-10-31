#!/usr/bin/env python
"""Phase 8 / step 27 — release packaging: dmg_atlas_v1.h5ad in CELLxGENE schema 7.1.0 layout + release checksum manifest.

Inputs : data/processed/atlas_core.h5ad (layers['counts'], obs, obsm, uns; read with h5py, counts streamed in row blocks),
         data/interim/gene_universe.csv (canonical symbol -> Ensembl ID), data/processed/metaprogram_scores.parquet +
         results/tables/metaprogram_families.csv (19_metaprogram_scores.py), data/processed/tme_subcluster_labels_clean.csv.gz +
         results/tables/tme_subclusters_clean.csv (22_tme_scvi_clean.py), GENCODE v48 primary-assembly GTF (schema 7.1.0 gene
         reference; downloaded into data/raw/gencode/ if missing).
Outputs (never overwritten if present, unless --force):
  data/processed/dmg_atlas_v1.h5ad        reference cells x GENCODE-v48 genes. X = raw UMI counts as float32 CSR (gzip); no raw / layers
                                          (CELLxGENE puts raw counts in X when no normalised matrix is given). var index = versionless
                                          Ensembl ID, var['gene_symbol'] = atlas symbol. obs = every atlas_core column (author_* verbatim)
                                          + CELLxGENE fields + metaprogram family + cleaned TME subcluster; atlas 'assay' / 'sex' are
                                          renamed assay_atlas / sex_atlas (reserved names). obsm and uns as atlas_core + CELLxGENE uns keys.
  results/tables/release_ontology_map.csv     field, atlas value, ontology term ID, term label, n_cells
  results/tables/release_genes_dropped.csv    atlas genes whose v32 Ensembl ID is absent from GENCODE v48: remapped v48 ID or dropped
  results/tables/release_manifest.csv         every release file: path, role, producing script, size, md5
  docs/release/release_v1.md                          generated release notes (field mapping, new obs fields, checksums, validation status)
  README.md                                   'Release v1' block appended (existing text untouched)
Run: nice -n 10 python scripts/27_package_release.py [--force] [--steps h5ad manifest docs]

Method (fixed before running)
  Ontology terms were checked against OLS4 (EBI) on 2026-10-18; labels are hard-coded below.
  assay: by chemistry; '10x (version not stated)' -> EFO:0030080 10x transcription profiling (3'/5' not stated in the deposits).
  suspension_type: scRNA -> cell; snRNA / multiome_RNA -> nucleus. tissue_type 'tissue'; tissue by site (pons, thalamus).
  disease: K27M -> MONDO:0957196; Okada 'K27-altered (series scope)' -> MONDO:1060171. self_reported_ethnicity 'unknown'.
  development_stage: age >= 1 y -> HsapDv '<n>-year-old stage' (floor); 0 < age < 1 -> '<n>-month-old stage' (floor of age x 12);
  age 0.0 (two Okada pons samples; recorded as 0 in the deposit, not plausible as newborn) and missing age -> 'pediatric stage' /
  'adult stage' from age_group, else 'unknown'.
  cell_type: cleaned TME subcluster (step 22, curated names from docs/phases/phase6_notes.md section 11) where available, else scanvi_label;
  scANVI 'vascular' without a subcluster (endothelial + pericyte mixture) -> 'unknown'. cell_type_atlas holds the label used.
  Audit P7-05 (2026-10-20): 'malignant' may only stand where consensus_call == 'malignant'. A cell that only scANVI called malignant
  (consensus 'ambiguous' or 'normal') is released as 'unspecified_cell' -> CL:0000000 'cell', flagged in cell_type_downgraded_p7_05,
  with malignant_call (three-tier consensus) and malignant_call_eventtrust alongside cnv_confidence / low_confidence_sample so either
  definition can be re-derived. qc_pass = False cells are kept and counted in the log and release docs, not silently dropped.
  is_primary_data False for every cell (meta-analysis re-using published deposits). donor_id = patient_id.
  Genes: Ensembl IDs of the atlas universe (GRCh38-2020-A, GENCODE v32) absent from GENCODE v48 are remapped to the v48 ID whose
  gene_name equals the atlas symbol when that match is unique and the ID is not already used (re-annotated genes, e.g. SOD2);
  otherwise dropped. Added after the first run dropped SOD2, PRRC2B, HOMEZ, PAXX. The cell count is unchanged (cells that would lose
  every count are reported; none are expected).
  Metaprogram fields (malignant cells selected in step 19 only; others 'not_selected' / NaN): mp_family (argmax family, step 19),
  mp_top, mp_max_z, mp_confident (max z >= 1, the cell2location reference subset), mp_z_<family>.
"""
import sys, gzip, json, time, hashlib, argparse, urllib.request
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import h5py, numpy as np, pandas as pd, anndata as ad  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import ROOT, RAW, INTERIM, PROCESSED, RESULTS, LOGS  # noqa: E402

log = U.StepLogger(LOGS / "phase8_release.log", context="27_package_release")
CORE = PROCESSED / "atlas_core.h5ad"
OUT = PROCESSED / "dmg_atlas_v1.h5ad"
TAB = RESULTS / "tables"
GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz"
GTF = RAW / "gencode" / "gencode.v48.primary_assembly.annotation.gtf.gz"
SCHEMA = "7.1.0"
TITLE = "Cross-cohort single-cell atlas of H3K27-altered diffuse midline glioma — integrated reference (v1)"
ROW_BLOCK = 20_000

ASSAY = {"10x 3' v2": ("EFO:0009899", "10x 3' v2"), "10x 3' v3": ("EFO:0009922", "10x 3' v3"),
         "10x Multiome (ARC) NEXT-GEM v1": ("EFO:0030059", "10x multiome"),
         "10x (version not stated)": ("EFO:0030080", "10x transcription profiling")}
SUSPENSION = {"scRNA": "cell", "snRNA": "nucleus", "multiome_RNA": "nucleus"}
TISSUE = {"pons": ("UBERON:0000988", "pons"), "thalamus": ("UBERON:0001897", "dorsal plus ventral thalamus")}
DISEASE = {"K27M": ("MONDO:0957196", "diffuse midline glioma, H3 K27M-mutant"),
           "K27-altered (series scope)": ("MONDO:1060171", "diffuse midline glioma, H3 K27-altered")}
SEX = {"F": ("PATO:0000383", "female"), "M": ("PATO:0000384", "male"), "unknown": ("unknown", "unknown")}
STAGE_GROUP = {"pediatric": ("HsapDv:0000264", "pediatric stage"), "adult": ("HsapDv:0000258", "adult stage")}
CELL_TYPE = {  # atlas label -> CL
    "malignant": ("CL:0001064", "malignant cell"), "myeloid": ("CL:0000766", "myeloid leukocyte"), "T_NK": ("CL:0000542", "lymphocyte"),
    "astrocyte": ("CL:0000127", "astrocyte"), "oligodendrocyte": ("CL:0000128", "oligodendrocyte"), "neuron": ("CL:0000540", "neuron"),
    "vascular": ("unknown", "unknown"),
    # audit P7-05: cells scANVI calls malignant but the CNV/marker consensus does not confirm get the unspecific term, not 'malignant cell'
    "unspecified_cell": ("CL:0000000", "cell"),
    "homeostatic_microglia": ("CL:0000129", "microglial cell"), "TGFBI_AOAH_macrophage": ("CL:0000235", "macrophage"),
    "hypoxic_glycolytic_TAM": ("CL:0000235", "macrophage"), "interferon_TAM": ("CL:0000235", "macrophage"),
    "ribosomal_residue_myeloid": ("CL:0000766", "myeloid leukocyte"), "proliferating_myeloid": ("CL:0000766", "myeloid leukocyte"),
    "THEMIS_CAMK4_T": ("CL:0000084", "T cell"), "CD8_T": ("CL:0000625", "CD8-positive, alpha-beta T cell"),
    "CD4_T": ("CL:0000624", "CD4-positive, alpha-beta T cell"), "NK": ("CL:0000623", "natural killer cell"),
    "IL7R_KLRB1_T": ("CL:0000084", "T cell"), "endothelial_neural_mixed": ("CL:0000115", "endothelial cell"),
    "endothelial_tip_angiogenic": ("CL:0000115", "endothelial cell"), "endothelial_capillary": ("CL:0000115", "endothelial cell"),
    "pericyte": ("CL:0000669", "pericyte"), "mature_oligodendrocyte": ("CL:0000128", "oligodendrocyte"),
    "committed_OPC": ("CL:0002453", "oligodendrocyte precursor cell")}
TME_NAME = {  # cleaned subcluster -> curated name (docs/phases/phase6_notes.md section 11; top genes in tme_subclusters_clean.csv)
    "myeloid_c0": "homeostatic_microglia", "myeloid_c1": "TGFBI_AOAH_macrophage", "myeloid_c2": "hypoxic_glycolytic_TAM",
    "myeloid_c3": "interferon_TAM", "myeloid_c4": "ribosomal_residue_myeloid", "myeloid_c5": "proliferating_myeloid",
    "T_NK_c0": "THEMIS_CAMK4_T", "T_NK_c1": "CD8_T", "T_NK_c2": "CD4_T", "T_NK_c3": "NK", "T_NK_c4": "IL7R_KLRB1_T",
    "vascular_c0": "endothelial_neural_mixed", "vascular_c1": "endothelial_tip_angiogenic", "vascular_c2": "endothelial_capillary",
    "vascular_c3": "endothelial_capillary", "vascular_c4": "pericyte", "oligodendrocyte_c0": "mature_oligodendrocyte",
    "oligodendrocyte_c1": "mature_oligodendrocyte", "oligodendrocyte_c2": "mature_oligodendrocyte",
    "oligodendrocyte_c3": "mature_oligodendrocyte", "oligodendrocyte_c4": "committed_OPC"}
RESERVED_RENAME = {"assay": "assay_atlas", "sex": "sex_atlas"}
CXG_RESERVED = {"assay", "cell_type", "development_stage", "disease", "experimental_condition", "observation_joinid",
                "perturbation_types", "self_reported_ethnicity", "sex", "tissue", "organism", "organism_ontology_term_id",
                "ethnicity", "ethnicity_ontology_term_id"}
MODEL_ROLES = {"scanvi": "reference model", "scvi_seed1": "reference model (chosen scVI, uns phase4/chosen_scvi)",
               "scvi_malignant": "reference model", "c2l_regression": "reference model", "c2l_regression_mp": "reference model",
               "scvi_seed0": "control", "scvi_seed2": "control", "scvi_seed0_e60": "control", "scvi_malignant_seed0": "control",
               "harmony": "control", "pca": "control"}


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def model_script(name: str) -> str:
    if name.startswith("c2l_regression"):
        return "modal/phase6_c2l.py"
    if name.startswith("scanvi_query_"):
        return "11_map_queries.py"
    if name.startswith("scvi_tme_"):
        return "22_tme_scvi_clean.py"
    return "modal/phase4_train.py"


def stage(age: float, group: str) -> tuple:
    if not np.isfinite(age) or age == 0:
        return STAGE_GROUP.get(group, ("unknown", "unknown"))
    if age < 1:
        m = int(age * 12)
        return (f"HsapDv:{173 + m:07d}", f"{m}-month-old stage") if m >= 1 else ("HsapDv:0000262", "newborn stage (0-28 days)")
    y = int(age)
    return ("HsapDv:0000246", "1-year-old stage") if y == 1 else (f"HsapDv:{94 + y:07d}", f"{y}-year-old stage")


def gencode_genes() -> pd.Series:  # versionless gene_id -> gene_name
    if not GTF.exists():  # resumable download into .part (the EBI FTP mirror drops long transfers); renamed when complete
        GTF.parent.mkdir(parents=True, exist_ok=True)
        part = GTF.with_suffix(".gz.part")
        total = int(urllib.request.urlopen(urllib.request.Request(GTF_URL, method="HEAD"), timeout=60).headers["Content-Length"])
        for attempt in range(20):
            have = part.stat().st_size if part.exists() else 0
            if have >= total:
                break
            log.write(f"downloading {GTF_URL} from byte {have:,} / {total:,} (attempt {attempt + 1})")
            try:
                req = urllib.request.Request(GTF_URL, headers={"Range": f"bytes={have}-"})
                with urllib.request.urlopen(req, timeout=120) as r, open(part, "ab") as fh:
                    for b in iter(lambda: r.read(1 << 20), b""):
                        fh.write(b)
            except OSError as e:
                log.write(f"download interrupted: {e}")
        assert part.stat().st_size == total, "GENCODE GTF download incomplete"
        part.rename(GTF)
    ids = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line[0] != "#":
                f = line.split("\t", 9)
                if f[2] == "gene":
                    gid = f[8].split('gene_id "', 1)[1].split('"', 1)[0].split(".", 1)[0]
                    ids[gid] = f[8].split('gene_name "', 1)[1].split('"', 1)[0] if 'gene_name "' in f[8] else gid
    return pd.Series(ids)


def build_obs(obs: pd.DataFrame) -> tuple:
    clash = CXG_RESERVED & (set(obs.columns) - set(RESERVED_RENAME))
    assert not clash, f"reserved CELLxGENE obs names present: {clash}"
    obs = obs.rename(columns=RESERVED_RENAME)
    unmapped = {c: sorted(set(obs[c].astype(str)) - set(d)) for c, d in
                [("chemistry", ASSAY), ("assay_atlas", SUSPENSION), ("site", TISSUE), ("k27m_status", DISEASE), ("sex_atlas", SEX)]}
    assert not any(unmapped.values()), f"unmapped atlas values: {unmapped}"
    new, maps = {}, []

    def put(field, src, mapping):
        s = obs[src].astype(str)
        new[field] = s.map(lambda v: mapping[v][0])
        for v, n in s.value_counts().items():
            maps.append(dict(field=field, atlas_column=src, atlas_value=v, term_id=mapping[v][0], term_label=mapping[v][1], n_cells=n))

    put("assay_ontology_term_id", "chemistry", ASSAY)
    put("tissue_ontology_term_id", "site", TISSUE)
    put("disease_ontology_term_id", "k27m_status", DISEASE)
    put("sex_ontology_term_id", "sex_atlas", SEX)
    new["suspension_type"] = obs["assay_atlas"].astype(str).map(SUSPENSION)
    new["tissue_type"] = "tissue"
    new["self_reported_ethnicity_ontology_term_id"] = "unknown"
    new["donor_id"] = obs["patient_id"].astype(str)
    new["is_primary_data"] = False
    st = [stage(a, g) for a, g in zip(obs["age_years"].astype(float), obs["age_group"].astype(str))]
    new["development_stage_ontology_term_id"] = [s[0] for s in st]
    sg = pd.DataFrame({"age": obs["age_years"].astype(float).fillna(-1).values, "term_id": new["development_stage_ontology_term_id"],
                       "term_label": [s[1] for s in st]}).value_counts().reset_index(name="n_cells")
    maps += [dict(field="development_stage_ontology_term_id", atlas_column="age_years / age_group",
                  atlas_value="missing" if r.age < 0 else f"{r.age:g}", term_id=r.term_id, term_label=r.term_label, n_cells=r.n_cells)
             for r in sg.itertuples()]

    # metaprogram families (step 19)
    ms = pd.read_parquet(PROCESSED / "metaprogram_scores.parquet")
    assert ms.index.isin(obs.index).all(), "metaprogram_scores cells missing from atlas_core"
    zc = [c for c in ms.columns if c.startswith("z_")]
    ms = ms.reindex(obs.index)
    new["malignant_selected"] = ms["metaprogram"].notna().values
    new["mp_family"] = ms["metaprogram"].fillna("not_selected").astype(str).values
    new["mp_top"] = ms["top_mp"].fillna("not_selected").astype(str).values
    new["mp_max_z"] = ms[zc].max(axis=1).astype("float32").values
    new["mp_confident"] = (ms[zc].max(axis=1) >= 1).values
    for c in zc:
        new[f"mp_{c}"] = ms[c].astype("float32").values

    # cleaned TME subclusters (step 22)
    t = pd.read_csv(PROCESSED / "tme_subcluster_labels_clean.csv.gz").set_index("cell_id")
    assert t.index.isin(obs.index).all() and set(t.subcluster) <= set(TME_NAME), "TME labels do not match atlas / curated names"
    t = t.reindex(obs.index)
    new["tme_compartment"] = t["compartment"].fillna("none").astype(str).values
    new["tme_subcluster"] = t["subcluster"].fillna("none").astype(str).values
    new["tme_marker_label"] = t["label"].fillna("none").astype(str).values
    new["tme_annotation"] = t["subcluster"].map(TME_NAME).fillna("none").astype(str).values
    base = pd.Series(np.where(new["tme_annotation"] != "none", new["tme_annotation"], obs["scanvi_label"].astype(str)), index=obs.index)
    # audit P7-05: 'malignant cell' (CL:0001064) may only stand where the three-way consensus confirms it. Cells that only scANVI
    # called malignant (consensus 'ambiguous' or 'normal') are released as CL:0000000 'cell'; the evidence travels in malignant_call,
    # malignant_call_eventtrust, cnv_confidence and low_confidence_sample so a user can re-derive either definition.
    cc = obs["consensus_call"].astype(str)
    downgraded = (base == "malignant") & (cc != "malignant")
    ct = base.mask(downgraded, "unspecified_cell")
    new["malignant_call"] = cc.values
    new["malignant_call_eventtrust"] = obs["consensus_call_eventtrust"].astype(str).values
    new["cell_type_downgraded_p7_05"] = downgraded.values
    assert set(ct) <= set(CELL_TYPE), f"unmapped cell types: {set(ct) - set(CELL_TYPE)}"
    new["cell_type_atlas"] = ct.values
    new["cell_type_ontology_term_id"] = ct.map(lambda v: CELL_TYPE[v][0]).values
    for v, n in ct.value_counts().items():
        maps.append(dict(field="cell_type_ontology_term_id", atlas_column="tme_annotation | scanvi_label (P7-05: gated on consensus_call)",
                         atlas_value=v, term_id=CELL_TYPE[v][0], term_label=CELL_TYPE[v][1], n_cells=n))
    log.write(f"P7-05 cell_type gating: scANVI-malignant cells not confirmed by consensus_call -> 'unspecified_cell' (CL:0000000): "
              f"{int(downgraded.sum()):,} cells (consensus mix: {cc[downgraded].value_counts().to_dict()}); "
              f"CL:0001064 retained for {int((ct == 'malignant').sum()):,} cells; "
              f"qc_pass = False cells kept and documented: {int((~obs['qc_pass'].astype(bool)).sum()):,}")
    xt = pd.crosstab(obs["scanvi_label"].astype(str), new["tme_compartment"])
    log.write("scanvi_label x tme_compartment:\n" + xt.to_string())

    new = pd.DataFrame(new, index=obs.index)
    for c in new.columns:
        if new[c].dtype == object:
            new[c] = new[c].astype("category")
    assert not set(new.columns) & set(obs.columns)
    return pd.concat([obs, new], axis=1), pd.DataFrame(maps)


def prune(d):  # CELLxGENE: uns values must not be empty
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = prune(v)
        if v is None or isinstance(v, (bool, np.bool_)) or (np.size(v) if not isinstance(v, (dict, str)) else len(v)) > 0:
            out[k] = v
    return out


def write_h5ad():
    with log.step("obs / var"):
        with h5py.File(CORE, "r") as f:
            obs = read_elem(f["obs"]); var = read_elem(f["var"]); uns = read_elem(f["uns"])
            obsm = {k: f["obsm"][k][:] for k in f["obsm"]}
            shape = tuple(f["layers/counts"].attrs["shape"])
        obs, maps = build_obs(obs)
        gu = pd.read_csv(INTERIM / "gene_universe.csv").set_index("canonical_symbol")["ensembl_id"]
        ens = gu.reindex(var.index)
        assert ens.notna().all() and not ens.duplicated().any(), "atlas symbols do not map 1:1 to Ensembl IDs"
        g48 = gencode_genes()
        by_name = g48.reset_index().groupby(0)["index"].agg(list)
        miss = ~ens.isin(g48.index)
        rescue = {s: by_name[s][0] for s in ens.index[miss] if s in by_name.index and len(by_name[s]) == 1}
        rescue = {s: e for s, e in rescue.items() if e not in set(ens.values)}  # never collide with an ID already in the atlas
        ens_v48 = ens.copy(); ens_v48[list(rescue)] = list(rescue.values())
        keep = ens_v48.isin(g48.index).values
        pd.DataFrame({"gene_symbol": var.index[miss.values], "ensembl_id_v32": ens[miss].values,
                      "ensembl_id_v48": [rescue.get(s, "") for s in var.index[miss.values]],
                      "action": ["remapped by unique v48 gene_name" if s in rescue else "dropped" for s in var.index[miss.values]]}
                     ).to_csv(TAB / "release_genes_dropped.csv", index=False)
        log.write(f"genes: {shape[1]:,} atlas; {int(miss.sum())} IDs absent from GENCODE v48: {len(rescue)} remapped by unique gene_name "
                  f"({', '.join(rescue)}), {(~keep).sum()} dropped -> {keep.sum():,} genes")
        ens = ens_v48
        v = var.loc[keep].copy()
        v.insert(0, "gene_symbol", v.index.astype(str))
        v.index = pd.Index(ens[keep].values, name="ensembl_id")
        v["feature_is_filtered"] = False
        maps.to_csv(TAB / "release_ontology_map.csv", index=False)
        uns = prune(uns)
        uns.update(title=TITLE, organism_ontology_term_id="NCBITaxon:9606", batch_condition=["sample_id"], default_embedding="X_umap",
                   release=dict(version="v1", date=time.strftime("%Y-%m-%d"), target_schema=f"CELLxGENE {SCHEMA}",
                                source=f"atlas_core.h5ad md5 {md5(CORE)}", gene_reference="GENCODE v48 (v32 IDs absent from v48 remapped by unique gene name, else dropped)",
                                X="raw UMI counts (float32); X_scvi / X_scanvi / X_umap from the Phase 4 reference models",
                                producer="scripts/27_package_release.py"))
    with log.step("write skeleton (obs, var, obsm, uns)"):
        adata = ad.AnnData(obs=obs, var=v, obsm=obsm, uns=uns)
        tmp = OUT.with_suffix(".tmp.h5ad")
        adata.write_h5ad(tmp, compression="gzip")
        del adata, obsm
    with log.step("stream X (counts -> float32 CSR, GENCODE v48 genes)"):
        newcol = np.full(shape[1], -1, np.int64); newcol[keep] = np.arange(keep.sum())
        with h5py.File(CORE, "r") as src, h5py.File(tmp, "a") as dst:
            c = src["layers/counts"]; indptr = c["indptr"][:]
            if "X" in dst:
                del dst["X"]
            g = dst.create_group("X")
            g.attrs.update({"encoding-type": "csr_matrix", "encoding-version": "0.1.0", "shape": np.array([shape[0], int(keep.sum())])})
            kw = dict(compression="gzip", compression_opts=4, chunks=(1 << 20,), maxshape=(None,))
            d = g.create_dataset("data", (0,), dtype="float32", **kw)
            ix = g.create_dataset("indices", (0,), dtype="int32", **kw)
            out_ptr = np.zeros(shape[0] + 1, np.int64); n, empty = 0, 0
            for r0 in range(0, shape[0], ROW_BLOCK):
                r1 = min(r0 + ROW_BLOCK, shape[0]); a, b = indptr[r0], indptr[r1]
                data = c["data"][a:b]; cols = newcol[c["indices"][a:b]]; m = cols >= 0
                assert (data > 0).all() and np.issubdtype(data.dtype, np.integer), "counts are not positive integers"
                rows = np.repeat(np.arange(r1 - r0), np.diff(indptr[r0:r1 + 1]))[m]
                per = np.bincount(rows, minlength=r1 - r0); empty += int((per == 0).sum())
                out_ptr[r0 + 1:r1 + 1] = n + np.cumsum(per)
                k = int(m.sum()); d.resize((n + k,)); ix.resize((n + k,))
                d[n:n + k] = data[m].astype("float32"); ix[n:n + k] = cols[m].astype("int32"); n += k
            g.create_dataset("indptr", data=out_ptr)
        log.write(f"X: {shape[0]:,} x {int(keep.sum()):,}, nnz {n:,} (source {int(indptr[-1]):,}); cells with no counts after gene drop: {empty}")
        assert empty == 0, "cells without counts after gene drop (CELLxGENE requires >= 1 non-zero per cell)"
        tmp.rename(OUT)
    with log.step("re-read check"):
        a = ad.read_h5ad(OUT, backed="r")
        x = a.X[:1000]
        log.write(f"re-read: {a.n_obs:,} x {a.n_vars:,}; X dtype {x.dtype}; obs cols {a.obs.shape[1]}; obsm {list(a.obsm)}; "
                  f"cell types {a.obs.cell_type_ontology_term_id.nunique()}; first-1000 row sum {float(x.sum()):.0f}")
        a.file.close()


def write_manifest() -> pd.DataFrame:
    rows = []

    # audit P7-13: results tables and figures previously carried no producing script. Attribute each file to the script whose
    # source mentions its filename literally; leave blank when zero or more than one script matches, so the column never guesses.
    # Packaging scripts are excluded as candidates: they LIST output filenames in their manifests without producing them
    # (e.g. 12_package_v1.py is the only script mentioning phase5_query_mapping.png, which it does not write).
    _src = {f.name: f.read_text(errors="ignore") for f in sorted((ROOT / "scripts").glob("*.py"))
            if ".pre_audit." not in f.name and "_package" not in f.name}

    def producing_script(p: Path) -> str:
        stem = p.name.replace(".pre_audit", "")
        hits = sorted(n for n, s in _src.items() if stem in s)
        return hits[0] if len(hits) == 1 else ""

    def add(p: Path, role: str, script: str):
        rows.append(dict(path=str(p.relative_to(ROOT)), role=role, script=script, size_mb=round(p.stat().st_size / 1e6, 3), md5=md5(p)))

    with log.step("release manifest (md5)"):
        for name, script, role in [("dmg_atlas_v1.h5ad", "27_package_release.py", "release atlas (CELLxGENE layout)"),
                                   ("atlas_core.h5ad", "09_assemble_atlas_core.py", "pipeline atlas (source of the release)"),
                                   ("query_mapping.h5ad", "11_map_queries.py", "query projections"),
                                   ("sample_manifest.csv", "02_build_sample_manifest.py", "sample metadata"),
                                   ("malignant_calls.csv", "06_cnv_aggregate.py", "malignant calls"),
                                   ("metaprogram_scores.parquet", "19_metaprogram_scores.py", "per-cell metaprogram scores"),
                                   ("tme_subcluster_labels_clean.csv.gz", "22_tme_scvi_clean.py", "per-cell TME subclusters"),
                                   ("foundation_embeddings.npz", "26_foundation_embeddings.py", "embeddings (methods note)")]:
            add(PROCESSED / name, role, script)
        for p in sorted((PROCESSED / "spatial").glob("spatial_deconvolution_*.h5ad")):
            add(p, "spatial deconvolution", "14_cell2location.py")
        for d in sorted((ROOT / "models").iterdir()):
            role = MODEL_ROLES.get(d.name, "query surgery model" if d.name.startswith("scanvi_query_")
                                   else "TME scVI model" if d.name.startswith("scvi_tme_") else None)
            if role is None:  # foundation (third-party weights cache), scanvi_nspl100_rejected
                continue
            for p in sorted(x for x in d.rglob("*") if x.is_file() and x.name != ".DS_Store"):
                add(p, role, model_script(d.name))
        for p in ["docs/design/decisions.md", "docs/design/build_plan.md", "docs/schema/atlas_schema.md", "docs/schema/label_ontology.csv",
                  "results/tables/release_ontology_map.csv", "results/tables/release_genes_dropped.csv"]:
            add(ROOT / p, "documentation", "27_package_release.py" if "release_" in p else producing_script(ROOT / p))
        for p in sorted((RESULTS / "tables").rglob("*")):
            if p.is_file() and p.name != "release_manifest.csv" and not p.name.startswith("release_"):
                add(p, "results table", producing_script(p))
        for p in sorted((RESULTS / "figures").glob("*.png")):
            add(p, "figure", producing_script(p))
        m = pd.DataFrame(rows)
        m.to_csv(TAB / "release_manifest.csv", index=False)
        log.write(f"release_manifest.csv: {len(m)} files, {m.size_mb.sum() / 1e3:.2f} GB")
    return m


def write_docs(m: pd.DataFrame):
    with log.step("docs/release/release_v1.md + README block"):
        maps = pd.read_csv(TAB / "release_ontology_map.csv"); dropped = pd.read_csv(TAB / "release_genes_dropped.csv")
        with h5py.File(OUT, "r") as f:
            n_obs, n_var = (int(x) for x in f["X"].attrs["shape"]); obs_cols = [c for c in f["obs"].keys() if c != "_index"]
            obsm = list(f["obsm"].keys())
        rel = m.set_index("path")
        lines = [f"# Release v1 — generated by scripts/27_package_release.py (regenerate rather than hand-edit)", "",
                 f"Generated {time.strftime('%Y-%m-%d %H:%M')}. Target: CELLxGENE schema {SCHEMA} (not yet validated with `cellxgene-schema`).", "",
                 "## dmg_atlas_v1.h5ad", "",
                 f"- {n_obs:,} cells / nuclei x {n_var:,} genes (GENCODE v48 versionless Ensembl IDs; of {len(dropped)} atlas IDs absent from v48, "
                 f"{(dropped.action != 'dropped').sum()} remapped by unique gene name and {(dropped.action == 'dropped').sum()} dropped: "
                 "`results/tables/release_genes_dropped.csv`).",
                 "- `X` = raw UMI counts, float32 CSR (no normalised matrix, no `raw`, no layers). `var['gene_symbol']` = atlas symbol, "
                 "`var['feature_is_filtered']` = False.",
                 f"- obsm: {', '.join(obsm)} (X_scvi_malignant is NaN outside the malignant subset). uns: title, organism_ontology_term_id, "
                 "batch_condition ['sample_id'], default_embedding X_umap, release, plus the atlas_core uns.",
                 f"- md5 `{rel.loc['data/processed/dmg_atlas_v1.h5ad', 'md5']}`; source atlas_core.h5ad md5 "
                 f"`{rel.loc['data/processed/atlas_core.h5ad', 'md5']}`.", "",
                 "## CELLxGENE field mapping", "",
                 "Fixed: tissue_type 'tissue'; self_reported_ethnicity_ontology_term_id 'unknown'; is_primary_data False (meta-analysis of "
                 "published deposits); donor_id = patient_id; suspension_type cell (scRNA) / nucleus (snRNA, multiome). Atlas columns "
                 "`assay` and `sex` are renamed `assay_atlas` / `sex_atlas` (reserved names).", "",
                 "| field | atlas value | term | label | cells |", "|---|---|---|---|---|"]
        lines += [f"| {r.field} | {r.atlas_value} | `{r.term_id}` | {r.term_label} | {r.n_cells:,} |" for r in maps.itertuples()]
        lines += ["", "## Added obs fields", "",
                  "| field | content |", "|---|---|",
                  "| cell_type_atlas | label behind cell_type_ontology_term_id: cleaned TME annotation where available, else scanvi_label, "
                  "with 'malignant' gated on consensus_call (audit P7-05) |",
                  "| malignant_call / malignant_call_eventtrust | three-tier malignancy consensus and its event-trust variant |",
                  "| cell_type_downgraded_p7_05 | scANVI called the cell malignant but the consensus did not: released as CL:0000000 'cell' |",
                  "| malignant_selected | cell in the step-19 malignant selection (metaprogram scoring) |",
                  "| mp_family / mp_top | assigned metaprogram family (argmax z; 'unassigned' if max z <= 0) / highest-scoring MP; 'not_selected' otherwise |",
                  "| mp_max_z / mp_confident / mp_z_<family> | family z-scores; confident = max z >= 1 (cell2location reference subset) |",
                  "| tme_compartment / tme_subcluster / tme_marker_label / tme_annotation | cleaned TME clustering (step 22); annotation = curated name (phase6_notes §11) |",
                  "", f"All {len(obs_cols)} obs columns: " + ", ".join(f"`{c}`" for c in obs_cols), "",
                  "## Release manifest", "",
                  f"`results/tables/release_manifest.csv`: {len(m)} files ({m.size_mb.sum() / 1e3:.2f} GB) with role, producing script and md5.", ""]
        lines += ["| role | files | MB |", "|---|---|---|"]
        lines += [f"| {r} | {len(g)} | {g.size_mb.sum():,.1f} |" for r, g in m.groupby("role", sort=False)]
        lines += ["", "Not in the release: models/foundation (third-party pretrained weights), models/scanvi_nspl100_rejected, "
                  "atlas_input*.h5ad (intermediate).", "",
                  "## Scope of the analytical layers", "",
                  "- Phases 1-6 (acquisition, QC, CNV, integration, query mapping, metaprograms and TME) are the planned atlas build.",
                  "- **Phase 7 (code) / build-plan Phase 8 - multimodal ATAC, bulk survival and foundation-model embeddings - is "
                  "exploratory and was not part of the approved plan** (audit P7-16). Its results are a methods note, not atlas claims.",
                  "- **The spatial layer is pending audit corrections** (track 'spatial', 12 actionable findings): Figure 5 and "
                  "Supplementary Figure 7 are computed from pre-audit spatial tables and are marked as such in their captions.", "",
                  "## Open before submission", "",
                  "- Run `cellxgene-schema validate data/processed/dmg_atlas_v1.h5ad` (package not in the multiverse env).",
                  "- Ensembl IDs come from GRCh38-2020-A (GENCODE v32) symbols; IDs present in v48 are kept as-is, but a symbol may differ "
                  "from the v48 gene name (the portal re-derives feature_name from the ID).",
                  "- '10x (version not stated)' studies (Okada, Collot, Zang) carry the generic EFO:0030080; refine if the chemistry is confirmed.",
                  "- Okada age 0.0 (two pons samples) is treated as unrecorded (pediatric stage).", ""]
        (ROOT / "docs/release/release_v1.md").write_text("\n".join(lines) + "\n")
        txt = (ROOT / "README.md").read_text()
        head = "## Release v1" if "## Release v1" not in txt else f"## Release v1 — update {time.strftime('%Y-%m-%d')}"
        block = ["", head, "",
                 f"`data/processed/dmg_atlas_v1.h5ad` — {n_obs:,} cells x {n_var:,} GENCODE-v48 genes in CELLxGENE schema {SCHEMA} layout (raw "
                 "counts in X; ontology fields, metaprogram families and cleaned TME subclusters added to obs). md5 "
                 f"`{rel.loc['data/processed/dmg_atlas_v1.h5ad', 'md5']}`. Field mapping and open items: `docs/release/release_v1.md`; every release "
                 "file with md5: `results/tables/release_manifest.csv`. Produced by `scripts/27_package_release.py`.", ""]
        with open(ROOT / "README.md", "a") as fh:
            fh.write("\n".join(block))
        log.write(f"docs/release/release_v1.md written; README block '{head}' appended")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--steps", nargs="+", default=["h5ad", "manifest", "docs"], choices=["h5ad", "manifest", "docs"])
    args = ap.parse_args()
    if "h5ad" in args.steps:
        if OUT.exists() and not args.force:
            log.write(f"{OUT.name} exists, skipping (use --force)")
        else:
            write_h5ad()
    man = write_manifest() if "manifest" in args.steps else pd.read_csv(TAB / "release_manifest.csv")
    if "docs" in args.steps:
        write_docs(man)
    log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
