"""Phase 8 / step 27b - make dmg_atlas_v1.h5ad fit to deposit.

The release built by 27_package_release.py inherited the working object's internal state:
  * `*_pre_audit` obs columns and obsm embeddings - rollback copies kept so the audit could compare the
    pre- and post-remediation labels. They are provenance for US, noise for a depositor.
  * `author_*` columns that are SAMPLE-level sequencing summaries broadcast onto every cell (raw reads,
    duplicate reads, ATAC fragment counts). A column that is constant within every sample is not cell
    metadata; it belongs in the sample manifest, which is released alongside.
Both are dropped here. Author ANNOTATIONS (cell types, malignant calls, cNMF programmes) and genuinely
per-cell QC (counts, genes, mitochondrial fraction) are kept - they are what makes the atlas auditable.

The pre-remediation release stays on disk as dmg_atlas_v1.pre_audit.h5ad, and atlas_core.h5ad keeps every
rollback array, so nothing is lost. Re-run scripts/29_validate_release.py afterwards.
"""
import hashlib, json, pathlib
import anndata as ad, numpy as np, pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "data/processed/dmg_atlas_v1.h5ad"
OUT = ROOT / "data/processed/dmg_atlas_v1.h5ad"
KEEP_AUTHOR_SUBSTRINGS = ("cell_type", "celltype", "cellType", "Malignant", "call", "cNMF", "Phase",
                          "Score", "score", "qc_pass", "nCount_RNA", "nFeature_RNA", "percent.mito",
                          "percent.mt", "percent.ribo", "leiden", "seurat_clusters", "snn_res",
                          "ID_paper", "Technology", "Study", "SampleID", "sample", "orig.ident",
                          "location", "treatement", "age", "sex")

a = ad.read_h5ad(SRC)
before = dict(n_obs_cols=a.obs.shape[1], obsm=list(a.obsm))

pre_obs = [c for c in a.obs.columns if c.endswith("_pre_audit")]
pre_obsm = [k for k in a.obsm if k.endswith("_pre_audit")]

# sample-level author columns: constant within every sample_id, and not an annotation we want to keep
author = [c for c in a.obs.columns if c.startswith("author_")]
const_in_sample = []
for c in author:
    if any(s in c for s in KEEP_AUTHOR_SUBSTRINGS):
        continue
    nun = a.obs.groupby("sample_id", observed=True)[c].nunique(dropna=False)
    if int(nun.max()) <= 1:
        const_in_sample.append(c)

a.obs = a.obs.drop(columns=pre_obs + const_in_sample)
for k in pre_obsm:
    del a.obsm[k]

a.uns["release_cleaning_27b"] = {
    "dropped_pre_audit_obs": pre_obs,
    "dropped_pre_audit_obsm": pre_obsm,
    "dropped_sample_level_author_obs": const_in_sample,
    "rationale": ("pre-audit arrays are internal rollback state (retained in data/processed/atlas_core.h5ad and in "
                  "dmg_atlas_v1.pre_audit.h5ad); the dropped author_* columns are constant within every sample and "
                  "are released per sample in sample_manifest.csv instead of per cell"),
}
a.write_h5ad(OUT, compression="gzip", compression_opts=4)

md5 = hashlib.md5(OUT.read_bytes()).hexdigest()
rep = dict(before=before, after=dict(n_obs_cols=a.obs.shape[1], obsm=list(a.obsm)),
           dropped_pre_audit_obs=len(pre_obs), dropped_pre_audit_obsm=len(pre_obsm),
           dropped_sample_level_author_obs=len(const_in_sample), md5=md5,
           size_bytes=OUT.stat().st_size)
(ROOT / "results/tables/release_cleaning_27b.json").write_text(json.dumps(rep, indent=1))
pd.DataFrame({"column": pre_obs + const_in_sample,
              "kind": ["pre_audit"] * len(pre_obs) + ["sample_level_author"] * len(const_in_sample)}
             ).to_csv(ROOT / "results/tables/release_dropped_columns.csv", index=False)
print(json.dumps(rep, indent=1))
