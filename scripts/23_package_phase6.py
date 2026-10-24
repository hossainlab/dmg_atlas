#!/usr/bin/env python
"""Phase 6 / step 23 — document every Phase-6 output: docs/phase6_outputs.md (generated; regenerate rather than hand-edit) and a
'Phase 6 outputs' block appended to README.md (existing text untouched; a dated update block if the section already exists).

docs/atlas_schema.md (12_package_v1.py) describes atlas_core.h5ad / query_mapping.h5ad, which Phase 6 does not modify, so it is not
regenerated. For each Phase-6 file: producing script, size, md5, and its structure (CSV: rows x columns + column names; parquet:
shape + columns; h5ad: shape + obs / var / obsm / layers / uns keys read with h5py; PNG / JSON: size only).
Run: nice -n 10 python scripts/23_package_phase6.py
"""
import sys, json, glob, hashlib, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import pandas as pd, h5py  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
from dmg_atlas import ROOT, LOGS  # noqa: E402

log = U.StepLogger(LOGS / "phase6_package.log", context="23_package_phase6")
GROUPS = [  # (producing script, globs relative to ROOT)
    ("13_cnmf_per_sample.py", ["results/tables/cnmf/k_selection.csv", "results/tables/cnmf/samples_analysed.csv", "data/processed/cnmf_usages_*.parquet"]),
    ("13_reference_signatures.py", ["results/tables/c2l_reference_classes*.csv", "results/tables/c2l_reference_state_scores.csv.gz", "results/tables/c2l_signatures*.csv"]),
    ("14_cell2location.py", ["results/tables/spatial_abundance_long.csv", "results/tables/c2l_deconvolution_runs.json", "data/processed/spatial/spatial_deconvolution_*.h5ad"]),
    ("15_spatial_tests.py", ["results/tables/spatial_colocalisation*.csv", "results/tables/spatial_niches*.csv", "results/tables/brainstem_gradient.csv",
                             "results/tables/brainstem_section_means.csv", "results/tables/ren_*.csv", "data/processed/spatial/spatial_niche_labels.csv.gz",
                             "results/figures/phase6_spatial_abundance_maps.png", "results/figures/phase6_colocalisation_matrix.png", "results/figures/phase6_brainstem_gradient.png"]),
    ("16_metaprograms.py", ["results/tables/metaprograms.csv", "results/tables/metaprogram_genes.csv", "results/tables/metaprogram_members.csv",
                            "results/tables/replication_matrix.csv", "results/tables/metaprogram_reference_overlap.csv", "results/tables/metaprogram_cut_sensitivity.csv",
                            "results/figures/phase6_metaprogram_jaccard.png", "results/figures/phase6_replication_matrix.png"]),
    ("17_tme_subclusters.py", ["results/tables/tme_subclusters.csv", "results/tables/tme_subcluster_markers.csv", "results/tables/tme_subcluster_marker_scores.csv",
                               "data/processed/tme_subcluster_labels.csv.gz", "results/figures/phase6_tme_umap.png", "results/figures/phase6_tme_marker_scores.png"]),
    ("18_composition.py", ["results/tables/composition_*.csv", "results/tables/malignant_state_assignment.csv.gz", "results/figures/phase6_composition_effects.png"]),
    ("19_metaprogram_scores.py", ["data/processed/metaprogram_scores.parquet", "results/tables/metaprogram_families.csv", "results/tables/metaprogram_family_composition.csv",
                                  "results/tables/jessa_cnmf_*.csv", "results/figures/phase6_jessa_cnmf_vs_metaprograms.png"]),
    ("20_assay_masked_metaprograms.py", ["results/tables/assay_biased_genes.csv", "results/tables/metaprograms_assay_masked.csv",
                                         "results/tables/metaprogram_members_assay_masked.csv", "results/tables/metaprogram_masking_summary.csv",
                                         "results/tables/metaprogram_cell_vs_nucleus.csv", "results/figures/phase6_metaprogram_masked_replication.png"]),
    ("21_brainstem_gradient_mp.py", ["results/tables/brainstem_gradient_mp.csv", "results/tables/brainstem_section_means_mp.csv", "results/tables/brainstem_mp_vs_fallback.csv",
                                     "results/tables/c2l_deconvolution_runs_mp.json", "results/figures/phase6_brainstem_gradient_mp.png"]),
    ("22_tme_scvi_clean.py", ["results/tables/tme_clean_exclusions.csv", "results/tables/tme_subclusters_clean.csv", "results/tables/tme_subcluster_markers_clean.csv",
                              "results/tables/tme_assay_purity.csv", "data/processed/tme_subcluster_labels_clean.csv.gz", "results/figures/phase6_tme_umap_clean.png"]),
]


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def describe(p: Path) -> str:
    n = p.name
    if n.endswith((".csv", ".csv.gz")):
        head = pd.read_csv(p, nrows=0).columns.tolist()
        rows = sum(len(c) for c in pd.read_csv(p, usecols=[0], chunksize=500000))
        cols = ", ".join(f"`{c}`" for c in head[:40]) + (f", … (+{len(head) - 40})" if len(head) > 40 else "")
        return f"{rows:,} × {len(head)}: {cols}"
    if n.endswith(".parquet"):
        md = pq.ParquetFile(p).metadata; names = pq.ParquetFile(p).schema_arrow.names
        return f"{md.num_rows:,} × {len(names)}: " + ", ".join(f"`{c}`" for c in names[:25]) + (f", … (+{len(names) - 25})" if len(names) > 25 else "")
    if n.endswith(".h5ad"):
        with h5py.File(p, "r") as f:
            nobs = len(f["obs"][f["obs"].attrs["_index"]]); nvar = len(f["var"][f["var"].attrs["_index"]])
            keys = {k: [x for x in f[k].keys() if x != "_index"] for k in ("obs", "var", "obsm", "layers", "uns") if k in f}
        return f"{nobs:,} × {nvar:,}; " + "; ".join(f"{k}: {', '.join(v)}" for k, v in keys.items() if v)
    if n.endswith(".json"):
        return "JSON: " + ", ".join(json.load(open(p)).keys())
    return "figure"


with log.step("describe outputs"):
    rows, missing = [], []
    for script, pats in GROUPS:
        for pat in pats:
            hits = sorted(glob.glob(str(ROOT / pat)))
            if not hits:
                missing.append(pat)
            for h in hits:
                p = Path(h)
                rows.append(dict(script=script, path=str(p.relative_to(ROOT)), size_mb=p.stat().st_size / 1e6, md5=md5(p), structure=describe(p)))
    out = pd.DataFrame(rows)
    log.write(f"{len(out)} files described; missing patterns: {missing}")

with log.step("write docs/phase6_outputs.md + README block"):
    lines = ["# Phase 6 outputs — generated by scripts/23_package_phase6.py (regenerate rather than hand-edit)", "",
             f"Generated {time.strftime('%Y-%m-%d %H:%M')}. Narrative and results: docs/phase6_notes.md. Per-sample cNMF programme tables "
             "(`results/tables/cnmf/<sample>_programs.csv`, `_kselection.csv`) are listed only via `k_selection.csv` / `samples_analysed.csv`.", ""]
    usage = out.path.str.contains("cnmf_usages_")
    for script, g in out[~usage].groupby("script", sort=False):
        lines += [f"## {script}", "", "| file | MB | md5 | structure |", "|---|---|---|---|"]
        lines += [f"| `{r.path}` | {r.size_mb:.2f} | `{r.md5}` | {r.structure} |" for r in g.itertuples()]
        lines.append("")
    u = out[usage]
    if len(u):
        lines += ["## 13_cnmf_per_sample.py — per-sample usages", "",
                  f"{len(u)} files `data/processed/cnmf_usages_<sample>.parquet` (cell × programme NNLS usages), {u.size_mb.sum():.1f} MB in total; "
                  "md5 per file in `results/tables/phase6_output_checksums.csv`.", ""]
    if missing:
        lines += ["## Not produced", ""] + [f"- `{m}`" for m in missing] + [""]
    (ROOT / "docs/phase6_outputs.md").write_text("\n".join(lines) + "\n")
    out.to_csv(ROOT / "results/tables/phase6_output_checksums.csv", index=False)
    readme = ROOT / "README.md"; txt = readme.read_text()
    head = "## Phase 6 outputs" if "## Phase 6 outputs" not in txt else f"## Phase 6 outputs — update {time.strftime('%Y-%m-%d')}"
    block = ["", head, "",
             "Malignant metaprograms (per-sample cNMF → 20 metaprograms; pre-registered replication rule), per-cell metaprogram scores and "
             "families, assay-masked re-analysis, microenvironment subclusters (raw and cleaned scVI), composition models, and cell2location "
             "deconvolution of 27 Visium sections (4 cohorts) with spatial tests. Results and caveats: `docs/phase6_notes.md`; every file with "
             "md5 and structure: `docs/phase6_outputs.md` (`results/tables/phase6_output_checksums.csv`).", "",
             "| script | key outputs |", "|---|---|"]
    for script, g in out[~usage].groupby("script", sort=False):
        block.append(f"| `{script}` | " + ", ".join(f"`{Path(p).name}`" for p in g.path[:8]) + (" …" if len(g) > 8 else "") + " |")
    block += ["", "Known limitations: cell- and nucleus-derived programmes and TME clusters separate by assay (only 6 cell samples, 3 publications); "
              "spatial state classes come from Filbin/Liu marker scoring (fallback) except the metaprogram-family brainstem re-run; the brainstem "
              "gradient is one autopsy; autopsy vs resection composition is not estimable (confounded with study).", ""]
    with open(readme, "a") as fh:
        fh.write("\n".join(block))
    log.write(f"docs/phase6_outputs.md written ({len(out)} files); README block '{head}' appended")
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
