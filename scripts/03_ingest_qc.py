#!/usr/bin/env python
"""Phase 2: ingest one study, compute per-sample QC + scrublet flags, harmonise genes, write data/interim/<study>.qc.h5ad.

Usage: python scripts/03_ingest_qc.py --study jessa_scrna [--no-scrublet] [--compression]
Outputs
  data/interim/<study>.qc.h5ad            included (atlas) samples, all cells, obs['qc_pass'] flag
  data/interim/<study>.outgroup.qc.h5ad   excluded integrated-role out-group samples (if any)
  results/qc/qc_per_sample_<study>.csv    per-sample QC table
  results/qc/gene_harmonisation.csv       appended per study (overwritten for the same study)
  logs/phase2_qc.log                      wall time + peak RSS per step
"""
import sys
import time
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from dmg_atlas import utils  # noqa: E402  (sets thread caps before numpy import)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import anndata as ad  # noqa: E402

from dmg_atlas import INTERIM, RESULTS, LOGS  # noqa: E402
from dmg_atlas.ingest import load_manifest, select_rows, LOADERS  # noqa: E402
from dmg_atlas.genes import GeneMapper  # noqa: E402
from dmg_atlas import qc as qcm  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True, choices=sorted(LOADERS))
    ap.add_argument("--no-scrublet", action="store_true")
    ap.add_argument("--compression", default=None, help="h5ad compression, e.g. lzf or gzip (default none)")
    args = ap.parse_args()
    study = args.study
    log = utils.StepLogger(LOGS / "phase2_qc.log", context=study)
    t_all = time.time()
    cfg = qcm.load_config()
    m = load_manifest()
    rows = select_rows(m, study)
    log.write(f"{len(rows)} manifest rows selected ({int(rows.include_in_atlas.sum())} atlas, {int((~rows.include_in_atlas).sum())} out-group)")
    mapper = GeneMapper()
    loader = LOADERS[study]
    kwargs = {"log": log.write} if study in ("collot", "okada", "liu", "labelle") else {}

    parts, summaries, harm_recs = [], [], []
    for row, a in loader(rows, **kwargs):
        with log.step(f"sample {row.sample_id}"):
            a = qcm.compute_metrics(a, cfg)
            a = qcm.apply_thresholds(a, cfg, row.assay, chemistry=row.chemistry)
            a = qcm.sample_gate(a, cfg)
            if not args.no_scrublet and row.counts_type == "umi":
                a = qcm.run_scrublet(a, cfg, log=log.write)
            else:
                a.obs["doublet_score"] = np.nan
                a.obs["predicted_doublet"] = pd.array([pd.NA] * a.n_obs, dtype="boolean")
                a.obs["predicted_doublet_fixed"] = pd.array([pd.NA] * a.n_obs, dtype="boolean")
            a, rec = mapper.harmonise(a, study)
            if "labelle_table_genes" in a.uns:  # audit DQ-09: which deposited table measured each gene (zero-filled elsewhere)
                for key, genes in a.uns["labelle_table_genes"].items():
                    canon, _ = mapper.map_symbols(pd.Index(genes))
                    a.var[f"measured_in_{key}"] = a.var_names.isin(set(c for c in canon if c is not None))
                a.uns["labelle_table_genes"] = {k: len(v) for k, v in a.uns["labelle_table_genes"].items()}
            rec["sample_id"] = row.sample_id
            harm_recs.append(rec)
            s = qcm.per_sample_summary(a)
            summaries.append(s)
            log.write(f"  {row.sample_id}: {s['cells_before_qc']} -> {s['cells_after_qc']} cells pass; median genes {s['median_genes']:.0f}; "
                      f"median mito {s['median_pct_mito']:.2f}%; doublets flagged {s['doublet_rate_flagged']:.3f}; genes {rec['n_input_genes']}->{rec['n_output_genes']}; "
                      f"mt_genes_present {s['mt_genes_present']}; mito rule {s['qc_mito_rule']}; sample_qc_pass {s['sample_qc_pass']}")
            # free per-cell uns not needed in the concatenated object
            a.uns.pop("scrublet", None)
            parts.append(a)

    with log.step("concatenate + write"):
        obs_str_cols = None
        for key, sel in (("", lambda a: bool(a.obs["include_in_atlas"].iloc[0])), ("outgroup.", lambda a: not bool(a.obs["include_in_atlas"].iloc[0]))):
            sub = [a for a in parts if sel(a)]
            if not sub:
                continue
            # within a study all samples share the same harmonised var except TPM cohorts: outer join, zero fill
            adata = sub[0] if len(sub) == 1 else ad.concat(sub, join="outer", merge="first", uns_merge="first", fill_value=0)
            adata.var = mapper.uni.set_index("canonical_symbol").reindex(adata.var_names)[["ensembl_id"]]
            adata.var.index.name = "symbol"
            # genes absent from a sample's deposited feature list are zero-filled by the outer join; record how many samples measured each gene
            cnt = pd.Series(0, index=adata.var_names, dtype=int)
            for a_ in sub:
                cnt.loc[a_.var_names] += 1
            adata.var["n_samples_measured"] = cnt.values
            adata.var["measured_in_all_samples"] = cnt.values == len(sub)
            for c in sorted({c for a_ in sub for c in a_.var.columns if c.startswith("measured_in_")}):  # LaBelle per-table flags (DQ-09)
                flag = pd.Series(False, index=adata.var_names)
                for a_ in sub:
                    if c in a_.var:
                        flag.loc[a_.var_names[a_.var[c].values.astype(bool)]] = True
                adata.var[c] = flag.values
            for L in list(adata.layers.keys()):
                adata.layers[L] = sp.csr_matrix(adata.layers[L])
            adata.X = sp.csr_matrix(adata.X)
            for c in adata.obs.columns:  # categorical/object hygiene for h5ad
                if adata.obs[c].dtype == object:
                    adata.obs[c] = adata.obs[c].astype(str).replace({"nan": np.nan, "None": np.nan, "<NA>": np.nan})
            adata.uns["phase2"] = {"study": study, "qc_config": cfg, "gene_universe": "GRCh38-2020-A symbols (see genes.py)",
                                   "n_samples": len(sub), "written": time.strftime("%Y-%m-%d %H:%M:%S")}
            if "value_scale" in sub[0].uns:
                adata.uns["value_scale"] = sub[0].uns["value_scale"]
            out = INTERIM / f"{study}.{key}qc.h5ad"
            adata.write_h5ad(out, compression=args.compression)
            log.write(f"wrote {out} : {adata.n_obs} cells x {adata.n_vars} genes; qc_pass {int(adata.obs.qc_pass.sum())}; layers {list(adata.layers)}")
        del parts

    (RESULTS / "qc").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(RESULTS / "qc" / f"qc_per_sample_{study}.csv", index=False)
    hf = RESULTS / "qc" / "gene_harmonisation.csv"
    h = pd.DataFrame(harm_recs)
    if hf.exists():
        old = pd.read_csv(hf)
        old = old[old.study != study]
        h = pd.concat([old, h], ignore_index=True)
    h.to_csv(hf, index=False)
    log.write(f"DONE {study} | total wall {(time.time() - t_all) / 60:.1f} min | peak RSS {utils.peak_rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
