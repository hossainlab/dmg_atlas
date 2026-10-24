#!/usr/bin/env python
"""Phase 6 / step 16a (audit MC-03, 2026-10-19) — write the curated Phase-6 gene sets from dmg_atlas.states.CURATED_PHASE6_SETS.

Output: results/tables/genesets/phase6_curated_gene_sets.csv (gene_set, source, gene, in_atlas) and
        results/tables/genesets/phase6_curated_gene_sets.sha256
The previous unscripted results/tables/genesets/phase6_gene_sets.csv is kept as phase6_gene_sets.pre_audit.csv; its Liu2022_* sets were
in-house DE sets (see reference_gene_sets_provenance.csv / phase1 notes) and are superseded by reference_gene_sets_scripted.csv.
Run: python scripts/16a_curated_gene_sets.py
"""
import sys, hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd  # noqa: E402
from dmg_atlas import RESULTS, PROCESSED  # noqa: E402
from dmg_atlas.states import CURATED_PHASE6_SETS  # noqa: E402

out = RESULTS / "tables" / "genesets"; out.mkdir(parents=True, exist_ok=True)
atlas_genes = set(pd.read_csv(PROCESSED / "gene_universe.csv").iloc[:, 0]) if (PROCESSED / "gene_universe.csv").exists() else None
if atlas_genes is None:
    import h5py
    from anndata.io import read_elem
    with h5py.File(PROCESSED / "atlas_core.h5ad", "r") as f:
        atlas_genes = set(read_elem(f["var"]).index)
rows = [dict(gene_set=k, source="curated_phase6_states.py", gene=g, in_atlas=g in atlas_genes) for k, v in CURATED_PHASE6_SETS.items() for g in v]
df = pd.DataFrame(rows)
old = out / "phase6_gene_sets.csv"
if old.exists() and not (out / "phase6_gene_sets.pre_audit.csv").exists():
    old.rename(out / "phase6_gene_sets.pre_audit.csv")
p = out / "phase6_curated_gene_sets.csv"; df.to_csv(p, index=False)
(out / "phase6_curated_gene_sets.sha256").write_text(hashlib.sha256(p.read_bytes()).hexdigest() + "  " + p.name + "\n")
print(df.groupby("gene_set").agg(n=("gene", "size"), in_atlas=("in_atlas", "sum")).to_string())
print("sha256", hashlib.sha256(p.read_bytes()).hexdigest())
