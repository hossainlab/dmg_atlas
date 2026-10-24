"""Regenerate results/tables/tier3_excluded_author_composition.csv (re-audit CNV-04).

The table written by 06_cnv_aggregate.py took each group's author_call from its first cell
(gg.author_call.iloc[0]) while grouping only on author_state, which mislabelled 666 Jessa-multiome
cells (author_state values there occur under both author calls). 06_cnv_aggregate.py is patched;
this script rebuilds the single affected table from data/processed/malignant_calls.csv without
re-running the whole aggregation. The old file is kept as *.pre_reaudit.csv.
"""
import numpy as np, pandas as pd, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "tables" / "tier3_excluded_author_composition.csv"
calls = pd.read_csv(ROOT / "data" / "processed" / "malignant_calls.csv",
                    usecols=["study", "cnv_call", "marker_call", "author_call", "author_state"])
labelled = calls[calls.author_call.notna()]
rows = []
for study, g in labelled.groupby("study"):
    mk_bin = np.where(g.marker_call == "normal", "normal",
                      np.where(g.marker_call == "malignant_candidate", "malignant", None))
    ex = g[g.cnv_call.where(mk_bin == g.cnv_call.values).isna()]          # tier-3 prediction undefined
    for (st_lab, acall, cnvc, mk), gg in ex.groupby(["author_state", "author_call", "cnv_call", "marker_call"]):
        rows.append(dict(study=study, author_state=st_lab, author_call=acall, cnv_call=cnvc, marker_call=mk,
                         n_cells=len(gg), frac_of_labelled_in_study=len(gg) / len(g)))
df = pd.DataFrame(rows).sort_values(["study", "n_cells"], ascending=[True, False])
if OUT.exists():
    shutil.copy2(OUT, OUT.with_suffix(".pre_reaudit.csv"))
df.to_csv(OUT, index=False)
m = df[(df.cnv_call == "malignant") & (df.marker_call == "normal")]
print(f"rows {len(df)}; cnv=malignant|marker=normal {m.n_cells.sum()} cells, "
      f"author-malignant {m.loc[m.author_call == 'malignant', 'n_cells'].sum()} "
      f"({100 * m.loc[m.author_call == 'malignant', 'n_cells'].sum() / m.n_cells.sum():.1f}%)")
for study, gg in m.groupby("study"):
    print(" ", study, gg.groupby("author_call").n_cells.sum().to_dict())
