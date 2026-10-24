"""Phase 6 step 1 — per-sample consensus NMF on malignant cells of the reference samples (>= 200 selected cells).

Outputs (never overwritten if present, unless --force):
    results/tables/cnmf<tag>/<sample>_programs.csv      top-50 genes per consensus programme (spectra score + loading)
    data/processed/cnmf_usages<tag>_<sample>.parquet     cell x programme NNLS usages
    results/tables/cnmf<tag>/k_selection.csv             stability / error per k and the chosen k per sample
    results/tables/cnmf<tag>/samples_analysed.csv        per-sample cell counts and selection sources
Run: nice -n 10 python scripts/13_cnmf_per_sample.py --selection primary [--only S1 S2 ...] [--force]

Audit correction (2026-10-19, MC-04 / MC-16): --selection primary|eventtrust reads the scripted selection table
data/processed/malignant_selection_audit.parquet (scripts/13a_malignant_selection.py; corrected consensus calls, qc_pass applied,
okada__1734_pons excluded). 'primary' writes to results/tables/cnmf/ (tag ''), 'eventtrust' to results/tables/cnmf_eventtrust/.
--selection legacy reproduces the pre-audit rule (consensus OR scANVI rescue; kept for the record only).
"""
import sys, time, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py, scipy.sparse as sp  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import PROCESSED, RESULTS, LOGS  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); ap.add_argument("--only", nargs="*")
ap.add_argument("--selection", default="primary", choices=["legacy", "primary", "eventtrust"])
ap.add_argument("--no-samples-table", action="store_true", help="do not rewrite samples_analysed.csv (parallel shards)")
args = ap.parse_args()
TAG = {"legacy": "", "primary": "", "eventtrust": "_eventtrust"}[args.selection]
CORE = PROCESSED / "atlas_core.h5ad"
OUT = RESULTS / "tables" / f"cnmf{TAG}"; OUT.mkdir(parents=True, exist_ok=True)
log = U.StepLogger(LOGS / "phase6_metaprograms.log", context=f"13_cnmf[{args.selection}]")

with h5py.File(CORE) as f:
    cols = ["sample_id", "study", "assay", "scanvi_label", "scanvi_maxprob", "consensus_call", "qc_pass"]
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in cols}, index=read_elem(f["obs"]["_index"])); var = read_elem(f["var"])
if args.selection == "legacy":
    obs["mal_sel_source"] = MP.select_malignant(obs)
    sel = obs["mal_sel_source"] != "not_selected"
else:
    sel_tab = pd.read_parquet(PROCESSED / "malignant_selection_audit.parquet", columns=[args.selection, "cnv_confidence"]).reindex(obs.index)
    sel = sel_tab[args.selection].fillna(False).astype(bool)
    obs["mal_sel_source"] = np.where(sel, "consensus_" + sel_tab.cnv_confidence.astype(str), "not_selected")
per = obs.groupby("sample_id", observed=True).agg(n_cells=("sample_id", "size"), study=("study", "first"), assay=("assay", "first"))
per["n_malignant_selected"] = sel.groupby(obs["sample_id"], observed=True).sum()
for s in sorted(set(obs["mal_sel_source"]) - {"not_selected"}):
    per[f"n_{s}"] = (obs["mal_sel_source"] == s).groupby(obs["sample_id"], observed=True).sum()
per["analysed"] = per["n_malignant_selected"] >= 200
per["selection"] = args.selection
if not args.no_samples_table:
    per.to_csv(OUT / "samples_analysed.csv")
log.write(f"selection={args.selection}: {per.analysed.sum()} of {len(per)} samples have >= 200 selected malignant cells; total selected {int(sel.sum())}")


def read_rows(idx):
    idx = np.asarray(idx); runs = np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1); parts = []
    with h5py.File(CORE) as f:
        g = f["layers/counts"]; indptr = g["indptr"]
        for r in runs:
            a, b = int(r[0]), int(r[-1]) + 1
            ip = indptr[a:b + 1]; s, e = int(ip[0]), int(ip[-1])
            parts.append(sp.csr_matrix((g["data"][s:e], g["indices"][s:e], ip - s), shape=(b - a, var.shape[0])))
    return sp.vstack(parts).tocsr()


ksel_all = []
todo = [s for s in per.index[per.analysed] if not args.only or s in args.only]
t_all = time.time()
for i, sid in enumerate(todo):
    p_csv = OUT / f"{sid}_programs.csv"; p_par = PROCESSED / f"cnmf_usages{TAG}_{sid}.parquet"; p_k = OUT / f"{sid}_kselection.csv"
    if p_csv.exists() and p_par.exists() and p_k.exists() and not args.force:
        ksel_all.append(pd.read_csv(p_k)); log.write(f"skip {sid} (exists)"); continue
    with log.step(f"cnmf {sid} ({i+1}/{len(todo)}) n={int(per.loc[sid,'n_malignant_selected'])}"):
        rows = np.flatnonzero(((obs.sample_id == sid) & sel).values)
        X = read_rows(rows)
        programs, usages, ksel = MP.run_cnmf_sample(X, var.index.to_numpy(), obs.index.to_numpy()[rows], sid)
        programs.to_csv(p_csv, index=False); usages.to_parquet(p_par); ksel.to_csv(p_k, index=False)
        ksel_all.append(ksel)
        kk = int(ksel.k[ksel.chosen].iloc[0]); log.write(f"  {sid}: chosen k={kk} ({ksel.rule.iloc[0]}); sil by k: " +
                                                        " ".join(f"{int(r.k)}:{r.sil:.2f}" for r in ksel.itertuples()))
if not args.only:
    pd.concat(ksel_all).to_csv(OUT / "k_selection.csv", index=False)
log.write(f"ALL DONE {len(todo)} samples | wall {(time.time()-t_all)/60:.1f} min | peak_rss {U.peak_rss_gb():.2f} GB")
