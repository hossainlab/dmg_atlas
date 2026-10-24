#!/usr/bin/env python
"""Phase 7 / step 26 — foundation-model embeddings (scGPT, Geneformer) vs scVI / PCA for malignant-state separability (methods note).

Inputs : data/processed/atlas_core.h5ad (layers['counts'] via h5py, obsm['X_scvi_malignant'], obs), data/processed/metaprogram_scores.parquet
         (family assignment + z scores, step 19), results/tables/atlas_malignant_hvg.csv, data/interim/gene_universe.csv (symbol -> Ensembl).
Pretrained weights (downloaded on first run, cached under models/foundation/): scGPT whole-human (Google Drive folder of the scGPT
authors, via gdown) and Geneformer V2-104M (Hugging Face ctheodoris/Geneformer; V1-10M if V2 fails).
Outputs (never overwritten if present, unless --force):
  results/tables/foundation_model_comparison.csv   per embedding: kNN family balanced accuracy (5-fold), cross-assay transfer accuracy
                                                   (train nuclei -> predict cells), silhouette, same-sample / same-assay neighbour fractions
  data/processed/foundation_embeddings.npz         cell ids + every embedding
Figure (results/figures/): phase7_foundation_models.png.
Run: python scripts/26_foundation_embeddings.py [--force]   (GPU; needs gdown + transformers: pip install gdown 'transformers==4.46.3' — newer releases need torch >= 2.6)

Method (fixed before running)
  Cells: malignant cells with a confident family (max z >= 1, step 19; 'unassigned' and cilia excluded), up to 2,500 per family drawn
  stratified by sample (seed 0). Labels are expression-derived, so every embedding is compared on the same (partly circular) footing.
  PCA: normalize_total(1e4) + log1p on the 3,000 malignant HVGs, scaled, 50 PCs. scVI: atlas X_scvi_malignant (batch = sample_id).
  scGPT: scgpt.tasks.embed_data on raw counts (gene symbols, max_length 1200, use_fast_transformer=False: no flash-attn in the image).
  Geneformer: rank-value encoding (counts / cell total x 1e4 / gene median from the model's dictionary, genes ranked descending,
  <cls> + top 4,095 tokens), cell embedding = mean of the penultimate layer over gene tokens (Geneformer default), fp16 inference.
  Metrics on L2-normalised embeddings: kNN (k = 15, cosine) balanced accuracy, stratified 5-fold; cross-assay = kNN trained on nuclei
  cells, scored on scRNA cells; silhouette (cosine, 5,000-cell subsample); mean fraction of 15 nearest neighbours from the same sample /
  same assay class (lower = more mixing; expected values from label frequencies are reported).
"""
import sys, argparse, pickle
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, h5py, scipy.sparse as sp  # noqa: E402
import anndata as ad, scanpy as sc  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors  # noqa: E402
from sklearn.metrics import balanced_accuracy_score, silhouette_score  # noqa: E402
from sklearn.preprocessing import normalize  # noqa: E402
from dmg_atlas import ROOT, PROCESSED, INTERIM, RESULTS, LOGS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
CORE = PROCESSED / "atlas_core.h5ad"; TABLES = RESULTS / "tables"; FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
FM = ROOT / "models" / "foundation"; FM.mkdir(parents=True, exist_ok=True)
log = U.StepLogger(LOGS / "phase7_multimodal.log", context="26_foundation_embeddings")
PER_FAMILY, MIN_Z, K, SEED, GF_LEN, GF_BATCH = 2500, 1.0, 15, 0, 4096, 8
SCGPT_FOLDER = "https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y"
ASSAY_CLASS = {"scRNA": "cells", "snRNA": "nuclei", "multiome_RNA": "nuclei"}
OUTS = [TABLES / "foundation_model_comparison.csv", PROCESSED / "foundation_embeddings.npz"]
if all(p.exists() for p in OUTS) and not args.force:
    log.write("all outputs exist; skip (use --force)"); sys.exit(0)
import torch  # noqa: E402
dev = "cuda" if torch.cuda.is_available() else "cpu"


def read_rows(g, rows: np.ndarray, n_genes: int, block: int = 20000) -> sp.csr_matrix:
    """Sorted row subset of an on-disk CSR group (layers/counts), read in contiguous row blocks via h5py."""
    indptr = g["indptr"][:]; parts = []
    for b0 in range(0, len(indptr) - 1, block):
        r = rows[(rows >= b0) & (rows < b0 + block)]
        if not len(r):
            continue
        b1 = min(b0 + block, len(indptr) - 1); s, e = int(indptr[b0]), int(indptr[b1])
        m = sp.csr_matrix((g["data"][s:e], g["indices"][s:e], indptr[b0:b1 + 1] - s), shape=(b1 - b0, n_genes))
        parts.append(m[r - b0])
    return sp.vstack(parts).tocsr()


with log.step("select cells"):
    f = h5py.File(CORE, "r")
    idx_all = pd.Index(read_elem(f["obs"]["_index"]))
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in ("sample_id", "assay")}, index=idx_all).astype(str)
    mp = pd.read_parquet(PROCESSED / "metaprogram_scores.parquet")
    zc = [c for c in mp.columns if c.startswith("z_")]
    mp = mp[(mp[zc].max(1) >= MIN_Z) & ~mp.metaprogram.isin(["unassigned", "novel_MP15"])]
    rng = np.random.default_rng(SEED); pick = []
    for fam, g in mp.groupby("metaprogram"):
        if len(g) <= PER_FAMILY:
            pick.extend(g.index); continue
        w = g.groupby("sample_id").sample_id.transform("size"); p = (1 / w) / (1 / w).sum()   # flatten over samples
        pick.extend(rng.choice(g.index.values, size=PER_FAMILY, replace=False, p=p.values))
    rows = np.sort(idx_all.get_indexer(pick)); cells = idx_all[rows]
    lab = mp.metaprogram.reindex(cells).values; samp = obs.sample_id.values[rows]; assay = obs.assay.map(ASSAY_CLASS).values[rows]
    var = read_elem(f["var"]); genes = var.index.to_numpy()
    X = read_rows(f["layers/counts"], rows, len(genes)).astype(np.float32)
    log.write(f"{len(cells)} cells; families {pd.Series(lab).value_counts().to_dict()}; assay {pd.Series(assay).value_counts().to_dict()}")

emb, status = {}, {}
with log.step("PCA + scVI"):
    a = ad.AnnData(X=X.copy(), var=pd.DataFrame(index=genes)); sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    hv = pd.read_csv(TABLES / "atlas_malignant_hvg.csv"); hv = [g for g in hv.symbol[hv.highly_variable] if g in a.var_names]
    a = a[:, hv].copy(); sc.pp.scale(a, max_value=10); sc.tl.pca(a, n_comps=50, random_state=SEED)
    emb["PCA"] = a.obsm["X_pca"]; emb["scVI"] = f["obsm"]["X_scvi_malignant"][rows]; status["PCA"] = status["scVI"] = "ok"
    del a

with log.step("scGPT"):
    try:
        import gdown
        from scgpt.tasks import embed_data
        mdir = FM / "scGPT_human"
        if not list(mdir.rglob("best_model.pt")):
            gdown.download_folder(SCGPT_FOLDER, output=str(mdir), quiet=True, use_cookies=False)
        mdir = next(mdir.rglob("best_model.pt")).parent   # gdown may nest the folder
        a = ad.AnnData(X=X.copy(), var=pd.DataFrame(index=genes)); a.var["gene_name"] = a.var_names
        out = embed_data(a, mdir, gene_col="gene_name", max_length=1200, batch_size=64, device=dev, use_fast_transformer=False, return_new_adata=True)
        emb["scGPT"] = np.asarray(out.X); status["scGPT"] = f"ok ({mdir.name})"
    except Exception as e:
        status["scGPT"] = f"failed: {type(e).__name__}: {str(e)[:200]}"
    log.write(f"scGPT: {status['scGPT']}")

with log.step("Geneformer"):
    try:
        from huggingface_hub import hf_hub_download
        from transformers import BertForMaskedLM
        model, used = None, None
        for sub, dicts, mlen in (("Geneformer-V2-104M", ("geneformer/token_dictionary_gc104M.pkl", "geneformer/gene_median_dictionary_gc104M.pkl"), 4096),
                                 ("Geneformer-V1-10M", ("geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl",
                                                        "geneformer/gene_dictionaries_30m/gene_median_dictionary_gc30M.pkl"), 2048)):
            try:
                tok = pickle.load(open(hf_hub_download("ctheodoris/Geneformer", dicts[0], cache_dir=FM), "rb"))
                med = pickle.load(open(hf_hub_download("ctheodoris/Geneformer", dicts[1], cache_dir=FM), "rb"))
                for fn in ("config.json", "model.safetensors"):
                    p = hf_hub_download("ctheodoris/Geneformer", f"{sub}/{fn}", cache_dir=FM)
                model = BertForMaskedLM.from_pretrained(Path(p).parent, output_hidden_states=True).to(dev).eval(); used = sub
                GF_MAX = min(mlen, model.config.max_position_embeddings); break
            except Exception as e:
                log.write(f"  {sub} unavailable: {type(e).__name__}: {repr(e)[:300]}")
        if model is None:
            raise RuntimeError("no Geneformer checkpoint could be loaded")
        uni = pd.read_csv(INTERIM / "gene_universe.csv").drop_duplicates("canonical_symbol").set_index("canonical_symbol").ensembl_id
        ens = uni.reindex(genes).values
        gi = np.array([i for i, e in enumerate(ens) if isinstance(e, str) and e in tok and e in med])
        gtok = np.array([tok[ens[i]] for i in gi]); gmed = np.array([med[ens[i]] for i in gi], dtype=np.float32)
        cls_id, pad_id = tok.get("<cls>"), tok.get("<pad>", 0)
        Xn = normalize(X, norm="l1", axis=1)[:, gi].tocsr()
        log.write(f"  {used}: {len(gi)} atlas genes in the Geneformer vocabulary; <cls> {'present' if cls_id is not None else 'absent'}")
        maxlen = GF_MAX - (1 if cls_id is not None else 0); out = []
        with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.float16, enabled=dev == "cuda"):
            for s in range(0, Xn.shape[0], GF_BATCH):
                seqs = []
                for i in range(s, min(s + GF_BATCH, Xn.shape[0])):
                    r = Xn[i]; v = r.data * 1e4 / gmed[r.indices]
                    t = gtok[r.indices[np.argsort(-v, kind="stable")]][:maxlen]
                    seqs.append(([cls_id] if cls_id is not None else []) + t.tolist())
                L = max(len(q) for q in seqs)
                ids_ = torch.full((len(seqs), L), pad_id, dtype=torch.long); att = torch.zeros((len(seqs), L), dtype=torch.long)
                for j, q in enumerate(seqs):
                    ids_[j, :len(q)] = torch.tensor(q); att[j, :len(q)] = 1
                hs = model(input_ids=ids_.to(dev), attention_mask=att.to(dev)).hidden_states[-2].float()
                m_ = att.to(dev).clone().float()
                if cls_id is not None:
                    m_[:, 0] = 0
                out.append(((hs * m_[..., None]).sum(1) / m_.sum(1, keepdim=True)).cpu().numpy())
        emb["Geneformer"] = np.vstack(out); status["Geneformer"] = f"ok ({used})"
    except Exception as e:
        status["Geneformer"] = f"failed: {type(e).__name__}: {str(e)[:200]}"
    log.write(f"Geneformer: {status['Geneformer']}")

with log.step("metrics"):
    rows_m = []
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    cells_mask, nuc_mask = assay == "cells", assay == "nuclei"
    sub = rng.choice(len(cells), size=min(5000, len(cells)), replace=False)
    exp_s = np.mean([(samp == s).mean() for s in samp]); exp_a = np.mean([(assay == x).mean() for x in assay])
    for name, E in emb.items():
        Z = normalize(np.nan_to_num(np.asarray(E, dtype=np.float32)))
        pred = np.empty(len(lab), dtype=object)
        for tr, te in skf.split(Z, lab):
            pred[te] = KNeighborsClassifier(K, metric="cosine").fit(Z[tr], lab[tr]).predict(Z[te])
        xa = balanced_accuracy_score(lab[cells_mask], KNeighborsClassifier(K, metric="cosine").fit(Z[nuc_mask], lab[nuc_mask]).predict(Z[cells_mask])) \
            if cells_mask.sum() and nuc_mask.sum() else np.nan
        nn = NearestNeighbors(n_neighbors=K + 1, metric="cosine").fit(Z).kneighbors(Z, return_distance=False)[:, 1:]
        rows_m.append(dict(embedding=name, status=status[name], dims=Z.shape[1], n_cells=len(Z), knn_balanced_accuracy=balanced_accuracy_score(lab, pred),
                           cross_assay_accuracy_nuclei_to_cells=xa, silhouette=silhouette_score(Z[sub], lab[sub], metric="cosine"),
                           same_sample_neighbours=(samp[nn] == samp[:, None]).mean(), same_sample_expected=exp_s,
                           same_assay_neighbours=(assay[nn] == assay[:, None]).mean(), same_assay_expected=exp_a))
    for name in status:
        if name not in emb:
            rows_m.append(dict(embedding=name, status=status[name]))
    res = pd.DataFrame(rows_m)
    log.write("\n" + res.drop(columns=["status"]).round(3).to_string(index=False) + "\n" + str(status))

with log.step("write + figure"):
    res.to_csv(OUTS[0], index=False)
    np.savez_compressed(OUTS[1], cell_ids=np.asarray(cells, dtype=str), family=lab.astype(str), **{k: np.asarray(v, dtype=np.float32) for k, v in emb.items()})
    apply_figure_style()
    r = res.dropna(subset=["knn_balanced_accuracy"])
    mets = [("knn_balanced_accuracy", "kNN family accuracy"), ("cross_assay_accuracy_nuclei_to_cells", "nuclei → cells transfer"),
            ("silhouette", "silhouette"), ("same_sample_neighbours", "same-sample neighbours")]
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.0))
    for ax, (m_, t) in zip(axes, mets):
        ax.bar(range(len(r)), r[m_], color="#4d4d4d"); ax.set_xticks(range(len(r))); ax.set_xticklabels(r.embedding, rotation=45, ha="right")
        ax.set_ylabel(t, fontsize=6)
        if m_ == "same_sample_neighbours":
            ax.axhline(r.same_sample_expected.iloc[0], color="grey", lw=0.5, ls=":")
    fig.savefig(FIGS / "phase7_foundation_models.png"); plt.close(fig)
log.write(f"done | peak RSS {U.peak_rss_gb():.2f} GB")
