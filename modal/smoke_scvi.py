
import time, json, os, numpy as np, scipy.sparse as sp, torch, anndata as ad, scvi, scanpy, cell2location
import importlib.metadata as md
out = {"torch": torch.__version__, "cuda": torch.version.cuda, "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
       "scvi": scvi.__version__, "scanpy": scanpy.__version__, "anndata": ad.__version__, "cell2location": md.version("cell2location"), "pyro": md.version("pyro-ppl")}
try:
    import scgpt; out["scgpt"] = md.version("scgpt")
except Exception as e: out["scgpt_import_error"] = repr(e)[:200]
n, g = 30000, 3000
rng = np.random.default_rng(0)
X = sp.random(n, g, density=0.10, random_state=0).tocsr(); X.data = (rng.negative_binomial(2, 0.3, X.nnz)+1).astype(np.float32)
a = ad.AnnData(X); a.obs["batch"] = rng.choice([f"s{i}" for i in range(40)], n)
scvi.model.SCVI.setup_anndata(a, batch_key="batch")
m = scvi.model.SCVI(a, n_layers=2, n_latent=30, gene_likelihood="nb")
t = time.time(); m.train(max_epochs=3, batch_size=256, enable_progress_bar=False, accelerator="gpu", devices=1); dt = (time.time()-t)/3
out["scvi_sec_per_epoch_30k_x_3k"] = round(dt, 2)
out["gpu_mem_peak_GB"] = round(torch.cuda.max_memory_allocated()/1e9, 2)
os.makedirs("out", exist_ok=True); json.dump(out, open("out/smoke.json","w"), indent=1)
print(json.dumps(out, indent=1))
