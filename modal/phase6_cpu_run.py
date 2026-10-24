#!/usr/bin/env python
"""Run CPU-only Phase 6 steps in a Modal sandbox (image im-gl8KES0WsH7PoFjf83hK6v) instead of on the Mac.

  python modal/phase6_cpu_run.py [--harvest cohort=<dir> ...] [--stage <rel> ...] [--link <rel>=<volume path> ...] -- "<cmd>" ...

Stages a slice of the repo (src/, scripts/, results/qc/spatial_sections.csv, results/tables/{cnmf,genesets}/,
results/tables/c2l_*.csv, logs/phase6_*.log, and each harvested cell2location out/ dir as /repo/harvest/<cohort>/)
into /repo (plus any --stage paths), mounts the Modal volume dmg-atlas-data at /vol and symlinks --link targets
(e.g. data/processed/atlas_core.h5ad=atlas_core.h5ad, uploaded with `modal volume put`), runs each command there in order
(stops at the first failure), then pulls back every regular file under results/, data/processed/ and logs/ that is new or
changed, unpacked over the local repo. Commands run with DMG_ATLAS_THREADS=8 on an 8-CPU, 64 GB sandbox (+ --gpu A10G when asked).
"""
import argparse, base64, io, shutil, sys, tarfile, tempfile
from pathlib import Path
import modal

ROOT = Path(__file__).resolve().parents[1]
IMAGE = "im-gl8KES0WsH7PoFjf83hK6v"
STAGE = ["src", "scripts", "results/qc/spatial_sections.csv", "results/tables/cnmf", "results/tables/genesets",
         "results/tables/c2l_signatures.csv", "results/tables/c2l_reference_classes.csv",
         "logs/phase6_spatial.log", "logs/phase6_metaprograms.log"]
PULL = "results data/processed logs models docs README.md"
VOLUME = "dmg-atlas-data"

ap = argparse.ArgumentParser()
ap.add_argument("--harvest", nargs="*", default=[], help="cohort=dir pairs")
ap.add_argument("--stage", nargs="*", default=[], help="extra repo-relative paths to stage")
ap.add_argument("--link", nargs="*", default=[], help="rel=volume_path symlinks into /vol")
ap.add_argument("--gpu", default=None, help="request a GPU (e.g. A10G) for GPU jobs such as modal/phase6_c2l.py or scVI")
ap.add_argument("--attach", help="reuse a live sandbox (sb-...) whose runner died: wait for its python jobs, then run cmds and pull")
ap.add_argument("cmds", nargs="*")
a = ap.parse_args()

stage = Path(tempfile.mkdtemp(prefix="dmg_stage_"))
if a.attach:
    a.harvest, a.stage, a.link = [], [], []
    STAGE = []
ign = shutil.ignore_patterns("__pycache__", "*.pyc")
for rel in STAGE + a.stage:
    s, d = ROOT / rel, stage / rel
    d.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(s, d, ignore=ign, dirs_exist_ok=True) if s.is_dir() else shutil.copy2(s, d)
for h in a.harvest:
    c, p = h.split("=", 1)
    shutil.copytree(p, stage / "harvest" / c)
if a.attach:
    sb = modal.Sandbox.from_id(a.attach)
else:
    app = modal.App.lookup("dmg-atlas-phase6-cpu", create_if_missing=True)
    img = modal.Image.from_id(IMAGE).add_local_dir(stage, "/repo", copy=True)
    sb = modal.Sandbox.create("sleep", "infinity", image=img, app=app, cpu=8, memory=65536, timeout=6 * 3600, workdir="/repo", gpu=a.gpu,
                              volumes={"/vol": modal.Volume.from_name(VOLUME)}, env={"DMG_ATLAS_THREADS": "8", "MPLBACKEND": "Agg"})
print("sandbox", sb.object_id, flush=True)
ok = True
try:
    links = " && ".join(f"mkdir -p $(dirname /repo/{r}) && ln -sf /vol/{v} /repo/{r}" for r, v in (x.split("=", 1) for x in a.link))
    if a.attach:  # wait until no python job is left in the sandbox (the stamp from the original run is kept)
        sb.exec("bash", "-c", "while grep -qs '[p]ython' /proc/[0-9]*/cmdline; do sleep 15; done").wait()
    else:
        sb.exec("bash", "-c", (links + " && " if links else "") + "touch /repo/.stamp && sleep 1").wait()  # image build resets mtimes; stamp marks job start
    for cmd in a.cmds:
        print(f"\n=== {cmd}", flush=True)
        p = sb.exec("bash", "-lc", cmd)
        for line in p.stdout:
            print(line, end="", flush=True)
        p.wait()
        err = p.stderr.read()
        if p.returncode:
            print(err[-4000:], file=sys.stderr); print(f"FAILED rc={p.returncode}", flush=True); ok = False; break
    p = sb.exec("bash", "-c", f"cd /repo && find {PULL} -type f -newer .stamp 2>/dev/null | tar czf - -T - | base64 -w0")
    blob = base64.b64decode(p.stdout.read()); p.wait()
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as t:
        names = t.getnames(); t.extractall(ROOT, filter="data")
    print(f"\npulled {len(names)} files:\n  " + "\n  ".join(sorted(names)))
finally:
    sb.terminate()
    shutil.rmtree(stage, ignore_errors=True)
sys.exit(0 if ok else 1)
