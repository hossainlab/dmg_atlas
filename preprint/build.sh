#!/usr/bin/env bash
# Build preprint/manuscript.pdf. Tectonic cannot fetch its TeX bundle in this sandbox (TLS verification
# fails on macOS), and the conda texlive-core package ships binaries only, so the bundle is downloaded
# with Python, extracted, and made readable (the tar stores every file mode 000, which Tectonic reports
# only as "Permission denied"). /tmp is cleared between sessions, so re-run this after a restart.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY=${PY:-/Users/jubayer/.claude-science/conda/envs/python/bin/python}
TECTONIC=${TECTONIC:-/Users/jubayer/.claude-science/conda/envs/multiverse/bin/tectonic}
BUNDLE=${BUNDLE:-/tmp/texbundle}

if [ ! -f "$BUNDLE/SHA256SUM" ]; then
  echo "fetching TeX bundle (2.9 GB, one-off per session)"
  "$PY" - <<'PYEOF'
import requests
u = "https://relay.fullyjustified.net/default_bundle_v33.tar"
with requests.get(u, stream=True, timeout=120) as r:
    r.raise_for_status()
    with open("/tmp/tectonic_bundle.tar", "wb") as f:
        for chunk in r.iter_content(1 << 22):
            f.write(chunk)
PYEOF
  mkdir -p "$BUNDLE" && tar -xf /tmp/tectonic_bundle.tar -C "$BUNDLE" && chmod -R u+rwX "$BUNDLE"
fi

"$PY" "$ROOT/preprint/make_tables.py"
cd "$ROOT/preprint"
HOME=/tmp/tecthome TECTONIC_CACHE_DIR=/tmp/tectonic-cache "$TECTONIC" --bundle "$BUNDLE" --outdir . manuscript.tex
