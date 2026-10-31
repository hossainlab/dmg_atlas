#!/usr/bin/env bash
# Rebuild the preprint: figure 7 + graphical abstract, tables, Supplementary Data, reference check, PDF.
#   bash preprint/scripts/build.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"; cd "$ROOT"
PY=${PY:-$HOME/.claude-science/conda/envs/multiverse/bin/python}
TECTONIC=${TECTONIC:-$HOME/.claude-science/conda/envs/latex/bin/tectonic}
"$PY" scripts/36_fig_robustness.py
"$PY" scripts/37_graphical_abstract.py
"$PY" preprint/scripts/make_tables.py
"$PY" preprint/scripts/make_supplementary_data.py
"$PY" preprint/scripts/verify_refs.py preprint/references.bib preprint/dmg_preprint.tex
cd preprint && "$TECTONIC" dmg_preprint.tex
