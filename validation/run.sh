#!/usr/bin/env bash
# dev_couple_sph_cfd regression harness - imposed-flow packed-bed seam.
#
#   ./validation/run.sh
#
# The gate delegates to examples/plume_surface/sweep.py, which runs the live
# coupled example, checks frozen local regression limits and fault sensitivity,
# and regenerates a diagnostic SVG. It does not validate a plume or crater.
set -euo pipefail
cd "$(dirname "$0")/.."

PY="${BENCH_PYTHON:-python3}"

echo "=== dev_couple_sph_cfd packed-bed seam regression ==="
"$PY" examples/plume_surface/sweep.py
echo "=== packed-bed seam regression passed ==="
