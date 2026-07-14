#!/usr/bin/env bash
# dev_couple_sph_cfd validation harness - SPH<->CFD plume/surface coupling.
#
#   ./validation/run.sh
#
# The gate delegates to examples/plume_surface/sweep.py, which runs the live
# coupled example, checks measured U_mf against DEM-CFD and Wen-Yu references,
# verifies the negative controls, and regenerates the committed SVG figure from
# the same parsed output.
set -euo pipefail
cd "$(dirname "$0")/.."

PY="${BENCH_PYTHON:-python3}"

echo "=== dev_couple_sph_cfd validation set ==="
"$PY" examples/psi_evidence/test_fail_closed.py
"$PY" examples/plume_surface/sweep.py
echo "=== all dev_couple_sph_cfd validation gates passed ==="
