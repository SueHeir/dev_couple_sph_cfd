#!/usr/bin/env bash
# dev_couple_sph_cfd executable smoke harness - imposed-flow packed-bed seam.
#
#   ./validation/run.sh
#
# The harness runs the live coupled example. Completion proves only that this
# configured cross-substrate path executed; it does not validate a plume or crater.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "=== dev_couple_sph_cfd packed-bed seam smoke ==="
cargo run --release --example packed_bed_seam -- examples/packed_bed_seam/config.toml
echo "=== packed-bed seam smoke completed (not physical validation) ==="
