#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

OUT="../outputs/test_utility_run_info_provenance_and_particles"
rm -rf "$OUT"

GALAXY_RUN_CONFIG="$ROOT/../configs/smoke_test.cfg" \
./galaxy_sim tpf_single_source_inspect \
  --output_dir="$OUT" \
  --physics_package=TPFCore \
  --save_run_info=true \
  --yes >/dev/null

[[ -f "$OUT/run_info.txt" ]]
# Utility-mode provenance uses explicit paths passed into write_run_info.
grep -q "^run_config\s\+$ROOT/../configs/smoke_test.cfg$" "$OUT/run_info.txt"
grep -q "^package_defaults\s\+physics/Newtonian/defaults.cfg$" "$OUT/run_info.txt"
# Particleless utility convention: resolver has no dynamic particle state; run_info records n_stars=0.
grep -q "^n_stars\s\+0$" "$OUT/run_info.txt"
