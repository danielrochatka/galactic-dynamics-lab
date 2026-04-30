#!/usr/bin/env bash
set -euo pipefail
ENGINE_ROOT="${ENGINE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
ROOT="$ENGINE_ROOT"
source "$ROOT/tests/integration/_env.sh"
OUT="$ROOT/../outputs/preflight_non_galaxy_ignored_smoke"
rm -rf "$OUT"
"$ROOT/galaxy_sim" bh_orbit_validation --output_dir="$OUT" --star_mass=1e40 --bh_mass=1 --n_stars=1000 --outer_radius=1 --softening=0.8 --save_run_info=true --yes

test -f "$OUT/run_info.txt"
! grep -q 'galaxy_preflight_warning_count' "$OUT/run_info.txt"
! grep -q '=== Galaxy preflight sanity check ===' "$OUT/run_info.txt"
