#!/usr/bin/env bash
set -euo pipefail

ENGINE_ROOT="${ENGINE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
ROOT="$ENGINE_ROOT"
source "$ROOT/tests/integration/_env.sh"

OUT="$ROOT/../outputs/preflight_warning_no_yes_fails_smoke"
rm -rf "$OUT"

"$ROOT/galaxy_sim" --output_dir="$OUT" --simulation_mode=galaxy --n_stars=100 --star_mass=2 --bh_mass=1 \
  --n_steps=1 --snapshot_every=1 --dt=0.01 --outer_radius=1 --softening=0.2 --save_run_info=true

test -f "$OUT/run_info.txt"
grep -q '=== Galaxy preflight sanity check ===' "$OUT/run_info.txt"
compgen -G "$OUT/snapshot_*.csv" > /dev/null
