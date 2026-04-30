#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
source "$ROOT/tests/integration/_env.sh"

OUT="$ROOT/../outputs/preflight_yes_bypass_smoke"
rm -rf "$OUT"

"$ROOT/galaxy_sim" --output_dir="$OUT" --simulation_mode=galaxy --n_stars=100 --star_mass=2 --bh_mass=1 \
  --n_steps=1 --dt=0.01 --outer_radius=1 --softening=0.2 --save_run_info=true --yes

test -f "$OUT/run_info.txt"
grep -q '=== Galaxy preflight sanity check ===' "$OUT/run_info.txt"
grep -q $'preflight_warning_count\t' "$OUT/run_info.txt"
