#!/usr/bin/env bash
set -euo pipefail
ENGINE_ROOT="${ENGINE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
ROOT="$ENGINE_ROOT"
source "$ROOT/tests/integration/_env.sh"
OUT="$ROOT/../outputs/preflight_compare_yes_bypass_smoke"
rm -rf "$OUT"
"$ROOT/galaxy_sim" --output_dir="$OUT" --simulation_mode=galaxy --physics_package=Newtonian --physics_package_compare=TPFCore \
  --n_stars=100 --star_mass=2 --bh_mass=1 --n_steps=1 --snapshot_every=1 --dt=0.01 --outer_radius=1 --softening=0.2 --yes

test -d "$OUT/left_Newtonian"
test -d "$OUT/right_TPFCore"
