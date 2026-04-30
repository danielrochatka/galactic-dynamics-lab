#!/usr/bin/env bash
set -euo pipefail
ENGINE_ROOT="${ENGINE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
ROOT="$ENGINE_ROOT"
source "$ROOT/tests/integration/_env.sh"
OUT="$ROOT/../outputs/preflight_compare_yes_bypass_smoke"
rm -rf "$OUT"
"$ROOT/galaxy_sim" galaxy --output_dir="$OUT" --physics_package=Newtonian --physics_package_compare=TPFCore \
  --tpfcore_enable_provisional_readout=true --n_stars=24 --n_steps=2 --snapshot_every=1 --yes

test -d "$OUT/left_Newtonian"
test -d "$OUT/right_TPFCore"
