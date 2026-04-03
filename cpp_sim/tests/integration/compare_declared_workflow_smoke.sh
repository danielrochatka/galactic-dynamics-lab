#!/usr/bin/env bash
set -euo pipefail
# shellcheck source=tests/integration/_env.sh
source "$(cd "$(dirname "$0")" && pwd)/_env.sh"

OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

./galaxy_sim galaxy --output_dir="$OUT" \
  --physics_package=Newtonian \
  --physics_package_compare=TPFCore \
  --tpfcore_enable_provisional_readout=true \
  --n_stars=24 --n_steps=8 --snapshot_every=4 --save_run_info=true

test -f "$OUT/compare_manifest.json"
test -f "$OUT/compare_manifest.txt"
test -d "$OUT/left_Newtonian"
test -d "$OUT/right_TPFCore"
test -f "$OUT/left_Newtonian/run_info.txt"
test -f "$OUT/right_TPFCore/run_info.txt"
test -f "$OUT/left_Newtonian/render_manifest.json"
test -f "$OUT/right_TPFCore/render_manifest.json"
test -f "$OUT/left_Newtonian/snapshot_00000.csv"
test -f "$OUT/right_TPFCore/snapshot_00000.csv"
cmp -s "$OUT/left_Newtonian/snapshot_00000.csv" "$OUT/right_TPFCore/snapshot_00000.csv"

