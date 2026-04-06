#!/usr/bin/env bash
set -euo pipefail
# shellcheck source=tests/integration/_env.sh
source "$(cd "$(dirname "$0")" && pwd)/_env.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
./galaxy_sim galaxy --output_dir="$OUT" --physics_package=Newtonian \
  --n_stars=30 --n_steps=8 --snapshot_every=4 --save_run_info=true
test -f "$OUT/run_info.txt"
test -n "$(find "$OUT" -maxdepth 1 -name 'snapshot_*.csv' -print -quit)"
grep -q $'\t' "$OUT/run_info.txt"
