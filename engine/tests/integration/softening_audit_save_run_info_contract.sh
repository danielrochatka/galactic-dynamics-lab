#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
ENGINE="$ROOT/engine"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT
cd "$ENGINE"
./galaxy_sim --simulation_mode=galaxy --physics_package=Newtonian --n_stars=4 --n_steps=1 --snapshot_every=1 \
  --softening_audit_enable=true --save_run_info=false --save_snapshots=false --output_dir="$OUT" >/dev/null
test -f "$OUT/softening_audit.txt"
test ! -f "$OUT/run_info.txt"
echo "softening_audit_save_run_info_contract: PASS"
