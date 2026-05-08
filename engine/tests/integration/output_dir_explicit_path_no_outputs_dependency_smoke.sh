#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

TMP_BASE="$(mktemp -d)"
trap 'rm -rf "$TMP_BASE"; rm -f ../outputs' EXIT
OUT="$TMP_BASE/custom/outdir"

rm -rf ../outputs
: > ../outputs

./galaxy_sim galaxy --output_dir="$OUT" --n_steps=0 --save_run_info=true --yes >/dev/null

[[ -d "$OUT" ]]
[[ -f "$OUT/run_info.txt" ]]
