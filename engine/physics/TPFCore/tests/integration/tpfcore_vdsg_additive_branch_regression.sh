#!/usr/bin/env bash
set -euo pipefail
if [[ -z "${ENGINE_ROOT:-}" ]]; then
  export ENGINE_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
fi
source "$ENGINE_ROOT/tests/integration/_env.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
set +e
./galaxy_sim galaxy --physics_package=TPFCore --tpf_dynamics_mode=direct_tpf \
  --tpfcore_enable_provisional_readout=true --n_stars=8 --n_steps=2 --snapshot_every=1 --save_run_info=true --yes \
  --output_dir="$OUT" >"$OUT/stdout.txt" 2>"$OUT/stderr.txt"
status=$?
set -e
if [[ $status -eq 0 ]]; then
  echo "expected provisional legacy gate to be rejected" >&2
  exit 1
fi
grep -q "tpfcore_enable_provisional_readout=true is deprecated and cannot be used with simulation_mode=galaxy" "$OUT/stderr.txt"
