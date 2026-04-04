#!/usr/bin/env bash
set -euo pipefail
if [[ -z "${CPP_SIM_ROOT:-}" ]]; then
  export CPP_SIM_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
fi
source "$CPP_SIM_ROOT/tests/integration/_env.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
./galaxy_sim galaxy --output_dir="$OUT" --physics_package=TPFCore \
  --tpfcore_enable_provisional_readout=true \
  --tpf_gdd_coupling=2.5e-42 \
  --n_stars=10 --n_steps=1 --snapshot_every=1 --save_run_info=true
grep -E '^tpf_vdsg_coupling[[:space:]]' "$OUT/run_info.txt" | grep -qE '2\.5e-42|2\.5e\-042'
