#!/usr/bin/env bash
set -euo pipefail
if [[ -z "${ENGINE_ROOT:-}" ]]; then
  export ENGINE_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
fi
source "$ENGINE_ROOT/tests/integration/_env.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
LOG="$OUT/run.log"
set +e
./galaxy_sim galaxy --output_dir="$OUT" --physics_package=TPFCore \
  --tpf_dynamics_mode=tpf_xi_theta_v1 --yes \
  --tpf_gdd_coupling=2.5e-42 \
  --n_stars=10 --n_steps=1 --snapshot_every=1 --save_run_info=true --yes >"$LOG" 2>&1
status=$?
set -e

if [[ $status -eq 0 ]]; then
  echo "Expected --tpf_gdd_coupling to be rejected, but command succeeded." >&2
  exit 1
fi
grep -q "Unknown CLI config key: tpf_gdd_coupling" "$LOG"
