#!/usr/bin/env bash
set -euo pipefail
if [[ -z "${ENGINE_ROOT:-}" ]]; then
  export ENGINE_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
fi
source "$ENGINE_ROOT/tests/integration/_env.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
./galaxy_sim galaxy --output_dir="$OUT" --physics_package=TPFCore \
  --tpf_dynamics_mode=tpf_xi_theta_v1 --yes \
  --tpf_vdsg_coupling=0.0 \
  --n_stars=15 --n_steps=3 --snapshot_every=1 --save_run_info=true --yes
test -f "$OUT/render_manifest.json"
grep -q '"active_dynamics_branch"' "$OUT/render_manifest.json"
grep -q 'tpf_xi_theta_v1' "$OUT/render_manifest.json"
! grep -q 'legacy_readout' "$OUT/render_manifest.json"
! grep -q 'provisional_readout' "$OUT/render_manifest.json"
! grep -q 'tpfcore_enable_provisional_readout": true' "$OUT/render_manifest.json"
