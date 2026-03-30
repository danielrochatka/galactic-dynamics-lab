#!/usr/bin/env bash
set -euo pipefail
# cpp_sim/physics/TPFCore/tests/integration -> four levels up to cpp_sim
if [[ -z "${CPP_SIM_ROOT:-}" ]]; then
  export CPP_SIM_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
fi
source "$CPP_SIM_ROOT/tests/integration/_env.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
./galaxy_sim galaxy --output_dir="$OUT" --physics_package=TPFCore \
  --tpfcore_enable_provisional_readout=true \
  --tpf_vdsg_coupling=1e-18 \
  --n_stars=15 --n_steps=3 --snapshot_every=1 --save_run_info=true
test -f "$OUT/render_manifest.json"
grep -q '"active_dynamics_branch"' "$OUT/render_manifest.json"
grep -q 'TPF_readout_acceleration' "$OUT/render_manifest.json"
grep -q 'accumulate_vdsg_velocity_modifier' "$OUT/render_manifest.json"
grep -q 'apply_global_accel_magnitude_shunt' "$OUT/render_manifest.json"
