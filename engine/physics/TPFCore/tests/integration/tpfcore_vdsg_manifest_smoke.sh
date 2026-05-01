#!/usr/bin/env bash
set -euo pipefail
# engine/physics/TPFCore/tests/integration -> four levels up to engine
if [[ -z "${ENGINE_ROOT:-}" ]]; then
  export ENGINE_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
fi
source "$ENGINE_ROOT/tests/integration/_env.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
./galaxy_sim galaxy --output_dir="$OUT" --physics_package=TPFCore \
  --tpf_dynamics_mode=legacy_readout --tpfcore_enable_provisional_readout=true --yes \
  --tpf_vdsg_coupling=1e-18 \
  --n_stars=15 --n_steps=3 --snapshot_every=1 --save_run_info=true --yes
test -f "$OUT/render_manifest.json"
grep -q '"active_dynamics_branch"' "$OUT/render_manifest.json"
grep -q 'tpf_dynamics_mode=legacy_readout; provisional readout;' "$OUT/render_manifest.json"
grep -q 'vdsg_coupling=1.000000e-18' "$OUT/render_manifest.json"
grep -q 'accumulate_vdsg_velocity_modifier' "$OUT/render_manifest.json"
grep -Fq 'shunt OFF' "$OUT/render_manifest.json"
