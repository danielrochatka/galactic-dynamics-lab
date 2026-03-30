#!/usr/bin/env bash
# Regression: nonzero tpf_vdsg_coupling must NOT switch active_dynamics_branch to a VDSG-only label.
# Baseline readout identity stays TPF_readout_acceleration:<mode>; acceleration_code_path gains +accumulate_vdsg_velocity_modifier.
set -euo pipefail
if [[ -z "${CPP_SIM_ROOT:-}" ]]; then
  export CPP_SIM_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi
source "$CPP_SIM_ROOT/tests/integration/_env.sh"
OUT0=$(mktemp -d)
OUT1=$(mktemp -d)
trap 'rm -rf "$OUT0" "$OUT1"' EXIT

common=(galaxy --physics_package=TPFCore --tpfcore_enable_provisional_readout=true
  --tpfcore_readout_mode=derived_tpf_radial_readout
  --n_stars=8 --n_steps=2 --snapshot_every=1 --save_run_info=true
  --galaxy_init_seed=424242 --galaxy_init_template=symmetric_disk)

./galaxy_sim "${common[@]}" --output_dir="$OUT0" --tpf_vdsg_coupling=0
./galaxy_sim "${common[@]}" --output_dir="$OUT1" --tpf_vdsg_coupling=1e-5

d0=$(grep -m1 '^active_dynamics_branch' "$OUT0/run_info.txt" | cut -f2)
d1=$(grep -m1 '^active_dynamics_branch' "$OUT1/run_info.txt" | cut -f2)
test "$d0" = "$d1"
echo "$d0" | grep -q '^TPF_readout_acceleration:derived_tpf_radial_readout$'
grep -q 'accumulate_vdsg_velocity_modifier' "$OUT1/run_info.txt"
test ! "$(grep -m1 '^acceleration_code_path' "$OUT0/run_info.txt" | cut -f2)" = "$(grep -m1 '^acceleration_code_path' "$OUT1/run_info.txt" | cut -f2)"
