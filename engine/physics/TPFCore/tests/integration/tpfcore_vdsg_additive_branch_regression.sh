#!/usr/bin/env bash
# Regression: nonzero tpf_vdsg_coupling must surface in active_dynamics_branch
# (legacy_readout branch string includes the configured readout mode and coupling).
# acceleration_code_path string is unchanged (same pipeline; VDSG is additive on legacy_readout).
set -euo pipefail
if [[ -z "${ENGINE_ROOT:-}" ]]; then
  export ENGINE_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
fi
source "$ENGINE_ROOT/tests/integration/_env.sh"
OUT0=$(mktemp -d)
OUT1=$(mktemp -d)
trap 'rm -rf "$OUT0" "$OUT1"' EXIT

common=(galaxy --physics_package=TPFCore --tpf_dynamics_mode=legacy_readout --tpfcore_enable_provisional_readout=true
  --tpfcore_readout_mode=derived_tpf_radial_readout
  --n_stars=8 --n_steps=2 --snapshot_every=1 --save_run_info=true --yes
  --tpf_global_accel_shunt_enable=true --tpf_global_accel_shunt_fraction=1e-3
  --galaxy_init_seed=424242 --galaxy_init_template=symmetric_disk)

./galaxy_sim "${common[@]}" --output_dir="$OUT0" --tpf_vdsg_coupling=0
./galaxy_sim "${common[@]}" --output_dir="$OUT1" --tpf_vdsg_coupling=1e-5

d0=$(grep -m1 '^active_dynamics_branch' "$OUT0/run_info.txt" | cut -f2)
d1=$(grep -m1 '^active_dynamics_branch' "$OUT1/run_info.txt" | cut -f2)
echo "$d0" | grep -q '^tpf_runtime_path_tier=deprecated_legacy; tpf_dynamics_mode=legacy_readout; provisional readout; mode=derived_tpf_radial_readout; vdsg_coupling=0.000000e+00$'
echo "$d1" | grep -q '^tpf_runtime_path_tier=deprecated_legacy; tpf_dynamics_mode=legacy_readout; provisional readout; mode=derived_tpf_radial_readout; vdsg_coupling=1.000000e-05$'
test "$d0" != "$d1"
grep -q 'accumulate_vdsg_velocity_modifier' "$OUT0/run_info.txt"
grep -q $'^tpf_global_accel_shunt_enable\t1$' "$OUT0/run_info.txt"
grep -q $'^tpf_last_global_accel_shunt_enabled\t1$' "$OUT0/run_info.txt"
grep -q $'^tpf_last_shunt_events\t' "$OUT0/run_info.txt"
grep -q $'^tpf_last_frac_capped\t' "$OUT0/run_info.txt"
grep -q '^tpf_runtime_path_tier\tdeprecated_legacy$' "$OUT0/run_info.txt"
! grep -q '^tpf_runtime_path_tier\tactive_supported$' "$OUT0/run_info.txt"
grep -q '^active_dynamics_branch\ttpf_runtime_path_tier=deprecated_legacy;' "$OUT0/run_info.txt"
! grep -q '^active_dynamics_branch\ttpf_runtime_path_tier=active_supported; tpf_dynamics_mode=xi_kernel_deformed' "$OUT0/run_info.txt"
grep -q '^active_dynamics_branch\ttpf_runtime_path_tier=deprecated_legacy;' "$OUT0/render_manifest.txt"
! grep -q '^active_dynamics_branch\ttpf_runtime_path_tier=active_supported; tpf_dynamics_mode=xi_kernel_deformed' "$OUT0/render_manifest.txt"
grep -q 'accumulate_vdsg_velocity_modifier' "$OUT1/run_info.txt"
grep -q $'^tpf_global_accel_shunt_enable\t1$' "$OUT1/run_info.txt"
grep -q $'^tpf_last_global_accel_shunt_enabled\t1$' "$OUT1/run_info.txt"
grep -q $'^tpf_last_shunt_events\t' "$OUT1/run_info.txt"
grep -q $'^tpf_last_frac_capped\t' "$OUT1/run_info.txt"
grep -q '^tpf_runtime_path_tier\tdeprecated_legacy$' "$OUT1/run_info.txt"
! grep -q '^tpf_runtime_path_tier\tactive_supported$' "$OUT1/run_info.txt"
grep -q '^active_dynamics_branch\ttpf_runtime_path_tier=deprecated_legacy;' "$OUT1/run_info.txt"
! grep -q '^active_dynamics_branch\ttpf_runtime_path_tier=active_supported; tpf_dynamics_mode=xi_kernel_deformed' "$OUT1/run_info.txt"
grep -q '^active_dynamics_branch\ttpf_runtime_path_tier=deprecated_legacy;' "$OUT1/render_manifest.txt"
! grep -q '^active_dynamics_branch\ttpf_runtime_path_tier=active_supported; tpf_dynamics_mode=xi_kernel_deformed' "$OUT1/render_manifest.txt"
a0=$(grep -m1 '^acceleration_code_path' "$OUT0/run_info.txt" | cut -f2)
a1=$(grep -m1 '^acceleration_code_path' "$OUT1/run_info.txt" | cut -f2)
test "$a0" = "$a1"
