#!/usr/bin/env bash
set -euo pipefail
# shellcheck source=tests/integration/_env.sh
source "$(cd "$(dirname "$0")" && pwd)/_env.sh"

OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

./galaxy_sim galaxy --output_dir="$OUT" \
  --physics_package=Newtonian \
  --physics_package_compare=TPFCore \
  --tpf_dynamics_mode=xi_kernel_deformed \
  --tpf_4d_xi_kernel_mode=scalar_beta \
  --tpf_4d_xi_kernel_coupling=1e-3 \
  --tpf_4d_xi_kernel_factor_mode=beta_power \
  --tpf_4d_xi_temporal_mode=off \
  --n_stars=24 --n_steps=8 --snapshot_every=4 --save_run_info=true

test -f "$OUT/compare_manifest.json"
test -f "$OUT/compare_manifest.txt"
test -d "$OUT/left_Newtonian"
test -d "$OUT/right_TPFCore"
test -f "$OUT/right_TPFCore/run_info.txt"
test -f "$OUT/right_TPFCore/render_manifest.txt"
test -f "$OUT/right_TPFCore/render_manifest.json"

grep -q $'tpf_dynamics_mode\txi_kernel_deformed' "$OUT/right_TPFCore/run_info.txt"
grep -q $'tpf_core_law_mode\txi_kernel_deformed' "$OUT/right_TPFCore/run_info.txt"
grep -q $'tpfcore_readout_mode_active_for_this_route\t0' "$OUT/right_TPFCore/run_info.txt"
grep -q $'tpf_vdsg_coupling_active_for_this_route\t0' "$OUT/right_TPFCore/run_info.txt"

grep -q '"tpf_core_law_mode": "xi_kernel_deformed"' "$OUT/right_TPFCore/render_manifest.json"
grep -q '"acceleration_formula": "a=-K_xi\*Xi_eff_spatial"' "$OUT/right_TPFCore/render_manifest.json"
grep -q '"provisional_readout_used": false' "$OUT/right_TPFCore/render_manifest.json"
grep -q '"principal_c_direct_tpf_used": false' "$OUT/right_TPFCore/render_manifest.json"
! grep -q 'legacy_readout_provisional' "$OUT/right_TPFCore/render_manifest.txt"
! grep -q 'direct_tpf_tensor_principal_part' "$OUT/right_TPFCore/render_manifest.txt"

if [[ -f "$OUT/right_TPFCore/galaxy_step0_accel_audit.csv" ]]; then
  grep -q 'ax_reference_direct_tpf' "$OUT/right_TPFCore/galaxy_step0_accel_audit.csv"
  ! grep -q 'ax_direct_tpf' "$OUT/right_TPFCore/galaxy_step0_accel_audit.csv"
  grep -q 'xi_kernel_deformed' "$OUT/right_TPFCore/galaxy_step0_accel_audit.csv"
fi
