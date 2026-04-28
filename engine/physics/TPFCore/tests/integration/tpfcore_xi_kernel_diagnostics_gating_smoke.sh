#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../../../.." && pwd)"
SIM="$ROOT/engine/galaxy_sim"
OUT0="$ROOT/engine/outputs/test_xi_diag_gate_off"
OUT1="$ROOT/engine/outputs/test_xi_diag_gate_on"

rm -rf "$OUT0" "$OUT1"

"$SIM" galaxy \
  --physics_package=TPFCore \
  --tpf_dynamics_mode=xi_kernel_deformed \
  --tpfcore_enable_provisional_readout=true \
  --tpf_4d_xi_kernel_mode=scalar_beta \
  --tpf_4d_xi_kernel_coupling=0.01 \
  --n_stars=8 --n_steps=2 --snapshot_every=1 \
  --save_run_info=true \
  --output_dir="$OUT0"

test ! -f "$OUT0/tpf_regime_diagnostics.txt"
test ! -f "$OUT0/tpf_readout_debug.csv"
grep -q $'^xi_runtime_theta_evaluations\t0$' "$OUT0/run_info.txt"
grep -q $'^xi_runtime_invariant_I_evaluations\t0$' "$OUT0/run_info.txt"
grep -q $'^xi_runtime_direct_tpf_evaluations\t0$' "$OUT0/run_info.txt"
grep -q $'^xi_runtime_provisional_readout_evaluations\t0$' "$OUT0/run_info.txt"

"$SIM" galaxy \
  --physics_package=TPFCore \
  --tpf_dynamics_mode=xi_kernel_deformed \
  --tpfcore_enable_provisional_readout=true \
  --tpf_4d_xi_kernel_mode=scalar_beta \
  --tpf_4d_xi_kernel_coupling=0.01 \
  --tpf_xi_kernel_dump_field_diagnostics=true \
  --n_stars=8 --n_steps=2 --snapshot_every=1 \
  --save_run_info=true \
  --output_dir="$OUT1"

test -f "$OUT1/tpf_regime_diagnostics.txt"
