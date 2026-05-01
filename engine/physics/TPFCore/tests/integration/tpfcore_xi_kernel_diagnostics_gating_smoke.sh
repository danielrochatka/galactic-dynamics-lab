#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../../../.." && pwd)"
ENGINE="$ROOT/engine"
SIM="$ENGINE/galaxy_sim"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

fail() {
  echo "[FAIL] $*" >&2
  exit 1
}

[[ -x "$SIM" ]] || fail "Missing executable: $SIM (run: make -C engine -j4)"

cd "$ENGINE"

OUT0="$TMPDIR/xi_diag_gate_off"
OUT1="$TMPDIR/xi_diag_gate_on"
mkdir -p "$OUT0" "$OUT1"

run_xi() {
  local out_dir="$1"
  shift
  "$SIM" galaxy \
    --physics_package=TPFCore \
    --tpf_dynamics_mode=xi_kernel_deformed \
    --tpf_4d_xi_kernel_mode=scalar_beta \
    --tpf_4d_xi_kernel_coupling=0.01 \
    --n_stars=8 --n_steps=2 --snapshot_every=1 \
    --save_run_info=true \
    --output_dir="$out_dir" \
    "$@"
}

echo "[INFO] Running default-off Xi diagnostics case..."
run_xi "$OUT0"

[[ ! -f "$OUT0/tpf_regime_diagnostics.txt" ]] || fail "default-off unexpectedly wrote tpf_regime_diagnostics.txt"

grep -q $'^xi_runtime_theta_evaluations\t0$' "$OUT0/run_info.txt" || fail "theta counter was not zero"
grep -q $'^xi_runtime_invariant_I_evaluations\t0$' "$OUT0/run_info.txt" || fail "invariant_I counter was not zero"
grep -q $'^xi_runtime_direct_tpf_evaluations\t0$' "$OUT0/run_info.txt" || fail "direct_tpf counter was not zero"
grep -Eq $'^xi_last_call_pair_evaluations\t[1-9][0-9]*$' "$OUT0/run_info.txt" || fail "xi_last_call_pair_evaluations was not > 0"
grep -Eq $'^xi_total_pair_evaluations\t[1-9][0-9]*$' "$OUT0/run_info.txt" || fail "xi_total_pair_evaluations was not > 0"
[[ "$(grep -c '^xi_runtime_theta_evaluations' "$OUT0/run_info.txt")" -eq 1 ]] || fail "theta counter duplicated in run_info"
[[ "$(grep -c '^xi_runtime_invariant_I_evaluations' "$OUT0/run_info.txt")" -eq 1 ]] || fail "invariant_I counter duplicated in run_info"
[[ "$(grep -c '^xi_runtime_direct_tpf_evaluations' "$OUT0/run_info.txt")" -eq 1 ]] || fail "direct_tpf counter duplicated in run_info"
[[ "$(grep -c '^xi_last_call_pair_evaluations' "$OUT0/run_info.txt")" -eq 1 ]] || fail "xi_last_call_pair counter duplicated in run_info"
[[ "$(grep -c '^xi_total_pair_evaluations' "$OUT0/run_info.txt")" -eq 1 ]] || fail "xi_total_pair counter duplicated in run_info"

echo "[INFO] Running opt-in Xi diagnostics case..."
run_xi "$OUT1" \
  --tpf_xi_kernel_dump_field_diagnostics=true

[[ ! -f "$OUT1/tpf_readout_debug.csv" ]] || fail "opt-in unexpectedly produced removed provisional readout debug csv"
[[ "$(grep -c '^xi_runtime_theta_evaluations' "$OUT1/run_info.txt")" -eq 1 ]] || fail "opt-in theta counter duplicated in run_info"
[[ "$(grep -c '^xi_runtime_invariant_I_evaluations' "$OUT1/run_info.txt")" -eq 1 ]] || fail "opt-in invariant_I counter duplicated in run_info"
[[ "$(grep -c '^xi_runtime_theta_evaluations' "$OUT1/run_info.txt")" -eq 1 ]] || fail "opt-in theta counter duplicated in run_info"
[[ "$(grep -c '^xi_runtime_invariant_I_evaluations' "$OUT1/run_info.txt")" -eq 1 ]] || fail "opt-in invariant_I counter duplicated in run_info"
[[ "$(grep -c '^xi_last_call_pair_evaluations' "$OUT1/run_info.txt")" -eq 1 ]] || fail "opt-in xi_last_call_pair counter duplicated in run_info"
[[ "$(grep -c '^xi_total_pair_evaluations' "$OUT1/run_info.txt")" -eq 1 ]] || fail "opt-in xi_total_pair counter duplicated in run_info"

echo "[PASS] xi_kernel_deformed diagnostics gating smoke test passed."
