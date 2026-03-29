#!/usr/bin/env bash
# Three 6,000-star galaxy runs: Newtonian control, pure VDSG (no cooling), VDSG + cooling.
# Requires: cpp_sim/galaxy_sim built; python3 + matplotlib + pandas + numpy for analyze_coherence.py
#
# Note: ACTIVE PHYSICS LEDGER (TPF branch auditor) prints to stderr only when physics_package=TPFCore.
#       Run 1 uses Newtonian — there is no TPF ledger line (expected).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIM="${ROOT}/cpp_sim/galaxy_sim"
OUT_BASE="${ROOT}/cpp_sim/outputs/pure_vdsg_audit"
LOG="${OUT_BASE}/run_pure_vdsg_tests.log"

N_STARS=6000
N_STEPS="${N_STEPS:-50000}"
SNAPSHOT_EVERY="${SNAPSHOT_EVERY:-500}"

if [[ ! -x "$SIM" ]]; then
  echo "error: build cpp_sim first: (cd cpp_sim && make)" >&2
  exit 1
fi

mkdir -p "$OUT_BASE"
: >"$LOG"

last_snapshot() {
  local dir="$1"
  local f
  f="$(ls -1 "$dir"/snapshot_*.csv 2>/dev/null | sort -V | tail -1 || true)"
  if [[ -z "$f" ]]; then
    echo "error: no snapshot_*.csv in $dir" >&2
    return 1
  fi
  echo "$f"
}

run_one() {
  local tag="$1"
  shift
  echo "" | tee -a "$LOG"
  echo "================================================================" | tee -a "$LOG"
  echo " $tag" | tee -a "$LOG"
  echo "================================================================" | tee -a "$LOG"
  (cd "${ROOT}/cpp_sim" && ./galaxy_sim galaxy "$@") 2>&1 | tee -a "$LOG"
}

# Writes coherence_diagnostic.png inside the run directory (data containment).
coherence_plot() {
  local run_dir="$1"
  local snap="$2"
  local png="${run_dir}/coherence_diagnostic.png"
  python3 "${ROOT}/analyze_coherence.py" --snapshot "$snap" --plot "$png"
  echo "Coherence diagnostic plot: ${png}" | tee -a "$LOG"
}

echo "Logging to $LOG"
echo "n_stars=$N_STARS n_steps=$N_STEPS snapshot_every=$SNAPSHOT_EVERY (override with env N_STEPS / SNAPSHOT_EVERY)"

# --- Run 1: Newtonian, no cooling, no VDSG ---
R1="${OUT_BASE}/run1_newtonian_control"
mkdir -p "$R1"
run_one "RUN 1: Newtonian control (no TDSG, no cooling; no TPF branch ledger on stderr)" \
  --output_dir="${R1}" \
  --physics_package=Newtonian \
  --n_stars="${N_STARS}" \
  --n_steps="${N_STEPS}" \
  --snapshot_every="${SNAPSHOT_EVERY}" \
  --tpf_vdsg_coupling=0 \
  --tpf_cooling_fraction=0

S1="$(last_snapshot "$R1")"
coherence_plot "$R1" "$S1"
P1="${R1}/coherence_diagnostic.png"

# --- Run 2: TPFCore + VDSG, no cooling ---
R2="${OUT_BASE}/run2_pure_vdsg_no_cooling"
mkdir -p "$R2"
run_one "RUN 2: Pure VDSG (tpf_vdsg_coupling=18497.1, tpf_cooling_fraction=0) — expect ACTIVE PHYSICS LEDGER on stderr" \
  --output_dir="${R2}" \
  --physics_package=TPFCore \
  --tpfcore_enable_provisional_readout=true \
  --n_stars="${N_STARS}" \
  --n_steps="${N_STEPS}" \
  --snapshot_every="${SNAPSHOT_EVERY}" \
  --tpf_vdsg_coupling=18497.1 \
  --tpf_cooling_fraction=0

S2="$(last_snapshot "$R2")"
coherence_plot "$R2" "$S2"
P2="${R2}/coherence_diagnostic.png"

# --- Run 3: TPFCore + VDSG + cooling ---
R3="${OUT_BASE}/run3_vdsg_with_cooling"
mkdir -p "$R3"
run_one "RUN 3: VDSG + cooling (reference-style; tpf_cooling_fraction=0.2) — expect ACTIVE PHYSICS LEDGER on stderr" \
  --output_dir="${R3}" \
  --physics_package=TPFCore \
  --tpfcore_enable_provisional_readout=true \
  --n_stars="${N_STARS}" \
  --n_steps="${N_STEPS}" \
  --snapshot_every="${SNAPSHOT_EVERY}" \
  --tpf_vdsg_coupling=18497.1 \
  --tpf_cooling_fraction=0.2

S3="$(last_snapshot "$R3")"
coherence_plot "$R3" "$S3"
P3="${R3}/coherence_diagnostic.png"

echo "" | tee -a "$LOG"
echo "Done (all artifacts under ${OUT_BASE}/)." | tee -a "$LOG"
echo "  Run 1 snapshot: $S1" | tee -a "$LOG"
echo "  Run 1 plot:     $P1" | tee -a "$LOG"
echo "  Run 2 snapshot: $S2" | tee -a "$LOG"
echo "  Run 2 plot:     $P2" | tee -a "$LOG"
echo "  Run 3 snapshot: $S3" | tee -a "$LOG"
echo "  Run 3 plot:     $P3" | tee -a "$LOG"
echo "  Full log:       $LOG" | tee -a "$LOG"
echo ""
echo "Tip: grep 'ACTIVE PHYSICS LEDGER' $LOG  # runs 2–3 only (TPFCore)"
