#!/usr/bin/env bash
set -euo pipefail
# shellcheck source=tests/integration/_env.sh
source "$(cd "$(dirname "$0")" && pwd)/_env.sh"

OUT_ROOT=$(mktemp -d)
DEV_PY="$ENGINE_ROOT/../dev/bin/python3"
DEV_PY_REAL="$ENGINE_ROOT/../dev/bin/python3.real"
RESTORE_DEV_PY=0

cleanup() {
  if [[ $RESTORE_DEV_PY -eq 1 && -f "$DEV_PY_REAL" ]]; then
    rm -f "$DEV_PY"
    mv "$DEV_PY_REAL" "$DEV_PY"
  fi
  rm -rf "$OUT_ROOT"
}
trap cleanup EXIT

NO_PNG_OUT="$OUT_ROOT/no_png_case"
LOG_NO_PNG="$OUT_ROOT/no_png.log"

if [[ -f "$DEV_PY" ]]; then
  mv "$DEV_PY" "$DEV_PY_REAL"
  RESTORE_DEV_PY=1
  cat > "$DEV_PY" <<'PY'
#!/usr/bin/env bash
exit 0
PY
  chmod +x "$DEV_PY"
fi

./galaxy_sim tpf_4d_static_residual_benchmark \
  --plot \
  --physics_package=TPFCore \
  --tpf_source_benchmark_shape=monopole \
  --tpf_source_benchmark_total_mass=1.0e12 \
  --tpf_4d_residual_grid_n=9 \
  --tpf_4d_residual_grid_half_extent=5 \
  --tpf_4d_residual_source_exclusion_radius=1.0 \
  --tpf_4d_residual_field_softening=0.1 \
  --output_dir="$NO_PNG_OUT" >"$LOG_NO_PNG" 2>&1

grep -q "plot script completed but no expected PNGs were found/generated" "$LOG_NO_PNG"
if grep -q "Generated optional PNGs in" "$LOG_NO_PNG"; then
  echo "unexpected success PNG log when no expected PNGs were generated" >&2
  exit 1
fi
test -f "$NO_PNG_OUT/tpf_4d_static_residual_slice_xy.csv"

if [[ $RESTORE_DEV_PY -eq 1 ]]; then
  rm -f "$DEV_PY"
  mv "$DEV_PY_REAL" "$DEV_PY"
  RESTORE_DEV_PY=0
fi

SPACE_OUT="$OUT_ROOT/with space"
LOG_SPACE="$OUT_ROOT/with_space.log"

./galaxy_sim tpf_4d_static_residual_benchmark \
  --plot \
  --physics_package=TPFCore \
  --tpf_source_benchmark_shape=monopole \
  --tpf_source_benchmark_total_mass=1.0e12 \
  --tpf_4d_residual_grid_n=9 \
  --tpf_4d_residual_grid_half_extent=5 \
  --tpf_4d_residual_source_exclusion_radius=1.0 \
  --tpf_4d_residual_field_softening=0.1 \
  --output_dir="$SPACE_OUT" >"$LOG_SPACE" 2>&1

grep -q "Generated optional PNGs in $SPACE_OUT" "$LOG_SPACE"
test -f "$SPACE_OUT/tpf_4d_static_residual_xy_normalized_residual.png"
