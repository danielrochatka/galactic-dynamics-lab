#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "$0")/../../.." && pwd)
OUT_DIR=$(mktemp -d)
trap 'rm -rf "$OUT_DIR"' EXIT

if ! python3 - <<'PY'
import importlib.util

mods = ("matplotlib", "numpy", "pandas")
missing = [m for m in mods if importlib.util.find_spec(m) is None]
raise SystemExit(1 if missing else 0)
PY
then
  echo "plotting dependencies unavailable; skipping" >&2
  exit 0
fi

cat >"$OUT_DIR/tpf_4d_static_residual_sources.csv" <<'CSV'
x,y,z
0,0,0
CSV

for plane in xy xz yz; do
  case "$plane" in
    xy)
      cat >"$OUT_DIR/tpf_4d_static_residual_slice_xy.csv" <<'CSV'
x,y,normalized_residual,residual_spatial_norm,xi_spatial_norm,invariant_I_4d,xi_x,xi_y,is_boundary,is_near_source,used_in_summary
-1,-1,1e-2,1e-2,10,0.5,1,0,0,0,1
0,-1,2e-2,2e-2,20,0.6,1,1,0,1,0
1,-1,4e-2,4e-2,30,0.7,0,1,1,0,0
-1,0,6e-2,6e-2,40,0.8,-1,1,0,0,1
0,0,8e-2,8e-2,50,0.9,-1,0,0,0,1
1,0,1e-1,1e-1,60,1.0,0,-1,0,0,1
-1,1,5e-3,5e-3,70,1.1,1,-1,0,0,1
0,1,2e-3,2e-3,80,1.2,1,0,0,0,1
1,1,1e-3,1e-3,90,1.3,0,1,0,0,1
CSV
      ;;
    xz)
      cat >"$OUT_DIR/tpf_4d_static_residual_slice_xz.csv" <<'CSV'
x,z,normalized_residual,residual_spatial_norm,xi_spatial_norm,invariant_I_4d,xi_x,xi_z,is_boundary,is_near_source,used_in_summary
-1,-1,1e-2,1e-2,10,0.5,1,0,0,0,1
0,-1,2e-2,2e-2,20,0.6,1,1,0,1,0
1,-1,4e-2,4e-2,30,0.7,0,1,1,0,0
-1,0,6e-2,6e-2,40,0.8,-1,1,0,0,1
0,0,8e-2,8e-2,50,0.9,-1,0,0,0,1
1,0,1e-1,1e-1,60,1.0,0,-1,0,0,1
-1,1,5e-3,5e-3,70,1.1,1,-1,0,0,1
0,1,2e-3,2e-3,80,1.2,1,0,0,0,1
1,1,1e-3,1e-3,90,1.3,0,1,0,0,1
CSV
      ;;
    yz)
      cat >"$OUT_DIR/tpf_4d_static_residual_slice_yz.csv" <<'CSV'
y,z,normalized_residual,residual_spatial_norm,xi_spatial_norm,invariant_I_4d,xi_y,xi_z,is_boundary,is_near_source,used_in_summary
-1,-1,1e-2,1e-2,10,0.5,1,0,0,0,1
0,-1,2e-2,2e-2,20,0.6,1,1,0,1,0
1,-1,4e-2,4e-2,30,0.7,0,1,1,0,0
-1,0,6e-2,6e-2,40,0.8,-1,1,0,0,1
0,0,8e-2,8e-2,50,0.9,-1,0,0,0,1
1,0,1e-1,1e-1,60,1.0,0,-1,0,0,1
-1,1,5e-3,5e-3,70,1.1,1,-1,0,0,1
0,1,2e-3,2e-3,80,1.2,1,0,0,0,1
1,1,1e-3,1e-3,90,1.3,0,1,0,0,1
CSV
      ;;
  esac
done

python3 "$ROOT_DIR/engine/physics/TPFCore/plots/plot_tpf_4d_static_residual.py" "$OUT_DIR" >/dev/null

test -f "$OUT_DIR/tpf_4d_static_residual_xy_normalized_residual.png"
test -f "$OUT_DIR/tpf_4d_static_residual_xz_xi_spatial_norm.png"
test -f "$OUT_DIR/tpf_4d_static_residual_yz_xi_quiver.png"
