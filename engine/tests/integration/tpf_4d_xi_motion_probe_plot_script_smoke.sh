#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "$0")/../../.." && pwd)
OUT_DIR=$(mktemp -d)
trap 'rm -rf "$OUT_DIR"' EXIT

cat >"$OUT_DIR/tpf_4d_xi_motion_probe_trajectories.csv" <<'CSV'
step,time,probe_id,x,y,z,vx,vy,vz,ax,ay,az,a_norm,xi_x,xi_y,xi_z,xi_spatial_norm,theta_trace_4d,invariant_I_4d,r_origin,radial_alignment_to_origin_inward,transverse_fraction_origin,is_near_source,escaped,valid
0,0.0,0,10.0,0.0,0.0,0.0,0.316,0.0,-0.01,0.0,0.0,0.01,1.0,0.0,0.0,1.0,0.1,0.2,10.0,1.0,0.0,0,0,1
1,0.1,0,9.999,0.0316,0.0,-0.001,0.316,0.0,-0.00999,-0.00003,0.0,0.00999,0.999,0.003,0.0,0.999,0.1,0.2,9.99905,1.0,0.0,0,0,1
2,0.2,0,9.996,0.0632,0.0,-0.002,0.3158,0.0,-0.00998,-0.00006,0.0,0.00998,0.998,0.006,0.0,0.998,0.1,0.2,9.99620,1.0,0.0,0,0,1
0,0.0,1,-10.0,0.0,0.0,0.0,-0.316,0.0,0.01,0.0,0.0,0.01,-1.0,0.0,0.0,1.0,0.1,0.2,10.0,1.0,0.0,0,0,1
1,0.1,1,-9.999,-0.0316,0.0,0.001,-0.316,0.0,0.00999,0.00003,0.0,0.00999,-0.999,-0.003,0.0,0.999,0.1,0.2,9.99905,1.0,0.0,0,0,1
2,0.2,1,-9.996,-0.0632,0.0,0.002,-0.3158,0.0,0.00998,0.00006,0.0,0.00998,-0.998,-0.006,0.0,0.998,0.1,0.2,9.99620,1.0,0.0,0,0,1
CSV

set +e
python3 "$ROOT_DIR/plot_tpf_4d_xi_motion_probe.py" "$OUT_DIR" >"$OUT_DIR/stdout.log" 2>"$OUT_DIR/stderr.log"
RET=$?
set -e

if [[ $RET -ne 0 ]]; then
  echo "plot script exited non-zero: $RET" >&2
  cat "$OUT_DIR/stderr.log" >&2
  exit 1
fi

if python3 - <<'PY'
import sys
try:
    import matplotlib.pyplot  # noqa: F401
    import numpy  # noqa: F401
    import pandas  # noqa: F401
except Exception:
    sys.exit(1)
sys.exit(0)
PY
then
  test -f "$OUT_DIR/tpf_4d_xi_motion_probe_xy_trajectories.png"
  test -f "$OUT_DIR/tpf_4d_xi_motion_probe_xz_trajectories.png"
  test -f "$OUT_DIR/tpf_4d_xi_motion_probe_yz_trajectories.png"
  test -f "$OUT_DIR/tpf_4d_xi_motion_probe_radius_vs_time.png"
  test -f "$OUT_DIR/tpf_4d_xi_motion_probe_acceleration_norm_vs_time.png"
else
  grep -q "plotting dependencies unavailable" "$OUT_DIR/stderr.log"
fi
