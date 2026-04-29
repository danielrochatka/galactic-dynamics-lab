#!/usr/bin/env bash
set -euo pipefail
# shellcheck source=tests/integration/_env.sh
source "$(cd "$(dirname "$0")" && pwd)/_env.sh"

OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

COMMON_ARGS=(
  --n_stars=48
  --star_mass=1.0
  --bh_mass=10000.0
  --galaxy_radius=120.0
  --dt=0.02
  --n_steps=16
  --snapshot_every=4
  --softening=0.25
  --enable_star_star_gravity=true
  --galaxy_init_template=preformed_spiral
  --galaxy_init_seed=20260429
  )

SINGLE="$OUT/single"
COMPARE="$OUT/compare"
./galaxy_sim galaxy --output_dir="$SINGLE" --physics_package=Newtonian "${COMMON_ARGS[@]}"
./galaxy_sim galaxy --output_dir="$COMPARE" --physics_package=Newtonian --physics_package_compare=TPFCore --tpfcore_enable_provisional_readout=true "${COMMON_ARGS[@]}"

python3 - <<'PY' "$SINGLE" "$COMPARE"
import csv, json, math, pathlib, sys
single = pathlib.Path(sys.argv[1])
compare = pathlib.Path(sys.argv[2])
left = compare / "left_Newtonian"
right = compare / "right_TPFCore"

def load_snap(path):
    rows=[]
    with open(path) as f:
        for line in f:
            parts=[p.strip() for p in line.strip().split(',')]
            if len(parts) < 6:
                continue
            try:
                vals=tuple(float(v) for v in parts[1:6])
            except ValueError:
                continue
            rows.append(vals)
    return rows

def max_abs_delta(a,b):
    m = 0.0
    for ra, rb in zip(a,b):
        for va,vb in zip(ra,rb):
            m=max(m,abs(va-vb))
    return m

single0 = load_snap(single / 'snapshot_00000.csv')
left0 = load_snap(left / 'snapshot_00000.csv')
right0 = load_snap(right / 'snapshot_00000.csv')
assert single0 == left0, 'single step0 != compare-left step0'
assert left0 == right0, 'compare left/right step0 mismatch'

singlef = load_snap(single / 'snapshot_00016.csv')
leftf = load_snap(left / 'snapshot_00016.csv')
assert len(singlef)==len(leftf)
assert max_abs_delta(singlef,leftf) < 1e-12, 'single final != compare-left final'

manifest = json.loads((compare/'compare_manifest.json').read_text())
assert manifest['ic_fingerprint_fnv1a64']

def meta(run):
    d={}
    for line in (run/'run_info.txt').read_text().splitlines():
        if '=' in line:
            k,v = line.split('=',1)
            d[k.strip()] = v.strip()
    return d
sm = meta(single)
lm = meta(left)
for k in ('n_steps','n_snapshots','effective_snapshot_every','configured_dt'):
    assert sm.get(k)==lm.get(k), f'metadata mismatch {k}: {sm.get(k)} vs {lm.get(k)}'
print('ok')
PY
