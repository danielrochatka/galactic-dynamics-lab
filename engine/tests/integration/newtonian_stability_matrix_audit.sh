#!/usr/bin/env bash
set -euo pipefail
# shellcheck source=tests/integration/_env.sh
source "$(cd "$(dirname "$0")" && pwd)/_env.sh"

OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

COMMON=(--physics_package=Newtonian --n_stars=64 --bh_mass=12000 --dt=0.02 --n_steps=20 --snapshot_every=5 --softening=0.3 --galaxy_init_seed=4242)

run_case() {
  local name="$1"; shift
  ./galaxy_sim galaxy --output_dir="$OUT/$name" "${COMMON[@]}" --yes "$@" --yes >/dev/null
}

run_case A_single_preformed_starstar_true --galaxy_init_template=preformed_spiral --enable_star_star_gravity=true
run_case B_single_symmetric_starstar_true --galaxy_init_template=symmetric_disk --enable_star_star_gravity=true
run_case C_single_preformed_starstar_false --galaxy_init_template=preformed_spiral --enable_star_star_gravity=false
./galaxy_sim galaxy --output_dir="$OUT/D_compare" "${COMMON[@]}" --galaxy_init_template=preformed_spiral --enable_star_star_gravity=true --physics_package_compare=TPFCore --tpf_dynamics_mode=xi_kernel_deformed --yes >/dev/null

python3 - <<'PY' "$OUT"
import csv, math, pathlib, statistics, sys
root = pathlib.Path(sys.argv[1])

def load_snap(p):
    lines=p.read_text().splitlines()
    hdr=None
    for i,line in enumerate(lines):
        cols=[c.strip().lower() for c in line.split(',')]
        if {'i','x','y','vx','vy','mass'}.issubset(set(cols)):
            hdr=i
            break
    if hdr is None:
        raise RuntimeError(f'missing header in {p}')
    rows=[]
    import io
    buf=io.StringIO('\n'.join(lines[hdr:]))
    r=csv.DictReader(buf)
    for row in r:
        rows.append({k: (int(v) if k=='i' else float(v)) for k,v in row.items()})
    return rows

def metrics(run_dir):
    files=sorted(run_dir.glob('snapshot_*.csv'))
    s0=load_snap(files[0]); sf=load_snap(files[-1])
    r0=[math.hypot(r['x'],r['y']) for r in s0]; rf=[math.hypot(r['x'],r['y']) for r in sf]
    m=[r['mass'] for r in sf]
    mt=sum(m)
    comx=sum(r['x']*r['mass'] for r in sf)/mt; comy=sum(r['y']*r['mass'] for r in sf)/mt
    comvx=sum(r['vx']*r['mass'] for r in sf)/mt; comvy=sum(r['vy']*r['mass'] for r in sf)/mt
    lz0=sum(r['mass']*(r['x']*r['vy']-r['y']*r['vx']) for r in s0)
    lzf=sum(r['mass']*(r['x']*r['vy']-r['y']*r['vx']) for r in sf)
    rv=[]
    for r in sf:
        rad=max(1e-30, math.hypot(r['x'],r['y']))
        rv.append((r['x']*r['vx']+r['y']*r['vy'])/rad)
    dr=[abs(a-b) for a,b in zip(rf,r0)]
    return dict(r0_p50=statistics.median(r0), rf_p50=statistics.median(rf), com=math.hypot(comx,comy), comv=math.hypot(comvx,comvy), lz_drift=(lzf-lz0)/(abs(lz0)+1e-30), med_rv=statistics.median(rv), med_dr=statistics.median(dr), max_dr=max(dr), snaps=len(files), final_step=files[-1].stem.split('_')[-1])

cases={
'A': root/'A_single_preformed_starstar_true',
'B': root/'B_single_symmetric_starstar_true',
'C': root/'C_single_preformed_starstar_false',
'D': root/'D_compare'/'left_Newtonian',
}
for k,p in cases.items():
    m=metrics(p)
    print(k, m)
PY
