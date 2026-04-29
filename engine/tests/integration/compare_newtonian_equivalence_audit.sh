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
./galaxy_sim galaxy --output_dir="$SINGLE" --physics_package=Newtonian "${COMMON_ARGS[@]}" >/dev/null
./galaxy_sim galaxy --output_dir="$COMPARE" --physics_package=Newtonian --physics_package_compare=TPFCore --tpfcore_enable_provisional_readout=true "${COMMON_ARGS[@]}" >/dev/null

IMMUNE_A="$OUT/newtonian_tpf_neutral"
IMMUNE_B="$OUT/newtonian_tpf_extreme"
./galaxy_sim galaxy --output_dir="$IMMUNE_A" --physics_package=Newtonian "${COMMON_ARGS[@]}" \
  --tpf_vdsg_coupling=0 --tpf_4d_xi_kernel_mode=off --tpf_4d_xi_kernel_coupling=0 \
  --tpf_global_accel_shunt_enable=false --tpf_cooling_fraction=0 >/dev/null
./galaxy_sim galaxy --output_dir="$IMMUNE_B" --physics_package=Newtonian "${COMMON_ARGS[@]}" \
  --tpf_vdsg_coupling=1e99 --tpf_4d_xi_kernel_mode=metric_transverse_wake --tpf_4d_xi_kernel_coupling=1e99 \
  --tpf_4d_xi_kernel_metric_min=1e-9 --tpf_4d_xi_kernel_metric_max=1e9 --tpf_4d_xi_temporal_mode=norm_scaled \
  --tpf_4d_xi_temporal_coupling=1e99 --tpf_global_accel_shunt_enable=true --tpf_global_accel_shunt_fraction=1e-9 \
  --tpf_cooling_fraction=0.99 --tpfcore_enable_provisional_readout=true --tpfcore_dump_readout_debug=true \
  --tpf_xi_kernel_dump_field_diagnostics=true >/dev/null

python3 - <<'PY' "$SINGLE" "$COMPARE" "$IMMUNE_A" "$IMMUNE_B"
import csv, json, pathlib, re, sys

single, compare, immune_a, immune_b = map(pathlib.Path, sys.argv[1:])
left = compare / "left_Newtonian"
right = compare / "right_TPFCore"

def parse_run_info(path: pathlib.Path) -> dict:
    out = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("===") or line.startswith("---"):
            continue
        if "\t" in line:
            key, value = line.split("\t", 1)
        elif "=" in line:
            key, value = line.split("=", 1)
        else:
            m = re.match(r"^([^:]+)\s:\s(.+)$", line)
            if not m:
                continue
            key, value = m.group(1), m.group(2)
        out[key.strip()] = value.strip()
    return out

def step_from_snapshot(path: pathlib.Path) -> int:
    return int(path.stem.split("_")[-1])

def list_snapshots(run_dir: pathlib.Path):
    snaps = sorted(run_dir.glob("snapshot_*.csv"), key=step_from_snapshot)
    if not snaps:
        raise AssertionError(f"no snapshots in {run_dir}")
    return snaps

def load_snapshot(path: pathlib.Path):
    lines = path.read_text().splitlines()
    hdr_idx = None
    for idx, line in enumerate(lines):
        cols = [c.strip().lower() for c in line.split(",")]
        if {"x", "y", "vx", "vy", "mass"}.issubset(set(cols)) and ("i" in cols or "particle_index" in cols):
            hdr_idx = idx
            break
    if hdr_idx is None:
        raise AssertionError(f"snapshot {path} missing required header row")
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for _ in range(hdr_idx):
            next(reader, None)
    with path.open(newline="") as f:
        for _ in range(hdr_idx):
            next(f)
        reader = csv.DictReader(f)
        field_map = {h.strip().lower(): h for h in (reader.fieldnames or []) if h is not None}
        req = ["x", "y", "vx", "vy", "mass"]
        missing = [k for k in req if k not in field_map]
        if missing:
            raise AssertionError(f"snapshot {path} missing columns {missing}; headers={reader.fieldnames}")
        if "i" not in field_map and "particle_index" not in field_map:
            raise AssertionError(f"snapshot {path} missing particle index column; headers={reader.fieldnames}")
        rows = []
        for row in reader:
            rows.append((float(row[field_map["x"]]), float(row[field_map["y"]]), float(row[field_map["vx"]]), float(row[field_map["vy"]]), float(row[field_map["mass"]])))
        return rows

def assert_snapshot_equal(a_path, b_path, tol=0.0):
    a = load_snapshot(a_path)
    b = load_snapshot(b_path)
    assert len(a) == len(b), f"row count mismatch: {a_path} vs {b_path}"
    for i, (ra, rb) in enumerate(zip(a, b)):
        for j, (va, vb) in enumerate(zip(ra, rb)):
            if abs(va - vb) > tol:
                raise AssertionError(f"snapshot mismatch @row {i} col {j}: {va} vs {vb} ({a_path} vs {b_path})")

# single vs compare-left equivalence + compare left/right init identity
single_snaps = list_snapshots(single)
left_snaps = list_snapshots(left)
right_snaps = list_snapshots(right)
assert step_from_snapshot(single_snaps[0]) == 0
assert step_from_snapshot(left_snaps[0]) == 0
assert step_from_snapshot(right_snaps[0]) == 0
assert_snapshot_equal(single_snaps[0], left_snaps[0], tol=0.0)
assert_snapshot_equal(left_snaps[0], right_snaps[0], tol=0.0)

single_final = single_snaps[-1]
left_final = left_snaps[-1]
assert step_from_snapshot(single_final) == step_from_snapshot(left_final), "single/left final step mismatch"
assert_snapshot_equal(single_final, left_final, tol=1e-12)

manifest = json.loads((compare / "compare_manifest.json").read_text())
if "shared_ic_fingerprint_fnv1a64" in manifest:
    fp = manifest["shared_ic_fingerprint_fnv1a64"]
elif "ic_fingerprint_fnv1a64" in manifest:
    fp = manifest["ic_fingerprint_fnv1a64"]
else:
    raise AssertionError(f"compare manifest missing fingerprint keys: {sorted(manifest.keys())}")
assert fp

smeta = parse_run_info(single / "run_info.txt")
lmeta = parse_run_info(left / "run_info.txt")
for key in ("effective_n_steps", "effective_number_of_snapshots", "effective_snapshot_every", "configured_dt", "effective_dt", "effective_n_steps_done"):
    assert smeta.get(key) == lmeta.get(key), f"run_info mismatch for {key}: {smeta.get(key)} vs {lmeta.get(key)}"

# Newtonian invariance under extreme TPF/VDSG config
na_snaps = list_snapshots(immune_a)
nb_snaps = list_snapshots(immune_b)
assert len(na_snaps) == len(nb_snaps), "immunity snapshot count mismatch"
assert step_from_snapshot(na_snaps[-1]) == step_from_snapshot(nb_snaps[-1]), "immunity final step mismatch"
assert_snapshot_equal(na_snaps[0], nb_snaps[0], tol=0.0)
assert_snapshot_equal(na_snaps[-1], nb_snaps[-1], tol=0.0)

nmeta_a = parse_run_info(immune_a / "run_info.txt")
nmeta_b = parse_run_info(immune_b / "run_info.txt")
for key in ("effective_number_of_snapshots", "effective_n_steps_done"):
    assert nmeta_a.get(key) == nmeta_b.get(key), f"immunity metadata mismatch {key}: {nmeta_a.get(key)} vs {nmeta_b.get(key)}"
for meta, tag in ((nmeta_a, "neutral"), (nmeta_b, "extreme")):
    assert meta.get("active_dynamics_branch") == "Newtonian_pairwise_G_SI", f"{tag}: unexpected branch {meta.get('active_dynamics_branch')}"
    assert meta.get("acceleration_code_path") == "NewtonianPackage::compute_accelerations", f"{tag}: unexpected path {meta.get('acceleration_code_path')}"

print("compare_newtonian_equivalence_audit: PASS")
PY
