import argparse
import csv
import glob
import io
import json
import math
import os
import tempfile

G = 6.67430e-11
REQUIRED_COLS = ("x", "y", "vx", "vy", "mass")

CASES = [
    {
        "case": "case1_single_source",
        "newtonian_dir": "outputs/val_case1_newt",
        "tpf_dir": "outputs/val_case1_tpf",
        "bh_mass": 1.98847e36,
        "softening": 0.0,
        "star_star": False,
    },
    {
        "case": "case2_galaxy_starstar_false",
        "newtonian_dir": "outputs/val_case2_newt",
        "tpf_dir": "outputs/val_case2_tpf",
        "bh_mass": 1.98847e36,
        "softening": 1e15,
        "star_star": False,
    },
    {
        "case": "case3_smallN_starstar_true",
        "newtonian_dir": "outputs/val_case3_newt",
        "tpf_dir": "outputs/val_case3_tpf",
        "bh_mass": 1.98847e35,
        "softening": 5e14,
        "star_star": True,
    },
]


def _clean_csv_lines(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"No non-comment CSV content found in snapshot: {path}")
    return lines


def read_snap(path):
    lines = _clean_csv_lines(path)
    header_idx = None
    for i, line in enumerate(lines):
        header = [c.strip() for c in line.split(",")]
        if all(col in header for col in REQUIRED_COLS):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"Snapshot missing required header columns {REQUIRED_COLS}: {path}")

    reader = csv.DictReader(io.StringIO("\n".join(lines[header_idx:])))
    fieldnames = reader.fieldnames or []
    missing = [c for c in REQUIRED_COLS if c not in fieldnames]
    if missing:
        raise ValueError(f"Snapshot header missing required columns {missing}: {path}")

    rows = []
    for row_idx, row in enumerate(reader, start=1):
        parsed = {}
        for col in REQUIRED_COLS:
            raw = row.get(col, None)
            if raw is None or str(raw).strip() == "":
                raise ValueError(f"Missing value for column '{col}' in {path} at data row {row_idx}")
            try:
                parsed[col] = float(raw)
            except ValueError as exc:
                raise ValueError(
                    f"Non-numeric value for column '{col}' in {path} at data row {row_idx}: {raw!r}"
                ) from exc
        rows.append(parsed)

    if not rows:
        raise ValueError(f"No data rows parsed from snapshot: {path}")
    return rows


def accel(state, bh_mass, soft, star_star):
    n = len(state)
    ax = [0.0] * n
    ay = [0.0] * n
    eps2 = soft * soft
    for i, a in enumerate(state):
        dx = -a["x"]
        dy = -a["y"]
        r2 = dx * dx + dy * dy + eps2
        inv = r2 ** -1.5
        ax[i] += G * bh_mass * dx * inv
        ay[i] += G * bh_mass * dy * inv
    if star_star:
        for i, a in enumerate(state):
            for j, b in enumerate(state):
                if i == j:
                    continue
                dx = b["x"] - a["x"]
                dy = b["y"] - a["y"]
                r2 = dx * dx + dy * dy + eps2
                inv = r2 ** -1.5
                ax[i] += G * b["mass"] * dx * inv
                ay[i] += G * b["mass"] * dy * inv
    return ax, ay


def compare_case(case_cfg):
    name = case_cfg["case"]
    nfiles = sorted(glob.glob(f"{case_cfg['newtonian_dir']}/snapshot_*.csv"))
    tfiles = sorted(glob.glob(f"{case_cfg['tpf_dir']}/snapshot_*.csv"))
    if not nfiles:
        raise FileNotFoundError(f"{name}: no Newtonian snapshots found in {case_cfg['newtonian_dir']}")
    if not tfiles:
        raise FileNotFoundError(f"{name}: no TPF snapshots found in {case_cfg['tpf_dir']}")
    if len(nfiles) != len(tfiles):
        raise RuntimeError(f"{name}: snapshot count mismatch Newtonian={len(nfiles)} TPF={len(tfiles)}")

    max_ad = max_ar = max_pd = max_vd = 0.0
    rows_per_snapshot = None
    total_rows = 0

    for k, (npath, tpath) in enumerate(zip(nfiles, tfiles)):
        ns = read_snap(npath)
        ts = read_snap(tpath)
        if len(ns) != len(ts):
            raise RuntimeError(
                f"{name}: particle row count mismatch at snapshot index {k}: Newtonian={len(ns)} TPF={len(ts)}"
            )
        if len(ns) == 0:
            raise RuntimeError(f"{name}: parsed zero rows at snapshot index {k}")
        if rows_per_snapshot is None:
            rows_per_snapshot = len(ns)

        total_rows += len(ns)
        nax, nay = accel(ns, case_cfg["bh_mass"], case_cfg["softening"], case_cfg["star_star"])
        tax, tay = accel(ts, case_cfg["bh_mass"], case_cfg["softening"], case_cfg["star_star"])

        for i in range(len(ns)):
            d = math.hypot(tax[i] - nax[i], tay[i] - nay[i])
            nmag = math.hypot(nax[i], nay[i])
            max_ad = max(max_ad, d)
            if nmag > 0:
                max_ar = max(max_ar, d / nmag)
            max_pd = max(max_pd, math.hypot(ts[i]["x"] - ns[i]["x"], ts[i]["y"] - ns[i]["y"]))
            max_vd = max(max_vd, math.hypot(ts[i]["vx"] - ns[i]["vx"], ts[i]["vy"] - ns[i]["vy"]))

    if total_rows == 0:
        raise RuntimeError(f"{name}: no parsed rows across all snapshots")

    return {
        "case": name,
        "newtonian_dir": case_cfg["newtonian_dir"],
        "tpf_dir": case_cfg["tpf_dir"],
        "newtonian_snapshot_count": len(nfiles),
        "tpf_snapshot_count": len(tfiles),
        "rows_per_snapshot": rows_per_snapshot,
        "snapshots_compared": len(nfiles),
        "max_abs_accel_diff": max_ad,
        "max_rel_accel_diff": max_ar,
        "max_pos_diff": max_pd,
        "max_vel_diff": max_vd,
        "pass": True,
    }


def run_self_test():
    with tempfile.TemporaryDirectory() as td:
        good = os.path.join(td, "good.csv")
        with open(good, "w", encoding="utf-8") as f:
            f.write("# step,0,time,0\n\n")
            f.write("i,x,y,vx,vy,mass\n")
            f.write("0,1.0,2.0,3.0,4.0,5.0\n")
        rows = read_snap(good)
        assert len(rows) == 1 and rows[0]["x"] == 1.0 and rows[0]["mass"] == 5.0

        empty_case = {
            "case": "selftest_empty",
            "newtonian_dir": os.path.join(td, "missing_newt"),
            "tpf_dir": os.path.join(td, "missing_tpf"),
            "bh_mass": 1.0,
            "softening": 0.0,
            "star_star": False,
        }
        try:
            compare_case(empty_case)
            raise AssertionError("compare_case should fail when snapshots are missing")
        except FileNotFoundError:
            pass

    print("Self-test PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true", help="run parser/compare sanity checks")
    args = parser.parse_args()

    if args.self_test:
        run_self_test()
        return

    results = [compare_case(case_cfg) for case_cfg in CASES]
    os.makedirs("validate/geodesic_correspondence_newtonian/artifacts", exist_ok=True)
    out = "validate/geodesic_correspondence_newtonian/artifacts/metrics.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
