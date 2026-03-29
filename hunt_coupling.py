#!/usr/bin/env python3
"""
Binary search for the largest stable tpf_vdsg_coupling in a bounded interval.

Stable = 100%% of stars have r <= 2e20 m in snapshot at n_steps (default 5000).

Requires: built galaxy_sim, run from repo root (or pass --cpp-sim). Uses your
configs/my.local.cfg (TPFCore + provisional readout, etc.) plus CLI overrides.

Example:
  python3 hunt_coupling.py
  python3 hunt_coupling.py --cpp-sim /path/to/cpp_sim
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path


R_ESCAPE_LIMIT = 2.0e20
STEP_TARGET = 5000


def repo_root() -> Path:
    return Path(__file__).resolve().parent


def find_galaxy_sim(cpp_sim: Path) -> Path:
    exe = cpp_sim / "galaxy_sim"
    if not exe.is_file():
        raise FileNotFoundError(f"Missing {exe} (run make -C cpp_sim first)")
    return exe


def newest_snapshot_at_step(cpp_sim: Path, step: int) -> Path | None:
    """Most recently modified snapshot_{step:05d}.csv under cpp_sim/outputs/."""
    outs = cpp_sim / "outputs"
    if not outs.is_dir():
        return None
    name = f"snapshot_{step:05d}.csv"
    best: Path | None = None
    best_mtime = 0.0
    for d in outs.iterdir():
        if not d.is_dir():
            continue
        p = d / name
        if p.is_file():
            mt = p.stat().st_mtime
            if mt > best_mtime:
                best_mtime = mt
                best = p
    return best


def snapshot_all_within_r(path: Path, r_max: float) -> tuple[bool, int, int, float]:
    """
    Return (ok, n_ok, n_total, r_worst) where ok means all r <= r_max.
    C++ format: # line, header i,x,y,vx,vy,mass, data rows.
    """
    text = path.read_text().strip().splitlines()
    if len(text) < 3:
        return False, 0, 0, float("nan")

    header = [c.strip().lower() for c in text[1].split(",")]
    rows: list[list[float]] = []
    for ln in text[2:]:
        ln = ln.strip()
        if not ln:
            continue
        parts = [p.strip() for p in ln.split(",")]
        try:
            rows.append([float(p) for p in parts])
        except ValueError:
            continue

    if not rows:
        return False, 0, 0, float("nan")

    try:
        ix = header.index("x")
        iy = header.index("y")
    except ValueError:
        if len(rows[0]) >= 5:
            xi, yi = 1, 2
        else:
            return False, 0, 0, float("nan")
    else:
        xi, yi = ix, iy

    r_worst = 0.0
    n_ok = 0
    n_total = 0
    for row in rows:
        if max(xi, yi) >= len(row):
            continue
        x = row[xi]
        y = row[yi]
        r = math.hypot(x, y)
        n_total += 1
        if r <= r_max:
            n_ok += 1
        if r > r_worst:
            r_worst = r
    return n_ok == n_total and n_total > 0, n_ok, n_total, r_worst


def run_once(
    cpp_sim: Path,
    exe: Path,
    coupling: float,
    n_steps: int,
    snapshot_every: int,
) -> tuple[bool, str]:
    """Run galaxy_sim; return (stable, status_message)."""
    cmd = [
        str(exe),
        "galaxy",
        f"--tpf_vdsg_coupling={coupling}",
        f"--n_steps={n_steps}",
        f"--snapshot_every={snapshot_every}",
    ]
    proc = subprocess.run(
        cmd,
        cwd=str(cpp_sim),
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        return False, f"galaxy_sim exit {proc.returncode}\n{proc.stderr[-2000:]}"

    snap = newest_snapshot_at_step(cpp_sim, n_steps)
    if snap is None:
        return False, f"No snapshot_{n_steps:05d}.csv under {cpp_sim / 'outputs'} (check snapshot_every divides n_steps)"

    ok, n_ok, n_tot, r_worst = snapshot_all_within_r(snap, R_ESCAPE_LIMIT)
    if not ok:
        return False, f"unstable: {n_ok}/{n_tot} within r<={R_ESCAPE_LIMIT:g} m, worst_r={r_worst:g} ({snap})"
    return True, f"stable {n_tot}/{n_tot}, worst_r={r_worst:g} ({snap})"


def main() -> int:
    parser = argparse.ArgumentParser(description="Binary search stable tpf_vdsg_coupling")
    parser.add_argument("--cpp-sim", type=Path, default=repo_root() / "cpp_sim", help="Path to cpp_sim directory")
    parser.add_argument("--low", type=float, default=17000.0)
    parser.add_argument("--high", type=float, default=18500.0)
    parser.add_argument("--tol", type=float, default=5.0, help="Stop when high - low < this")
    parser.add_argument("--n-steps", type=int, default=STEP_TARGET)
    parser.add_argument(
        "--snapshot-every",
        type=int,
        default=50,
        help="Must divide n_steps so snapshot at final step exists",
    )
    args = parser.parse_args()

    cpp_sim = args.cpp_sim.resolve()
    try:
        exe = find_galaxy_sim(cpp_sim)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        return 1

    n_steps = args.n_steps
    if n_steps % args.snapshot_every != 0:
        print(
            f"Error: n_steps={n_steps} not divisible by snapshot_every={args.snapshot_every}",
            file=sys.stderr,
        )
        return 1

    low, high = float(args.low), float(args.high)
    if high <= low:
        print("Error: need high > low", file=sys.stderr)
        return 1

    print(f"Checking low={low:g} (sanity)...")
    ok_low, msg_low = run_once(cpp_sim, exe, low, n_steps, args.snapshot_every)
    print(f"  {msg_low}")
    if not ok_low:
        print("Error: lower bound is not stable; widen range or fix physics.", file=sys.stderr)
        return 1

    best_stable = low
    tol = float(args.tol)

    while high - low >= tol:
        mid = 0.5 * (low + high)
        print(f"Trying mid={mid:.6g} (bracket [{low:g}, {high:g}])...")
        ok, msg = run_once(cpp_sim, exe, mid, n_steps, args.snapshot_every)
        print(f"  {msg}")
        if ok:
            best_stable = mid
            low = mid
        else:
            high = mid

    print()
    print(f"Highest stable alpha found (within tol={tol:g}): {best_stable:.6g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
