#!/usr/bin/env python3
"""
Post-process C++ galaxy simulation outputs and generate visualizations.

Reads snapshot CSV files and run_info.txt from a cpp_sim run directory,
then produces the same kinds of plots as the Python pipeline:
  - initial and final scatter plots (galaxy view)
  - optional MP4/GIF animation
  - optional diagnostic time-series plots

Usage:
  python plot_cpp_run.py <run_dir> [--no-animation] [--no-diagnostics]
  python plot_cpp_run.py cpp_sim/outputs/20260308_175421

Outputs are written into the same run_dir (or run_dir/plots if you prefer;
  currently we write into run_dir to match "same run folder").

Requires: numpy, pandas, matplotlib. Optional: ffmpeg or Pillow for animation.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Import from existing project (run from repo root or with PYTHONPATH)
from render import save_static_plot, create_animation, has_ffmpeg
from diagnostics import compute_diagnostics, plot_and_save_all


# -----------------------------------------------------------------------------
# C++ output format
# -----------------------------------------------------------------------------
# run_info.txt: one key-value per line, tab-separated (e.g. "dt\t0.01").
# snapshot_00000.csv: first line "# step,0,time,0.0e+00"; second "i,x,y,vx,vy,mass"; then data rows.


@dataclass
class Snapshot:
    """One snapshot: step, time, positions (n,2), velocities (n,2)."""
    step: int
    time: float
    positions: np.ndarray
    velocities: np.ndarray


def load_run_info(run_dir: Path) -> dict[str, str | int | float]:
    """Parse run_info.txt into a dict of config values."""
    path = run_dir / "run_info.txt"
    if not path.exists():
        return {}
    info = {}
    for line in path.read_text().strip().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        key, val = parts[0].strip(), parts[1].strip()
        try:
            if "." in val:
                info[key] = float(val)
            else:
                info[key] = int(val)
        except ValueError:
            info[key] = val
    return info


def load_snapshot_csv(path: Path) -> Snapshot | None:
    """
    Load one snapshot_XXXXX.csv file.
    First line: # step,<step>,time,<time>
    Second line: i,x,y,vx,vy,mass
    Then data rows: i,x,y,vx,vy,mass
    """
    if not path.exists():
        return None
    lines = path.read_text().strip().splitlines()
    if len(lines) < 3:
        return None

    # Parse comment line: "# step,0,time,0.000000e+00"
    comment = lines[0]
    step = 0
    time = 0.0
    m = re.search(r"step\s*,\s*(\d+)", comment, re.IGNORECASE)
    if m:
        step = int(m.group(1))
    m = re.search(r"time\s*,\s*([\d.eE+-]+)", comment, re.IGNORECASE)
    if m:
        time = float(m.group(1))

    # Data: skip header (i,x,y,vx,vy,mass), parse rest as CSV
    data_lines = [ln for ln in lines[2:] if ln.strip()]
    if not data_lines:
        return None

    rows = []
    for ln in data_lines:
        parts = ln.split(",")
        if len(parts) >= 5:
            rows.append([float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])])
    arr = np.array(rows)
    positions = arr[:, :2]
    velocities = arr[:, 2:4]

    return Snapshot(step=step, time=time, positions=positions, velocities=velocities)


def load_all_snapshots(run_dir: Path) -> list[Snapshot]:
    """Find all snapshot_*.csv in run_dir, sort by step, load and return."""
    pattern = "snapshot_*.csv"
    files = sorted(run_dir.glob(pattern), key=lambda p: p.stem)
    snapshots = []
    for p in files:
        snap = load_snapshot_csv(p)
        if snap is not None:
            snapshots.append(snap)
    return snapshots


def galaxy_zoom_target_limit(positions: np.ndarray, fallback: float = 1e20) -> float:
    """
    Half-axis range from the active disk: 90th-percentile outer radius, then mean r
    of stars inside that shell, times 1.10 cushion; floor 1e18. (Static plots: no smoothing.)
    """
    if positions.size == 0 or len(positions) == 0:
        return fallback
    df = pd.DataFrame({"x": positions[:, 0], "y": positions[:, 1]})
    r = np.sqrt(df["x"] ** 2 + df["y"] ** 2)
    if len(r) > 0:
        active_radius_limit = np.percentile(r, 90)
        inner_stars = r[r <= active_radius_limit]
        trimmed_mean = float(np.mean(inner_stars)) if len(inner_stars) > 0 else 1e20
        target_limit = trimmed_mean * 1.10
        target_limit = max(float(target_limit), 1e18)
    else:
        target_limit = 1e20
    if not np.isfinite(target_limit):
        return fallback
    return float(target_limit)


def galaxy_zoom_smoothed_limit(positions: np.ndarray, fallback: float = 1e20) -> float:
    """
    Same target as galaxy_zoom_target_limit, then exponential smoothing on plt.current_zoom_limit
    (alpha=0.1) for fluid camera motion between animation frames.
    """
    if not hasattr(plt, "current_zoom_limit"):
        plt.current_zoom_limit = 1e20
    target_limit = galaxy_zoom_target_limit(positions, fallback=fallback)
    alpha = 0.1
    plt.current_zoom_limit = (1.0 - alpha) * plt.current_zoom_limit + alpha * target_limit
    return float(plt.current_zoom_limit)


def get_masses_from_snapshots(snapshots: list[Snapshot], run_dir: Path) -> np.ndarray:
    """Read masses from the first snapshot CSV (mass column)."""
    if not snapshots:
        return np.array([])
    path = run_dir / f"snapshot_{snapshots[0].step:05d}.csv"
    if not path.exists():
        return np.array([])
    lines = path.read_text().strip().splitlines()
    masses = []
    for ln in lines[2:]:
        ln = ln.strip()
        if not ln:
            continue
        parts = ln.split(",")
        if len(parts) >= 6:
            masses.append(float(parts[5]))
    return np.array(masses) if masses else np.array([])


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot C++ galaxy simulation outputs (initial/final, animation, diagnostics)."
    )
    parser.add_argument(
        "run_dir",
        type=Path,
        help="Path to the C++ run directory (e.g. cpp_sim/outputs/YYYYMMDD_HHMMSS)",
    )
    parser.add_argument(
        "--no-animation",
        action="store_true",
        help="Skip MP4/GIF animation generation",
    )
    parser.add_argument(
        "--no-diagnostics",
        action="store_true",
        help="Skip diagnostic time-series plots",
    )
    parser.add_argument(
        "--render-radius",
        type=float,
        default=150.0,
        help="Fallback half-axis range if snapshot extents are degenerate (default 150)",
    )
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"Not a directory: {run_dir}")

    # Load run info and snapshots
    run_info = load_run_info(run_dir)
    snapshots = load_all_snapshots(run_dir)
    if not snapshots:
        raise SystemExit(f"No snapshots found in {run_dir} (look for snapshot_*.csv)")

    n_stars = len(snapshots[0].positions)
    print(f"Run dir: {run_dir}")
    print(f"Snapshots: {len(snapshots)}, particles: {n_stars}")
    if run_info:
        print(f"  n_steps: {run_info.get('n_steps', '?')}, dt: {run_info.get('dt', '?')}")

    fallback_radius = float(args.render_radius)
    initial = snapshots[0]
    final = snapshots[-1]
    render_radius_initial = galaxy_zoom_target_limit(initial.positions, fallback=fallback_radius)
    render_radius_final = galaxy_zoom_target_limit(final.positions, fallback=fallback_radius)

    # Initial and final scatter plots
    save_static_plot(
        initial.positions,
        run_dir / "galaxy_initial.png",
        title="Galaxy – Initial (C++ run)",
        render_radius=render_radius_initial,
    )
    print(f"Saved: {run_dir / 'galaxy_initial.png'}")
    save_static_plot(
        final.positions,
        run_dir / "galaxy_final.png",
        title="Galaxy – Final (C++ run)",
        render_radius=render_radius_final,
    )
    print(f"Saved: {run_dir / 'galaxy_final.png'}")

    # Optional animation
    if not args.no_animation:
        print("Creating animation...")
        plt.current_zoom_limit = 1e20
        if has_ffmpeg():
            print("  Using ffmpeg for MP4")
        else:
            print("  ffmpeg not found, using GIF if available")
        ok = create_animation(
            snapshots,
            run_dir / "galaxy",
            render_radius=lambda pos: galaxy_zoom_smoothed_limit(pos, fallback=fallback_radius),
            interval=50,
            progress_interval=max(1, len(snapshots) // 20),
        )
        if ok:
            if (run_dir / "galaxy.mp4").exists():
                print(f"Saved: {run_dir / 'galaxy.mp4'}")
            else:
                print(f"Saved: {run_dir / 'galaxy.gif'}")
        else:
            print("  Animation failed (install ffmpeg or Pillow)")
    else:
        print("Skipping animation (--no-animation).")

    # Optional diagnostics (same as Python pipeline)
    if not args.no_diagnostics and len(snapshots) > 1:
        masses = get_masses_from_snapshots(snapshots, run_dir)
        if len(masses) != n_stars:
            masses = np.ones(n_stars) * run_info.get("star_mass", 1.0)
        cutoff = 50.0  # C++ run_info does not store this; use default
        diag = compute_diagnostics(snapshots, masses, cutoff)
        plot_and_save_all(diag, run_dir, cutoff)
        print(f"Saved diagnostic plots in {run_dir}")
        print(f"  Final median_r: {diag['median_r'][-1]:.2f}, L_z: {diag['L_z'][-1]:.2f}")
    else:
        print("Skipping diagnostics (--no-diagnostics or single snapshot).")

    print(f"Done. Outputs in {run_dir}")


if __name__ == "__main__":
    main()
