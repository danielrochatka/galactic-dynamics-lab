#!/usr/bin/env python3
"""
Rotation curve from N-body snapshot CSV: compare simulated speeds to a Newtonian Keplerian baseline.

Primary glob: output_*.csv (highest step inferred from filename). Falls back to snapshot_*.csv
(C++ galaxy_sim format) if no output_*.csv is found.

`save_rotation_curve_png` is also used by plot_cpp_run.py so rotation_curve.png is written with other outputs.
The Newtonian overlay uses fixed NEWTONIAN_REFERENCE_M_BH (1e41 kg), not run_info bh_mass.

Usage:
  python plot_rotation_curve.py
  python plot_rotation_curve.py cpp_sim/outputs/20260327_120000
  python plot_rotation_curve.py --search-dir . --output rotation_curve.png
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

G_SI = 6.6743e-11
# Fixed Newtonian benchmark for rotation-curve overlay (do not use simulation bh_mass / effective mass).
NEWTONIAN_REFERENCE_M_BH = 1.0e41
DEFAULT_X_MAX = 1.2e20


def _step_from_stem(stem: str) -> int:
    m = re.search(r"(?:snapshot_|output_)(\d+)", stem, re.IGNORECASE)
    if m:
        return int(m.group(1))
    nums = re.findall(r"\d+", stem)
    return int(max(nums, key=len)) if nums else -1


def find_latest_csv(search_dir: Path) -> Path:
    """
    Prefer output_*.csv with largest step number; if none exist, use snapshot_*.csv (e.g. C++ runs).
    """
    search_dir = search_dir.resolve()
    output_files = [Path(p) for p in glob.glob(str(search_dir / "output_*.csv"))]
    if output_files:
        return max(output_files, key=lambda p: _step_from_stem(p.stem))
    snapshot_files = [Path(p) for p in glob.glob(str(search_dir / "snapshot_*.csv"))]
    if snapshot_files:
        return max(snapshot_files, key=lambda p: _step_from_stem(p.stem))
    raise FileNotFoundError(
        f"No output_*.csv or snapshot_*.csv under {search_dir}. "
        "Pass --search-dir to your run folder (e.g. cpp_sim/outputs/...)."
    )


def load_particle_kinematics(path: Path) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Load x,y,(z) and vx,vy,(vz). C++ snapshots: comment line, then i,x,y,vx,vy,mass.
    Returns r_3d, v_3d, title_suffix (filename + step if parsed).
    """
    lines = path.read_text().strip().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Not enough lines in {path}")

    step_note = ""
    m = re.search(r"step\s*,\s*(\d+)", lines[0], re.IGNORECASE)
    if m:
        step_note = f" (step {m.group(1)})"

    header = [c.strip().lower() for c in lines[1].split(",")]
    rows: list[list[float]] = []
    for ln in lines[2:]:
        ln = ln.strip()
        if not ln:
            continue
        parts = [p.strip() for p in ln.split(",")]
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        rows.append(vals)

    if not rows:
        raise ValueError(f"No numeric data rows in {path}")

    arr = np.array(rows)

    def idx(name: str) -> int | None:
        try:
            return header.index(name)
        except ValueError:
            return None

    ix, iy = idx("x"), idx("y")
    iz, ivx, ivy, ivz = idx("z"), idx("vx"), idx("vy"), idx("vz")

    if ix is None or iy is None or ivx is None or ivy is None:
        if arr.shape[1] >= 5:
            x = arr[:, 1]
            y = arr[:, 2]
            z = np.zeros_like(x)
            vx = arr[:, 3]
            vy = arr[:, 4]
            vz = np.zeros_like(vx)
        else:
            raise ValueError(f"Unrecognized columns in {path}: {header}")
    else:
        x = arr[:, ix]
        y = arr[:, iy]
        z = arr[:, iz] if iz is not None and iz < arr.shape[1] else np.zeros_like(x)
        vx = arr[:, ivx]
        vy = arr[:, ivy]
        vz = arr[:, ivz] if ivz is not None and ivz < arr.shape[1] else np.zeros_like(vx)

    r = np.sqrt(x * x + y * y + z * z)
    v = np.sqrt(vx * vx + vy * vy + vz * vz)
    title_suffix = f"{path.name}{step_note}"
    return r, v, title_suffix


def rv_from_plane_arrays(positions: np.ndarray, velocities: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    3D r and speed from 2D (n,2) or 3D (n,3) positions/velocities (C++ snapshot arrays).
    """
    x = positions[:, 0]
    y = positions[:, 1]
    z = positions[:, 2] if positions.shape[1] >= 3 else np.zeros_like(x)
    vx = velocities[:, 0]
    vy = velocities[:, 1]
    vz = velocities[:, 2] if velocities.shape[1] >= 3 else np.zeros_like(vx)
    r = np.sqrt(x * x + y * y + z * z)
    v = np.sqrt(vx * vx + vy * vy + vz * vz)
    return r, v


def save_rotation_curve_png(
    r: np.ndarray,
    v: np.ndarray,
    output_path: Path,
    title_suffix: str,
    *,
    M_bh: float = NEWTONIAN_REFERENCE_M_BH,
    x_max: float = DEFAULT_X_MAX,
) -> None:
    """
    Scatter r vs v and overlay Newtonian v = sqrt(G*M_bh/r). Arrays need not be pre-sorted.
    """
    order = np.argsort(r)
    r_s = r[order]
    v_s = v[order]

    r_line = np.linspace(max(float(r_s.min()), 1.0), x_max, 2000)
    v_newton = np.sqrt(np.maximum(0.0, (G_SI * M_bh) / r_line))

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(r_s, v_s, s=5, alpha=0.6, c="blue", edgecolors="none", label="TPF Simulated Stars")
    ax.plot(
        r_line,
        v_newton,
        "r--",
        linewidth=2,
        label="Newtonian Prediction (1/sqrt(r))",
    )
    ax.set_xlim(0.0, x_max)
    ax.set_xlabel("Distance from Galactic Center (m)")
    ax.set_ylabel("Orbital Velocity (m/s)")
    ax.set_title(f"Rotation curve — {title_suffix}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot rotation curve from latest snapshot/output CSV.")
    parser.add_argument(
        "search_dir",
        type=Path,
        nargs="?",
        default=Path("."),
        help="Directory to glob for output_*.csv / snapshot_*.csv (default: current directory)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("rotation_curve.png"),
        help="Output image path (default: rotation_curve.png)",
    )
    parser.add_argument(
        "--M-bh",
        type=float,
        default=NEWTONIAN_REFERENCE_M_BH,
        help="Black hole mass for Newtonian baseline (kg); default fixed benchmark 1e41",
    )
    args = parser.parse_args()

    csv_path = find_latest_csv(args.search_dir)
    r, v, title_suffix = load_particle_kinematics(csv_path)
    save_rotation_curve_png(r, v, args.output, title_suffix, M_bh=args.M_bh)
    print(f"Loaded: {csv_path}")
    print(f"Saved: {args.output.resolve()}")


if __name__ == "__main__":
    main()
