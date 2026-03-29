#!/usr/bin/env python3
"""
Post-process C++ galaxy simulation outputs and generate visualizations.

Reads snapshot CSV files and run_info.txt from a cpp_sim run directory,
then produces the same kinds of plots as the Python pipeline:
  - initial and final scatter plots (galaxy view)
  - rotation_curve.png (final snapshot vs Keplerian √(GM/r); M_bh from run_info when present)
  - optional MP4/GIF animation
  - optional diagnostic time-series plots

Animation viewport: default uses temporally smoothed limits (dynamic zoom). When disabled
  (run_info plot_animation_dynamic_zoom=0 or --no-dynamic-zoom), the axis half-range is held
  constant at the max velocity-gated limit over all snapshots (no per-frame zoom pumping).

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
from plot_rotation_curve import (
    load_run_info,
    newtonian_m_bh_from_run_info,
    rotation_curve_x_max,
    rv_from_plane_arrays,
    save_rotation_curve_png,
    scatter_label_from_run_info,
)


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


def resolve_galaxy_radius_meters(
    run_info: dict[str, str | int | float],
    snapshots: list[Snapshot],
    render_radius_arg: float,
) -> float:
    """Prefer run_info galaxy_radius; else max r on the initial snapshot; else CLI fallback."""
    raw = run_info.get("galaxy_radius")
    if isinstance(raw, (int, float)) and np.isfinite(raw) and float(raw) > 0:
        return float(raw)
    if snapshots and len(snapshots[0].positions) > 0:
        p0 = snapshots[0].positions
        r0 = np.sqrt(p0[:, 0] ** 2 + p0[:, 1] ** 2)
        if len(r0) > 0:
            return float(np.max(r0))
    return max(float(render_radius_arg), 1.0)


def galaxy_velocity_gated_target_limit(
    positions: np.ndarray,
    velocities: np.ndarray,
    galaxy_radius_m: float,
    degenerate_fallback: float,
) -> float:
    """
    Velocity-gated viewport: mean r of 'stable' stars (mostly tangential motion) + 30% cushion;
    if too few stable stars, use 1.2 * galaxy_radius_m. Static plots use this directly.
    """
    if positions.size == 0 or len(positions) == 0:
        return degenerate_fallback
    if velocities is None or len(velocities) != len(positions):
        return float(max(1.2 * galaxy_radius_m, degenerate_fallback))
    df = pd.DataFrame(
        {
            "x": positions[:, 0],
            "y": positions[:, 1],
            "vx": velocities[:, 0],
            "vy": velocities[:, 1],
        }
    )
    r = np.sqrt(df["x"] ** 2 + df["y"] ** 2).to_numpy()
    if len(r) == 0:
        return degenerate_fallback
    r_safe = np.maximum(r, 1e-30)
    v_rad = (df["x"] * df["vx"] + df["y"] * df["vy"]).to_numpy() / r_safe
    v_total = np.sqrt(df["vx"] ** 2 + df["vy"] ** 2).to_numpy()
    stable_mask = np.abs(v_rad) < (0.3 * v_total)
    stable_stars = r[stable_mask]
    n = len(r)
    if len(stable_stars) > (0.1 * n):
        target_limit = float(np.mean(stable_stars) * 1.30)
    else:
        target_limit = float(1.2 * galaxy_radius_m)
    if not np.isfinite(target_limit) or target_limit <= 0:
        return degenerate_fallback
    return float(target_limit)


def plot_animation_dynamic_zoom_from_run_info(
    run_info: dict[str, str | int | float],
) -> bool:
    """
    Read plot_animation_dynamic_zoom from run_info (written by cpp_sim).
    Default True if missing (preserve legacy animation smoothing).
    """
    raw = run_info.get("plot_animation_dynamic_zoom")
    if raw is None:
        return True
    if isinstance(raw, (int, float)):
        return int(raw) != 0
    if isinstance(raw, str):
        s = raw.strip().lower()
        if s in ("0", "false", "no"):
            return False
        if s in ("1", "true", "yes"):
            return True
    return True


def galaxy_velocity_gated_smoothed_limit(
    positions: np.ndarray,
    velocities: np.ndarray,
    galaxy_radius_m: float,
    degenerate_fallback: float,
) -> float:
    """Same target as galaxy_velocity_gated_target_limit, smoothed on plt.current_zoom_limit (alpha=0.1)."""
    if not hasattr(plt, "current_zoom_limit"):
        plt.current_zoom_limit = 1e20
    target_limit = galaxy_velocity_gated_target_limit(
        positions, velocities, galaxy_radius_m, degenerate_fallback
    )
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
    zoom_group = parser.add_mutually_exclusive_group()
    zoom_group.add_argument(
        "--no-dynamic-zoom",
        action="store_true",
        help="Animation: constant axis half-range = max velocity-gated limit over all snapshots "
        "(no temporal smoothing, no per-frame zoom). Overrides run_info plot_animation_dynamic_zoom.",
    )
    zoom_group.add_argument(
        "--dynamic-zoom",
        action="store_true",
        help="Animation: force smoothed viewport (overrides run_info).",
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
    galaxy_radius_m = resolve_galaxy_radius_meters(run_info, snapshots, fallback_radius)
    if args.dynamic_zoom:
        use_animation_dynamic_zoom = True
    elif args.no_dynamic_zoom:
        use_animation_dynamic_zoom = False
    else:
        use_animation_dynamic_zoom = plot_animation_dynamic_zoom_from_run_info(run_info)

    # region agent log
    _pdz_raw = run_info.get("plot_animation_dynamic_zoom")
    _rip = run_dir / "run_info.txt"
    _zoom_line_count = 0
    if _rip.exists():
        _zoom_line_count = sum(
            1
            for ln in _rip.read_text().splitlines()
            if ln.strip().startswith("plot_animation_dynamic_zoom")
        )
    _agent_log(
        "H2",
        "plot_cpp_run.py:zoom_resolve",
        "resolved use_animation_dynamic_zoom and run_info",
        {
            "raw_key_value": repr(_pdz_raw),
            "raw_key_type": type(_pdz_raw).__name__,
            "use_animation_dynamic_zoom": use_animation_dynamic_zoom,
            "cli_dynamic_zoom": bool(args.dynamic_zoom),
            "cli_no_dynamic_zoom": bool(args.no_dynamic_zoom),
            "run_info_plot_zoom_line_count": _zoom_line_count,
        },
    )
    # endregion

    initial = snapshots[0]
    final = snapshots[-1]
    render_radius_initial = galaxy_velocity_gated_target_limit(
        initial.positions, initial.velocities, galaxy_radius_m, fallback_radius
    )
    render_radius_final = galaxy_velocity_gated_target_limit(
        final.positions, final.velocities, galaxy_radius_m, fallback_radius
    )

    # Initial and final scatter plots
    save_static_plot(
        initial.positions,
        run_dir / "galaxy_initial.png",
        title="Galaxy – Initial (C++ run)",
        render_radius=render_radius_initial,
        velocities=initial.velocities,
    )
    print(f"Saved: {run_dir / 'galaxy_initial.png'}")
    save_static_plot(
        final.positions,
        run_dir / "galaxy_final.png",
        title="Galaxy – Final (C++ run)",
        render_radius=render_radius_final,
        velocities=final.velocities,
    )
    print(f"Saved: {run_dir / 'galaxy_final.png'}")

    try:
        r_rc, v_rc = rv_from_plane_arrays(final.positions, final.velocities)
        title_rc = f"snapshot_{final.step:05d}.csv (step {final.step}, t={final.time:g})"
        M_bh_rc = newtonian_m_bh_from_run_info(run_info)
        x_max_rc = rotation_curve_x_max(run_info, r_rc)
        scatter_lbl = scatter_label_from_run_info(run_info)
        save_rotation_curve_png(
            r_rc,
            v_rc,
            run_dir / "rotation_curve.png",
            title_rc,
            M_bh=M_bh_rc,
            x_max=x_max_rc,
            scatter_label=scatter_lbl,
        )
        print(
            f"Saved: {run_dir / 'rotation_curve.png'} "
            f"(Keplerian overlay M_bh={M_bh_rc:g} kg; x_max={x_max_rc:g} m; "
            "red curve is κ-independent — only blue scatter changes with TPF)"
        )
    except Exception as exc:
        print(f"Warning: could not write rotation_curve.png: {exc}")

    # Optional animation
    if not args.no_animation:
        print("Creating animation...")
        if use_animation_dynamic_zoom:
            plt.current_zoom_limit = 1e20
            render_radius_cb = lambda pos, vel: galaxy_velocity_gated_smoothed_limit(
                pos, vel, galaxy_radius_m, fallback_radius
            )
        else:
            # Constant axis limits for the whole clip: per-frame velocity-gated limits still change
            # with each snapshot (H1); max over frames avoids false "zoom" when dynamic zoom is off.
            _fixed_anim_r = max(
                galaxy_velocity_gated_target_limit(
                    s.positions, s.velocities, galaxy_radius_m, fallback_radius
                )
                for s in snapshots
            )
            render_radius_cb = lambda pos, vel, rfix=_fixed_anim_r: rfix

        # region agent log
        _anim_cb_calls = [0]

        def _wrap_anim_render_radius(inner):
            def _wrapped(pos, vel):
                lim = inner(pos, vel)
                i = _anim_cb_calls[0]
                ntot = len(snapshots)
                if i in (0, ntot // 2, max(0, ntot - 1)):
                    _agent_log(
                        "H1",
                        "plot_cpp_run.py:anim_render_cb",
                        "per-frame axis half-range from callback",
                        {
                            "callback_invocation": i,
                            "n_snapshots": ntot,
                            "render_radius_half_axis": float(lim),
                            "use_smoothed_dynamic_zoom": use_animation_dynamic_zoom,
                        },
                    )
                _anim_cb_calls[0] += 1
                return lim

            return _wrapped

        render_radius_cb = _wrap_anim_render_radius(render_radius_cb)
        # endregion

        print(
            f"  Animation viewport: {'dynamic zoom (smoothed)' if use_animation_dynamic_zoom else 'constant (max half-axis over snapshots; no per-frame zoom)'}"
        )
        if has_ffmpeg():
            print("  Using ffmpeg for MP4")
        else:
            print("  ffmpeg not found, using GIF if available")
        ok = create_animation(
            snapshots,
            run_dir / "galaxy",
            render_radius=render_radius_cb,
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
