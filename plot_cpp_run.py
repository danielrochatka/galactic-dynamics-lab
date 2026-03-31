#!/usr/bin/env python3
"""
Post-process C++ galaxy simulation outputs and generate visualizations.

Reads snapshot CSV files and run_info.txt from a cpp_sim run directory,
then produces the same kinds of plots as the Python pipeline:
  - initial and final scatter plots (galaxy view)
  - rotation_curve.png (final snapshot vs Keplerian √(GM/r); M_bh from run_info when present)
  - optional MP4/GIF animation
  - optional diagnostic time-series plots

Animation viewport: default uses **smart framing** — one fixed axis half-range from the **median**
  radial distance of all stars across all snapshots, times 1.2 (see `calculate_smart_bounds`).
  Set **`plot_animation_dynamic_zoom = true`** in the run config (or **`--dynamic-zoom`** on the
  command line) for the legacy **per-frame** velocity-gated smoothed zoom. **`--no-dynamic-zoom`**
  forces smart framing even if run_info requests dynamic zoom.
Burn-in filter (plotting only): `--skip-initial-steps` and/or `--skip-initial-snapshots` ignore
early snapshots for plotting products (animation/PNGs/diagnostics). You can also set
`plot_skip_initial_steps` / `plot_skip_initial_snapshots` in the run config (written to
`run_info.txt`). CLI flags override run_info; raw snapshot files are unchanged on disk.

Usage:
  python plot_cpp_run.py <run_dir> [--no-animation] [--no-diagnostics]
  python plot_cpp_run.py cpp_sim/outputs/20260308_175421

Outputs are written into the same run_dir (or run_dir/plots if you prefer;
  currently we write into run_dir to match "same run folder").

Requires: numpy, pandas, matplotlib. Optional: ffmpeg or Pillow for animation.

Galaxy runs: cpp_sim writes render_manifest.json and render_manifest.txt (with active_dynamics_branch,
active_metrics_branch, acceleration_code_path). plot_cpp_run draws optional text overlay on galaxy_*.png
and animation frames (run_info render_overlay_mode or --render-overlay-mode).
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
from render_overlay import build_overlay_spec, resolve_overlay_mode
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
    records = load_all_snapshot_records(run_dir)
    return [snap for _, snap in records]


def load_all_snapshot_records(run_dir: Path) -> list[tuple[Path, Snapshot]]:
    """Find all snapshot_*.csv in run_dir, sort by step, load and return (path, snapshot)."""
    pattern = "snapshot_*.csv"
    files = sorted(run_dir.glob(pattern), key=lambda p: p.stem)
    records: list[tuple[Path, Snapshot]] = []
    for p in files:
        snap = load_snapshot_csv(p)
        if snap is not None:
            records.append((p, snap))
    return records


def filter_snapshots_for_plotting(
    records: list[tuple[Path, Snapshot]],
    skip_initial_steps: int = 0,
    skip_initial_snapshots: int = 0,
) -> list[tuple[Path, Snapshot]]:
    """
    Burn-in filter for plotting products only.
    Applies both filters conservatively: (step >= skip_initial_steps) then drop first N remaining.
    """
    filtered = records
    if skip_initial_steps > 0:
        filtered = [(p, s) for (p, s) in filtered if s.step >= skip_initial_steps]
    if skip_initial_snapshots > 0:
        filtered = filtered[skip_initial_snapshots:]
    return filtered


def resolve_burnin_skip_settings(
    run_info: dict[str, str | int | float],
    cli_skip_initial_steps: int | None,
    cli_skip_initial_snapshots: int | None,
) -> tuple[int, int]:
    """Resolve burn-in skip settings with precedence: CLI overrides run_info; defaults to 0."""

    def ri_int(key: str, default: int) -> int:
        if not run_info or key not in run_info:
            return default
        try:
            return int(run_info[key])
        except Exception:
            return default

    effective_steps = (
        cli_skip_initial_steps
        if cli_skip_initial_steps is not None
        else ri_int("plot_skip_initial_steps", 0)
    )
    effective_snaps = (
        cli_skip_initial_snapshots
        if cli_skip_initial_snapshots is not None
        else ri_int("plot_skip_initial_snapshots", 0)
    )
    return effective_steps, effective_snaps


def resolve_diagnostic_cutoff_radius(
    run_info: dict[str, str | int | float],
    cli_cutoff_radius: float | None,
) -> tuple[float, str]:
    """
    Resolve diagnostic cutoff precedence:
      1) explicit CLI cutoff
      2) run_info diagnostic_cutoff_radius
      3) run_info galaxy_radius
      4) fail (no silent fallback).
    Returns (cutoff_radius, source_label).
    """
    if cli_cutoff_radius is not None:
        if not np.isfinite(cli_cutoff_radius) or cli_cutoff_radius <= 0:
            raise SystemExit("--diagnostic-cutoff-radius must be > 0")
        return float(cli_cutoff_radius), "cli(--diagnostic-cutoff-radius)"

    raw_cfg = run_info.get("diagnostic_cutoff_radius")
    if isinstance(raw_cfg, (int, float)) and np.isfinite(raw_cfg) and float(raw_cfg) > 0:
        return float(raw_cfg), "run_info(diagnostic_cutoff_radius)"

    raw_gr = run_info.get("galaxy_radius")
    if isinstance(raw_gr, (int, float)) and np.isfinite(raw_gr) and float(raw_gr) > 0:
        return float(raw_gr), "run_info(galaxy_radius)"

    raise SystemExit(
        "Diagnostics cutoff radius is undefined. Provide --diagnostic-cutoff-radius, or ensure "
        "run_info.txt contains diagnostic_cutoff_radius or galaxy_radius."
    )


def resolve_cooling_audit_flags(
    run_info: dict[str, str | int | float],
    snapshots: list[Snapshot],
) -> tuple[bool, int, int, float]:
    """
    Cooling audit metadata from run_info (with fallback to loaded snapshot list).
    Returns (cooling_active, cooling_steps, first_saved_step, first_saved_time).
    """

    def ri_int(key: str, default: int) -> int:
        raw = run_info.get(key)
        if isinstance(raw, (int, float)):
            return int(raw)
        try:
            return int(str(raw)) if raw is not None else default
        except Exception:
            return default

    def ri_float(key: str, default: float) -> float:
        raw = run_info.get(key)
        if isinstance(raw, (int, float)):
            return float(raw)
        try:
            return float(str(raw)) if raw is not None else default
        except Exception:
            return default

    cooling_active = ri_int("cooling_active", 0) != 0
    cooling_steps = ri_int("cooling_steps", 0)
    fs_step = ri_int("first_saved_snapshot_step", snapshots[0].step if snapshots else 0)
    fs_time = ri_float("first_saved_snapshot_time", snapshots[0].time if snapshots else 0.0)
    return cooling_active, cooling_steps, fs_step, fs_time


def initial_snapshot_plot_title(
    cooling_active: bool,
    initial_step: int,
    *,
    burn_in_plotting: bool = False,
) -> str:
    """Truthful first-frame title: cooling, burn-in plotting filter, or true initial."""
    if cooling_active and initial_step > 0:
        return "Galaxy – First saved snapshot after cooling (C++ run)"
    if burn_in_plotting:
        return "Galaxy – First plotted snapshot (C++ run)"
    return "Galaxy – Initial (C++ run)"


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


def fraction_stars_inside_square_viewport(
    positions: np.ndarray,
    half_axis: float,
) -> float:
    """Fraction of finite (x,y) inside the square |x|,|y| <= half_axis (matches matplotlib axis limits)."""
    if positions.size == 0 or len(positions) == 0:
        return 0.0
    x = positions[:, 0]
    y = positions[:, 1]
    m = np.isfinite(x) & np.isfinite(y)
    if not np.any(m):
        return 0.0
    ins = (np.abs(x[m]) <= half_axis) & (np.abs(y[m]) <= half_axis)
    return float(np.mean(ins))


def bulk_chebyshev_extent_half_axis(
    positions: np.ndarray,
    percentile: float = 95.0,
    cushion: float = 1.18,
) -> float:
    """
    Robust half-axis for square limits: percentile of max(|x|,|y|) among finite points, times cushion.
    Resists a handful of escapers vs raw max(r).
    """
    if positions.size == 0 or len(positions) == 0:
        return 1.0
    x = positions[:, 0]
    y = positions[:, 1]
    m = np.isfinite(x) & np.isfinite(y)
    if not np.any(m):
        return 1.0
    d = np.maximum(np.abs(x[m]), np.abs(y[m]))
    if d.size == 0:
        return 1.0
    q = float(np.percentile(d, percentile))
    if not np.isfinite(q) or q <= 0:
        return 1.0
    return float(q * cushion)


def static_viewport_radius_validated(
    positions: np.ndarray,
    velocities: np.ndarray | None,
    galaxy_radius_m: float,
    degenerate_fallback: float,
    min_inside_fraction: float = 0.50,
) -> float:
    """
    Static galaxy frame: start from velocity-gated candidate; if too few stars would appear in the
    square viewport, fall back to bulk (percentile) extent. Prevents blank plots when escapers or
    degenerate velocity masks make galaxy_velocity_gated_target_limit far smaller than the data extent.
    """
    cand = galaxy_velocity_gated_target_limit(
        positions, velocities, galaxy_radius_m, degenerate_fallback
    )
    if not np.isfinite(cand) or cand <= 0:
        cand = float(degenerate_fallback)
    frac = fraction_stars_inside_square_viewport(positions, cand)
    if frac >= min_inside_fraction:
        return float(cand)
    bulk = bulk_chebyshev_extent_half_axis(positions, 95.0, 1.18)
    bulk_wide = bulk_chebyshev_extent_half_axis(positions, 99.5, 1.12)
    out = max(cand, bulk, bulk_wide, float(degenerate_fallback))
    if not np.isfinite(out) or out <= 0:
        return float(degenerate_fallback)
    # Second check: if still almost empty, prefer widest bulk
    frac2 = fraction_stars_inside_square_viewport(positions, out)
    if frac2 < min_inside_fraction:
        out = max(out, bulk_chebyshev_extent_half_axis(positions, 99.9, 1.08))
    return float(out)


def calculate_smart_bounds(
    output_dir: Path,
    fallback: float = 150.0,
    snapshot_paths: list[Path] | None = None,
) -> float:
    """
    Global animation axis half-range: median of r = sqrt(x^2 + y^2) over every star in every
    snapshot_*.csv, times 1.20. Only the x and y columns are read from each CSV.
    """
    paths = snapshot_paths if snapshot_paths is not None else sorted(
        output_dir.glob("snapshot_*.csv"), key=lambda p: p.stem
    )
    if not paths:
        return float(max(fallback, 1.0))
    radii_chunks: list[np.ndarray] = []
    for path in paths:
        try:
            df = pd.read_csv(path, skiprows=1, usecols=["x", "y"])
        except (ValueError, KeyError):
            try:
                df = pd.read_csv(path, skiprows=2, header=None, usecols=[1, 2], names=["x", "y"])
            except Exception:
                continue
        x = df["x"].to_numpy(dtype=np.float64, copy=False)
        y = df["y"].to_numpy(dtype=np.float64, copy=False)
        radii_chunks.append(np.sqrt(x * x + y * y))
    if not radii_chunks:
        return float(max(fallback, 1.0))
    r_all = np.concatenate(radii_chunks)
    if r_all.size == 0:
        return float(max(fallback, 1.0))
    med = float(np.median(r_all))
    static_bound = med * 1.20
    if not np.isfinite(static_bound) or static_bound <= 0:
        return float(max(fallback, 1.0))
    return float(static_bound)


def plot_animation_dynamic_zoom_from_run_info(
    run_info: dict[str, str | int | float],
) -> bool:
    """
    Read plot_animation_dynamic_zoom from run_info (written by cpp_sim).
    Default False if missing: smart framing (fixed bounds from calculate_smart_bounds).
    """
    raw = run_info.get("plot_animation_dynamic_zoom")
    if raw is None:
        return False
    if isinstance(raw, (int, float)):
        return int(raw) != 0
    if isinstance(raw, str):
        s = raw.strip().lower()
        if s in ("0", "false", "no"):
            return False
        if s in ("1", "true", "yes"):
            return True
    return False


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
        help="Animation: smart framing (median r × 1.2 over all snapshots). "
        "Overrides run_info plot_animation_dynamic_zoom.",
    )
    zoom_group.add_argument(
        "--dynamic-zoom",
        action="store_true",
        help="Animation: legacy per-frame velocity-gated smoothed zoom (overrides run_info).",
    )
    parser.add_argument(
        "--render-overlay-mode",
        type=str,
        choices=("none", "minimal", "audit_full"),
        default=None,
        help="On-frame audit overlay for galaxy PNG/animation (default: from run_info render_overlay_mode; "
        "if missing, none for backward compatibility).",
    )
    parser.add_argument(
        "--skip-initial-steps",
        type=int,
        default=None,
        help="Burn-in filter (plotting only): ignore snapshots with step < this value. "
        "Overrides run_info plot_skip_initial_steps when set.",
    )
    parser.add_argument(
        "--skip-initial-snapshots",
        type=int,
        default=None,
        help="Burn-in filter (plotting only): after step filtering, drop first N snapshots. "
        "Overrides run_info plot_skip_initial_snapshots when set.",
    )
    parser.add_argument(
        "--diagnostic-cutoff-radius",
        type=float,
        default=None,
        help="Diagnostics-only radial cutoff (meters). Overrides run_info diagnostic_cutoff_radius "
        "and galaxy_radius when set.",
    )
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"Not a directory: {run_dir}")

    # Load run info and snapshots
    run_info = load_run_info(run_dir)
    skip_initial_steps, skip_initial_snapshots = resolve_burnin_skip_settings(
        run_info,
        cli_skip_initial_steps=args.skip_initial_steps,
        cli_skip_initial_snapshots=args.skip_initial_snapshots,
    )
    if skip_initial_steps < 0:
        raise SystemExit("--skip-initial-steps must be >= 0 (effective)")
    if skip_initial_snapshots < 0:
        raise SystemExit("--skip-initial-snapshots must be >= 0 (effective)")
    snapshot_records = load_all_snapshot_records(run_dir)
    if not snapshot_records:
        raise SystemExit(f"No snapshots found in {run_dir} (look for snapshot_*.csv)")
    filtered_records = filter_snapshots_for_plotting(
        snapshot_records,
        skip_initial_steps=skip_initial_steps,
        skip_initial_snapshots=skip_initial_snapshots,
    )
    if not filtered_records:
        raise SystemExit(
            "Burn-in filter removed all snapshots. "
            f"Found={len(snapshot_records)}, skip_initial_steps={skip_initial_steps}, "
            f"skip_initial_snapshots={skip_initial_snapshots}."
        )
    snapshots = [s for _, s in filtered_records]
    filtered_paths = [p for p, _ in filtered_records]

    n_stars = len(snapshots[0].positions)
    print(f"Run dir: {run_dir}")
    print(f"Snapshots found (raw): {len(snapshot_records)}")
    print(
        "Burn-in filter: "
        f"skip_initial_steps={skip_initial_steps}, "
        f"skip_initial_snapshots={skip_initial_snapshots}"
    )
    print(f"Snapshots used for plotting: {len(snapshots)}, particles: {n_stars}")
    print(
        f"Plotted range: first step/time={snapshots[0].step}/{float(snapshots[0].time):g}, "
        f"last step/time={snapshots[-1].step}/{float(snapshots[-1].time):g}"
    )
    if run_info:
        print(f"  n_steps: {run_info.get('n_steps', '?')}, dt: {run_info.get('dt', '?')}")
    cooling_active, cooling_steps, first_saved_step, first_saved_time = resolve_cooling_audit_flags(
        run_info, snapshots
    )
    if cooling_active:
        print(
            "WARNING: cooling phase was active; early intermediate snapshots were suppressed "
            "during cooling for audit/performance."
        )
        print(
            "  cooling_audit: "
            f"cooling_steps={cooling_steps}, first_saved_snapshot_step={first_saved_step}, "
            f"first_saved_snapshot_time={first_saved_time:g}"
        )

    fallback_radius = float(args.render_radius)
    galaxy_radius_m = resolve_galaxy_radius_meters(run_info, snapshots, fallback_radius)
    if args.dynamic_zoom:
        use_animation_dynamic_zoom = True
    elif args.no_dynamic_zoom:
        use_animation_dynamic_zoom = False
    else:
        use_animation_dynamic_zoom = plot_animation_dynamic_zoom_from_run_info(run_info)

    overlay_mode = resolve_overlay_mode(run_info, args.render_overlay_mode)
    overlay_spec = build_overlay_spec(run_dir, run_info)
    overlay_kw = {}
    if overlay_mode != "none":
        overlay_kw = {
            "overlay_mode": overlay_mode,
            "overlay_spec": overlay_spec,
            "run_info": run_info,
        }

    initial = snapshots[0]
    final = snapshots[-1]
    render_radius_initial = static_viewport_radius_validated(
        initial.positions, initial.velocities, galaxy_radius_m, fallback_radius
    )
    render_radius_final = static_viewport_radius_validated(
        final.positions, final.velocities, galaxy_radius_m, fallback_radius
    )

    burn_in_plotting = skip_initial_steps > 0 or skip_initial_snapshots > 0
    # Initial and final scatter plots
    initial_title = initial_snapshot_plot_title(
        cooling_active, initial.step, burn_in_plotting=burn_in_plotting
    )
    save_static_plot(
        initial.positions,
        run_dir / "galaxy_initial.png",
        title=initial_title,
        render_radius=render_radius_initial,
        velocities=initial.velocities,
        overlay_step=initial.step,
        overlay_time=float(initial.time),
        **overlay_kw,
    )
    print(f"Saved: {run_dir / 'galaxy_initial.png'}")
    save_static_plot(
        final.positions,
        run_dir / "galaxy_final.png",
        title="Galaxy – Final (C++ run)",
        render_radius=render_radius_final,
        velocities=final.velocities,
        overlay_step=final.step,
        overlay_time=float(final.time),
        **overlay_kw,
    )
    print(f"Saved: {run_dir / 'galaxy_final.png'}")

    try:
        r_rc, v_rc = rv_from_plane_arrays(final.positions, final.velocities)
        title_rc = f"snapshot_{final.step:05d}.csv (step {final.step}, t={final.time:g})"
        M_bh_rc = newtonian_m_bh_from_run_info(run_info)
        x_max_rc = rotation_curve_x_max(run_info, r_rc)
        scatter_lbl = scatter_label_from_run_info(run_info)
        rc_prov = run_info.get("code_version_label") if run_info else None
        rc_prov_s = (
            str(rc_prov).strip()
            if isinstance(rc_prov, str) and rc_prov.strip()
            else None
        )
        save_rotation_curve_png(
            r_rc,
            v_rc,
            run_dir / "rotation_curve.png",
            title_rc,
            M_bh=M_bh_rc,
            x_max=x_max_rc,
            scatter_label=scatter_lbl,
            provenance_label=rc_prov_s,
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
            static_bound = calculate_smart_bounds(
                run_dir,
                fallback=fallback_radius,
                snapshot_paths=filtered_paths,
            )
            render_radius_cb = static_bound

        if use_animation_dynamic_zoom:
            print("  Animation viewport: dynamic zoom (velocity-gated, smoothed per frame)")
        else:
            print(
                f"  Animation viewport: smart framing (fixed half-axis = {float(render_radius_cb):.6g} m, "
                "median r × 1.2 over all snapshots)"
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
            **overlay_kw,
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
        cutoff, cutoff_source = resolve_diagnostic_cutoff_radius(
            run_info, cli_cutoff_radius=args.diagnostic_cutoff_radius
        )
        print(f"Diagnostics cutoff radius: {cutoff:g} m (source={cutoff_source})")
        diag = compute_diagnostics(snapshots, masses, cutoff)
        plot_and_save_all(diag, run_dir, cutoff)
        print(f"Saved diagnostic plots in {run_dir}")
        print(f"  Final median_r: {diag['median_r'][-1]:.2f}, L_z: {diag['L_z'][-1]:.2f}")
    else:
        print("Skipping diagnostics (--no-diagnostics or single snapshot).")

    print(f"Done. Outputs in {run_dir}")


if __name__ == "__main__":
    main()
