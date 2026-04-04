#!/usr/bin/env python3
"""
Post-process C++ galaxy simulation outputs and generate visualizations.

Reads snapshot CSV files and run_info.txt from a cpp_sim run directory,
then produces the same kinds of plots as the Python pipeline:
  - initial and final scatter plots (galaxy view)
  - rotation_curve.png (final snapshot vs Keplerian √(GM/r); M_bh from run_info when present)
  - optional MP4/GIF animation
  - optional diagnostic time-series plots
  - **bh_orbit_validation**: primary pair diagnostics plus `bh_orbit_trajectory_xy.png`, `bh_orbit_trajectory_xy_zoom.png`,
    `bh_orbit_separation_extrema.png` (footnotes include package and coupling labels from run_info)

Animation viewport: default **smart framing** — one square viewport from **x/y quantile bounds**
  over all plotted particles (all snapshots), plus optional origin marker (BH-present runs only), margin, and a square
  half-axis (see `framing.global_viewport_from_snapshots`). Same rule for **galaxy** and non-galaxy
  modes; no median-radius proxy.
  Set **`plot_animation_dynamic_zoom = true`** (or **`--dynamic-zoom`**) for **windowed** geometric
  framing per frame with exponential smoothing on center and log(half-axis). **`--no-dynamic-zoom`**
  forces fixed global smart framing even if run_info requests dynamic zoom.
Burn-in filter (plotting only): `--skip-initial-steps` and/or `--skip-initial-snapshots` ignore
early snapshots for plotting products (animation/PNGs/diagnostics). You can also set
`plot_skip_initial_steps` / `plot_skip_initial_snapshots` in the run config (written to
`run_info.txt`). CLI flags override run_info; raw snapshot files are unchanged on disk.

Usage:
  python plot_cpp_run.py <run_dir> [--no-animation] [--no-diagnostics]
  python plot_cpp_run.py cpp_sim/outputs/20260308_175421

Outputs are written into the same run_dir (or run_dir/plots if you prefer;
  currently we write into run_dir to match "same run folder").

Requires: numpy, matplotlib. Optional: ffmpeg or Pillow for animation.

Galaxy runs: cpp_sim writes render_manifest.json and render_manifest.txt (with active_dynamics_branch,
active_metrics_branch, acceleration_code_path). plot_cpp_run draws optional text overlay on galaxy_*.png
and animation frames (run_info render_overlay_mode or --render-overlay-mode).
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from dataclasses import dataclass

import numpy as np

# Import from existing project (run from repo root or with PYTHONPATH)
from framing import (
    DynamicViewportSmoother,
    SquareViewport,
    global_viewport_from_snapshots,
    viewport_window_snapshots,
)
from render import save_static_plot, create_animation, has_ffmpeg
from render_overlay import build_overlay_spec, resolve_overlay_mode
from diagnostics import (
    bh_orbit_validation_plot_footnote,
    compute_diagnostics,
    compute_two_body_pair_diagnostics,
    plot_and_save_all,
    plot_bh_orbit_validation_extras,
    plot_two_body_pair_diagnostics,
    save_two_body_timeseries_csv,
    write_two_body_diagnostics_readme,
)
from display_units import (
    display_unit_config_from_run_info,
    series_display_generic_validation,
    series_display_for_galaxy_diagnostics,
    series_display_for_two_body,
    spatial_display_for_xy_plot,
)
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


# Matches cpp_sim/config.hpp enum class SimulationMode (declaration order).
# Index 1 is two_body_orbit in C++; postprocess uses the same diagnostics as earth_moon_benchmark.
_SIMULATION_MODE_INT_TO_NAME: dict[int, str] = {
    0: "galaxy",
    1: "earth_moon_benchmark",
    2: "symmetric_pair",
    3: "small_n_conservation",
    4: "timestep_convergence",
    5: "tpf_single_source_inspect",
    6: "tpf_symmetric_pair_inspect",
    7: "tpf_two_body_sweep",
    8: "tpf_weak_field_calibration",
    9: "tpf_newtonian_force_compare",
    10: "tpf_diagnostic_consistency_audit",
    11: "tpf_bound_orbit_sweep",
    12: "tpf_v11_weak_field_correspondence",
    13: "earth_moon_benchmark",
    14: "bh_orbit_validation",
}


def simulation_mode_name_from_run_info(run_info: dict[str, str | int | float]) -> str:
    """
    Resolved simulation mode name from run_info.txt.
    The file may list simulation_mode twice (string then int); the last line wins (typically int).
    """
    if not run_info:
        return "galaxy"
    raw = run_info.get("simulation_mode")
    if isinstance(raw, int):
        return _SIMULATION_MODE_INT_TO_NAME.get(raw, "galaxy")
    if isinstance(raw, str) and raw.strip():
        return raw.strip()
    return "galaxy"


def _slug(raw: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9]+", "_", str(raw).strip().lower()).strip("_")
    return s or "unknown"


def physics_label_from_run_info(
    run_info: dict[str, str | int | float],
    mode_name: str,
) -> str:
    if mode_name == "tpf_v11_weak_field_correspondence":
        return "paper_mode"
    pkg = str(run_info.get("physics_package", "Newtonian") or "Newtonian")
    if pkg == "Newtonian":
        return "newtonian"
    if pkg != "TPFCore":
        return _slug(pkg)
    dyn = str(run_info.get("tpf_dynamics_mode", "legacy_readout") or "legacy_readout")
    if dyn == "direct_tpf":
        return "tpfcore_direct_tpf"
    readout = str(run_info.get("tpfcore_readout_mode", "legacy_readout") or "legacy_readout")
    vdsg_raw = run_info.get("tpf_vdsg_coupling", 0.0)
    vdsg = float(vdsg_raw) if isinstance(vdsg_raw, (int, float)) else 0.0
    vdsg_lbl = "vdsg_off" if vdsg == 0.0 else "vdsg_on"
    return f"tpfcore_{_slug(readout)}_{vdsg_lbl}"


def mode_aware_name(
    mode: str,
    physics: str,
    scope: str,
    quantity: str,
    stage: str,
    ext: str,
) -> str:
    return f"{mode}__{physics}__{scope}__{quantity}__{stage}.{ext}"


def write_compat_alias(src: Path, dst: Path) -> None:
    if src == dst or not src.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


# Numeric step in snapshot_<step>.csv stem (not lex sort: e.g. snapshot_100000 must follow snapshot_99999).
_SNAPSHOT_STEM_NUMERIC_STEP = re.compile(r"^snapshot_(\d+)$", re.IGNORECASE)


def snapshot_csv_numeric_sort_key(path: Path) -> tuple[int, str]:
    m = _SNAPSHOT_STEM_NUMERIC_STEP.match(path.stem)
    if m:
        return (int(m.group(1)), path.name)
    return (10**18, path.name)


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
    """Find all snapshot_*.csv in run_dir, sort by numeric step (then name), load and return (path, snapshot)."""
    pattern = "snapshot_*.csv"
    files = sorted(run_dir.glob(pattern), key=snapshot_csv_numeric_sort_key)
    records: list[tuple[Path, Snapshot]] = []
    for p in files:
        snap = load_snapshot_csv(p)
        if snap is not None:
            records.append((p, snap))
    # Tie-break by CSV header (filename width can differ from step for edge cases).
    records.sort(key=lambda ps: (ps[1].step, ps[1].time, ps[0].name))
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
    """
    Prefer run_info galaxy_radius for **galaxy** mode only.

    For validation / test modes (earth_moon_benchmark, bh_orbit_validation, symmetric_pair, …), config galaxy_radius is often a
    leftover disk scale and must not override the true extent; we use max r on the first snapshot
    (then CLI fallback) instead.
    """
    mode = simulation_mode_name_from_run_info(run_info)
    raw = run_info.get("galaxy_radius")
    if (
        mode == "galaxy"
        and isinstance(raw, (int, float))
        and np.isfinite(raw)
        and float(raw) > 0
    ):
        return float(raw)
    if snapshots and len(snapshots[0].positions) > 0:
        p0 = snapshots[0].positions
        r0 = np.sqrt(p0[:, 0] ** 2 + p0[:, 1] ** 2)
        if len(r0) > 0:
            return float(np.max(r0))
    return max(float(render_radius_arg), 1.0)


# Smart framing defaults (display layer only; SI in/out)
_SMART_FRAMING_TRIM = 0.01
_SMART_FRAMING_MARGIN = 1.15
_DYNAMIC_FRAMING_WINDOW = 7
_DYNAMIC_FRAMING_ALPHA = 0.12
_BH_MARKER_XY = np.array([[0.0, 0.0]], dtype=np.float64)


def plot_animation_dynamic_zoom_from_run_info(
    run_info: dict[str, str | int | float],
) -> bool:
    """
    Read plot_animation_dynamic_zoom from run_info (written by cpp_sim).
    Default False if missing: fixed global geometric smart framing.
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


def _run_info_float(run_info: dict[str, str | int | float], key: str, default: float = 0.0) -> float:
    raw = run_info.get(key, default)
    if isinstance(raw, (int, float)):
        return float(raw)
    try:
        return float(str(raw))
    except Exception:
        return float(default)


def _run_info_bool(run_info: dict[str, str | int | float], key: str, default: bool = False) -> bool:
    raw = run_info.get(key)
    if raw is None:
        return default
    if isinstance(raw, bool):
        return raw
    if isinstance(raw, (int, float)):
        return int(raw) != 0
    if isinstance(raw, str):
        s = raw.strip().lower()
        if s in ("1", "true", "yes", "on"):
            return True
        if s in ("0", "false", "no", "off"):
            return False
    return default


def should_draw_central_bh_marker(
    run_info: dict[str, str | int | float],
    mode_name: str,
) -> bool:
    """
    Render BH/origin marker only when semantically part of the run interpretation.
    Prefer suppressing ambiguous marker for non-BH runs.
    """
    bh_mass = _run_info_float(run_info, "bh_mass", 0.0)
    if bh_mass <= 0.0:
        return False
    if mode_name in ("earth_moon_benchmark", "two_body_orbit"):
        return False
    if mode_name == "symmetric_pair":
        return _run_info_bool(run_info, "validation_symmetric_include_bh", default=True)
    return mode_name in ("galaxy", "bh_orbit_validation")


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
        help="Animation: fixed global geometric smart framing (quantile x/y bounds + margin). "
        "Overrides run_info plot_animation_dynamic_zoom.",
    )
    zoom_group.add_argument(
        "--dynamic-zoom",
        action="store_true",
        help="Animation: windowed geometric framing with smoothed center and log(half-axis) (overrides run_info).",
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
    mode_name = simulation_mode_name_from_run_info(run_info)
    mode_label = "tpf_v11_correspondence" if mode_name == "tpf_v11_weak_field_correspondence" else mode_name
    physics_label = physics_label_from_run_info(run_info, mode_name)
    title_context = f"{mode_label} / {physics_label}"
    show_bh_marker = should_draw_central_bh_marker(run_info, mode_name)

    global_viewport = global_viewport_from_snapshots(
        snapshots,
        extra_xy=_BH_MARKER_XY if show_bh_marker else None,
        trim_fraction=_SMART_FRAMING_TRIM,
        margin=_SMART_FRAMING_MARGIN,
        fallback_half_axis=fallback_radius,
    )
    render_radius_initial = global_viewport
    render_radius_final = global_viewport
    spatial_half_axis_m = float(global_viewport.half_axis)
    unit_cfg = display_unit_config_from_run_info(run_info)
    spatial_display = spatial_display_for_xy_plot(
        mode_name,
        spatial_half_axis_m,
        preferred_unit=unit_cfg.distance_unit,
    )
    max_time_s = float(max(s.time for s in snapshots)) if snapshots else 0.0
    max_speed_m_s = 0.0
    if snapshots:
        max_speed_m_s = float(
            max(np.nanmax(np.linalg.norm(s.velocities, axis=1)) for s in snapshots)
        )
    overlay_spec["active_display_distance_unit"] = spatial_display.unit
    if mode_name in ("earth_moon_benchmark", "bh_orbit_validation"):
        series_preview = series_display_for_two_body(
            mode_name,
            max_distance_m=max(spatial_half_axis_m, 1.0),
            max_time_s=max_time_s,
            max_speed_m_s=max_speed_m_s,
            preferred_distance_unit=unit_cfg.distance_unit,
            preferred_time_unit=unit_cfg.time_unit,
            preferred_velocity_unit=unit_cfg.velocity_unit,
        )
    elif mode_name == "galaxy":
        series_preview = series_display_for_galaxy_diagnostics(
            max(spatial_half_axis_m, 1.0),
            max_time_s,
            max_speed_m_s=max_speed_m_s,
            preferred_distance_unit=unit_cfg.distance_unit,
            preferred_time_unit=unit_cfg.time_unit,
            preferred_velocity_unit=unit_cfg.velocity_unit,
        )
    else:
        series_preview = series_display_generic_validation(
            max(spatial_half_axis_m, 1.0),
            max_time_s=max_time_s,
            max_speed_m_s=max_speed_m_s,
            preferred_distance_unit=unit_cfg.distance_unit,
            preferred_time_unit=unit_cfg.time_unit,
            preferred_velocity_unit=unit_cfg.velocity_unit,
        )
    overlay_spec["active_display_time_unit"] = (
        series_preview.time_unit
    )
    overlay_spec["active_display_velocity_unit"] = (
        series_preview.speed_unit
    )
    unit_reference_text = None
    if unit_cfg.show_unit_reference:
        unit_reference_text = (
            f"distance display = {overlay_spec['active_display_distance_unit']}\\n"
            f"time display = {overlay_spec['active_display_time_unit']}\\n"
            f"velocity display = {overlay_spec['active_display_velocity_unit']}"
        )
    print(
        f"  Smart framing (simulation_mode={mode_name}): center=({global_viewport.center_x:.6g}, "
        f"{global_viewport.center_y:.6g}) m, half_axis={spatial_half_axis_m:.6g} m "
        f"(geometric quantile trim={_SMART_FRAMING_TRIM:g}, margin={_SMART_FRAMING_MARGIN:g}, "
        f"origin_in_cloud={show_bh_marker})"
    )
    print(f"  Central BH marker: {'shown' if show_bh_marker else 'hidden'}")

    burn_in_plotting = skip_initial_steps > 0 or skip_initial_snapshots > 0
    # Initial and final scatter plots
    legacy_initial_title = initial_snapshot_plot_title(
        cooling_active, initial.step, burn_in_plotting=burn_in_plotting
    )
    initial_title = f"{title_context} — primary trajectory x-y (initial)"
    if "first plotted snapshot" in legacy_initial_title.lower():
        initial_title = f"{title_context} — primary trajectory x-y (first plotted snapshot)"
    final_title = f"{title_context} — primary trajectory x-y (final)"

    initial_legacy = run_dir / "galaxy_initial.png"
    initial_mode_aware = run_dir / mode_aware_name(
        mode_label, physics_label, "primary", "trajectory_xy", "initial", "png"
    )
    final_legacy = run_dir / "galaxy_final.png"
    final_mode_aware = run_dir / mode_aware_name(
        mode_label, physics_label, "primary", "trajectory_xy", "final", "png"
    )
    save_static_plot(
        initial.positions,
        initial_mode_aware,
        title=initial_title,
        render_radius=render_radius_initial,
        velocities=initial.velocities,
        overlay_step=initial.step,
        overlay_time=float(initial.time),
        show_bh_marker=show_bh_marker,
        spatial_display=spatial_display,
        unit_reference_text=unit_reference_text,
        **overlay_kw,
    )
    write_compat_alias(initial_mode_aware, initial_legacy)
    print(f"Saved: {initial_mode_aware} (legacy alias: {initial_legacy.name})")
    save_static_plot(
        final.positions,
        final_mode_aware,
        title=final_title,
        render_radius=render_radius_final,
        velocities=final.velocities,
        overlay_step=final.step,
        overlay_time=float(final.time),
        show_bh_marker=show_bh_marker,
        spatial_display=spatial_display,
        unit_reference_text=unit_reference_text,
        **overlay_kw,
    )
    write_compat_alias(final_mode_aware, final_legacy)
    print(f"Saved: {final_mode_aware} (legacy alias: {final_legacy.name})")

    if mode_name in ("earth_moon_benchmark", "bh_orbit_validation"):
        print(
            "Skipping rotation_curve.png for two-body mode (Keplerian disk overlay is not the primary "
            "interpretation tool; use diagnostic_pair_*.png and diagnostic_two_body_timeseries.csv)."
        )
    else:
        try:
            r_rc, v_rc = rv_from_plane_arrays(final.positions, final.velocities)
            title_rc = (
                f"{title_context} — rotation curve from snapshot_{final.step:05d}.csv "
                f"(step {final.step}, t={final.time:g})"
            )
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
                run_dir / mode_aware_name(
                    mode_label, physics_label, "primary", "rotation_curve", "final", "png"
                ),
                title_rc,
                M_bh=M_bh_rc,
                x_max=x_max_rc,
                scatter_label=scatter_lbl,
                provenance_label=rc_prov_s,
                preferred_distance_unit=unit_cfg.distance_unit,
                preferred_velocity_unit=unit_cfg.velocity_unit,
            )
            write_compat_alias(
                run_dir / mode_aware_name(mode_label, physics_label, "primary", "rotation_curve", "final", "png"),
                run_dir / "rotation_curve.png",
            )
            print(
                f"Saved: {run_dir / mode_aware_name(mode_label, physics_label, 'primary', 'rotation_curve', 'final', 'png')} "
                f"(Keplerian overlay M_bh={M_bh_rc:g} kg; x_max={x_max_rc:g} m; "
                "red curve is κ-independent; blue scatter varies with configured package parameters)"
            )
        except Exception as exc:
            print(f"Warning: could not write rotation_curve.png: {exc}")

    # Optional animation
    if not args.no_animation:
        print("Creating animation...")
        anim_mutable_idx: dict[str, int] = {"i": 0}
        if use_animation_dynamic_zoom:
            smoother = DynamicViewportSmoother(
                alpha=_DYNAMIC_FRAMING_ALPHA,
                min_half_axis=max(1.0, fallback_radius * 1e-18),
            )

            def render_radius_cb(
                pos: np.ndarray, vel: np.ndarray | None
            ) -> SquareViewport:
                _ = pos, vel
                win = viewport_window_snapshots(
                    snapshots, anim_mutable_idx["i"], _DYNAMIC_FRAMING_WINDOW
                )
                raw = global_viewport_from_snapshots(
                    win,
                    extra_xy=_BH_MARKER_XY if show_bh_marker else None,
                    trim_fraction=_SMART_FRAMING_TRIM,
                    margin=_SMART_FRAMING_MARGIN,
                    fallback_half_axis=fallback_radius,
                )
                return smoother.step(raw)

            render_radius_anim = render_radius_cb
            print(
                f"  Animation viewport: dynamic geometric (window={_DYNAMIC_FRAMING_WINDOW} frames, "
                f"EMA alpha={_DYNAMIC_FRAMING_ALPHA:g} on center and log half-axis)"
            )
        else:
            render_radius_anim = global_viewport
            print(
                f"  Animation viewport: fixed global smart framing (half_axis={global_viewport.half_axis:.6g} m)"
            )
        if has_ffmpeg():
            print("  Using ffmpeg for MP4")
        else:
            print("  ffmpeg not found, using GIF if available")
        ok = create_animation(
            snapshots,
            run_dir / "galaxy",
            render_radius=render_radius_anim,
            interval=50,
            progress_interval=max(1, len(snapshots) // 20),
            spatial_display=spatial_display,
            simulation_mode=mode_name,
            preferred_time_unit=unit_cfg.time_unit,
            active_time_unit=series_preview.time_unit,
            unit_reference_text=unit_reference_text,
            mutable_frame_index=anim_mutable_idx,
            show_bh_marker=show_bh_marker,
            **overlay_kw,
        )
        if ok:
            if (run_dir / "galaxy.mp4").exists():
                write_compat_alias(
                    run_dir / "galaxy.mp4",
                    run_dir / mode_aware_name(mode_label, physics_label, "primary", "trajectory_xy", "animation", "mp4"),
                )
                print(f"Saved: {run_dir / 'galaxy.mp4'}")
            else:
                write_compat_alias(
                    run_dir / "galaxy.gif",
                    run_dir / mode_aware_name(mode_label, physics_label, "primary", "trajectory_xy", "animation", "gif"),
                )
                print(f"Saved: {run_dir / 'galaxy.gif'}")
        else:
            print("  Animation failed (install ffmpeg or Pillow)")
    else:
        print("Skipping animation (--no-animation).")

    # Optional diagnostics (same as Python pipeline)
    if not args.no_diagnostics and len(snapshots) > 1:
        masses = get_masses_from_snapshots(snapshots, run_dir)
        if len(masses) != n_stars:
            masses = np.ones(n_stars) * run_info.get("star_mass", 1.98847e30)
        cutoff, cutoff_source = resolve_diagnostic_cutoff_radius(
            run_info, cli_cutoff_radius=args.diagnostic_cutoff_radius
        )
        print(f"Diagnostics cutoff radius: {cutoff:g} m (source={cutoff_source})")
        physics_pkg = str(run_info.get("physics_package", "") or "Newtonian")
        if mode_name in ("earth_moon_benchmark", "bh_orbit_validation"):
            bh_m = run_info.get("bh_mass", 0.0)
            bh_m_f = float(bh_m) if isinstance(bh_m, (int, float)) else 0.0
            try:
                pair_diag = compute_two_body_pair_diagnostics(
                    snapshots, masses, mode_name, bh_m_f
                )
                save_two_body_timeseries_csv(pair_diag, run_dir)
                write_two_body_diagnostics_readme(run_dir, mode_name, physics_pkg)
                max_pair_d = max(
                    float(np.max(pair_diag["pair_separation"])),
                    float(np.max(pair_diag["center_of_mass_radius"])),
                )
                pair_series = series_display_for_two_body(
                    mode_name,
                    max_distance_m=max_pair_d,
                    max_time_s=float(np.max(pair_diag["time"])),
                    max_speed_m_s=float(np.max(pair_diag["pair_relative_speed"])),
                    preferred_distance_unit=unit_cfg.distance_unit,
                    preferred_time_unit=unit_cfg.time_unit,
                    preferred_velocity_unit=unit_cfg.velocity_unit,
                )
                if mode_name == "bh_orbit_validation":
                    vdsg_raw = run_info.get("tpf_vdsg_coupling", 0.0)
                    vdsg_f = float(vdsg_raw) if isinstance(vdsg_raw, (int, float)) else 0.0
                    foot = bh_orbit_validation_plot_footnote(physics_pkg, vdsg_f)
                    title_suf = "(bh_orbit_validation · experimental)"
                    plot_two_body_pair_diagnostics(
                        pair_diag,
                        run_dir,
                        physics_pkg,
                        title_suffix=title_suf,
                        footnote=foot,
                        context_label=title_context,
                        series=pair_series,
                    )
                    plot_bh_orbit_validation_extras(
                        snapshots,
                        pair_diag,
                        run_dir,
                        physics_pkg,
                        vdsg_f,
                        context_label=title_context,
                        unit_config=unit_cfg,
                    )
                else:
                    plot_two_body_pair_diagnostics(
                        pair_diag,
                        run_dir,
                        physics_pkg,
                        context_label=title_context,
                        series=pair_series,
                    )
                write_compat_alias(
                    run_dir / "diagnostic_two_body_timeseries.csv",
                    run_dir
                    / mode_aware_name(
                        mode_label, physics_label, "primary", "two_body_timeseries", "final", "csv"
                    ),
                )
                write_compat_alias(
                    run_dir / "two_body_diagnostics_README.txt",
                    run_dir
                    / mode_aware_name(
                        mode_label, physics_label, "primary", "two_body_summary", "final", "txt"
                    ),
                )
                primary_aliases = {
                    "diagnostic_pair_separation.png": "pair_separation",
                    "diagnostic_pair_relative_speed.png": "pair_relative_speed",
                    "diagnostic_com_radius.png": "center_of_mass_radius",
                    "diagnostic_relative_angular_momentum_z.png": "relative_angular_momentum_z",
                    "diagnostic_relative_energy.png": "newtonian_reference_specific_energy",
                }
                for old, quantity in primary_aliases.items():
                    write_compat_alias(
                        run_dir / old,
                        run_dir / mode_aware_name(mode_label, physics_label, "primary", quantity, "timeseries", "png"),
                    )
                bh_aliases = {
                    "bh_orbit_trajectory_xy.png": "trajectory_xy",
                    "bh_orbit_trajectory_xy_zoom.png": "trajectory_xy_zoom",
                    "bh_orbit_separation_extrema.png": "pair_separation_extrema",
                }
                for old, quantity in bh_aliases.items():
                    write_compat_alias(
                        run_dir / old,
                        run_dir
                        / mode_aware_name(mode_label, physics_label, "primary", quantity, "final", "png"),
                    )
                print(
                    f"Primary two-body diagnostics in {run_dir}: "
                    "diagnostic_pair_separation.png, diagnostic_pair_relative_speed.png, "
                    "diagnostic_com_radius.png, diagnostic_relative_angular_momentum_z.png, "
                    "diagnostic_relative_energy.png (if finite), diagnostic_two_body_timeseries.csv, "
                    "two_body_diagnostics_README.txt"
                )
                if mode_name == "bh_orbit_validation":
                    print(
                        "  bh_orbit_validation extras: bh_orbit_trajectory_xy.png, "
                        "bh_orbit_trajectory_xy_zoom.png, bh_orbit_separation_extrema.png "
                        "(footnotes: package and coupling labels from run_info; "
                        "configuration labels only)"
                    )
                print(
                    f"  Final pair separation: {float(pair_diag['pair_separation'][-1]):.6g} m; "
                    f"pair relative speed: {float(pair_diag['pair_relative_speed'][-1]):.6g} m/s; "
                    f"|R_COM| (lab): {float(pair_diag['center_of_mass_radius'][-1]):.6g} m"
                )
            except ValueError as exc:
                print(f"Warning: primary two-body diagnostics skipped: {exc}")
            diag = compute_diagnostics(snapshots, masses, cutoff)
            plot_and_save_all(
                diag,
                run_dir,
                cutoff,
                lab_frame_secondary=True,
                context_label=title_context,
                simulation_mode=mode_name,
                two_body_secondary=True,
                unit_config=unit_cfg,
            )
            secondary_aliases = {
                "diagnostic_median_radius.png": "median_radius_origin",
                "diagnostic_mean_radius.png": "mean_radius_origin",
                "diagnostic_std_radius.png": "std_radius_origin",
                "diagnostic_max_radius.png": "max_radius_origin",
                "diagnostic_frac_vr_positive.png": "frac_vr_positive_origin",
                "diagnostic_frac_vr_negative.png": "frac_vr_negative_origin",
                "diagnostic_frac_beyond_cutoff.png": "frac_beyond_cutoff_origin",
                "diagnostic_angular_momentum_z.png": "angular_momentum_z_origin",
            }
            for old, quantity in secondary_aliases.items():
                write_compat_alias(
                    run_dir / old,
                    run_dir / mode_aware_name(mode_label, physics_label, "secondary", quantity, "timeseries", "png"),
                )
            print(
                f"Secondary lab-frame / origin-radial plots in {run_dir} (titles prefixed); "
                f"final median_r (secondary): {diag['median_r'][-1]:.6g}, "
                f"L_z about origin (secondary): {diag['L_z'][-1]:.6g}"
            )
        else:
            diag = compute_diagnostics(snapshots, masses, cutoff)
            plot_and_save_all(
                diag,
                run_dir,
                cutoff,
                context_label=title_context,
                simulation_mode=mode_name,
                two_body_secondary=False,
                unit_config=unit_cfg,
            )
            galaxy_diag_aliases = {
                "diagnostic_median_radius.png": "median_radius_origin",
                "diagnostic_mean_radius.png": "mean_radius_origin",
                "diagnostic_std_radius.png": "std_radius_origin",
                "diagnostic_max_radius.png": "max_radius_origin",
                "diagnostic_frac_vr_positive.png": "frac_vr_positive_origin",
                "diagnostic_frac_vr_negative.png": "frac_vr_negative_origin",
                "diagnostic_frac_beyond_cutoff.png": "frac_beyond_cutoff_origin",
                "diagnostic_angular_momentum_z.png": "angular_momentum_z_origin",
            }
            for old, quantity in galaxy_diag_aliases.items():
                write_compat_alias(
                    run_dir / old,
                    run_dir / mode_aware_name(mode_label, physics_label, "primary", quantity, "timeseries", "png"),
                )
            print(f"Saved diagnostic plots in {run_dir}")
            print(f"  Final median_r: {diag['median_r'][-1]:.2f}, L_z: {diag['L_z'][-1]:.2f}")
    else:
        print("Skipping diagnostics (--no-diagnostics or single snapshot).")

    print(f"Done. Outputs in {run_dir}")


if __name__ == "__main__":
    main()
