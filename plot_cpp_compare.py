#!/usr/bin/env python3
"""
Render declared side-by-side compare outputs from a compare parent directory.
"""

from __future__ import annotations

import argparse
import bisect
import json
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from framing import SquareViewport
from display_units import (
    DisplayUnitConfig,
    display_unit_config_from_run_info,
    format_animation_time_caption,
    series_display_generic_validation,
    spatial_display_for_xy_plot,
)
from render_overlay import build_overlay_spec, draw_galaxy_render_overlay, resolve_overlay_mode


@dataclass
class SideData:
    label: str
    run_dir: Path
    run_info: dict
    overlay_mode: str
    overlay_spec: dict
    snapshots_by_step: dict[int, object]


@dataclass(frozen=True)
class CompareDisplaySelection:
    config: DisplayUnitConfig
    active_distance_unit: str
    active_time_unit: str
    active_velocity_unit: str


def _run_info_configured_value(run_info: dict, key: str, default=None):
    """Configured setting precedence: configured_<key> then legacy <key>."""
    if not run_info:
        return default
    cfg_key = f"configured_{key}"
    if cfg_key in run_info:
        return run_info[cfg_key]
    if key in run_info:
        return run_info[key]
    return default


def _run_info_effective_value(run_info: dict, key: str, default=None):
    """Effective runtime precedence: effective_<key> then configured_<key> then legacy <key>."""
    if not run_info:
        return default
    eff_key = f"effective_{key}"
    if eff_key in run_info:
        return run_info[eff_key]
    cfg = _run_info_configured_value(run_info, key, None)
    if cfg is not None:
        return cfg
    return default


def _slug(raw: str) -> str:
    import re

    s = re.sub(r"[^a-zA-Z0-9]+", "_", str(raw).strip().lower()).strip("_")
    return s or "unknown"


def _simulation_mode_name(run_info: dict) -> str:
    from plot_cpp_run import simulation_mode_name_from_run_info

    return simulation_mode_name_from_run_info(run_info)


def _physics_label(run_info: dict) -> str:
    from plot_cpp_run import physics_label_from_run_info

    mode = _simulation_mode_name(run_info)
    return physics_label_from_run_info(run_info, mode)


def _panel_physics_text(run_info: dict) -> str:
    pkg = str(_run_info_effective_value(run_info, "physics_package", "?") or "?").strip()
    if not pkg:
        return "?"
    if pkg == "Newtonian":
        return "Newtonian baseline package"
    if pkg == "TPFCore":
        dyn = str(_run_info_effective_value(run_info, "tpf_dynamics_mode", "") or "").strip()
        if dyn == "xi_kernel_deformed":
            kernel_mode = str(_run_info_effective_value(run_info, "tpf_4d_xi_kernel_mode", "unknown_kernel") or "unknown_kernel")
            return f"TPFCore xi_kernel_deformed / {kernel_mode}"
        if dyn:
            return f"{pkg} ({dyn})"
    return pkg


def _mode_aware_compare_name(stage: str, left_run_info: dict, right_run_info: dict, *, ext: str) -> str:
    left_lbl = _slug(_physics_label(left_run_info))
    right_lbl = _slug(_physics_label(right_run_info))
    return f"galaxy_compare__{left_lbl}_vs_{right_lbl}__compare__{stage}_side_by_side.{ext}"


def _load_compare_manifest(parent: Path) -> dict:
    p = parent / "compare_manifest.json"
    if not p.exists():
        raise SystemExit(f"Missing compare_manifest.json in {parent}")
    return json.loads(p.read_text(encoding="utf-8"))


def _resolve_side_run_dir(manifest_path_str: str, compare_parent: Path) -> Path:
    """
    C++ writes left_dir/right_dir relative to engine cwd (e.g. outputs/RUN/left_TPFCore).
    compare_parent is the resolved compare folder (.../outputs/RUN).
    Resolving Path(manifest) against repo-root cwd breaks; prefer compare_parent / basename,
    then engine / manifest path.
    """
    p = Path(manifest_path_str)
    if p.is_absolute():
        return p
    direct = compare_parent / p.name
    if direct.is_dir():
        return direct.resolve()
    # Manifest path is relative to engine: outputs/RUN/left_*
    engine = compare_parent.parent.parent
    under_cpp = (engine / p).resolve()
    if under_cpp.is_dir():
        return under_cpp
    cwd_rel = (Path.cwd() / p).resolve()
    if cwd_rel.is_dir():
        return cwd_rel
    return direct


def _load_side_data(run_dir: Path, label: str, overlay_mode_override: str | None) -> SideData:
    from plot_cpp_run import load_all_snapshot_records, load_run_info

    if not run_dir.is_dir():
        raise SystemExit(f"Missing compare child directory: {run_dir}")
    run_info = load_run_info(run_dir)
    overlay_mode = resolve_overlay_mode(run_info, overlay_mode_override)
    overlay_spec = build_overlay_spec(run_dir, run_info)
    records = load_all_snapshot_records(run_dir)
    if not records:
        raise SystemExit(f"No snapshots found in {run_dir}")
    by_step = {snap.step: snap for _, snap in records}
    return SideData(
        label=label,
        run_dir=run_dir,
        run_info=run_info,
        overlay_mode=overlay_mode,
        overlay_spec=overlay_spec,
        snapshots_by_step=by_step,
    )


def matched_steps_strict(left_steps: set[int], right_steps: set[int]) -> list[int]:
    if left_steps != right_steps:
        missing_left = sorted(right_steps - left_steps)
        missing_right = sorted(left_steps - right_steps)
        raise ValueError(
            "Strict step matching failed: "
            f"missing in left={missing_left[:10]}, missing in right={missing_right[:10]}"
        )
    return sorted(left_steps)


def _fallback_render_radius_m(run_info: dict) -> float:
    """CLI-style fallback (m) when viewport extent is degenerate; matches plot_cpp_run default."""
    raw = run_info.get("render_radius", 150.0)
    try:
        return max(1.0, float(raw))
    except (TypeError, ValueError):
        return 150.0


def calculate_compare_smart_viewport(
    left_snaps: list[object],
    right_snaps: list[object],
    fallback: float,
    *,
    coverage: float = 0.95,
) -> SquareViewport:
    """
    Compare-only fixed shared viewport over matched frames.

    For each matched frame, pool left+right points (+origin), compute a minimum-area
    axis-aligned rectangle containing the requested point fraction, convert it to a
    square viewport, and apply a small fixed margin. The final compare viewport is
    the union of those per-frame squares across time.
    """
    MARGIN = 1.10
    fb = float(max(1.0, fallback))
    cov = _validate_compare_coverage(coverage)

    union_xmin: float | None = None
    union_xmax: float | None = None
    union_ymin: float | None = None
    union_ymax: float | None = None

    for left_snap, right_snap in zip(left_snaps, right_snaps):
        left_xy = np.asarray(getattr(left_snap, "positions", []), dtype=np.float64).reshape(-1, 2)
        right_xy = np.asarray(getattr(right_snap, "positions", []), dtype=np.float64).reshape(-1, 2)
        pooled_xy = np.vstack(
            [left_xy, right_xy, np.array([[0.0, 0.0]], dtype=np.float64)]
        )
        x = pooled_xy[:, 0]
        y = pooled_xy[:, 1]
        finite = np.isfinite(x) & np.isfinite(y)
        if not np.any(finite):
            continue
        points = np.column_stack([x[finite], y[finite]])
        xmin, xmax, ymin, ymax = minimum_axis_aligned_box_covering_fraction(points, cov)
        cx = 0.5 * (xmin + xmax)
        cy = 0.5 * (ymin + ymax)
        half_axis = 0.5 * max(xmax - xmin, ymax - ymin) * MARGIN
        if not np.isfinite(half_axis) or half_axis <= 0.0:
            half_axis = fb
        sq_xmin = cx - half_axis
        sq_xmax = cx + half_axis
        sq_ymin = cy - half_axis
        sq_ymax = cy + half_axis
        union_xmin = sq_xmin if union_xmin is None else min(union_xmin, sq_xmin)
        union_xmax = sq_xmax if union_xmax is None else max(union_xmax, sq_xmax)
        union_ymin = sq_ymin if union_ymin is None else min(union_ymin, sq_ymin)
        union_ymax = sq_ymax if union_ymax is None else max(union_ymax, sq_ymax)

    if union_xmin is None or union_xmax is None or union_ymin is None or union_ymax is None:
        return SquareViewport(0.0, 0.0, fb)

    shared_cx = 0.5 * (union_xmin + union_xmax)
    shared_cy = 0.5 * (union_ymin + union_ymax)
    shared_half_axis = 0.5 * max(union_xmax - union_xmin, union_ymax - union_ymin)
    if not np.isfinite(shared_half_axis) or shared_half_axis <= 0.0:
        shared_half_axis = fb
    return SquareViewport(shared_cx, shared_cy, max(shared_half_axis, 1e-30))


def _validate_compare_coverage(raw: object) -> float:
    try:
        cov = float(raw)
    except (TypeError, ValueError) as e:
        raise ValueError(f"plot_compare_smart_zoom_coverage must be numeric, got {raw!r}") from e
    if not np.isfinite(cov):
        raise ValueError("plot_compare_smart_zoom_coverage must be finite")
    if cov <= 0.0 or cov > 1.0:
        raise ValueError("plot_compare_smart_zoom_coverage must be in (0, 1]")
    return cov


def minimum_axis_aligned_box_covering_fraction(points_xy: np.ndarray, coverage: float) -> tuple[float, float, float, float]:
    """
    Return (xmin, xmax, ymin, ymax) for the minimum-area axis-aligned rectangle
    that contains at least coverage * N points.
    """
    arr = np.asarray(points_xy, dtype=np.float64).reshape(-1, 2)
    if arr.size == 0:
        return (0.0, 0.0, 0.0, 0.0)
    finite = np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 1])
    pts = arr[finite]
    if pts.shape[0] == 0:
        return (0.0, 0.0, 0.0, 0.0)

    cov = _validate_compare_coverage(coverage)
    n = int(pts.shape[0])
    k = max(1, int(np.ceil(cov * n)))

    x = pts[:, 0]
    y = pts[:, 1]
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    ys = y[order]

    best_area = np.inf
    best: tuple[float, float, float, float] | None = None
    for i in range(n):
        y_sorted_window: list[float] = []
        for j in range(i, n):
            bisect.insort(y_sorted_window, float(ys[j]))
            m = j - i + 1
            if m < k:
                continue

            x_min = float(xs[i])
            x_max = float(xs[j])
            width = x_max - x_min
            best_y_span = np.inf
            best_y_min = float(y_sorted_window[0])
            best_y_max = float(y_sorted_window[-1])
            for t in range(0, m - k + 1):
                y_min = y_sorted_window[t]
                y_max = y_sorted_window[t + k - 1]
                y_span = y_max - y_min
                if y_span < best_y_span:
                    best_y_span = y_span
                    best_y_min = float(y_min)
                    best_y_max = float(y_max)
            area = width * best_y_span
            if area < best_area:
                best_area = area
                best = (x_min, x_max, best_y_min, best_y_max)

    if best is not None:
        return best

    return (float(np.min(x)), float(np.max(x)), float(np.min(y)), float(np.max(y)))


def _resolve_compare_coverage(left_run_info: dict, right_run_info: dict) -> float:
    raw = _run_info_effective_value(left_run_info, "plot_compare_smart_zoom_coverage", None)
    if raw is None:
        raw = _run_info_effective_value(right_run_info, "plot_compare_smart_zoom_coverage", None)
    if raw is None:
        raw = 0.95
    return _validate_compare_coverage(raw)


def _resolve_compare_preferred_unit(left: str, right: str) -> str:
    l = (left or "auto").strip()
    r = (right or "auto").strip()
    if l == "auto" and r == "auto":
        return "auto"
    if l == "auto":
        return r
    if r == "auto":
        return l
    if l == r:
        return l
    return "auto"


def resolve_compare_display_selection(
    left_cfg: DisplayUnitConfig,
    right_cfg: DisplayUnitConfig,
    shared_half_axis_m: float,
    max_time_s: float,
    max_speed_m_s: float,
) -> CompareDisplaySelection:
    cfg = DisplayUnitConfig(
        distance_unit=_resolve_compare_preferred_unit(left_cfg.distance_unit, right_cfg.distance_unit),
        time_unit=_resolve_compare_preferred_unit(left_cfg.time_unit, right_cfg.time_unit),
        velocity_unit=_resolve_compare_preferred_unit(left_cfg.velocity_unit, right_cfg.velocity_unit),
        units_in_overlay=left_cfg.units_in_overlay and right_cfg.units_in_overlay,
        show_unit_reference=left_cfg.show_unit_reference and right_cfg.show_unit_reference,
    )
    spatial = spatial_display_for_xy_plot(
        "galaxy_compare",
        float(max(shared_half_axis_m, 1.0)),
        preferred_unit=cfg.distance_unit,
    )
    series = series_display_generic_validation(
        float(max(shared_half_axis_m, 1.0)),
        max_time_s=float(max_time_s),
        max_speed_m_s=float(max_speed_m_s),
        preferred_distance_unit=cfg.distance_unit,
        preferred_time_unit=cfg.time_unit,
        preferred_velocity_unit=cfg.velocity_unit,
    )
    return CompareDisplaySelection(
        config=cfg,
        active_distance_unit=spatial.unit,
        active_time_unit=series.time_unit,
        active_velocity_unit=series.speed_unit,
    )




def radial_kinematics(snapshot):
    pos = np.asarray(getattr(snapshot, "positions", []), dtype=np.float64).reshape(-1, 2)
    vel = np.asarray(getattr(snapshot, "velocities", []), dtype=np.float64).reshape(-1, 2)
    n = min(pos.shape[0], vel.shape[0])
    pos = pos[:n]
    vel = vel[:n]
    r = np.linalg.norm(pos, axis=1)
    safe = np.where(r > 0.0, r, 1.0)
    r_hat = pos / safe[:, None]
    r_hat[r <= 0.0] = 0.0
    t_hat = np.column_stack([-r_hat[:, 1], r_hat[:, 0]])
    v_radial = np.sum(vel * r_hat, axis=1)
    v_tangential = np.sum(vel * t_hat, axis=1)
    return {
        "radius": r,
        "r_hat": r_hat,
        "t_hat": t_hat,
        "v_radial": v_radial,
        "v_tangential": v_tangential,
        "v_t": np.abs(v_tangential),
    }


def estimate_snapshot_accelerations(snaps):
    if not snaps:
        return []
    out = []
    n = len(snaps)
    for i in range(n):
        cur_v = np.asarray(getattr(snaps[i], "velocities", []), dtype=np.float64)
        if n == 1:
            out.append(None)
            continue
        if i == 0:
            v0 = cur_v
            v1 = np.asarray(getattr(snaps[1], "velocities", []), dtype=np.float64)
            dt = float(snaps[1].time) - float(snaps[0].time)
        elif i == n - 1:
            v0 = np.asarray(getattr(snaps[n - 2], "velocities", []), dtype=np.float64)
            v1 = cur_v
            dt = float(snaps[n - 1].time) - float(snaps[n - 2].time)
        else:
            v0 = np.asarray(getattr(snaps[i - 1], "velocities", []), dtype=np.float64)
            v1 = np.asarray(getattr(snaps[i + 1], "velocities", []), dtype=np.float64)
            dt = float(snaps[i + 1].time) - float(snaps[i - 1].time)
        if (not np.isfinite(dt)) or dt == 0.0:
            out.append(None)
            continue
        m = min(v0.shape[0], v1.shape[0], cur_v.shape[0])
        out.append((v1[:m] - v0[:m]) / dt)
    return out


def binned_profile(x, y, bins=30):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    if x.size == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    if xmax <= xmin:
        return np.array([xmin]), np.array([np.median(y)]), np.array([np.percentile(y, 25)]), np.array([np.percentile(y, 75)])
    edges = np.linspace(xmin, xmax, int(max(2, bins)) + 1)
    mids=[]; med=[]; p25=[]; p75=[]
    for i in range(edges.size - 1):
        left = edges[i]
        right = edges[i + 1]
        mask = (x >= left) & (x < right if i < edges.size - 2 else x <= right)
        if not np.any(mask):
            continue
        yy = y[mask]
        mids.append(0.5 * (left + right))
        med.append(float(np.median(yy)))
        p25.append(float(np.percentile(yy, 25)))
        p75.append(float(np.percentile(yy, 75)))
    return np.asarray(mids), np.asarray(med), np.asarray(p25), np.asarray(p75)


def save_compare_profile_plot(parent_dir, mode_name, alias_name, left_xy, right_xy, *, left_label, right_label, xlabel, ylabel, title, add_ref1=False):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(9, 6), facecolor="white")
    style_compare_diagnostic_axes(ax)

    def _draw(side_xy, label, color):
        if side_xy is None:
            return
        x, y = side_xy
        bx, bm, bp25, bp75 = binned_profile(x, y, bins=30)
        if bx.size == 0:
            return
        ax.plot(bx, bm, color=color, lw=2, label=label)
        ax.fill_between(bx, bp25, bp75, color=color, alpha=0.2)

    _draw(left_xy, left_label, "tab:cyan")
    _draw(right_xy, right_label, "tab:orange")
    if add_ref1:
        ax.axhline(1.0, color="0.2", lw=1, ls="--", alpha=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    style_compare_diagnostic_legend(ax)
    fig.tight_layout()
    for nm in (mode_name, alias_name):
        fig.savefig(parent_dir / nm, dpi=150, facecolor="white", edgecolor="none")
    plt.close(fig)


def style_compare_diagnostic_axes(ax) -> None:
    ax.set_facecolor("white")
    ax.tick_params(colors="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.title.set_color("black")
    for spine in ax.spines.values():
        spine.set_color("0.35")
    ax.grid(True, alpha=0.25, color="0.75")


def style_compare_diagnostic_legend(ax) -> None:
    handles, _ = ax.get_legend_handles_labels()
    if not handles:
        return
    leg = ax.legend(facecolor="white", edgecolor="0.5")
    for text in leg.get_texts():
        text.set_color("black")


def _time_display_factor(unit: str) -> float:
    u = str(unit or "s")
    seconds_per_day = 86400.0
    seconds_per_julian_year = 365.25 * seconds_per_day
    seconds_per_kyr = 1.0e3 * seconds_per_julian_year
    seconds_per_myr = 1.0e6 * seconds_per_julian_year
    if u == "Myr":
        return 1.0 / seconds_per_myr
    if u == "kyr":
        return 1.0 / seconds_per_kyr
    if u == "yr":
        return 1.0 / seconds_per_julian_year
    if u == "day":
        return 1.0 / seconds_per_day
    if u == "hr":
        return 1.0 / 3600.0
    if u == "min":
        return 1.0 / 60.0
    return 1.0


def _velocity_display_factor(unit: str) -> float:
    return 1.0e-3 if str(unit or "m/s") == "km/s" else 1.0


def _draw_panel(ax, side: SideData, snap, viewport: SquareViewport | float, *, spatial_display) -> None:
    from render import scatter_frame

    scatter_frame(
        ax,
        snap.positions,
        velocities=getattr(snap, "velocities", None),
        render_radius=viewport,
        spatial_display=spatial_display,
    )
    ax.set_title("")
    panel_lines = _compare_panel_metadata_lines(side)
    if panel_lines:
        ax.text(
            0.02,
            0.98,
            "\n".join(panel_lines),
            transform=ax.transAxes,
            ha="left",
            va="top",
            color="#66ff66",
            fontsize=5.5,
            bbox={
                "facecolor": (0.0, 0.0, 0.0, 0.72),
                "edgecolor": (0.2, 0.8, 0.2, 0.85),
                "linewidth": 0.6,
                "boxstyle": "round,pad=0.20",
            },
            zorder=30,
        )
    if side.overlay_mode != "none":
        draw_galaxy_render_overlay(
            ax,
            side.overlay_mode,
            side.overlay_spec,
            run_info=side.run_info,
            step=int(snap.step),
            time_s=float(snap.time),
        )


def _compare_panel_metadata_lines(side: SideData) -> list[str]:
    run_info = side.run_info or {}
    pkg = str(_run_info_effective_value(run_info, "physics_package", "?") or "?").strip() or "?"
    role = "left / baseline" if side.label == "left_primary" else "right / compare"
    lines = [pkg, f"package: {pkg}", f"role: {role}"]

    dyn = str(_run_info_effective_value(run_info, "tpf_dynamics_mode", "") or "").strip()
    kernel = str(_run_info_effective_value(run_info, "tpf_4d_xi_kernel_mode", "") or "").strip()
    if pkg == "TPFCore":
        if dyn:
            lines.append(f"dynamics: {dyn}")
        if dyn == "xi_kernel_deformed" and kernel:
            lines.append(f"kernel: {kernel}")
        coupling = _run_info_effective_value(run_info, "tpf_4d_xi_kernel_coupling", None)
        if coupling is None:
            coupling = _run_info_effective_value(run_info, "tpf_vdsg_coupling", None)
        if coupling is not None and str(coupling).strip() != "":
            lines.append(f"coupling: {coupling}")
        factor_mode = _run_info_effective_value(run_info, "tpf_factor_mode", None)
        if factor_mode is not None and str(factor_mode).strip() != "":
            lines.append(f"factor_mode: {factor_mode}")
    return lines


def render_compare(
    parent_dir: Path,
    no_animation: bool = False,
    overlay_mode: str | None = None,
) -> None:
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt

    mf = _load_compare_manifest(parent_dir)
    left_dir = _resolve_side_run_dir(str(mf["left_dir"]), parent_dir)
    right_dir = _resolve_side_run_dir(str(mf["right_dir"]), parent_dir)
    compare_run_id = str(mf.get("compare_run_id", parent_dir.name))

    left = _load_side_data(left_dir, "left_primary", overlay_mode)
    right = _load_side_data(right_dir, "right_compare", overlay_mode)
    steps = matched_steps_strict(set(left.snapshots_by_step.keys()), set(right.snapshots_by_step.keys()))
    left_snaps = [left.snapshots_by_step[s] for s in steps]
    right_snaps = [right.snapshots_by_step[s] for s in steps]

    fb_l = _fallback_render_radius_m(left.run_info)
    fb_r = _fallback_render_radius_m(right.run_info)
    fallback = float(max(fb_l, fb_r, 1.0))
    compare_coverage = _resolve_compare_coverage(left.run_info, right.run_info)
    shared_vp = calculate_compare_smart_viewport(
        left_snaps, right_snaps, fallback, coverage=compare_coverage
    )
    max_time_s = max(float(left_snaps[-1].time), float(right_snaps[-1].time))
    max_speed_m_s = max(
        float(np.max(np.linalg.norm(s.velocities, axis=1))) for s in (left_snaps + right_snaps)
    )
    left_disp_cfg = display_unit_config_from_run_info(left.run_info)
    right_disp_cfg = display_unit_config_from_run_info(right.run_info)
    shared_display = resolve_compare_display_selection(
        left_disp_cfg, right_disp_cfg, shared_vp.half_axis, max_time_s, max_speed_m_s
    )
    spatial_display = spatial_display_for_xy_plot(
        "galaxy_compare",
        shared_vp.half_axis,
        preferred_unit=shared_display.config.distance_unit,
    )
    for side in (left, right):
        side.overlay_spec["active_display_distance_unit"] = shared_display.active_distance_unit
        side.overlay_spec["active_display_time_unit"] = shared_display.active_time_unit
        side.overlay_spec["active_display_velocity_unit"] = shared_display.active_velocity_unit
        side.overlay_spec["display_units_in_overlay"] = shared_display.config.units_in_overlay
        side.overlay_spec["display_show_unit_reference"] = shared_display.config.show_unit_reference
    print(
        f"Compare smart framing: center=({shared_vp.center_x:.6g}, {shared_vp.center_y:.6g}) m, "
        f"half_axis={shared_vp.half_axis:.6g} m, coverage={compare_coverage:.6g}, "
        f"display_distance_unit={shared_display.active_distance_unit}, "
        f"display_time_unit={shared_display.active_time_unit}, display_velocity_unit={shared_display.active_velocity_unit}"
    )
    print(
        f"Compare matched-step span: first_step={steps[0]}, last_step={steps[-1]}, "
        f"shared_half_axis_m={shared_vp.half_axis:.6g}"
    )

    def save_static(step_idx: int, out_name: str) -> None:
        ls = left_snaps[step_idx]
        rs = right_snaps[step_idx]
        fig, axes = plt.subplots(1, 2, figsize=(16, 8), facecolor="black")
        for ax in axes:
            ax.set_facecolor("black")
            ax.tick_params(colors="gray")
            for s in ax.spines.values():
                s.set_color("gray")
        _draw_panel(axes[0], left, ls, shared_vp, spatial_display=spatial_display)
        _draw_panel(axes[1], right, rs, shared_vp, spatial_display=spatial_display)
        step = int(steps[step_idx])
        tc = format_animation_time_caption(
            float(ls.time),
            "galaxy_compare",
            preferred_time_unit=shared_display.config.time_unit,
            active_time_unit=shared_display.active_time_unit,
        )
        fig.suptitle(
            f"Compare {compare_run_id} | step={step} | {tc} "
            f"| display units: d={shared_display.active_distance_unit}, t={shared_display.active_time_unit}, "
            f"v={shared_display.active_velocity_unit}",
            color="white",
            fontsize=11,
        )
        if shared_display.config.show_unit_reference:
            fig.text(
                0.995,
                0.01,
                f"distance display = {shared_display.active_distance_unit}\n"
                f"time display = {shared_display.active_time_unit}\n"
                f"velocity display = {shared_display.active_velocity_unit}",
                ha="right",
                va="bottom",
                fontsize=8,
                color="white",
            )
        fig.tight_layout()
        fig.savefig(parent_dir / out_name, dpi=150, facecolor="black", edgecolor="none")
        plt.close(fig)

    left_panel_label = _panel_physics_text(left.run_info)
    right_panel_label = _panel_physics_text(right.run_info)
    mode_aware_initial = _mode_aware_compare_name("initial", left.run_info, right.run_info, ext="png")
    mode_aware_final = _mode_aware_compare_name("final", left.run_info, right.run_info, ext="png")
    save_static(0, mode_aware_initial)
    save_static(len(steps) - 1, mode_aware_final)
    # Compatibility aliases: keep historical filenames as copies of the mode-aware primary artifacts.
    shutil.copy2(parent_dir / mode_aware_initial, parent_dir / "galaxy_initial_compare.png")
    shutil.copy2(parent_dir / mode_aware_final, parent_dir / "galaxy_final_compare.png")
    required = [
        parent_dir / mode_aware_initial,
        parent_dir / mode_aware_final,
        parent_dir / "galaxy_initial_compare.png",
        parent_dir / "galaxy_final_compare.png",
    ]
    missing = [str(p.name) for p in required if not p.exists()]
    if missing:
        raise SystemExit(
            "Compare render failed: missing expected PNG outputs in "
            f"{parent_dir}: {', '.join(missing)}"
        )
    print(
        "Compare static renders generated: "
        + ", ".join([mode_aware_initial, mode_aware_final, "galaxy_initial_compare.png", "galaxy_final_compare.png"])
    )

    try:
        dist_scale = spatial_display.factor
        dist_unit = getattr(spatial_display, "unit", shared_display.active_distance_unit)
        time_scale = _time_display_factor(shared_display.active_time_unit)
        vel_scale = _velocity_display_factor(shared_display.active_velocity_unit)
        vel_unit = shared_display.active_velocity_unit
        lk = radial_kinematics(left_snaps[-1])
        rk = radial_kinematics(right_snaps[-1])
        la_all = estimate_snapshot_accelerations(left_snaps)
        ra_all = estimate_snapshot_accelerations(right_snaps)
        la = la_all[-1] if la_all else None
        ra = ra_all[-1] if ra_all else None

        def inward_accel(kin, acc):
            if acc is None:
                return None
            m = min(acc.shape[0], kin["r_hat"].shape[0])
            ar = np.sum(acc[:m] * kin["r_hat"][:m], axis=1)
            return -ar

        l_in = inward_accel(lk, la)
        r_in = inward_accel(rk, ra)
        left_name = _slug(_physics_label(left.run_info))
        right_name = _slug(_physics_label(right.run_info))
        pref = f"galaxy_compare__{left_name}_vs_{right_name}__compare__"
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(9, 6), facecolor="white")
        style_compare_diagnostic_axes(ax)
        if l_in is not None:
            ax.scatter(lk["radius"] * dist_scale, l_in, s=8, alpha=0.2, color="tab:cyan", label=left_panel_label)
        if r_in is not None:
            ax.scatter(rk["radius"] * dist_scale, r_in, s=8, alpha=0.2, color="tab:orange", label=right_panel_label)
        ax.set_xlabel(f"radius [{dist_unit}]")
        ax.set_ylabel("inward radial acceleration [m/s^2]")
        ax.set_title("Final snapshot inward radial acceleration vs radius (raw)")
        leg = ax.legend(facecolor="white", edgecolor="0.5")
        for text in leg.get_texts():
            text.set_color("black")
        fig.tight_layout()
        fig.savefig(parent_dir / (pref + "radial_acceleration_vs_radius_final.png"), dpi=150, facecolor="white", edgecolor="none")
        fig.savefig(parent_dir / "compare_radial_acceleration_vs_radius_final.png", dpi=150, facecolor="white", edgecolor="none")
        plt.close(fig)
        save_compare_profile_plot(parent_dir, pref+"binned_inward_acceleration_vs_radius_final.png", "compare_binned_inward_acceleration_vs_radius_final.png", (lk["radius"]*dist_scale, l_in) if l_in is not None else None, (rk["radius"]*dist_scale, r_in) if r_in is not None else None, left_label=left_panel_label, right_label=right_panel_label, xlabel=f"radius [{dist_unit}]", ylabel="inward radial acceleration [m/s^2]", title="Binned inward radial acceleration vs radius")
        save_compare_profile_plot(parent_dir, pref+"rotation_curve_final.png", "compare_rotation_curve_final.png", (lk["radius"]*dist_scale, lk["v_t"] * vel_scale), (rk["radius"]*dist_scale, rk["v_t"] * vel_scale), left_label=left_panel_label, right_label=right_panel_label, xlabel=f"radius [{dist_unit}]", ylabel=f"tangential speed [{vel_unit}]", title="Rotation curve (final snapshot)")
        l_cent = np.divide(lk["v_t"]**2, lk["radius"], out=np.full_like(lk["radius"], np.nan), where=lk["radius"]>0)
        r_cent = np.divide(rk["v_t"]**2, rk["radius"], out=np.full_like(rk["radius"], np.nan), where=rk["radius"]>0)
        save_compare_profile_plot(parent_dir, pref+"centripetal_profile_final.png", "compare_centripetal_profile_final.png", (lk["radius"]*dist_scale, l_cent), (rk["radius"]*dist_scale, r_cent), left_label=left_panel_label, right_label=right_panel_label, xlabel=f"radius [{dist_unit}]", ylabel="v_t^2 / r [m/s^2]", title="Centripetal profile (final snapshot)")
        save_compare_profile_plot(parent_dir, pref+"radial_velocity_vs_radius_final.png", "compare_radial_velocity_vs_radius_final.png", (lk["radius"]*dist_scale, lk["v_radial"] * vel_scale), (rk["radius"]*dist_scale, rk["v_radial"] * vel_scale), left_label=left_panel_label, right_label=right_panel_label, xlabel=f"radius [{dist_unit}]", ylabel=f"radial velocity [{vel_unit}]", title="Radial velocity vs radius (final snapshot)")

        lmed=[];lp25=[];lp75=[];rmed=[];rp25=[];rp75=[];tvals=[]
        for ls, rs in zip(left_snaps, right_snaps):
            lr = radial_kinematics(ls)["radius"]*dist_scale
            rr = radial_kinematics(rs)["radius"]*dist_scale
            tvals.append(float(ls.time))
            lmed.append(np.median(lr)); lp25.append(np.percentile(lr,25)); lp75.append(np.percentile(lr,75))
            rmed.append(np.median(rr)); rp25.append(np.percentile(rr,25)); rp75.append(np.percentile(rr,75))
        fig, ax = plt.subplots(1,1,figsize=(9,6),facecolor='white')
        style_compare_diagnostic_axes(ax)
        t=np.asarray(tvals) * time_scale; ax.plot(t,lmed,color='tab:cyan',label=left_panel_label); ax.fill_between(t,lp25,lp75,color='tab:cyan',alpha=0.2)
        ax.plot(t,rmed,color='tab:orange',label=right_panel_label); ax.fill_between(t,rp25,rp75,color='tab:orange',alpha=0.2)
        ax.set_xlabel(f"time [{shared_display.active_time_unit}]"); ax.set_ylabel(f"radius [{dist_unit}]"); ax.set_title('Radius percentiles over time')
        style_compare_diagnostic_legend(ax)
        fig.tight_layout()
        fig.savefig(parent_dir/(pref+"radius_percentiles_over_time.png"),dpi=150,facecolor='white',edgecolor='none'); fig.savefig(parent_dir/"compare_radius_percentiles_over_time.png",dpi=150,facecolor='white',edgecolor='none'); plt.close(fig)

        if l_in is not None and r_in is not None:
            l_x, l_med, _, _ = binned_profile(lk["radius"] * dist_scale, l_in, bins=30)
            r_x, r_med, _, _ = binned_profile(rk["radius"] * dist_scale, r_in, bins=30)
            common_x = np.union1d(l_x, r_x)
            if common_x.size > 0:
                l_interp = np.interp(common_x, l_x, l_med, left=np.nan, right=np.nan)
                r_interp = np.interp(common_x, r_x, r_med, left=np.nan, right=np.nan)
                valid = np.isfinite(l_interp) & np.isfinite(r_interp) & (np.abs(l_interp) > 0.0)
                if np.any(~valid):
                    print(f"Warning: skipped {int(np.sum(~valid))} invalid acceleration-ratio bins")
                ratio = np.divide(r_interp[valid], l_interp[valid])
                save_compare_profile_plot(parent_dir, pref+"acceleration_ratio_vs_radius_final.png", "compare_acceleration_ratio_vs_radius_final.png", None, (common_x[valid], ratio), left_label=left_panel_label, right_label=f"{right_panel_label}/{left_panel_label}", xlabel=f"radius [{dist_unit}]", ylabel="inward_accel_right / inward_accel_left", title="Acceleration ratio vs radius (final snapshot)", add_ref1=True)
    except Exception as e:
        print(f"Warning: compare diagnostic overlay plots failed: {e}")

    if no_animation:
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 8), facecolor="black")
    for ax in axes:
        ax.set_facecolor("black")
        ax.tick_params(colors="gray")
        for s in ax.spines.values():
            s.set_color("gray")
    if shared_display.config.show_unit_reference:
        fig.text(
            0.995,
            0.01,
            f"distance display = {shared_display.active_distance_unit}\n"
            f"time display = {shared_display.active_time_unit}\n"
            f"velocity display = {shared_display.active_velocity_unit}",
            ha="right",
            va="bottom",
            fontsize=8,
            color="white",
        )

    def animate(i: int):
        _draw_panel(axes[0], left, left_snaps[i], shared_vp, spatial_display=spatial_display)
        _draw_panel(axes[1], right, right_snaps[i], shared_vp, spatial_display=spatial_display)
        t = float(left_snaps[i].time)
        tc = format_animation_time_caption(
            t,
            "galaxy_compare",
            preferred_time_unit=shared_display.config.time_unit,
            active_time_unit=shared_display.active_time_unit,
        )
        fig.suptitle(
            f"Compare {compare_run_id} | step={steps[i]} | {tc} "
            f"| display units: d={shared_display.active_distance_unit}, t={shared_display.active_time_unit}, "
            f"v={shared_display.active_velocity_unit}",
            color="white",
            fontsize=11,
        )
        return []

    anim = animation.FuncAnimation(fig, animate, frames=len(steps), interval=50, blit=False)
    mode_aware_anim_mp4 = parent_dir / _mode_aware_compare_name(
        "animation", left.run_info, right.run_info, ext="mp4"
    )
    mode_aware_anim_gif = parent_dir / _mode_aware_compare_name(
        "animation", left.run_info, right.run_info, ext="gif"
    )
    try:
        anim.save(str(mode_aware_anim_mp4), writer="ffmpeg", fps=20, dpi=100)
    except Exception:
        anim.save(str(mode_aware_anim_gif), writer="pillow", fps=15, dpi=100)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Render side-by-side compare outputs for engine compare runs.")
    ap.add_argument("compare_parent_dir", type=Path)
    ap.add_argument("--no-animation", action="store_true")
    ap.add_argument(
        "--render-overlay-mode",
        choices=("none", "minimal", "audit_full"),
        default=None,
        help="Override per-side run_info overlay mode.",
    )
    args = ap.parse_args()
    render_compare(args.compare_parent_dir.resolve(), args.no_animation, args.render_overlay_mode)


if __name__ == "__main__":
    main()
