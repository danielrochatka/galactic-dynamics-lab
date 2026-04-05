#!/usr/bin/env python3
"""
Render declared side-by-side compare outputs from a compare parent directory.
"""

from __future__ import annotations

import argparse
import json
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from framing import SquareViewport, global_viewport_from_snapshots
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
    if pkg == "TPFCore":
        dyn = str(_run_info_effective_value(run_info, "tpf_dynamics_mode", "") or "").strip()
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
    C++ writes left_dir/right_dir relative to cpp_sim cwd (e.g. outputs/RUN/left_TPFCore).
    compare_parent is the resolved compare folder (.../cpp_sim/outputs/RUN).
    Resolving Path(manifest) against repo-root cwd breaks; prefer compare_parent / basename,
    then cpp_sim / manifest path.
    """
    p = Path(manifest_path_str)
    if p.is_absolute():
        return p
    direct = compare_parent / p.name
    if direct.is_dir():
        return direct.resolve()
    # Manifest path is relative to cpp_sim: outputs/RUN/left_*
    cpp_sim = compare_parent.parent.parent
    under_cpp = (cpp_sim / p).resolve()
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
    left_snaps: list[object], right_snaps: list[object], fallback: float
) -> SquareViewport:
    """One square viewport containing all particles from both sides (all matched frames) plus origin."""
    combo = list(left_snaps) + list(right_snaps)
    return global_viewport_from_snapshots(
        combo,
        extra_xy=np.array([[0.0, 0.0]], dtype=np.float64),
        trim_fraction=0.01,
        margin=1.15,
        fallback_half_axis=float(max(1.0, fallback)),
    )


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


def _draw_panel(ax, side: SideData, snap, viewport: SquareViewport | float, *, spatial_display) -> None:
    from render import scatter_frame

    scatter_frame(
        ax,
        snap.positions,
        velocities=getattr(snap, "velocities", None),
        render_radius=viewport,
        spatial_display=spatial_display,
    )
    ax.set_title(
        f"{side.label} / {_panel_physics_text(side.run_info)} / primary trajectory x-y / step={int(snap.step)}",
        color="white",
        fontsize=11,
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
    shared_vp = calculate_compare_smart_viewport(left_snaps, right_snaps, fallback)
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
        f"half_axis={shared_vp.half_axis:.6g} m, display_distance_unit={shared_display.active_distance_unit}, "
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
        fig.suptitle(
            f"Compare {compare_run_id} | left={left_panel_label} "
            f"| right={right_panel_label} "
            f"| step={step} "
            f"| rev={left_rev} "
            f"| display: d={shared_display.active_distance_unit}, t={shared_display.active_time_unit}, "
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
    left_rev = str(_run_info_effective_value(left.run_info, "code_version_label", "unknown") or "unknown")
    mode_aware_initial = _mode_aware_compare_name("initial", left.run_info, right.run_info, ext="png")
    mode_aware_final = _mode_aware_compare_name("final", left.run_info, right.run_info, ext="png")
    save_static(0, mode_aware_initial)
    save_static(len(steps) - 1, mode_aware_final)
    # Compatibility aliases: keep historical filenames as copies of the mode-aware primary artifacts.
    shutil.copy2(parent_dir / mode_aware_initial, parent_dir / "galaxy_initial_compare.png")
    shutil.copy2(parent_dir / mode_aware_final, parent_dir / "galaxy_final_compare.png")

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
            f"| left={left_panel_label} "
            f"| right={right_panel_label} "
            f"| display: d={shared_display.active_distance_unit}, t={shared_display.active_time_unit}, "
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
    out_mp4 = parent_dir / f"{compare_run_id}.mp4"
    out_gif = parent_dir / f"{compare_run_id}.gif"
    try:
        anim.save(str(mode_aware_anim_mp4), writer="ffmpeg", fps=20, dpi=100)
        shutil.copy2(mode_aware_anim_mp4, out_mp4)
    except Exception:
        anim.save(str(mode_aware_anim_gif), writer="pillow", fps=15, dpi=100)
        shutil.copy2(mode_aware_anim_gif, out_gif)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Render side-by-side compare outputs for cpp_sim compare runs.")
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
