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
from render_overlay import build_overlay_spec, draw_galaxy_render_overlay, resolve_overlay_mode


@dataclass
class SideData:
    label: str
    run_dir: Path
    run_info: dict
    overlay_mode: str
    overlay_spec: dict
    snapshots_by_step: dict[int, object]


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


def _draw_panel(ax, side: SideData, snap, viewport: SquareViewport | float) -> None:
    from render import scatter_frame

    scatter_frame(
        ax,
        snap.positions,
        velocities=getattr(snap, "velocities", None),
        render_radius=viewport,
    )
    ax.set_title(
        f"{side.label} / {str(side.run_info.get('physics_package', '?')).lower()} / primary trajectory x-y",
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
    print(
        f"Compare smart framing: center=({shared_vp.center_x:.6g}, {shared_vp.center_y:.6g}) m, "
        f"half_axis={shared_vp.half_axis:.6g} m"
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
        _draw_panel(axes[0], left, ls, shared_vp)
        _draw_panel(axes[1], right, rs, shared_vp)
        fig.suptitle(
            f"Compare {compare_run_id} | left={left.run_info.get('physics_package', '?')} "
            f"| right={right.run_info.get('physics_package', '?')} "
            f"| rev={left.run_info.get('code_version_label', 'unknown')}",
            color="white",
            fontsize=11,
        )
        fig.tight_layout()
        fig.savefig(parent_dir / out_name, dpi=150, facecolor="black", edgecolor="none")
        plt.close(fig)

    mode_aware_initial = (
        f"galaxy_compare__{str(left.run_info.get('physics_package', '?')).lower()}_vs_"
        f"{str(right.run_info.get('physics_package', '?')).lower()}__compare__initial_side_by_side.png"
    )
    mode_aware_final = (
        f"galaxy_compare__{str(left.run_info.get('physics_package', '?')).lower()}_vs_"
        f"{str(right.run_info.get('physics_package', '?')).lower()}__compare__final_side_by_side.png"
    )
    save_static(0, mode_aware_initial)
    save_static(len(steps) - 1, mode_aware_final)
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

    def animate(i: int):
        _draw_panel(axes[0], left, left_snaps[i], shared_vp)
        _draw_panel(axes[1], right, right_snaps[i], shared_vp)
        t = float(left_snaps[i].time)
        fig.suptitle(
            f"Compare {compare_run_id} | step={steps[i]} | t={t:.6g} "
            f"| left={left.run_info.get('physics_package', '?')} "
            f"| right={right.run_info.get('physics_package', '?')}",
            color="white",
            fontsize=11,
        )
        return []

    anim = animation.FuncAnimation(fig, animate, frames=len(steps), interval=50, blit=False)
    mode_aware_anim_mp4 = (
        parent_dir
        / f"galaxy_compare__{str(left.run_info.get('physics_package', '?')).lower()}_vs_"
          f"{str(right.run_info.get('physics_package', '?')).lower()}__compare__animation_side_by_side.mp4"
    )
    mode_aware_anim_gif = (
        parent_dir
        / f"galaxy_compare__{str(left.run_info.get('physics_package', '?')).lower()}_vs_"
          f"{str(right.run_info.get('physics_package', '?')).lower()}__compare__animation_side_by_side.gif"
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
