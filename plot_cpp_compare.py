#!/usr/bin/env python3
"""
Render declared side-by-side compare outputs from a compare parent directory.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from render_overlay import build_overlay_spec, draw_galaxy_render_overlay, resolve_overlay_mode
from plot_cpp_run import resolve_galaxy_radius_meters, static_viewport_radius_validated


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
    """CLI-style fallback (m) when velocity gating needs a floor; matches plot_cpp_run default."""
    raw = run_info.get("render_radius", 150.0)
    try:
        return max(1.0, float(raw))
    except (TypeError, ValueError):
        return 150.0


def _validated_radii_pair(
    snap_l,
    snap_r,
    galaxy_radius_m_l: float,
    galaxy_radius_m_r: float,
    fallback_l: float,
    fallback_r: float,
) -> tuple[float, float]:
    """Per-side validated static viewport half-axes (plot_cpp_run)."""
    r_l = static_viewport_radius_validated(
        snap_l.positions,
        snap_l.velocities,
        galaxy_radius_m_l,
        fallback_l,
    )
    r_r = static_viewport_radius_validated(
        snap_r.positions,
        snap_r.velocities,
        galaxy_radius_m_r,
        fallback_r,
    )
    return r_l, r_r


def _compare_shared_radius_for_pair(
    snap_l,
    snap_r,
    galaxy_radius_m_l: float,
    galaxy_radius_m_r: float,
    fallback_l: float,
    fallback_r: float,
) -> float:
    """
    Single half-axis for both panels: max(left, right) validated radii.
    Use only with --shared-viewport: identical scale, but when one side is much more extended
    than the other, the compact side can look empty and motion can be invisible.
    """
    r_l, r_r = _validated_radii_pair(
        snap_l, snap_r, galaxy_radius_m_l, galaxy_radius_m_r, fallback_l, fallback_r
    )
    return max(r_l, r_r)


def _draw_panel(ax, side: SideData, snap, radius: float) -> None:
    from render import scatter_frame

    scatter_frame(ax, snap.positions, velocities=getattr(snap, "velocities", None), render_radius=radius)
    ax.set_title(f"{side.label} ({side.run_info.get('physics_package', '?')})", color="white", fontsize=11)
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
    *,
    shared_viewport: bool = False,
) -> None:
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt

    mf = _load_compare_manifest(parent_dir)
    left_dir = _resolve_side_run_dir(str(mf["left_dir"]), parent_dir)
    right_dir = _resolve_side_run_dir(str(mf["right_dir"]), parent_dir)
    compare_run_id = str(mf.get("compare_run_id", parent_dir.name))

    left = _load_side_data(left_dir, "left", overlay_mode)
    right = _load_side_data(right_dir, "right", overlay_mode)
    steps = matched_steps_strict(set(left.snapshots_by_step.keys()), set(right.snapshots_by_step.keys()))
    left_snaps = [left.snapshots_by_step[s] for s in steps]
    right_snaps = [right.snapshots_by_step[s] for s in steps]

    fb_l = _fallback_render_radius_m(left.run_info)
    fb_r = _fallback_render_radius_m(right.run_info)
    galaxy_radius_m_l = resolve_galaxy_radius_meters(left.run_info, left_snaps, fb_l)
    galaxy_radius_m_r = resolve_galaxy_radius_meters(right.run_info, right_snaps, fb_r)

    print(
        "Compare viewport: "
        + (
            "shared max(left, right) — same axis scale on both panels"
            if shared_viewport
            else "per-panel — each side uses its own validated radius (recommended when scales differ a lot)"
        )
    )

    def save_static(step_idx: int, out_name: str) -> None:
        ls = left_snaps[step_idx]
        rs = right_snaps[step_idx]
        r_l, r_r = _validated_radii_pair(
            ls, rs, galaxy_radius_m_l, galaxy_radius_m_r, fb_l, fb_r
        )
        if shared_viewport:
            r_l = r_r = max(r_l, r_r)
        fig, axes = plt.subplots(1, 2, figsize=(16, 8), facecolor="black")
        for ax in axes:
            ax.set_facecolor("black")
            ax.tick_params(colors="gray")
            for s in ax.spines.values():
                s.set_color("gray")
        _draw_panel(axes[0], left, ls, r_l)
        _draw_panel(axes[1], right, rs, r_r)
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

    save_static(0, "galaxy_initial_compare.png")
    save_static(len(steps) - 1, "galaxy_final_compare.png")

    if no_animation:
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 8), facecolor="black")
    for ax in axes:
        ax.set_facecolor("black")
        ax.tick_params(colors="gray")
        for s in ax.spines.values():
            s.set_color("gray")

    raw_radii = [
        _validated_radii_pair(
            left_snaps[i],
            right_snaps[i],
            galaxy_radius_m_l,
            galaxy_radius_m_r,
            fb_l,
            fb_r,
        )
        for i in range(len(steps))
    ]
    if shared_viewport:
        fixed_shared = max(max(r_l, r_r) for r_l, r_r in raw_radii)
        fixed_l = fixed_r = float(fixed_shared)
    else:
        fixed_l = float(max(r_l for r_l, _ in raw_radii))
        fixed_r = float(max(r_r for _, r_r in raw_radii))
    print(
        "Animation viewport lock: "
        f"left={fixed_l:.6g} m, right={fixed_r:.6g} m"
    )

    def animate(i: int):
        # Use fixed viewport(s) across animation so true expansion/contraction is visible.
        # Per-frame autoscaling can keep normalized positions nearly constant and look static.
        r_l, r_r = fixed_l, fixed_r
        if shared_viewport:
            r_l = r_r = fixed_l
        _draw_panel(axes[0], left, left_snaps[i], r_l)
        _draw_panel(axes[1], right, right_snaps[i], r_r)
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
    out_mp4 = parent_dir / "galaxy_compare.mp4"
    out_gif = parent_dir / "galaxy_compare.gif"
    try:
        anim.save(str(out_mp4), writer="ffmpeg", fps=20, dpi=100)
    except Exception:
        anim.save(str(out_gif), writer="pillow", fps=15, dpi=100)
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
    ap.add_argument(
        "--shared-viewport",
        action="store_true",
        help="Use the same axis half-range on both panels (max of validated left/right radii). "
        "Same physical scale on both sides, but if one run is much more extended, the compact side "
        "can look empty and stellar motion can be invisible.",
    )
    args = ap.parse_args()
    render_compare(
        args.compare_parent_dir.resolve(),
        args.no_animation,
        args.render_overlay_mode,
        shared_viewport=args.shared_viewport,
    )


if __name__ == "__main__":
    main()

