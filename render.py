"""
Rendering: static scatter plots and animation.
Top-down 2D view of the galaxy.
"""

import subprocess
import time
from pathlib import Path
from typing import Any, Callable, Optional, Union

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

from text_layout import set_fitted_title

from display_units import (
    SpatialDisplay,
    apply_suppress_tick_offset,
    format_animation_time_caption,
)
from framing import SquareViewport

# float: centered origin; SquareViewport: explicit center + half-axis; callable: from current positions
ViewportSpec = Union[float, SquareViewport, Callable[[np.ndarray, Optional[np.ndarray]], SquareViewport]]


def _resolve_square_viewport(
    spec: ViewportSpec,
    positions: np.ndarray,
    velocities: Optional[np.ndarray] = None,
) -> SquareViewport:
    if isinstance(spec, SquareViewport):
        return spec
    if callable(spec):
        if (
            velocities is not None
            and velocities.size > 0
            and len(velocities) == len(positions)
        ):
            return spec(positions, velocities)
        return spec(positions, None)
    return SquareViewport.centered_origin(float(spec))


def _setup_axes_square(
    ax: plt.Axes,
    center_x: float,
    center_y: float,
    display_half_extent: float,
    xy_unit: str,
) -> None:
    """Equal-aspect square |x-cx|,|y-cy| <= half (display length units)."""
    ax.set_aspect("equal")
    ax.set_xlim(center_x - display_half_extent, center_x + display_half_extent)
    ax.set_ylim(center_y - display_half_extent, center_y + display_half_extent)
    ax.set_xlabel(f"x ({xy_unit})")
    ax.set_ylabel(f"y ({xy_unit})")


def scatter_frame(
    ax: plt.Axes,
    positions: np.ndarray,
    bh_position: tuple[float, float] = (0.0, 0.0),
    render_radius: ViewportSpec = 150.0,
    velocities: Optional[np.ndarray] = None,
    star_size: float = 2.0,
    bh_size: float = 80.0,
    spatial_display: Optional[SpatialDisplay] = None,
) -> None:
    """Draw one frame: black hole + stars. Positions are SI (m); axes use display units when set."""
    ax.clear()
    vp = _resolve_square_viewport(render_radius, positions, velocities)
    sd = spatial_display if spatial_display is not None else SpatialDisplay(1.0, "m")
    f = sd.factor
    cx = vp.center_x * f
    cy = vp.center_y * f
    h_disp = vp.half_axis * f
    _setup_axes_square(ax, cx, cy, h_disp, sd.unit)

    # Stars
    ax.scatter(
        positions[:, 0] * f,
        positions[:, 1] * f,
        s=star_size,
        c="white",
        edgecolors="none",
        alpha=0.8,
        rasterized=True,
    )
    # Black hole
    ax.scatter(
        [bh_position[0] * f],
        [bh_position[1] * f],
        s=bh_size,
        c="yellow",
        marker="*",
        edgecolors="orange",
        linewidths=0.5,
        zorder=10,
    )
    apply_suppress_tick_offset(ax)


def _radial_velocity(positions: np.ndarray, velocities: np.ndarray) -> np.ndarray:
    """Radial velocity v_r = (x*vx + y*vy) / r relative to origin; safe for small r."""
    r = np.sqrt(np.sum(positions**2, axis=1))
    r_safe = np.maximum(r, 1e-10)
    return (positions[:, 0] * velocities[:, 0] + positions[:, 1] * velocities[:, 1]) / r_safe


def save_radial_velocity_plot(
    positions: np.ndarray,
    velocities: np.ndarray,
    output_path: Path,
    title: str = "Radial velocity",
    render_radius: float = 150.0,
) -> None:
    """Save a scatter plot with stars colored by radial velocity (sign and magnitude)."""
    vr = _radial_velocity(positions, velocities)
    fig, ax = plt.subplots(figsize=(10, 10), facecolor="black")
    ax.set_facecolor("black")
    ax.tick_params(colors="gray")
    for spine in ax.spines.values():
        spine.set_color("gray")
    sd = SpatialDisplay(1.0, "m")
    h = float(render_radius) * sd.factor
    _setup_axes_square(ax, 0.0, 0.0, h, sd.unit)

    # Diverging colormap: blue = inward (v_r < 0), red = outward (v_r > 0)
    v_abs_max = float(np.nanmax(np.abs(vr))) if np.any(np.isfinite(vr)) else 1.0
    if v_abs_max <= 0:
        v_abs_max = 1.0
    scatter = ax.scatter(
        positions[:, 0],
        positions[:, 1],
        c=vr,
        s=2.0,
        cmap="RdBu_r",
        vmin=-v_abs_max,
        vmax=v_abs_max,
        edgecolors="none",
        alpha=0.9,
        rasterized=True,
    )
    cbar = fig.colorbar(scatter, ax=ax, label="Radial velocity v_r")
    cbar.ax.yaxis.set_tick_params(color="gray")
    cbar.ax.yaxis.label.set_color("gray")
    plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color="gray")

    ax.scatter(
        [0.0], [0.0],
        s=80, c="yellow", marker="*", edgecolors="orange", linewidths=0.5, zorder=10,
    )
    set_fitted_title(ax, title, color="white", fontsize=14, min_fontsize=6)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, facecolor="black", edgecolor="none")
    plt.close(fig)


def save_static_plot(
    positions: np.ndarray,
    output_path: Path,
    title: str = "Galaxy",
    render_radius: ViewportSpec = 150.0,
    velocities: Optional[np.ndarray] = None,
    *,
    overlay_mode: str = "none",
    overlay_spec: Optional[dict[str, Any]] = None,
    run_info: Optional[dict[str, Any]] = None,
    overlay_step: int = 0,
    overlay_time: float = 0.0,
    spatial_display: Optional[SpatialDisplay] = None,
    unit_reference_text: str | None = None,
) -> None:
    """Save a single static scatter plot. Positions remain SI; axes use display units when provided."""
    fig, ax = plt.subplots(figsize=(10, 10), facecolor="black")
    ax.set_facecolor("black")
    ax.tick_params(colors="gray")
    ax.spines["bottom"].set_color("gray")
    ax.spines["top"].set_color("gray")
    ax.spines["left"].set_color("gray")
    ax.spines["right"].set_color("gray")

    scatter_frame(
        ax,
        positions,
        render_radius=render_radius,
        velocities=velocities,
        spatial_display=spatial_display,
    )
    set_fitted_title(ax, title, color="white", fontsize=14, min_fontsize=6)
    if unit_reference_text:
        ax.text(
            0.99,
            0.01,
            unit_reference_text,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7,
            color="white",
            bbox=dict(boxstyle="round,pad=0.25", fc="black", ec="white", alpha=0.5),
        )
    if (
        overlay_mode != "none"
        and overlay_spec is not None
        and run_info is not None
    ):
        from render_overlay import draw_galaxy_render_overlay

        draw_galaxy_render_overlay(
            ax,
            overlay_mode,
            overlay_spec,
            run_info=run_info,
            step=overlay_step,
            time_s=overlay_time,
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, facecolor="black", edgecolor="none")
    plt.close(fig)


def create_animation(
    snapshots: list,
    output_path: Path,
    render_radius: ViewportSpec = 150.0,
    interval: int = 50,
    progress_interval: int = 50,
    *,
    overlay_mode: str = "none",
    overlay_spec: Optional[dict[str, Any]] = None,
    run_info: Optional[dict[str, Any]] = None,
    spatial_display: Optional[SpatialDisplay] = None,
    simulation_mode: str = "galaxy",
    preferred_time_unit: str = "auto",
    unit_reference_text: str | None = None,
    mutable_frame_index: Optional[dict[str, int]] = None,
) -> bool:
    """
    Create animation (MP4 or GIF) from snapshots.
    Tries ffmpeg for MP4 first; falls back to GIF if unavailable.
    Prints encoding progress every progress_interval frames.
    Returns True if successful.
    """
    if not snapshots:
        return False

    n_frames = len(snapshots)
    fig, ax = plt.subplots(figsize=(10, 10), facecolor="black")
    ax.set_facecolor("black")
    ax.tick_params(colors="gray")
    ax.spines["bottom"].set_color("gray")
    ax.spines["top"].set_color("gray")
    ax.spines["left"].set_color("gray")
    ax.spines["right"].set_color("gray")

    def animate(frame_idx: int) -> list:
        if mutable_frame_index is not None:
            mutable_frame_index["i"] = int(frame_idx)
        snap = snapshots[frame_idx]
        vel = getattr(snap, "velocities", None)
        scatter_frame(
            ax,
            snap.positions,
            render_radius=render_radius,
            velocities=vel,
            spatial_display=spatial_display,
        )
        tc = format_animation_time_caption(
            float(snap.time), simulation_mode, preferred_time_unit=preferred_time_unit
        )
        ax.set_title(f"Step {snap.step}  |  {tc}", color="white")
        if unit_reference_text:
            ax.text(
                0.99,
                0.01,
                unit_reference_text,
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=7,
                color="white",
                bbox=dict(boxstyle="round,pad=0.25", fc="black", ec="white", alpha=0.5),
            )
        if (
            overlay_mode != "none"
            and overlay_spec is not None
            and run_info is not None
        ):
            from render_overlay import draw_galaxy_render_overlay

            draw_galaxy_render_overlay(
                ax,
                overlay_mode,
                overlay_spec,
                run_info=run_info,
                step=int(snap.step),
                time_s=float(snap.time),
            )
        if progress_interval and (frame_idx % progress_interval == 0 or frame_idx == n_frames - 1):
            pct = 100 * (frame_idx + 1) / n_frames
            print(f"    Frame {frame_idx + 1}/{n_frames} ({pct:.0f}%)")
        return []

    anim = animation.FuncAnimation(
        fig,
        animate,
        frames=n_frames,
        interval=interval,
        blit=False,
    )

    mp4_path = output_path.with_suffix(".mp4")
    gif_path = output_path.with_suffix(".gif")

    print(f"  Total frames: {n_frames}")
    print("  Encoding started...")
    t0 = time.perf_counter()

    try:
        anim.save(str(mp4_path), writer="ffmpeg", fps=20, dpi=100)
        elapsed = time.perf_counter() - t0
        plt.close(fig)
        print(f"  Encoding finished in {elapsed:.1f} s")
        print(f"  Saved: {mp4_path.resolve()}")
        return True
    except Exception:
        pass

    try:
        anim.save(str(gif_path), writer="pillow", fps=15, dpi=100)
        elapsed = time.perf_counter() - t0
        plt.close(fig)
        print(f"  Encoding finished in {elapsed:.1f} s")
        print(f"  Saved: {gif_path.resolve()}")
        return True
    except Exception as e:
        plt.close(fig)
        print(f"  Animation save failed: {e}")
        return False


def has_ffmpeg() -> bool:
    """Check if ffmpeg is available."""
    try:
        subprocess.run(
            ["ffmpeg", "-version"],
            capture_output=True,
            check=False,
        )
        return True
    except (FileNotFoundError, subprocess.SubprocessError):
        return False
