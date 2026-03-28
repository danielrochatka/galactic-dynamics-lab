"""
Rendering: static scatter plots and animation.
Top-down 2D view of the galaxy.
"""

import subprocess
import time
from collections.abc import Callable
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


RenderRadius = Union[float, Callable[[np.ndarray], float]]


def _resolve_render_radius(render_radius: RenderRadius, positions: np.ndarray) -> float:
    if callable(render_radius):
        return float(render_radius(positions))
    return float(render_radius)


def _setup_axes(ax: plt.Axes, render_radius: float) -> None:
    """Configure axes: equal aspect, fixed range."""
    ax.set_aspect("equal")
    ax.set_xlim(-render_radius, render_radius)
    ax.set_ylim(-render_radius, render_radius)
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def scatter_frame(
    ax: plt.Axes,
    positions: np.ndarray,
    bh_position: tuple[float, float] = (0.0, 0.0),
    render_radius: RenderRadius = 150.0,
    star_size: float = 2.0,
    bh_size: float = 80.0,
) -> None:
    """Draw one frame: black hole + stars."""
    ax.clear()
    r = _resolve_render_radius(render_radius, positions)
    _setup_axes(ax, r)

    # Stars
    ax.scatter(
        positions[:, 0],
        positions[:, 1],
        s=star_size,
        c="white",
        edgecolors="none",
        alpha=0.8,
        rasterized=True,
    )
    # Black hole
    ax.scatter(
        [bh_position[0]],
        [bh_position[1]],
        s=bh_size,
        c="yellow",
        marker="*",
        edgecolors="orange",
        linewidths=0.5,
        zorder=10,
    )


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
    _setup_axes(ax, render_radius)

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
    ax.set_title(title, color="white", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, facecolor="black", edgecolor="none")
    plt.close(fig)


def save_static_plot(
    positions: np.ndarray,
    output_path: Path,
    title: str = "Galaxy",
    render_radius: RenderRadius = 150.0,
) -> None:
    """Save a single static scatter plot."""
    fig, ax = plt.subplots(figsize=(10, 10), facecolor="black")
    ax.set_facecolor("black")
    ax.tick_params(colors="gray")
    ax.spines["bottom"].set_color("gray")
    ax.spines["top"].set_color("gray")
    ax.spines["left"].set_color("gray")
    ax.spines["right"].set_color("gray")

    scatter_frame(ax, positions, render_radius=render_radius)
    ax.set_title(title, color="white", fontsize=14)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, facecolor="black", edgecolor="none")
    plt.close(fig)


def create_animation(
    snapshots: list,
    output_path: Path,
    render_radius: RenderRadius = 150.0,
    interval: int = 50,
    progress_interval: int = 50,
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
        snap = snapshots[frame_idx]
        scatter_frame(ax, snap.positions, render_radius=render_radius)
        ax.set_title(f"Step {snap.step}  |  t = {snap.time:.1f}", color="white")
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
