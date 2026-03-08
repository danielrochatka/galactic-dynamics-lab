"""
Diagnostics from simulation snapshots: time series and plots.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _radii(positions: np.ndarray) -> np.ndarray:
    """Star radii from origin."""
    return np.sqrt(np.sum(positions**2, axis=1))


def _radial_velocities(positions: np.ndarray, velocities: np.ndarray) -> np.ndarray:
    """Radial velocity v_r = (x*vx + y*vy) / r; safe for small r."""
    r = _radii(positions)
    r_safe = np.maximum(r, 1e-10)
    return (positions[:, 0] * velocities[:, 0] + positions[:, 1] * velocities[:, 1]) / r_safe


def _angular_momentum_z(positions: np.ndarray, velocities: np.ndarray, masses: np.ndarray) -> float:
    """Total z-component of angular momentum: sum over m * (x*vy - y*vx)."""
    L_z_per_star = masses * (positions[:, 0] * velocities[:, 1] - positions[:, 1] * velocities[:, 0])
    return float(np.sum(L_z_per_star))


def compute_diagnostics(
    snapshots: list,
    masses: np.ndarray,
    cutoff_radius: float,
) -> dict[str, np.ndarray]:
    """
    Compute diagnostic time series from every snapshot.
    Time is simulation time (snapshot.time). Returns dict with keys:
    time, median_r, mean_r, std_r, max_r, frac_vr_pos, frac_vr_neg, frac_beyond_cutoff, L_z.
    """
    n_snap = len(snapshots)
    time = np.array([s.time for s in snapshots])
    median_r = np.zeros(n_snap)
    mean_r = np.zeros(n_snap)
    std_r = np.zeros(n_snap)
    max_r = np.zeros(n_snap)
    frac_vr_pos = np.zeros(n_snap)
    frac_vr_neg = np.zeros(n_snap)
    frac_beyond_cutoff = np.zeros(n_snap)
    L_z = np.zeros(n_snap)

    for i, snap in enumerate(snapshots):
        r = _radii(snap.positions)
        median_r[i] = np.median(r)
        mean_r[i] = np.mean(r)
        std_r[i] = np.std(r)
        max_r[i] = np.max(r)
        vr = _radial_velocities(snap.positions, snap.velocities)
        n = len(r)
        frac_vr_pos[i] = np.sum(vr > 0) / n
        frac_vr_neg[i] = np.sum(vr < 0) / n
        frac_beyond_cutoff[i] = np.sum(r > cutoff_radius) / n
        L_z[i] = _angular_momentum_z(snap.positions, snap.velocities, masses)

    return {
        "time": time,
        "median_r": median_r,
        "mean_r": mean_r,
        "std_r": std_r,
        "max_r": max_r,
        "frac_vr_pos": frac_vr_pos,
        "frac_vr_neg": frac_vr_neg,
        "frac_beyond_cutoff": frac_beyond_cutoff,
        "L_z": L_z,
    }


def _save_plot(
    output_path: Path,
    time: np.ndarray,
    y: np.ndarray,
    ylabel: str,
    title: str,
) -> None:
    """Time-series plot: x-axis is simulation time, one point per snapshot."""
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(time, y, color="steelblue", linewidth=1)
    ax.scatter(time, y, s=8, c="steelblue", zorder=5)  # show every snapshot point
    ax.set_xlabel("Time (simulation units)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_and_save_all(
    diagnostics: dict[str, np.ndarray],
    output_dir: Path,
    cutoff_radius: float,
) -> None:
    """Save one PNG per diagnostic in output_dir (one point per snapshot, time = simulation time)."""
    t = diagnostics["time"]
    _save_plot(
        output_dir / "diagnostic_median_radius.png",
        t, diagnostics["median_r"],
        "Median radius",
        "Median radius vs time",
    )
    _save_plot(
        output_dir / "diagnostic_mean_radius.png",
        t, diagnostics["mean_r"],
        "Mean radius",
        "Mean radius vs time",
    )
    _save_plot(
        output_dir / "diagnostic_std_radius.png",
        t, diagnostics["std_r"],
        "Std radius",
        "Radius standard deviation vs time",
    )
    _save_plot(
        output_dir / "diagnostic_max_radius.png",
        t, diagnostics["max_r"],
        "Max radius",
        "Max radius vs time",
    )
    _save_plot(
        output_dir / "diagnostic_frac_vr_positive.png",
        t, diagnostics["frac_vr_pos"],
        "Fraction with v_r > 0",
        "Fraction of stars with positive radial velocity vs time",
    )
    _save_plot(
        output_dir / "diagnostic_frac_vr_negative.png",
        t, diagnostics["frac_vr_neg"],
        "Fraction with v_r < 0",
        "Fraction of stars with negative radial velocity vs time",
    )
    _save_plot(
        output_dir / "diagnostic_frac_beyond_cutoff.png",
        t, diagnostics["frac_beyond_cutoff"],
        f"Fraction beyond r = {cutoff_radius}",
        f"Fraction of stars beyond r = {cutoff_radius} vs time",
    )
    _save_plot(
        output_dir / "diagnostic_angular_momentum_z.png",
        t, diagnostics["L_z"],
        "Total L_z",
        "Total angular momentum (z) vs time",
    )
