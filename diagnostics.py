"""
Diagnostics from simulation snapshots: time series and plots.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# SI gravitational constant (matches cpp_sim/init_conditions.cpp)
G_SI = 6.6743e-11


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


def compute_two_body_pair_diagnostics(
    snapshots: list,
    masses: np.ndarray,
    mode: str,
    bh_mass: float,
) -> dict[str, np.ndarray | str]:
    """
    Pair- and COM-based time series for earth_moon_benchmark (n=2) or bh_orbit_validation (n=1 + BH at lab origin).

    Frames:
    - Lab frame: snapshot positions/velocities as stored (origin = simulation origin; BH treated at (0,0) for bh mode).
    - Relative vectors: r_rel = r1 - r0 (two masses) or star position vs fixed BH at origin (n=1).
    - COM position R_com is computed in the lab frame from particle masses (and M_bh at origin for bh mode).
    - relative_Lz: z-component of total angular momentum about COM (same in lab as in COM frame for L_tot).
    - newtonian_specific_energy: point-mass Newtonian formula in SI; for TPF runs this is a *post-hoc Newtonian
      equivalent*, not a conserved quantity of the TPF integrator.

    Returns dict with numpy arrays plus key "variant" in {"two_mass", "star_bh"}.
    """
    if not snapshots:
        raise ValueError("compute_two_body_pair_diagnostics: empty snapshots")

    n_snap = len(snapshots)
    time = np.array([s.time for s in snapshots], dtype=np.float64)
    pair_sep = np.zeros(n_snap)
    pair_rel_speed = np.zeros(n_snap)
    com_x = np.zeros(n_snap)
    com_y = np.zeros(n_snap)
    com_radius = np.zeros(n_snap)
    rel_Lz = np.zeros(n_snap)
    E_spec = np.full(n_snap, np.nan)

    if mode == "earth_moon_benchmark":
        if masses.shape[0] != 2:
            raise ValueError("earth_moon_benchmark expects 2 particles in snapshots/masses")
        m0, m1 = float(masses[0]), float(masses[1])
        Mtot = m0 + m1
        mu = (m0 * m1) / Mtot if Mtot > 0 else 0.0
        variant = "two_mass"
        for i, snap in enumerate(snapshots):
            pos = snap.positions
            vel = snap.velocities
            if pos.shape[0] != 2:
                raise ValueError("earth_moon_benchmark: expected 2 positions per snapshot")
            r0 = pos[0].astype(np.float64)
            r1 = pos[1].astype(np.float64)
            v0 = vel[0].astype(np.float64)
            v1 = vel[1].astype(np.float64)
            r_rel = r1 - r0
            v_rel = v1 - v0
            sep = float(np.linalg.norm(r_rel))
            pair_sep[i] = sep
            pair_rel_speed[i] = float(np.linalg.norm(v_rel))
            R = (m0 * r0 + m1 * r1) / Mtot
            com_x[i], com_y[i] = R[0], R[1]
            com_radius[i] = float(np.linalg.norm(R))
            rel_Lz[i] = mu * (r_rel[0] * v_rel[1] - r_rel[1] * v_rel[0])
            if sep > 1e-30 and mu > 0:
                E_spec[i] = 0.5 * pair_rel_speed[i] ** 2 - G_SI * (m0 + m1) / sep
    elif mode == "bh_orbit_validation":
        Mbh = float(bh_mass)
        if masses.shape[0] != 1:
            raise ValueError("bh_orbit_validation expects 1 particle in snapshots/masses")
        ms = float(masses[0])
        Mtot = ms + Mbh
        variant = "star_bh"
        for i, snap in enumerate(snapshots):
            pos = snap.positions
            vel = snap.velocities
            r_s = pos[0].astype(np.float64)
            v_s = vel[0].astype(np.float64)
            r_bh = np.zeros(2, dtype=np.float64)
            v_bh = np.zeros(2, dtype=np.float64)
            pair_sep[i] = float(np.linalg.norm(r_s - r_bh))
            pair_rel_speed[i] = float(np.linalg.norm(v_s - v_bh))
            R = (ms * r_s + Mbh * r_bh) / Mtot if Mtot > 0 else np.zeros(2)
            com_x[i], com_y[i] = R[0], R[1]
            com_radius[i] = float(np.linalg.norm(R))
            V = (ms * v_s + Mbh * v_bh) / Mtot if Mtot > 0 else np.zeros(2)
            r_s_rel = r_s - R
            v_s_rel = v_s - V
            r_bh_rel = r_bh - R
            v_bh_rel = v_bh - V
            Lz_star = ms * (r_s_rel[0] * v_s_rel[1] - r_s_rel[1] * v_s_rel[0])
            Lz_bh = Mbh * (r_bh_rel[0] * v_bh_rel[1] - r_bh_rel[1] * v_bh_rel[0])
            rel_Lz[i] = Lz_star + Lz_bh
            r_mag = pair_sep[i]
            if r_mag > 1e-30 and Mbh > 0:
                E_spec[i] = 0.5 * float(np.dot(v_s, v_s)) - G_SI * Mbh / r_mag
    else:
        raise ValueError(f"unknown two-body diagnostics mode: {mode}")

    return {
        "time": time,
        "pair_separation": pair_sep,
        "pair_relative_speed": pair_rel_speed,
        "center_of_mass_x": com_x,
        "center_of_mass_y": com_y,
        "center_of_mass_radius": com_radius,
        "relative_angular_momentum_z": rel_Lz,
        "newtonian_specific_energy": E_spec,
        "variant": variant,
    }


def write_two_body_diagnostics_readme(output_dir: Path, mode: str, physics_package: str) -> None:
    """Short text file: frames, formulas, Newtonian vs TPF caveat."""
    tpf = physics_package.strip() == "TPFCore"
    lines = [
        "Two-body primary diagnostics (pair / COM frame)",
        "==============================================",
        f"simulation_mode (postprocess): {mode}",
        f"physics_package: {physics_package}",
        "",
        "Lab frame: snapshot x,y,vx,vy are in the simulation inertial frame (origin at code origin).",
        "For bh_orbit_validation the central mass is fixed at the origin and is NOT a CSV particle row.",
        "",
        "Columns in diagnostic_two_body_timeseries.csv:",
        "  pair_separation: |r1 - r0| (two-mass) or |r_star| (star vs origin/BH).",
        "  pair_relative_speed: |v1 - v0| or |v_star| (BH velocity taken as 0 in lab).",
        "  center_of_mass_x, center_of_mass_y, center_of_mass_radius: R_com from lab-frame masses.",
        "  relative_angular_momentum_z: L_z about COM (sum of r'_i x (m_i v'_i) in 2D).",
        "  newtonian_specific_energy_J_per_kg:",
        "    earth_moon: 0.5*|v_rel|^2 - G*(m0+m1)/|r_rel|  (per unit reduced mass; SI).",
        "    bh_orbit:   0.5*|v_star|^2 - G*M_bh/|r_star|  (per unit star mass; SI).",
        "    Uses point-mass Newtonian U; simulator softening is NOT subtracted here.",
    ]
    if tpf:
        lines += [
            "",
            "TPFCore: the plotted energy is a Newtonian post-processing equivalent only.",
            "It is NOT claimed equal to a TPF Hamiltonian or conserved TPF scalar.",
        ]
    (output_dir / "two_body_diagnostics_README.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_two_body_timeseries_csv(diag: dict[str, np.ndarray | str], output_dir: Path) -> None:
    """Write numeric time series for spreadsheets."""
    t = diag["time"]
    rows = np.column_stack(
        [
            t,
            diag["pair_separation"],
            diag["pair_relative_speed"],
            diag["center_of_mass_x"],
            diag["center_of_mass_y"],
            diag["center_of_mass_radius"],
            diag["relative_angular_momentum_z"],
            diag["newtonian_specific_energy"],
        ]
    )
    header = (
        "time,pair_separation_m,pair_relative_speed_m_s,center_of_mass_x_m,center_of_mass_y_m,"
        "center_of_mass_radius_m,relative_angular_momentum_z_kg_m2_s,newtonian_specific_energy_J_per_kg"
    )
    out = output_dir / "diagnostic_two_body_timeseries.csv"
    np.savetxt(out, rows, delimiter=",", header=header, comments="")


def _save_plot_with_note(
    output_path: Path,
    time: np.ndarray,
    y: np.ndarray,
    ylabel: str,
    title: str,
    footnote: str | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(time, y, color="steelblue", linewidth=1)
    ax.scatter(time, y, s=8, c="steelblue", zorder=5)
    ax.set_xlabel("Time (simulation units)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    if footnote:
        fig.text(0.5, 0.02, footnote, ha="center", fontsize=8)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18 if footnote else 0.12)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_two_body_pair_diagnostics(
    diag: dict[str, np.ndarray | str],
    output_dir: Path,
    physics_package: str,
) -> None:
    """Primary PNG set for two-body modes."""
    tpf_note = (
        "Newtonian equivalent only for TPFCore — not a TPF conserved scalar."
        if physics_package.strip() == "TPFCore"
        else None
    )
    t = diag["time"]
    _save_plot(
        output_dir / "diagnostic_pair_separation.png",
        t,
        diag["pair_separation"],
        "Pair separation (m)",
        "Primary: pair separation |Δr| vs time (lab-frame positions)",
    )
    _save_plot(
        output_dir / "diagnostic_pair_relative_speed.png",
        t,
        diag["pair_relative_speed"],
        "Pair relative speed (m/s)",
        "Primary: |Δv| vs time (lab-frame velocities)",
    )
    _save_plot(
        output_dir / "diagnostic_com_radius.png",
        t,
        diag["center_of_mass_radius"],
        "|R_COM| from origin (m)",
        "COM radius in lab frame (drift here is not 'escape' of the pair)",
    )
    _save_plot(
        output_dir / "diagnostic_relative_angular_momentum_z.png",
        t,
        diag["relative_angular_momentum_z"],
        "L_z about COM (SI)",
        "Primary: angular momentum about center of mass (z)",
    )
    E = diag["newtonian_specific_energy"]
    if np.any(np.isfinite(E)):
        _save_plot_with_note(
            output_dir / "diagnostic_relative_energy.png",
            t,
            E,
            "Newtonian specific energy (J/kg)",
            "Primary: Newtonian specific mechanical energy (see README)",
            footnote=tpf_note,
        )


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
    *,
    lab_frame_secondary: bool = False,
) -> None:
    """
    Save one PNG per diagnostic in output_dir (one point per snapshot, time = simulation time).

    If lab_frame_secondary is True, titles state that these metrics are origin-radial / lab-frame
    summaries (for N-body disk-style interpretation), not primary two-body pair diagnostics.
    """
    prefix = (
        "Secondary (lab / origin radial): "
        if lab_frame_secondary
        else ""
    )
    t = diagnostics["time"]
    _save_plot(
        output_dir / "diagnostic_median_radius.png",
        t, diagnostics["median_r"],
        "Median radius from origin (m)",
        prefix + "Median |r| from origin vs time",
    )
    _save_plot(
        output_dir / "diagnostic_mean_radius.png",
        t, diagnostics["mean_r"],
        "Mean radius from origin (m)",
        prefix + "Mean |r| from origin vs time",
    )
    _save_plot(
        output_dir / "diagnostic_std_radius.png",
        t, diagnostics["std_r"],
        "Std of |r| from origin (m)",
        prefix + "Std of star radii from origin vs time",
    )
    _save_plot(
        output_dir / "diagnostic_max_radius.png",
        t, diagnostics["max_r"],
        "Max radius from origin (m)",
        prefix + "Max |r| from origin vs time",
    )
    _save_plot(
        output_dir / "diagnostic_frac_vr_positive.png",
        t, diagnostics["frac_vr_pos"],
        "Fraction with v_r > 0 (origin radial)",
        prefix + "Fraction of bodies with v_r > 0 vs time",
    )
    _save_plot(
        output_dir / "diagnostic_frac_vr_negative.png",
        t, diagnostics["frac_vr_neg"],
        "Fraction with v_r < 0 (origin radial)",
        prefix + "Fraction of bodies with v_r < 0 vs time",
    )
    _save_plot(
        output_dir / "diagnostic_frac_beyond_cutoff.png",
        t, diagnostics["frac_beyond_cutoff"],
        f"Fraction beyond r = {cutoff_radius} (origin)",
        prefix + f"Fraction of bodies beyond r = {cutoff_radius} vs time",
    )
    _save_plot(
        output_dir / "diagnostic_angular_momentum_z.png",
        t, diagnostics["L_z"],
        "Total L_z about origin (SI)",
        prefix + "Total L_z about simulation origin vs time",
    )
