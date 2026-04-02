"""
Diagnostics from simulation snapshots: time series and plots.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from text_layout import add_fitted_footer, set_fitted_title

from display_units import (
    SeriesDisplay,
    SpatialDisplay,
    apply_suppress_tick_offset,
    scale_angular_momentum_display,
    scale_angular_momentum_origin_display,
    scale_energy_display,
    series_display_for_galaxy_diagnostics,
    series_display_for_two_body,
    series_display_generic_validation,
    spatial_display_for_xy_plot,
)

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
    if mode == "bh_orbit_validation":
        lines += [
            "",
            "bh_orbit_validation extras (plot_cpp_run.py): bh_orbit_trajectory_xy.png, bh_orbit_trajectory_xy_zoom.png,",
            "bh_orbit_separation_extrema.png — experimental visuals; footnotes state Newtonian baseline vs TPFCore legacy_readout,",
            "VDSG coupling from run_info, and that this is not paper correspondence mode and not direct_tpf.",
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
    time_si: np.ndarray,
    y: np.ndarray,
    ylabel: str,
    title: str,
    *,
    series: SeriesDisplay,
    footnote: str | None = None,
) -> None:
    t_plot = np.asarray(time_si, dtype=np.float64) * series.time_factor
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(t_plot, y, color="steelblue", linewidth=1)
    ax.scatter(t_plot, y, s=8, c="steelblue", zorder=5)
    ax.set_xlabel(f"Time ({series.time_unit})")
    ax.set_ylabel(ylabel)
    set_fitted_title(ax, title, fontsize=12, min_fontsize=6)
    ax.grid(True, alpha=0.3)
    if series.suppress_tick_offset:
        apply_suppress_tick_offset(ax)
    if footnote:
        add_fitted_footer(fig, footnote, fontsize=8, min_fontsize=5)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18 if footnote else 0.12)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def _save_plot(
    output_path: Path,
    time_si: np.ndarray,
    y_si: np.ndarray,
    ylabel: str,
    title: str,
    *,
    series: SeriesDisplay,
    y_scale: float,
    footnote: str | None = None,
) -> None:
    """Time-series plot: SI arrays converted to display units via series + y_scale."""
    t_plot = np.asarray(time_si, dtype=np.float64) * series.time_factor
    y_plot = np.asarray(y_si, dtype=np.float64) * float(y_scale)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_plot, y_plot, color="steelblue", linewidth=1)
    ax.scatter(t_plot, y_plot, s=8, c="steelblue", zorder=5)  # show every snapshot point
    ax.set_xlabel(f"Time ({series.time_unit})")
    ax.set_ylabel(ylabel)
    set_fitted_title(ax, title, fontsize=12, min_fontsize=6)
    ax.grid(True, alpha=0.3)
    if series.suppress_tick_offset:
        apply_suppress_tick_offset(ax)
    if footnote:
        add_fitted_footer(fig, footnote, fontsize=8, min_fontsize=5)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18 if footnote else 0.12)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def _save_plot_fraction(
    output_path: Path,
    time_si: np.ndarray,
    y: np.ndarray,
    ylabel: str,
    title: str,
    *,
    series: SeriesDisplay,
) -> None:
    """Dimensionless fraction vs time (y in [0,1]); only time axis uses display units."""
    t_plot = np.asarray(time_si, dtype=np.float64) * series.time_factor
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_plot, y, color="steelblue", linewidth=1)
    ax.scatter(t_plot, y, s=8, c="steelblue", zorder=5)
    ax.set_xlabel(f"Time ({series.time_unit})")
    ax.set_ylabel(ylabel)
    set_fitted_title(ax, title, fontsize=12, min_fontsize=6)
    ax.grid(True, alpha=0.3)
    if series.suppress_tick_offset:
        apply_suppress_tick_offset(ax)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.12)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_two_body_pair_diagnostics(
    diag: dict[str, np.ndarray | str],
    output_dir: Path,
    physics_package: str,
    *,
    title_suffix: str = "",
    footnote: str | None = None,
    context_label: str = "",
    series: SeriesDisplay,
) -> None:
    """Primary PNG set for two-body modes (display units from `series`; data arrays remain SI in memory)."""
    suf = (" " + title_suffix.strip()) if title_suffix.strip() else ""
    tpf_note = (
        "Newtonian equivalent only for TPFCore — not a TPF conserved scalar."
        if physics_package.strip() == "TPFCore"
        else None
    )
    t = diag["time"]
    ctx = (context_label.strip() + " — ") if context_label.strip() else ""
    disp_note = "Axes use display units; CSV remains SI."
    foot_disp = (footnote + " | " + disp_note) if footnote else disp_note
    _save_plot(
        output_dir / "diagnostic_pair_separation.png",
        t,
        diag["pair_separation"],
        f"Pair separation ({series.distance_unit})",
        ctx + "Primary: pair separation |Δr| vs time (lab-frame positions)" + suf,
        series=series,
        y_scale=series.distance_factor,
        footnote=foot_disp,
    )
    _save_plot(
        output_dir / "diagnostic_pair_relative_speed.png",
        t,
        diag["pair_relative_speed"],
        f"Pair relative speed ({series.speed_unit})",
        ctx + "Primary: |Δv| vs time (lab-frame velocities)" + suf,
        series=series,
        y_scale=series.speed_factor,
        footnote=foot_disp,
    )
    _save_plot(
        output_dir / "diagnostic_com_radius.png",
        t,
        diag["center_of_mass_radius"],
        f"|R_COM| from origin ({series.distance_unit})",
        ctx + "Primary: COM radius in lab frame (drift here is not 'escape' of the pair)" + suf,
        series=series,
        y_scale=series.distance_factor,
        footnote=foot_disp,
    )
    Lz_plot, Lz_ylab = scale_angular_momentum_display(np.asarray(diag["relative_angular_momentum_z"]))
    _save_plot_with_note(
        output_dir / "diagnostic_relative_angular_momentum_z.png",
        t,
        Lz_plot,
        Lz_ylab,
        ctx + "Primary: angular momentum about center of mass (z)" + suf,
        series=series,
        footnote=foot_disp,
    )
    E = diag["newtonian_specific_energy"]
    if np.any(np.isfinite(E)):
        e_foot = tpf_note
        if footnote:
            e_foot = (tpf_note + " | " + footnote) if tpf_note else footnote
        e_foot = (e_foot + " | " + disp_note) if e_foot else disp_note
        E_plot, E_ylab = scale_energy_display(np.asarray(E))
        _save_plot_with_note(
            output_dir / "diagnostic_relative_energy.png",
            t,
            E_plot,
            E_ylab,
            ctx + "Primary: Newtonian specific mechanical energy (see README)" + suf,
            series=series,
            footnote=e_foot,
        )


def bh_orbit_validation_plot_footnote(physics_package: str, tpf_vdsg_coupling: float) -> str:
    """Honest labels for experimental bh_orbit_validation postprocess (not paper / not direct_tpf)."""
    pkg = physics_package.strip()
    vdsg = float(tpf_vdsg_coupling)
    vdsg_part = (
        "VDSG off (tpf_vdsg_coupling = 0) — baseline for legacy_readout compare."
        if vdsg == 0.0
        else f"tpf_vdsg_coupling = {vdsg:g} (nonzero; not the recommended clean baseline)."
    )
    if pkg == "Newtonian":
        return f"Newtonian baseline · {vdsg_part} Not paper correspondence mode."
    if pkg == "TPFCore":
        return (
            f"TPFCore legacy_readout experimental · not direct_tpf · not paper mode · {vdsg_part}"
        )
    return f"{physics_package} · {vdsg_part}"


def _local_extrema_indices(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Interior local minima and maxima (for periapsis / apoapsis markers on r(t))."""
    n = len(y)
    if n < 3:
        return np.array([], dtype=int), np.array([], dtype=int)
    mins: list[int] = []
    maxs: list[int] = []
    for i in range(1, n - 1):
        if y[i] < y[i - 1] and y[i] <= y[i + 1]:
            mins.append(i)
        elif y[i] > y[i - 1] and y[i] >= y[i + 1]:
            maxs.append(i)
    return np.array(mins, dtype=int), np.array(maxs, dtype=int)


def plot_bh_orbit_validation_extras(
    snapshots: list,
    diag: dict[str, np.ndarray | str],
    output_dir: Path,
    physics_package: str,
    tpf_vdsg_coupling: float,
    *,
    context_label: str = "",
) -> None:
    """
    bh_orbit_validation only: orbit in x–y, separation with periapsis/apoapsis markers, zoomed orbit shape.
    """
    foot = bh_orbit_validation_plot_footnote(physics_package, tpf_vdsg_coupling)
    foot = foot + " | Axes: display units; CSV SI."
    ctx = (context_label.strip() + " — ") if context_label.strip() else ""
    t = diag["time"]
    sep = np.asarray(diag["pair_separation"], dtype=np.float64)
    mins_i, maxs_i = _local_extrema_indices(sep)
    series = series_display_for_two_body(
        "bh_orbit_validation", max_distance_m=float(np.max(sep)) if sep.size else 1.0
    )

    fig, ax = plt.subplots(figsize=(8, 8))
    xs: list[float] = []
    ys: list[float] = []
    for snap in snapshots:
        p = snap.positions
        if p.shape[0] >= 1:
            xs.append(float(p[0, 0]))
            ys.append(float(p[0, 1]))
    xs_a = np.array(xs, dtype=np.float64)
    ys_a = np.array(ys, dtype=np.float64)
    r_arr = np.hypot(xs_a, ys_a)
    half_extent = float(np.nanmax(r_arr)) * 1.2 if r_arr.size else 1.0
    spatial: SpatialDisplay = spatial_display_for_xy_plot("bh_orbit_validation", half_extent)
    xf = spatial.factor
    u = spatial.unit
    ax.plot(xs_a * xf, ys_a * xf, "-", color="steelblue", lw=1.2, label="Star path (lab frame)")
    ax.scatter(xs_a * xf, ys_a * xf, s=10, c="steelblue", zorder=5)
    ax.scatter([0.0], [0.0], s=120, c="black", marker="*", zorder=6, label="BH at origin")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(f"x ({u})")
    ax.set_ylabel(f"y ({u})")
    set_fitted_title(
        ax,
        ctx + "Primary: trajectory x-y (one star, fixed BH at origin) (experimental)",
        fontsize=12,
        min_fontsize=6,
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    apply_suppress_tick_offset(ax)
    add_fitted_footer(fig, foot, fontsize=7, min_fontsize=5)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.14)
    fig.savefig(output_dir / "bh_orbit_trajectory_xy.png", dpi=150)
    plt.close(fig)

    r_max = float(np.nanmax(r_arr)) if r_arr.size else 1.0
    zoom = max(r_max * 1.25, 1.0)
    spatial_z = spatial_display_for_xy_plot("bh_orbit_validation", zoom)
    xfz = spatial_z.factor
    uz = spatial_z.unit
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(xs_a * xfz, ys_a * xfz, "-", color="#1f77b4", lw=1.2)
    ax.scatter(xs_a * xfz, ys_a * xfz, s=8, c="#1f77b4", zorder=5)
    ax.scatter([0.0], [0.0], s=100, c="black", marker="*", zorder=6)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlim(-zoom * xfz, zoom * xfz)
    ax.set_ylim(-zoom * xfz, zoom * xfz)
    ax.set_xlabel(f"x ({uz})")
    ax.set_ylabel(f"y ({uz})")
    set_fitted_title(
        ax,
        ctx + "Primary: trajectory x-y zoom to extent (experimental)",
        fontsize=12,
        min_fontsize=6,
    )
    ax.grid(True, alpha=0.3)
    apply_suppress_tick_offset(ax)
    add_fitted_footer(fig, foot, fontsize=7, min_fontsize=5)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.12)
    fig.savefig(output_dir / "bh_orbit_trajectory_xy_zoom.png", dpi=150)
    plt.close(fig)

    t_plot = np.asarray(t, dtype=np.float64) * series.time_factor
    sep_plot = sep * series.distance_factor
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(t_plot, sep_plot, "-", color="steelblue", lw=1.2, label=f"|r_star| ({series.distance_unit})")
    ax.scatter(t_plot, sep_plot, s=8, c="steelblue", zorder=5)
    for idx in mins_i[:12]:
        ax.axvline(float(t_plot[idx]), color="#2ca02c", ls="--", lw=0.8, alpha=0.7)
    for idx in maxs_i[:12]:
        ax.axvline(float(t_plot[idx]), color="#d62728", ls=":", lw=0.8, alpha=0.7)
    ax.set_xlabel(f"Time ({series.time_unit})")
    ax.set_ylabel(f"Separation ({series.distance_unit})")
    set_fitted_title(
        ax,
        ctx + "Primary: star–BH separation vs time (green dashed ≈ periapsis-like minima; "
        "red dotted ≈ apoapsis-like maxima; sampling-limited)",
        fontsize=12,
        min_fontsize=6,
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    apply_suppress_tick_offset(ax)
    add_fitted_footer(fig, foot, fontsize=7, min_fontsize=5)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18)
    fig.savefig(output_dir / "bh_orbit_separation_extrema.png", dpi=150)
    plt.close(fig)


def _series_for_bulk_diagnostics(
    simulation_mode: str,
    two_body_secondary: bool,
    diagnostics: dict[str, np.ndarray],
    cutoff_radius: float,
) -> SeriesDisplay:
    max_r = max(
        float(np.max(diagnostics["median_r"])),
        float(np.max(diagnostics["mean_r"])),
        float(np.max(diagnostics["std_r"])),
        float(np.max(diagnostics["max_r"])),
        float(cutoff_radius),
    )
    max_t = float(np.max(diagnostics["time"]))
    if two_body_secondary and simulation_mode == "earth_moon_benchmark":
        return series_display_for_two_body("earth_moon_benchmark", max_distance_m=max_r)
    if simulation_mode == "galaxy":
        return series_display_for_galaxy_diagnostics(max_r, max_t)
    return series_display_generic_validation(max_r)


def plot_and_save_all(
    diagnostics: dict[str, np.ndarray],
    output_dir: Path,
    cutoff_radius: float,
    *,
    lab_frame_secondary: bool = False,
    context_label: str = "",
    simulation_mode: str = "galaxy",
    two_body_secondary: bool = False,
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
    ctx = (context_label.strip() + " — ") if context_label.strip() else ""
    series = _series_for_bulk_diagnostics(
        simulation_mode, two_body_secondary, diagnostics, cutoff_radius
    )
    disp_note = "Axes use display units; CSV remains SI."
    t = diagnostics["time"]
    du = series.distance_unit
    _save_plot(
        output_dir / "diagnostic_median_radius.png",
        t,
        diagnostics["median_r"],
        f"Median radius from origin ({du})",
        ctx + prefix + "Median |r| from origin vs time | " + disp_note,
        series=series,
        y_scale=series.distance_factor,
    )
    _save_plot(
        output_dir / "diagnostic_mean_radius.png",
        t,
        diagnostics["mean_r"],
        f"Mean radius from origin ({du})",
        ctx + prefix + "Mean |r| from origin vs time | " + disp_note,
        series=series,
        y_scale=series.distance_factor,
    )
    _save_plot(
        output_dir / "diagnostic_std_radius.png",
        t,
        diagnostics["std_r"],
        f"Std of |r| from origin ({du})",
        ctx + prefix + "Std of star radii from origin vs time | " + disp_note,
        series=series,
        y_scale=series.distance_factor,
    )
    _save_plot(
        output_dir / "diagnostic_max_radius.png",
        t,
        diagnostics["max_r"],
        f"Max radius from origin ({du})",
        ctx + prefix + "Max |r| from origin vs time | " + disp_note,
        series=series,
        y_scale=series.distance_factor,
    )
    _save_plot_fraction(
        output_dir / "diagnostic_frac_vr_positive.png",
        t,
        diagnostics["frac_vr_pos"],
        "Fraction with v_r > 0 (origin radial)",
        ctx + prefix + "Fraction of bodies with v_r > 0 vs time | " + disp_note,
        series=series,
    )
    _save_plot_fraction(
        output_dir / "diagnostic_frac_vr_negative.png",
        t,
        diagnostics["frac_vr_neg"],
        "Fraction with v_r < 0 (origin radial)",
        ctx + prefix + "Fraction of bodies with v_r < 0 vs time | " + disp_note,
        series=series,
    )
    cut_disp = float(cutoff_radius) * series.distance_factor
    _save_plot_fraction(
        output_dir / "diagnostic_frac_beyond_cutoff.png",
        t,
        diagnostics["frac_beyond_cutoff"],
        "Fraction beyond cutoff (origin radial)",
        ctx
        + prefix
        + f"Fraction beyond r = {cut_disp:.6g} {du} (diagnostic cutoff; display) vs time | "
        + disp_note,
        series=series,
    )
    Lz_plot, Lz_ylab = scale_angular_momentum_origin_display(np.asarray(diagnostics["L_z"]))
    _save_plot_with_note(
        output_dir / "diagnostic_angular_momentum_z.png",
        t,
        Lz_plot,
        Lz_ylab,
        ctx + prefix + "Total L_z about simulation origin vs time | " + disp_note,
        series=series,
    )
