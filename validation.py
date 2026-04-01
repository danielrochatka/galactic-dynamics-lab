"""
Validation modes for testing force calculation and integrator correctness.
Reuses physics and simulation modules; saves outputs to config.output_dir.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from config import VALIDATION_MODES, SimulationConfig
from physics import (
    compute_accelerations_newtonian,
    compute_kinetic_energy,
    compute_potential_energy,
)
from simulation import run_simulation


def _L_z_total(positions: np.ndarray, velocities: np.ndarray, masses: np.ndarray) -> float:
    return float(np.sum(masses * (positions[:, 0] * velocities[:, 1] - positions[:, 1] * velocities[:, 0])))


def _P_total(velocities: np.ndarray, masses: np.ndarray) -> np.ndarray:
    return np.sum(masses[:, np.newaxis] * velocities, axis=0)


def _run_and_report(
    config: SimulationConfig,
    positions: np.ndarray,
    velocities: np.ndarray,
    masses: np.ndarray,
    bh_mass: float,
    star_star: bool,
    n_steps: int,
    snapshot_every: int,
) -> list:
    """Run simulation with given IC and return snapshots."""
    def acc(pos, m, bh, soft):
        return compute_accelerations_newtonian(pos, m, bh, soft, star_star=star_star)
    return run_simulation(
        positions, velocities, masses,
        bh_mass=bh_mass,
        softening=config.softening,
        dt=config.dt,
        n_steps=n_steps,
        snapshot_every=snapshot_every,
        compute_acc=acc,
    )


def run_bh_orbit_validation(config: SimulationConfig) -> None:
    """One star orbiting fixed BH; same IC family as C++ bh_orbit_validation / timestep_convergence."""
    r0 = config.validation_two_body_radius
    v_circ = np.sqrt(config.bh_mass / r0)
    v0 = config.validation_two_body_speed_ratio * v_circ
    positions = np.array([[r0, 0.0]])
    velocities = np.array([[0.0, v0]])
    masses = np.array([config.star_mass])
    n_steps = config.validation_n_steps
    snap_every = config.validation_snapshot_every

    print("BH orbit validation: 1 star, r0={}, v0={} (speed_ratio={})".format(
        r0, v0, config.validation_two_body_speed_ratio))
    snapshots = _run_and_report(
        config, positions, velocities, masses,
        bh_mass=config.bh_mass, star_star=False,
        n_steps=n_steps, snapshot_every=snap_every,
    )

    # Trajectory plot
    xs = np.array([s.positions[0, 0] for s in snapshots])
    ys = np.array([s.positions[0, 1] for s in snapshots])
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    ax.plot(xs, ys, "b-", alpha=0.7, label="Star")
    ax.scatter([0], [0], s=200, c="yellow", marker="*", edgecolors="orange", zorder=10, label="BH")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Two-body orbit trajectory")
    ax.legend()
    ax.grid(True, alpha=0.3)
    out = config.output_dir / "validation_two_body_trajectory.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    # Radii and conservation
    radii = np.array([np.sqrt(s.positions[0, 0]**2 + s.positions[0, 1]**2) for s in snapshots])
    L_z_arr = np.array([_L_z_total(s.positions, s.velocities, masses) for s in snapshots])
    E_arr = np.array([
        compute_kinetic_energy(s.velocities, masses) + compute_potential_energy(s.positions, masses, config.bh_mass, config.softening)
        for s in snapshots
    ])

    print("\n--- Two-body orbit report ---")
    print(f"  Initial radius:  {radii[0]:.4f}")
    print(f"  Final radius:    {radii[-1]:.4f}")
    print(f"  Min radius:      {radii.min():.4f}")
    print(f"  Max radius:      {radii.max():.4f}")
    print(f"  Angular momentum L_z: {L_z_arr[0]:.4f} (initial), {L_z_arr[-1]:.4f} (final)")
    print(f"  Total energy:    {E_arr[0]:.4f} (initial), {E_arr[-1]:.4f} (final), drift = {abs(E_arr[-1]-E_arr[0]):.2e}")


def run_symmetric_pair(config: SimulationConfig) -> None:
    """Two stars mirrored about origin; verify symmetry and momentum conservation."""
    a = config.validation_symmetric_separation
    v = config.validation_symmetric_speed
    positions = np.array([[-a, 0.0], [a, 0.0]])
    velocities = np.array([[0.0, v], [0.0, -v]])
    masses = np.array([config.star_mass, config.star_mass])
    bh_mass = config.bh_mass if config.validation_symmetric_include_bh else 0.0
    n_steps = config.validation_n_steps
    snap_every = config.validation_snapshot_every

    print("Symmetric pair: 2 stars at ±({}, 0), v=(0, ±{}), include_bh={}".format(
        a, v, config.validation_symmetric_include_bh))
    snapshots = _run_and_report(
        config, positions, velocities, masses,
        bh_mass=bh_mass, star_star=True,
        n_steps=n_steps, snapshot_every=snap_every,
    )

    # Symmetry error: expect pos[1] = -pos[0]
    sym_errors = np.array([
        np.sqrt((s.positions[0, 0] + s.positions[1, 0])**2 + (s.positions[0, 1] + s.positions[1, 1])**2)
        for s in snapshots
    ])
    P_mags = np.array([np.linalg.norm(_P_total(s.velocities, masses)) for s in snapshots])
    L_z_arr = np.array([_L_z_total(s.positions, s.velocities, masses) for s in snapshots])
    E_arr = np.array([
        compute_kinetic_energy(s.velocities, masses) + compute_potential_energy(s.positions, masses, bh_mass, config.softening)
        for s in snapshots
    ])

    # Trajectory plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    x0 = [s.positions[0, 0] for s in snapshots]
    y0 = [s.positions[0, 1] for s in snapshots]
    x1 = [s.positions[1, 0] for s in snapshots]
    y1 = [s.positions[1, 1] for s in snapshots]
    ax.plot(x0, y0, "b-", alpha=0.7, label="Star 1")
    ax.plot(x1, y1, "r-", alpha=0.7, label="Star 2")
    if bh_mass > 0:
        ax.scatter([0], [0], s=200, c="yellow", marker="*", edgecolors="orange", zorder=10, label="BH")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Symmetric pair trajectory")
    ax.legend()
    ax.grid(True, alpha=0.3)
    out = config.output_dir / "validation_symmetric_pair_trajectory.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    print("\n--- Symmetric pair report ---")
    print(f"  Symmetry error (|pos0+pos1|): max = {sym_errors.max():.2e}, final = {sym_errors[-1]:.2e}")
    print(f"  Total linear momentum |P|: initial = {P_mags[0]:.4f}, final = {P_mags[-1]:.4f}")
    print(f"  Total L_z: initial = {L_z_arr[0]:.4f}, final = {L_z_arr[-1]:.4f}")
    print(f"  Total energy: initial = {E_arr[0]:.4f}, final = {E_arr[-1]:.4f}, drift = {abs(E_arr[-1]-E_arr[0]):.2e}")


def run_small_n_conservation(config: SimulationConfig) -> None:
    """Small N-body with full pairwise gravity; check E, L_z, P conservation."""
    n = max(3, min(10, config.validation_small_n))
    rng = np.random.default_rng(42)
    # Deterministic: ring of radius ~25
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False) + 0.1
    r = 20 + rng.uniform(-2, 2, n)
    positions = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    v_circ = np.sqrt(config.bh_mass / r)
    velocities = np.column_stack([-v_circ * np.sin(theta), v_circ * np.cos(theta)])
    masses = np.full(n, config.star_mass)
    n_steps = config.validation_n_steps
    snap_every = config.validation_snapshot_every

    print("Small-N conservation: n={}, full pairwise gravity".format(n))
    snapshots = _run_and_report(
        config, positions, velocities, masses,
        bh_mass=config.bh_mass, star_star=True,
        n_steps=n_steps, snapshot_every=snap_every,
    )

    time_arr = np.array([s.time for s in snapshots])
    E_arr = np.array([
        compute_kinetic_energy(s.velocities, masses) + compute_potential_energy(s.positions, masses, config.bh_mass, config.softening)
        for s in snapshots
    ])
    L_z_arr = np.array([_L_z_total(s.positions, s.velocities, masses) for s in snapshots])
    P_mag_arr = np.array([np.linalg.norm(_P_total(s.velocities, masses)) for s in snapshots])

    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    axes[0].plot(time_arr, E_arr, "b-")
    axes[0].set_ylabel("Total energy")
    axes[0].set_title("Small-N conservation")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(time_arr, L_z_arr, "g-")
    axes[1].set_ylabel("Total L_z")
    axes[1].grid(True, alpha=0.3)
    axes[2].plot(time_arr, P_mag_arr, "r-")
    axes[2].set_ylabel("|Total P|")
    axes[2].set_xlabel("Time")
    axes[2].grid(True, alpha=0.3)
    fig.tight_layout()
    out = config.output_dir / "validation_small_n_conservation.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    print("\n--- Small-N conservation report ---")
    print(f"  Total energy:   initial = {E_arr[0]:.4f}, final = {E_arr[-1]:.4f}, rel drift = {abs(E_arr[-1]-E_arr[0])/abs(E_arr[0]) if E_arr[0]!=0 else 0:.2e}")
    print(f"  Total L_z:      initial = {L_z_arr[0]:.4f}, final = {L_z_arr[-1]:.4f}, rel drift = {abs(L_z_arr[-1]-L_z_arr[0])/abs(L_z_arr[0]) if L_z_arr[0]!=0 else 0:.2e}")
    print(f"  |Total P|:      initial = {P_mag_arr[0]:.4f}, final = {P_mag_arr[-1]:.4f}")


def run_timestep_convergence(config: SimulationConfig) -> None:
    """Run star-around-BH IC at dt, dt/2, dt/4 (same total time) and compare outputs."""
    r0 = config.validation_two_body_radius
    v_circ = np.sqrt(config.bh_mass / r0)
    v0 = config.validation_two_body_speed_ratio * v_circ
    positions0 = np.array([[r0, 0.0]])
    velocities0 = np.array([[0.0, v0]])
    masses = np.array([config.star_mass])
    n_steps_base = config.validation_n_steps
    snap_every = config.validation_snapshot_every
    total_time = n_steps_base * config.dt

    dts = [config.dt, config.dt / 2, config.dt / 4]
    results = []
    for dt in dts:
        n_steps = int(round(total_time / dt))
        def acc(pos, m, bh, soft):
            return compute_accelerations_newtonian(pos, m, bh, soft, star_star=False)
        snapshots = run_simulation(
            positions0.copy(), velocities0.copy(), masses,
            bh_mass=config.bh_mass, softening=config.softening,
            dt=dt, n_steps=n_steps, snapshot_every=max(1, snap_every),
            compute_acc=acc,
        )
        s = snapshots[-1]
        r_final = np.sqrt(s.positions[0, 0]**2 + s.positions[0, 1]**2)
        L_z = _L_z_total(s.positions, s.velocities, masses)
        E0 = compute_kinetic_energy(velocities0, masses) + compute_potential_energy(positions0, masses, config.bh_mass, config.softening)
        E_final = compute_kinetic_energy(s.velocities, masses) + compute_potential_energy(s.positions, masses, config.bh_mass, config.softening)
        results.append({
            "dt": dt,
            "final_pos": s.positions[0].copy(),
            "final_r": r_final,
            "L_z": L_z,
            "E_drift": abs(E_final - E0),
        })
    base = results[0]
    print("\n--- Timestep convergence (bh_orbit_validation IC) ---")
    print("  dt      final_x   final_y   final_r   L_z       E_drift")
    for r in results:
        print("  {:.6f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.2e}".format(
            r["dt"], r["final_pos"][0], r["final_pos"][1], r["final_r"], r["L_z"], r["E_drift"]))
    print("\n  Change from dt to dt/4:")
    print("    final_pos delta: ", np.linalg.norm(results[2]["final_pos"] - base["final_pos"]))
    print("    final_r  delta: ", abs(results[2]["final_r"] - base["final_r"]))
    print("    L_z      delta: ", abs(results[2]["L_z"] - base["L_z"]))

    out = config.output_dir / "validation_timestep_convergence.txt"
    lines = [
        "Timestep convergence (bh_orbit_validation IC)",
        "dt\tfinal_x\tfinal_y\tfinal_r\tL_z\tE_drift",
    ]
    for r in results:
        lines.append("{}\t{}\t{}\t{}\t{}\t{}".format(r["dt"], r["final_pos"][0], r["final_pos"][1], r["final_r"], r["L_z"], r["E_drift"]))
    out.write_text("\n".join(lines))
    print(f"Saved: {out}")


def run_earth_moon_benchmark_notice(_config: SimulationConfig) -> None:
    """Earth–Moon SI benchmark exists only in cpp_sim (galaxy_sim); Python pipeline does not implement it."""
    print(
        "simulation_mode=earth_moon_benchmark is implemented by cpp_sim/galaxy_sim (Earth–Moon SI), "
        "not this Python pipeline. Use: ./cpp_sim/galaxy_sim earth_moon_benchmark"
    )


def run_validation(config: SimulationConfig) -> None:
    """Dispatch to the requested validation mode."""
    mode = config.simulation_mode
    if mode not in VALIDATION_MODES or mode == "galaxy":
        raise ValueError("simulation_mode must be one of: " + ", ".join(VALIDATION_MODES) + " (and not 'galaxy')")
    print("Validation mode:", mode)
    if mode == "bh_orbit_validation":
        run_bh_orbit_validation(config)
    elif mode == "earth_moon_benchmark":
        run_earth_moon_benchmark_notice(config)
    elif mode == "two_body_orbit":
        print(
            'Warning: simulation_mode "two_body_orbit" is deprecated in the Python pipeline. '
            "It previously meant one star + BH; use bh_orbit_validation. "
            "(In C++, the string two_body_orbit is deprecated: use earth_moon_benchmark for Earth–Moon.)"
        )
        run_bh_orbit_validation(config)
    elif mode == "symmetric_pair":
        run_symmetric_pair(config)
    elif mode == "small_n_conservation":
        run_small_n_conservation(config)
    elif mode == "timestep_convergence":
        run_timestep_convergence(config)
