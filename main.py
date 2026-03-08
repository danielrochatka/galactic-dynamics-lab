"""
2D Galaxy N-Body Simulation

Run from project root:
    python main.py   # or: python3 main.py

Requires: numpy, matplotlib (pip install numpy matplotlib).
Optional: ffmpeg (for MP4 animation) or Pillow (for GIF fallback).

Outputs (in outputs/<run_id>/ per run, run_id = YYYYMMDD_HHMMSS):
    galaxy_initial.png, galaxy_final.png, galaxy.mp4/.gif, diagnostics, etc.

Gravity is pluggable: see physics.compute_accelerations_newtonian.
Swap for compute_accelerations_tpf (or similar) when implementing TPF gravity.
"""

from datetime import datetime
from pathlib import Path

import numpy as np

from config import SimulationConfig
from init_conditions import generate_disk
from physics import (
    compute_accelerations_newtonian,
    compute_kinetic_energy,
    compute_potential_energy,
)
from render import create_animation, has_ffmpeg, save_static_plot, save_radial_velocity_plot
from simulation import run_simulation
from diagnostics import compute_diagnostics, plot_and_save_all


def main() -> None:
    config = SimulationConfig()
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    config.output_dir = Path("outputs") / run_id
    config.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 50)
    print("2D Galaxy N-Body Simulation")
    print("=" * 50)
    print(f"Output directory: {config.output_dir.resolve()}")
    print(f"Stars: {config.n_stars}")
    print(f"Steps: {config.n_steps}, dt = {config.dt}")
    print(f"Snapshot every: {config.snapshot_every} steps")
    print(f"Rendering: initial_final={config.render_initial_final}, diagnostics={config.render_diagnostics}, animation={config.render_animation}")
    print()

    # Initial conditions
    print("Generating initial conditions (rotating disk)...")
    positions, velocities, masses = generate_disk(
        n_stars=config.n_stars,
        inner_radius=config.inner_radius,
        outer_radius=config.outer_radius,
        bh_mass=config.bh_mass,
        star_mass=config.star_mass,
        velocity_noise=config.velocity_noise,
    )
    print(f"  Inner radius: {config.inner_radius}, Outer: {config.outer_radius}")
    print(f"  Velocity noise: {config.velocity_noise}")
    print()

    # Static initial plot
    if config.render_initial_final:
        init_path = config.output_dir / "galaxy_initial.png"
        save_static_plot(
            positions,
            init_path,
            title="Galaxy – Initial",
            render_radius=config.render_radius,
        )
        print(f"Saved: {init_path}")

    # Run simulation
    print("\nRunning simulation (velocity Verlet)...")

    def progress(step: int, total: int) -> None:
        pct = 100 * step / total
        print(f"  Progress: {step}/{total} ({pct:.0f}%)")

    def gravity_acc(positions, masses, bh_mass, softening):
        return compute_accelerations_newtonian(
            positions, masses, bh_mass, softening,
            star_star=config.enable_star_star_gravity,
        )

    snapshots = run_simulation(
        positions=positions,
        velocities=velocities,
        masses=masses,
        bh_mass=config.bh_mass,
        softening=config.softening,
        dt=config.dt,
        n_steps=config.n_steps,
        snapshot_every=config.snapshot_every,
        compute_acc=gravity_acc,
        progress_callback=progress,
    )

    print(f"  Done. Collected {len(snapshots)} snapshots.")
    print()

    # Run metadata
    total_time = config.n_steps * config.dt
    run_info_path = config.output_dir / "run_info.txt"
    run_info_path.write_text(
        f"dt\t{config.dt}\n"
        f"n_steps\t{config.n_steps}\n"
        f"snapshot_every\t{config.snapshot_every}\n"
        f"softening\t{config.softening}\n"
        f"star_mass\t{config.star_mass}\n"
        f"bh_mass\t{config.bh_mass}\n"
        f"enable_star_star_gravity\t{config.enable_star_star_gravity}\n"
        f"total_simulated_time\t{total_time}\n"
        f"number_of_snapshots\t{len(snapshots)}\n"
    )
    print(f"Saved: {run_info_path}")

    # Static final plot
    final_snap = snapshots[-1]
    if config.render_initial_final:
        final_path = config.output_dir / "galaxy_final.png"
        save_static_plot(
            final_snap.positions,
            final_path,
            title="Galaxy – Final",
            render_radius=config.render_radius,
        )
        print(f"Saved: {final_path}")

    # Animation
    if config.render_animation:
        anim_base = config.output_dir / "galaxy"
        print("\nRendering animation...")
        if has_ffmpeg():
            print("  Using ffmpeg for MP4")
        else:
            print("  ffmpeg not found, using GIF (may require Pillow)")
        success = create_animation(
            snapshots,
            anim_base,
            render_radius=config.render_radius,
            interval=50,
            progress_interval=50,
        )
        if not success:
            print("  Animation save failed (install ffmpeg or Pillow)")
    else:
        print("\nSkipping animation (render_animation=False).")

    # Diagnostics: compute for summary; optionally save plots
    print("\nComputing diagnostics...")
    diag = compute_diagnostics(
        snapshots, masses, config.diagnostic_cutoff_radius
    )
    if config.render_diagnostics:
        plot_and_save_all(diag, config.output_dir, config.diagnostic_cutoff_radius)
        for name in [
            "diagnostic_median_radius.png",
            "diagnostic_mean_radius.png",
            "diagnostic_std_radius.png",
            "diagnostic_max_radius.png",
            "diagnostic_frac_vr_positive.png",
            "diagnostic_frac_vr_negative.png",
            "diagnostic_frac_beyond_cutoff.png",
            "diagnostic_angular_momentum_z.png",
        ]:
            print(f"Saved: {config.output_dir / name}")

        # Radial-velocity colored diagnostic: early, middle, late snapshot
        indices = [0, len(snapshots) // 2, len(snapshots) - 1]
        labels = ["early", "middle", "late"]
        for idx, label in zip(indices, labels):
            snap = snapshots[idx]
            path = config.output_dir / f"diagnostic_vr_{label}.png"
            save_radial_velocity_plot(
                snap.positions,
                snap.velocities,
                path,
                title=f"Radial velocity (v_r) — {label}, step {snap.step}, t = {snap.time:.1f}",
                render_radius=config.render_radius,
            )
            print(f"Saved: {path}")
    else:
        print("  Skipping diagnostic plots (render_diagnostics=False).")

    # Summary
    print("\n" + "=" * 50)
    print("Diagnostics")
    print("=" * 50)
    print(f"Number of stars:     {config.n_stars}")
    print(f"Number of snapshots: {len(snapshots)}")

    r_initial = np.sqrt(np.sum(positions**2, axis=1))
    r_final = np.sqrt(np.sum(final_snap.positions**2, axis=1))
    print(f"Initial radius: min = {r_initial.min():.2f}, max = {r_initial.max():.2f}")
    print(f"Final radius:   min = {r_final.min():.2f}, max = {r_final.max():.2f}")

    ke_init = compute_kinetic_energy(velocities, masses)
    pe_init = compute_potential_energy(positions, masses, config.bh_mass, config.softening)
    ke_final = compute_kinetic_energy(final_snap.velocities, masses)
    pe_final = compute_potential_energy(
        final_snap.positions, masses, config.bh_mass, config.softening
    )
    print(f"\nInitial KE: {ke_init:.2f}, PE: {pe_init:.2f}, Total: {ke_init + pe_init:.2f}")
    print(f"Final   KE: {ke_final:.2f}, PE: {pe_final:.2f}, Total: {ke_final + pe_final:.2f}")
    drift = abs((ke_final + pe_final) - (ke_init + pe_init))
    print(f"Energy drift: {drift:.4f} (should be small for Verlet)")

    print("\nDiagnostic summary (final snapshot):")
    print(f"  Median radius:        {diag['median_r'][-1]:.2f}")
    print(f"  Mean radius:          {diag['mean_r'][-1]:.2f}")
    print(f"  Std radius:           {diag['std_r'][-1]:.2f}")
    print(f"  Max radius:           {diag['max_r'][-1]:.2f}")
    print(f"  Fraction v_r > 0:     {diag['frac_vr_pos'][-1]:.4f}")
    print(f"  Fraction v_r < 0:     {diag['frac_vr_neg'][-1]:.4f}")
    print(f"  Fraction beyond r={config.diagnostic_cutoff_radius}: {diag['frac_beyond_cutoff'][-1]:.4f}")
    print(f"  Total L_z:            {diag['L_z'][-1]:.2f}")
    print(f"\nOutput directory: {config.output_dir.resolve()}")
    print()


if __name__ == "__main__":
    main()
