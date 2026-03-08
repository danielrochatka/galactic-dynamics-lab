"""
N-body simulation with velocity Verlet integrator.
Accepts a pluggable gravity function for computing accelerations.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np

from physics import compute_accelerations_newtonian


@dataclass
class Snapshot:
    """A single snapshot of the simulation state at one timestep."""

    step: int
    time: float
    positions: np.ndarray
    velocities: np.ndarray


def velocity_verlet_step(
    positions: np.ndarray,
    velocities: np.ndarray,
    masses: np.ndarray,
    bh_mass: float,
    softening: float,
    dt: float,
    compute_acc: Callable[
        [np.ndarray, np.ndarray, float, float],
        np.ndarray,
    ],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Single velocity Verlet integration step.
    x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
    v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    """
    acc = compute_acc(positions, masses, bh_mass, softening)

    # Position half-step
    positions_new = positions + velocities * dt + 0.5 * acc * dt**2

    # New acceleration
    acc_new = compute_acc(positions_new, masses, bh_mass, softening)

    # Velocity update
    velocities_new = velocities + 0.5 * (acc + acc_new) * dt

    return positions_new, velocities_new


def run_simulation(
    positions: np.ndarray,
    velocities: np.ndarray,
    masses: np.ndarray,
    bh_mass: float,
    softening: float,
    dt: float,
    n_steps: int,
    snapshot_every: int,
    compute_acc: Callable[
        [np.ndarray, np.ndarray, float, float],
        np.ndarray,
    ] | None = None,
    progress_callback: Callable[[int, int], None] | None = None,
) -> list[Snapshot]:
    """
    Run the N-body simulation and collect snapshots.

    Args:
        positions: Initial positions (n_stars, 2)
        velocities: Initial velocities (n_stars, 2)
        masses: Star masses (n_stars,)
        bh_mass: Central black hole mass
        softening: Gravitational softening
        dt: Timestep
        n_steps: Total number of integration steps
        snapshot_every: Save snapshot every N steps
        compute_acc: Gravity function. Defaults to Newtonian.
        progress_callback: Optional callback(step, n_steps) for progress

    Returns:
        List of Snapshot objects
    """
    if compute_acc is None:
        compute_acc = compute_accelerations_newtonian

    pos = positions.copy()
    vel = velocities.copy()
    snapshots: list[Snapshot] = []

    # Save initial snapshot
    snapshots.append(
        Snapshot(step=0, time=0.0, positions=pos.copy(), velocities=vel.copy())
    )

    for step in range(1, n_steps + 1):
        pos, vel = velocity_verlet_step(
            pos, vel, masses, bh_mass, softening, dt, compute_acc
        )

        if step % snapshot_every == 0:
            snapshots.append(
                Snapshot(
                    step=step,
                    time=step * dt,
                    positions=pos.copy(),
                    velocities=vel.copy(),
                )
            )

        if progress_callback and step % 500 == 0:
            progress_callback(step, n_steps)

    return snapshots
