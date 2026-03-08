"""
Gravity calculations for the N-body simulation.
Designed with a pluggable interface: swap compute_accelerations_newtonian
for compute_accelerations_tpf (or similar) when implementing alternative gravity.
"""

from typing import Protocol

import numpy as np


def compute_accelerations_newtonian(
    positions: np.ndarray,
    masses: np.ndarray,
    bh_mass: float,
    softening: float,
    star_star: bool = True,
) -> np.ndarray:
    """
    Compute gravitational accelerations using Newtonian gravity.
    Central black hole always included; star-star gravity optional.

    Args:
        positions: (n_stars, 2) array of (x, y) positions
        masses: (n_stars,) array of star masses
        bh_mass: Mass of central black hole (fixed at origin)
        softening: Softening length (epsilon)
        star_star: If True, include pairwise star-star gravity; if False, BH only.

    Returns:
        accelerations: (n_stars, 2) array of (ax, ay) accelerations
    """
    n = len(positions)
    accelerations = np.zeros_like(positions)

    # Central black hole contribution (fixed at origin)
    r_sq_bh = np.sum(positions**2, axis=1) + softening**2
    r_mag_bh = np.sqrt(r_sq_bh)
    acc_mag_bh = bh_mass / (r_sq_bh * r_mag_bh)
    accelerations -= acc_mag_bh[:, np.newaxis] * positions

    if star_star:
        # Pairwise star-star interactions (fully vectorized O(n^2))
        dx = positions[:, 0, np.newaxis] - positions[:, 0]  # (n, n)
        dy = positions[:, 1, np.newaxis] - positions[:, 1]  # (n, n)
        r_sq = dx**2 + dy**2 + softening**2
        np.fill_diagonal(r_sq, np.inf)
        r_mag = np.sqrt(r_sq)
        acc_mag = masses[np.newaxis, :] / (r_sq * r_mag)
        accelerations[:, 0] += np.sum(acc_mag * dx, axis=1)
        accelerations[:, 1] += np.sum(acc_mag * dy, axis=1)

    return accelerations


# --- Pluggable gravity interface ---
# Define a protocol so any acceleration function can be swapped in.
# Later: implement compute_accelerations_tpf(...) with same signature
# and pass it to the integrator instead of compute_accelerations_newtonian.


class GravityComputer(Protocol):
    """Protocol for gravity/acceleration computation functions."""

    def __call__(
        self,
        positions: np.ndarray,
        masses: np.ndarray,
        bh_mass: float,
        softening: float,
    ) -> np.ndarray:
        """Compute accelerations. Same signature as compute_accelerations_newtonian."""
        ...


def compute_potential_energy(
    positions: np.ndarray,
    masses: np.ndarray,
    bh_mass: float,
    softening: float,
) -> float:
    """
    Compute total gravitational potential energy (Newtonian) for diagnostics.
    PE = -G * sum over pairs of (m_i * m_j / r_ij) + black hole contributions.
    Uses softening in denominators.
    """
    n = len(positions)
    pe = 0.0

    # Black hole - star potential
    for i in range(n):
        r = np.sqrt(np.sum(positions[i] ** 2) + softening**2)
        pe -= bh_mass * masses[i] / r

    # Star-star potential (each pair counted once)
    for i in range(n):
        for j in range(i + 1, n):
            dr = positions[j] - positions[i]
            r = np.sqrt(np.sum(dr**2) + softening**2)
            pe -= masses[i] * masses[j] / r

    return pe


def compute_kinetic_energy(velocities: np.ndarray, masses: np.ndarray) -> float:
    """Compute total kinetic energy: 0.5 * sum(m_i * v_i^2)."""
    v_sq = np.sum(velocities**2, axis=1)
    return 0.5 * np.sum(masses * v_sq)
