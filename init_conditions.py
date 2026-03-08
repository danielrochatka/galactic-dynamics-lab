"""
Initial conditions for the rotating disk galaxy.
Stars are distributed radially and given mostly tangential velocities.
Velocities are set from enclosed mass (BH + stars inside r) for self-consistency.
"""

import numpy as np

# Minimum radius for v_circ to avoid division by zero
_MIN_RADIUS = 1e-8


def generate_disk(
    n_stars: int,
    inner_radius: float,
    outer_radius: float,
    bh_mass: float,
    star_mass: float,
    velocity_noise: float,
    seed: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate initial positions and velocities for a rotating disk.
    Stars are distributed in radius (inner to outer) with random polar angle.
    Tangential velocities are set from enclosed mass (BH + stars inside r)
    so the disk is self-consistent with the full gravity model.
    Optional small random velocity noise.

    Returns:
        positions: (n_stars, 2) array of (x, y)
        velocities: (n_stars, 2) array of (vx, vy)
        masses: (n_stars,) array, all equal to star_mass
    """
    rng = np.random.default_rng(seed)

    # Radial distribution: uniform in area (so more stars at larger r)
    u = rng.uniform(0, 1, n_stars)
    r_sq = inner_radius**2 + u * (outer_radius**2 - inner_radius**2)
    radii = np.sqrt(r_sq)

    # Random polar angle
    theta = rng.uniform(0, 2 * np.pi, n_stars)

    # Positions
    x = radii * np.cos(theta)
    y = radii * np.sin(theta)
    positions = np.column_stack([x, y])

    # Enclosed mass: for each star, count stars with strictly smaller radius
    # Sort by radius; rank[i] = number of stars with radius < radii[i]
    order = np.argsort(radii)
    n_inside = np.empty(n_stars, dtype=np.intp)
    n_inside[order] = np.arange(n_stars)
    enclosed_stellar_mass = n_inside * star_mass
    enclosed_mass = bh_mass + enclosed_stellar_mass

    # Circular velocity: v_circ = sqrt(G * M_enclosed / r), G=1 in normalized units
    r_safe = np.maximum(radii, _MIN_RADIUS)
    v_circ = np.sqrt(enclosed_mass / r_safe)

    # Tangential direction: (-sin(theta), cos(theta)) for CCW rotation
    vx = -v_circ * np.sin(theta)
    vy = v_circ * np.cos(theta)

    # Add velocity noise (small random perturbation)
    if velocity_noise > 0:
        noise_scale = velocity_noise * v_circ
        vx += rng.normal(0, noise_scale, n_stars)
        vy += rng.normal(0, noise_scale, n_stars)

    velocities = np.column_stack([vx, vy])
    masses = np.full(n_stars, star_mass, dtype=float)

    return positions, velocities, masses
