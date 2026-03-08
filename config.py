"""
Configuration parameters for the 2D galaxy N-body simulation.
All values use normalized units.
"""

from dataclasses import dataclass
from pathlib import Path


VALIDATION_MODES = ("galaxy", "two_body_orbit", "symmetric_pair", "small_n_conservation", "timestep_convergence")


@dataclass
class SimulationConfig:
    """Configuration for the galaxy simulation."""

    # Mode: galaxy (default) or validation mode
    simulation_mode: str = "timestep_convergence"

    # Star distribution
    n_stars: int = 1000
    star_mass: float = 0.5

    # Central black hole (fixed at origin)
    bh_mass: float = 1000.0

    # Disk geometry
    inner_radius: float = 5.0
    outer_radius: float = 50.0

    # Integration
    dt: float = 0.01
    n_steps: int = 120000
    snapshot_every: int = 10

    # Physics
    softening: float = 1.0
    enable_star_star_gravity: bool = True  # False = stars feel only central BH

    # Initial conditions
    velocity_noise: float = 0.05  # Fractional random perturbation to tangential velocity

    # Output
    output_dir: Path = Path("outputs")
    animation_format: str = "mp4"  # "mp4" or "gif"

    # Rendering
    render_radius: float = 150.0  # Axis limit: plot from -render_radius to +render_radius
    render_animation: bool = False  # Skip MP4/GIF to speed up test runs
    render_diagnostics: bool = True  # Diagnostic time-series and vr plots
    render_initial_final: bool = True  # Static initial/final scatter plots

    # Diagnostics
    diagnostic_cutoff_radius: float = 50.0  # Fraction of stars beyond this radius (diagnostic 5)

    # Validation mode parameters
    validation_two_body_radius: float = 20.0
    validation_two_body_speed_ratio: float = 1.0  # 1.0 = circular, <1 or >1 = elliptical
    validation_symmetric_include_bh: bool = True
    validation_symmetric_separation: float = 30.0  # distance of each star from origin along x
    validation_symmetric_speed: float = 4.0  # tangential speed (y) for each star
    validation_small_n: int = 5  # n=3..10 for small_n_conservation
    validation_n_steps: int = 5000
    validation_snapshot_every: int = 5

    def __post_init__(self) -> None:
        """Create output directory if it doesn't exist."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
