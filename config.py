"""
Configuration parameters for the 2D galaxy N-body simulation.
All values use normalized units.
"""

from dataclasses import dataclass
from pathlib import Path


@dataclass
class SimulationConfig:
    """Configuration for the galaxy simulation."""

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

    def __post_init__(self) -> None:
        """Create output directory if it doesn't exist."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
