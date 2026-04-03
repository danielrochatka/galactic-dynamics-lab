"""
Configuration parameters for the 2D galaxy N-body simulation.
Masses and lengths use SI (kg, m) with astronomical-scale defaults (solar masses, kpc-order disk).

Local overrides: if configs/my.local.cfg or any configs/local/*.cfg exists,
it is loaded after defaults and applied to SimulationConfig. Keys match
attribute names (e.g. n_stars, dt). Use # for comments.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any


VALIDATION_MODES = (
    "galaxy",
    "earth_moon_benchmark",
    "bh_orbit_validation",
    "two_body_orbit",
    "symmetric_pair",
    "small_n_conservation",
    "timestep_convergence",
)

# Project root (directory containing config.py)
_PROJECT_ROOT = Path(__file__).resolve().parent

# SI defaults (match cpp_sim/config.hpp): one solar mass; ~10^6 M_sun SMBH; ~10 kpc disk scale.
SOLAR_MASS_KG = 1.98847e30
DEFAULT_BH_MASS_KG = 1.0e6 * SOLAR_MASS_KG
_DEFAULT_GALAXY_RADIUS_M = 3.0e20
_DEFAULT_INNER_RADIUS_M = 3.0e19
_DEFAULT_SOFTENING_M = 1.0e16

# Paths to check for local config (relative to project root)
CONFIG_SEARCH_PATHS = [
    _PROJECT_ROOT / "configs" / "my.local.cfg",
    _PROJECT_ROOT / "configs" / "local" / "my.local.cfg",
]
CONFIG_DIR_LOCAL = _PROJECT_ROOT / "configs" / "local"


@dataclass
class SimulationConfig:
    """Configuration for the galaxy simulation."""

    # Mode: galaxy (default) or validation mode
    simulation_mode: str = "timestep_convergence"

    # Star distribution
    n_stars: int = 1000
    star_mass: float = SOLAR_MASS_KG

    # Central black hole (fixed at origin)
    bh_mass: float = DEFAULT_BH_MASS_KG

    # Disk geometry (meters; ~1 kpc inner, ~10 kpc outer)
    inner_radius: float = _DEFAULT_INNER_RADIUS_M
    outer_radius: float = _DEFAULT_GALAXY_RADIUS_M

    # Integration
    dt: float = 0.01
    n_steps: int = 120000
    snapshot_every: int = 10

    # Physics
    softening: float = _DEFAULT_SOFTENING_M
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
    diagnostic_cutoff_radius: float = _DEFAULT_GALAXY_RADIUS_M

    # Validation mode parameters
    validation_two_body_radius: float = 5.0e18  # ~0.5 pc (works with default SMBH masses)
    validation_two_body_speed_ratio: float = 1.0  # 1.0 = circular, <1 or >1 = elliptical
    validation_symmetric_include_bh: bool = True
    validation_symmetric_separation: float = 7.48e10  # half-separation along x; full sep = 2x (~1 AU)
    validation_symmetric_speed: float = 3.0e4  # ~30 km/s (order of AU-scale binary)
    validation_small_n: int = 5  # n=3..10 for small_n_conservation
    validation_n_steps: int = 5000
    validation_snapshot_every: int = 5

    def __post_init__(self) -> None:
        """Create output directory if it doesn't exist."""
        self.output_dir.mkdir(parents=True, exist_ok=True)


def _parse_cfg_value(raw: str) -> Any:
    """Parse a single config value (bool, int, float, or str)."""
    raw = raw.strip()
    if raw.lower() in ("true", "yes", "1"):
        return True
    if raw.lower() in ("false", "no", "0"):
        return False
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        pass
    return raw


def load_cfg_file(path: Path) -> dict[str, Any]:
    """Load key = value pairs from a .cfg file. Skips comments and empty lines."""
    out: dict[str, Any] = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        val = val.strip()
        if key:
            out[key] = _parse_cfg_value(val)
    return out


def _config_file_paths() -> list[Path]:
    """Return paths to load: single-file overrides first, then configs/local/*.cfg."""
    paths: list[Path] = []
    for p in CONFIG_SEARCH_PATHS:
        if p.exists():
            paths.append(p)
    if CONFIG_DIR_LOCAL.exists():
        for p in sorted(CONFIG_DIR_LOCAL.glob("*.cfg")):
            paths.append(p)
    return paths


def load_local_config() -> dict[str, Any]:
    """Load and merge all local config files (later files override earlier)."""
    merged: dict[str, Any] = {}
    for path in _config_file_paths():
        merged.update(load_cfg_file(path))
    return merged


def apply_cfg_to_config(cfg: dict[str, Any], config: SimulationConfig) -> None:
    """Apply a parsed cfg dict onto a SimulationConfig instance. Ignores unknown keys."""
    for key, value in cfg.items():
        if not hasattr(config, key):
            continue
        if key == "output_dir":
            config.output_dir = Path(value) if isinstance(value, str) else value
        else:
            setattr(config, key, value)


def load_config() -> SimulationConfig:
    """Build SimulationConfig from defaults, then apply any configs/my.local.cfg and configs/local/*.cfg."""
    config = SimulationConfig()
    local = load_local_config()
    if local:
        apply_cfg_to_config(local, config)
    return config
