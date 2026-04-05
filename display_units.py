"""
Display-only unit selection for plots and renders.

All simulation state, snapshots, and CSV exports remain in SI. This module only
converts values at draw time and chooses axis labels / tick formatting.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# IAU / CODATA-friendly constants (display only)
METERS_PER_AU = 1.495978707e11
METERS_PER_LY = 9.460731e15
METERS_PER_PC = 3.085677581491367e16
METERS_PER_KPC = 1.0e3 * METERS_PER_PC
SECONDS_PER_MINUTE = 60.0
SECONDS_PER_HOUR = 3600.0
SECONDS_PER_DAY = 86400.0
SECONDS_PER_JULIAN_YEAR = 365.25 * SECONDS_PER_DAY
SECONDS_PER_KYR = 1.0e3 * SECONDS_PER_JULIAN_YEAR
SECONDS_PER_MYR = 1.0e6 * SECONDS_PER_JULIAN_YEAR


@dataclass(frozen=True)
class SpatialDisplay:
    """Multiply SI position (meters) by `factor` to get axis coordinates."""

    factor: float
    unit: str


@dataclass(frozen=True)
class SeriesDisplay:
    """Scale SI series to human-readable plot axes."""

    distance_factor: float
    distance_unit: str
    speed_factor: float
    speed_unit: str
    time_factor: float
    time_unit: str
    suppress_tick_offset: bool


@dataclass(frozen=True)
class DisplayUnitConfig:
    distance_unit: str = "auto"
    time_unit: str = "auto"
    velocity_unit: str = "auto"
    units_in_overlay: bool = True
    show_unit_reference: bool = True


def _boolish(raw: object, default: bool) -> bool:
    if raw is None:
        return default
    if isinstance(raw, bool):
        return raw
    if isinstance(raw, (int, float)):
        return int(raw) != 0
    s = str(raw).strip().lower()
    if s in ("1", "true", "yes", "on"):
        return True
    if s in ("0", "false", "no", "off"):
        return False
    return default


def display_unit_config_from_run_info(run_info: dict[str, str | int | float]) -> DisplayUnitConfig:
    return DisplayUnitConfig(
        distance_unit=str(run_info.get("display_distance_unit", "auto") or "auto").strip(),
        time_unit=str(run_info.get("display_time_unit", "auto") or "auto").strip(),
        velocity_unit=str(run_info.get("display_velocity_unit", "auto") or "auto").strip(),
        units_in_overlay=_boolish(run_info.get("display_units_in_overlay"), True),
        show_unit_reference=_boolish(run_info.get("display_show_unit_reference"), True),
    )


def _distance_from_unit(unit: str) -> SpatialDisplay:
    mapping = {
        "m": SpatialDisplay(1.0, "m"),
        "km": SpatialDisplay(1e-3, "km"),
        "AU": SpatialDisplay(1.0 / METERS_PER_AU, "AU"),
        "ly": SpatialDisplay(1.0 / METERS_PER_LY, "ly"),
        "pc": SpatialDisplay(1.0 / METERS_PER_PC, "pc"),
        "kpc": SpatialDisplay(1.0 / METERS_PER_KPC, "kpc"),
    }
    if unit not in mapping:
        raise ValueError(f"Unsupported distance unit: {unit}")
    return mapping[unit]


def _time_from_unit(unit: str) -> tuple[float, str]:
    mapping = {
        "s": (1.0, "s"),
        "min": (1.0 / SECONDS_PER_MINUTE, "min"),
        "hr": (1.0 / SECONDS_PER_HOUR, "hr"),
        "day": (1.0 / SECONDS_PER_DAY, "day"),
        "yr": (1.0 / SECONDS_PER_JULIAN_YEAR, "yr"),
        "kyr": (1.0 / SECONDS_PER_KYR, "kyr"),
        "Myr": (1.0 / SECONDS_PER_MYR, "Myr"),
    }
    if unit not in mapping:
        raise ValueError(f"Unsupported time unit: {unit}")
    return mapping[unit]


def _velocity_from_unit(unit: str) -> tuple[float, str]:
    mapping = {
        "m/s": (1.0, "m/s"),
        "km/s": (1e-3, "km/s"),
    }
    if unit not in mapping:
        raise ValueError(f"Unsupported velocity unit: {unit}")
    return mapping[unit]


def _auto_distance_unit(scale_m: float) -> str:
    x = float(scale_m)
    if not np.isfinite(x) or x <= 0:
        return "m"
    if x >= 0.2 * METERS_PER_KPC:
        return "kpc"
    if x >= 0.1 * METERS_PER_PC:
        return "pc"
    if x >= 0.5 * METERS_PER_LY:
        return "ly"
    if x >= 0.5 * METERS_PER_AU:
        return "AU"
    if x >= 1e3:
        return "km"
    return "m"


def _auto_time_unit(scale_s: float) -> str:
    x = float(scale_s)
    if not np.isfinite(x) or x < 0:
        return "s"
    if x >= 0.5 * SECONDS_PER_MYR:
        return "Myr"
    if x >= 0.5 * SECONDS_PER_KYR:
        return "kyr"
    if x >= 2.0 * SECONDS_PER_JULIAN_YEAR:
        return "yr"
    if x >= 2.0 * SECONDS_PER_DAY:
        return "day"
    if x >= 2.0 * SECONDS_PER_HOUR:
        return "hr"
    if x >= 2.0 * SECONDS_PER_MINUTE:
        return "min"
    return "s"


def _auto_velocity_unit(scale_m_s: float) -> str:
    x = float(scale_m_s)
    if not np.isfinite(x) or x < 0:
        return "m/s"
    if x >= 1e3:
        return "km/s"
    return "m/s"


def _resolve_distance_unit(preferred_unit: str, scale_m: float) -> SpatialDisplay:
    unit = preferred_unit if preferred_unit != "auto" else _auto_distance_unit(scale_m)
    return _distance_from_unit(unit)


def _resolve_time_unit(preferred_unit: str, scale_s: float) -> tuple[float, str]:
    unit = preferred_unit if preferred_unit != "auto" else _auto_time_unit(scale_s)
    return _time_from_unit(unit)


def _resolve_velocity_unit(preferred_unit: str, scale_m_s: float) -> tuple[float, str]:
    unit = preferred_unit if preferred_unit != "auto" else _auto_velocity_unit(scale_m_s)
    return _velocity_from_unit(unit)


def spatial_display_for_xy_plot(
    mode: str,
    half_axis_m: float,
    *,
    preferred_unit: str = "auto",
) -> SpatialDisplay:
    _ = mode
    return _resolve_distance_unit(preferred_unit, half_axis_m)


def series_display_for_two_body(
    mode: str,
    *,
    max_distance_m: float,
    max_time_s: float | None = None,
    max_speed_m_s: float | None = None,
    preferred_distance_unit: str = "auto",
    preferred_time_unit: str = "auto",
    preferred_velocity_unit: str = "auto",
) -> SeriesDisplay:
    if mode not in ("earth_moon_benchmark", "bh_orbit_validation"):
        raise ValueError(f"series_display_for_two_body: unsupported mode {mode!r}")
    md = float(max_distance_m)
    mt = float(max_time_s) if max_time_s is not None else 0.0
    mv = float(max_speed_m_s) if max_speed_m_s is not None else 0.0
    dist = _resolve_distance_unit(preferred_distance_unit, md)
    tf, tu = _resolve_time_unit(preferred_time_unit, mt)
    vf, vu = _resolve_velocity_unit(preferred_velocity_unit, mv)
    return SeriesDisplay(
        distance_factor=dist.factor,
        distance_unit=dist.unit,
        speed_factor=vf,
        speed_unit=vu,
        time_factor=tf,
        time_unit=tu,
        suppress_tick_offset=True,
    )


def series_display_for_galaxy_diagnostics(
    max_radius_m: float,
    max_time_s: float,
    *,
    max_speed_m_s: float = 0.0,
    preferred_distance_unit: str = "auto",
    preferred_time_unit: str = "auto",
    preferred_velocity_unit: str = "auto",
) -> SeriesDisplay:
    dist = _resolve_distance_unit(preferred_distance_unit, max_radius_m)
    tf, tu = _resolve_time_unit(preferred_time_unit, max_time_s)
    vf, vu = _resolve_velocity_unit(preferred_velocity_unit, max_speed_m_s)
    return SeriesDisplay(
        distance_factor=dist.factor,
        distance_unit=dist.unit,
        speed_factor=vf,
        speed_unit=vu,
        time_factor=tf,
        time_unit=tu,
        suppress_tick_offset=True,
    )


def series_display_generic_validation(
    max_radius_m: float,
    *,
    max_time_s: float = 0.0,
    max_speed_m_s: float = 0.0,
    preferred_distance_unit: str = "auto",
    preferred_time_unit: str = "auto",
    preferred_velocity_unit: str = "auto",
) -> SeriesDisplay:
    dist = _resolve_distance_unit(preferred_distance_unit, max_radius_m)
    tf, tu = _resolve_time_unit(preferred_time_unit, max_time_s)
    vf, vu = _resolve_velocity_unit(preferred_velocity_unit, max_speed_m_s)
    return SeriesDisplay(
        distance_factor=dist.factor,
        distance_unit=dist.unit,
        speed_factor=vf,
        speed_unit=vu,
        time_factor=tf,
        time_unit=tu,
        suppress_tick_offset=True,
    )


def rotation_curve_display(
    x_max_m: float,
    *,
    preferred_distance_unit: str = "auto",
    preferred_velocity_unit: str = "auto",
) -> tuple[float, float, str, str]:
    d = _resolve_distance_unit(preferred_distance_unit, x_max_m)
    v = _resolve_velocity_unit(preferred_velocity_unit, x_max_m / max(SECONDS_PER_DAY, 1.0))
    return (
        d.factor,
        v[0],
        f"Distance from galactic center ({d.unit})",
        f"Orbital speed ({v[1]})",
    )


def apply_suppress_tick_offset(ax) -> None:
    """Disable matplotlib's +1eN offset text on both axes (display-layer readability)."""
    import matplotlib.ticker as mticker

    for axis in (ax.xaxis, ax.yaxis):
        axis.get_major_formatter().set_useOffset(False)
        axis.set_major_formatter(mticker.ScalarFormatter(useOffset=False))


def scale_angular_momentum_display(y_si: np.ndarray) -> tuple[np.ndarray, str]:
    y = np.asarray(y_si, dtype=np.float64)
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return y, "L_z (SI)"
    m = float(np.nanmax(np.abs(finite)))
    if m <= 0 or not np.isfinite(m):
        return y, "L_z (SI)"
    exp = int(np.floor(np.log10(m)))
    if abs(exp) >= 6:
        exp3 = (exp // 3) * 3
        s = 10.0 ** (-exp3)
        return y * s, f"L_z about COM (×10^{exp3} kg·m²/s)"
    return y, "L_z about COM (kg·m²/s)"


def scale_angular_momentum_origin_display(y_si: np.ndarray) -> tuple[np.ndarray, str]:
    y_plot, lab = scale_angular_momentum_display(y_si)
    if "about COM" in lab:
        return y_plot, lab.replace("L_z about COM", "Total L_z about origin")
    return y_plot, "Total L_z about origin (kg·m²/s)"


def scale_energy_display(y_si: np.ndarray) -> tuple[np.ndarray, str]:
    y = np.asarray(y_si, dtype=np.float64)
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return y, "Newtonian specific energy (J/kg)"
    m = float(np.nanmax(np.abs(finite)))
    if m >= 1e6:
        return y / 1e6, "Newtonian specific energy (MJ/kg)"
    return y, "Newtonian specific energy (J/kg)"


def format_animation_time_caption(
    time_s: float,
    mode: str,
    *,
    preferred_time_unit: str = "auto",
    active_time_unit: str | None = None,
) -> str:
    t = float(time_s)
    if not np.isfinite(t):
        return "t = ?"
    _ = mode
    if active_time_unit is not None:
        tf, tu = _time_from_unit(str(active_time_unit))
    else:
        tf, tu = _resolve_time_unit(preferred_time_unit, t)
    return f"t = {t * tf:.2f} {tu}"
