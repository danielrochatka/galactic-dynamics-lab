"""
Display-only unit selection for plots and renders.

All simulation state, snapshots, and CSV exports remain in SI. This module only
converts values at draw time and chooses axis labels / tick formatting.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# IAU / CODATA-friendly constants (display only)
METERS_PER_LY = 9.460731e15
SECONDS_PER_JULIAN_YEAR = 365.25 * 86400.0


@dataclass(frozen=True)
class SpatialDisplay:
    """Multiply SI position (meters) by `factor` to get axis coordinates."""

    factor: float
    unit: str  # short label: "m", "km", "ly", "kly", "Mly"


def spatial_display_for_xy_plot(mode: str, half_axis_m: float) -> SpatialDisplay:
    """
    Choose x/y axis units for top-down scatter/animation from viewport half-axis (meters).
    """
    ha = float(half_axis_m)
    if not np.isfinite(ha) or ha <= 0:
        ha = 1.0

    if mode == "earth_moon_benchmark":
        return SpatialDisplay(1e-3, "km")

    if mode == "bh_orbit_validation":
        if ha >= 5.0e5:
            return SpatialDisplay(1e-3, "km")
        return SpatialDisplay(1.0, "m")

    if mode == "galaxy":
        return _pick_astronomical_distance_scale(ha)

    # Other modes: heuristic from extent
    if ha >= 0.5 * 1e6 * METERS_PER_LY:
        return SpatialDisplay(1.0 / (1e6 * METERS_PER_LY), "Mly")
    if ha >= 0.5 * 1e3 * METERS_PER_LY:
        return SpatialDisplay(1.0 / (1e3 * METERS_PER_LY), "kly")
    if ha >= 0.5 * METERS_PER_LY:
        return SpatialDisplay(1.0 / METERS_PER_LY, "ly")
    if ha >= 5.0e5:
        return SpatialDisplay(1e-3, "km")
    return SpatialDisplay(1.0, "m")


def _pick_astronomical_distance_scale(half_axis_m: float) -> SpatialDisplay:
    """Galaxy-scale runs: ly / kly / Mly from viewport size; fall back to km/m for compact toy runs."""
    if half_axis_m >= 0.5 * 1e6 * METERS_PER_LY:
        return SpatialDisplay(1.0 / (1e6 * METERS_PER_LY), "Mly")
    if half_axis_m >= 0.5 * 1e3 * METERS_PER_LY:
        return SpatialDisplay(1.0 / (1e3 * METERS_PER_LY), "kly")
    if half_axis_m >= 0.5 * METERS_PER_LY:
        return SpatialDisplay(1.0 / METERS_PER_LY, "ly")
    if half_axis_m >= 5.0e5:
        return SpatialDisplay(1e-3, "km")
    return SpatialDisplay(1.0, "m")


@dataclass(frozen=True)
class SeriesDisplay:
    """Scale SI series to human-readable plot axes (two-body + generic diagnostics)."""

    distance_factor: float
    distance_unit: str
    speed_factor: float
    speed_unit: str
    time_factor: float
    time_unit: str
    suppress_tick_offset: bool


def series_display_for_two_body(mode: str, *, max_distance_m: float) -> SeriesDisplay:
    """Pair separation, COM radius, speeds vs time for earth_moon / bh_orbit."""
    md = float(max_distance_m)
    if not np.isfinite(md) or md <= 0:
        md = 1.0

    if mode == "earth_moon_benchmark":
        return SeriesDisplay(
            distance_factor=1e-3,
            distance_unit="km",
            speed_factor=1e-3,
            speed_unit="km/s",
            time_factor=1.0,
            time_unit="s",
            suppress_tick_offset=True,
        )

    if mode == "bh_orbit_validation":
        if md >= 5.0e5:
            return SeriesDisplay(
                distance_factor=1e-3,
                distance_unit="km",
                speed_factor=1.0,
                speed_unit="m/s",
                time_factor=1.0,
                time_unit="s",
                suppress_tick_offset=True,
            )
        return SeriesDisplay(
            distance_factor=1.0,
            distance_unit="m",
            speed_factor=1.0,
            speed_unit="m/s",
            time_factor=1.0,
            time_unit="s",
            suppress_tick_offset=True,
        )

    raise ValueError(f"series_display_for_two_body: unsupported mode {mode!r}")


def series_display_for_galaxy_diagnostics(
    max_radius_m: float, max_time_s: float
) -> SeriesDisplay:
    """
    Secondary / galaxy diagnostic time series: radii in ly/kly/Mly, time usually seconds.
    """
    mr = float(max_radius_m)
    if not np.isfinite(mr) or mr <= 0:
        mr = 1.0
    mt = float(max_time_s)
    if not np.isfinite(mt) or mt < 0:
        mt = 0.0

    dist = _pick_astronomical_distance_scale(mr)

    # Time: keep seconds unless the run spans many years (galaxy orbits)
    if mt >= 1.0e9 * SECONDS_PER_JULIAN_YEAR:
        tf = 1.0 / (1e6 * SECONDS_PER_JULIAN_YEAR)
        tu = "Myr"
    elif mt >= 1.0e3 * SECONDS_PER_JULIAN_YEAR:
        tf = 1.0 / (1e3 * SECONDS_PER_JULIAN_YEAR)
        tu = "kyr"
    elif mt >= 3600.0 * 24 * 365.25 * 10:
        tf = 1.0 / SECONDS_PER_JULIAN_YEAR
        tu = "yr"
    else:
        tf = 1.0
        tu = "s"

    return SeriesDisplay(
        distance_factor=dist.factor,
        distance_unit=dist.unit,
        speed_factor=1e-3,
        speed_unit="km/s",
        time_factor=tf,
        time_unit=tu,
        suppress_tick_offset=True,
    )


def series_display_generic_validation(max_radius_m: float) -> SeriesDisplay:
    """symmetric_pair, small_n, etc.: m/km/ly heuristics from spatial extent."""
    mr = float(max_radius_m)
    if not np.isfinite(mr) or mr <= 0:
        mr = 1.0
    if mr >= 0.5 * METERS_PER_LY:
        d = _pick_astronomical_distance_scale(mr)
        return SeriesDisplay(
            distance_factor=d.factor,
            distance_unit=d.unit,
            speed_factor=1e-3,
            speed_unit="km/s",
            time_factor=1.0,
            time_unit="s",
            suppress_tick_offset=True,
        )
    if mr >= 5.0e5:
        return SeriesDisplay(
            distance_factor=1e-3,
            distance_unit="km",
            speed_factor=1e-3,
            speed_unit="km/s",
            time_factor=1.0,
            time_unit="s",
            suppress_tick_offset=True,
        )
    return SeriesDisplay(
        distance_factor=1.0,
        distance_unit="m",
        speed_factor=1.0,
        speed_unit="m/s",
        time_factor=1.0,
        time_unit="s",
        suppress_tick_offset=True,
    )


def rotation_curve_display(x_max_m: float) -> tuple[float, float, str, str]:
    """
    Returns (r_factor, v_factor, x_label, y_label) for r (m) and v (m/s) inputs.
    Plotted values are r_plot = r_m * r_factor, v_plot = v_m_s * v_factor.
    """
    xm = float(x_max_m)
    if not np.isfinite(xm) or xm <= 0:
        xm = 1.0
    if xm >= 0.5 * 1e6 * METERS_PER_LY:
        rf = 1.0 / (1e6 * METERS_PER_LY)
        xl = "Distance from galactic center (Mly)"
    elif xm >= 0.5 * 1e3 * METERS_PER_LY:
        rf = 1.0 / (1e3 * METERS_PER_LY)
        xl = "Distance from galactic center (kly)"
    elif xm >= 0.5 * METERS_PER_LY:
        rf = 1.0 / METERS_PER_LY
        xl = "Distance from galactic center (ly)"
    elif xm >= 5.0e8:
        rf = 1e-3
        xl = "Distance from galactic center (km)"
    else:
        rf = 1.0
        xl = "Distance from galactic center (m)"

    # Speed: km/s for galaxy-scale; m/s when distances are small (tests / toy)
    if xm >= 1e9:
        vf = 1e-3
        yl = "Orbital speed (km/s)"
    else:
        vf = 1.0
        yl = "Orbital speed (m/s)"
    return rf, vf, xl, yl


def apply_suppress_tick_offset(ax) -> None:
    """Disable matplotlib's +1eN offset text on both axes (display-layer readability)."""
    import matplotlib.ticker as mticker

    for axis in (ax.xaxis, ax.yaxis):
        axis.get_major_formatter().set_useOffset(False)
        axis.set_major_formatter(mticker.ScalarFormatter(useOffset=False))


def scale_angular_momentum_display(y_si: np.ndarray) -> tuple[np.ndarray, str]:
    """
    Scale L_z (SI) to a compact magnitude with honest y-label suffix.
    """
    y = np.asarray(y_si, dtype=np.float64)
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return y, "L_z (SI)"
    m = float(np.nanmax(np.abs(finite)))
    if m <= 0 or not np.isfinite(m):
        return y, "L_z (SI)"
    exp = int(np.floor(np.log10(m)))
    # Keep exponent a multiple of 3 when large
    if abs(exp) >= 6:
        exp3 = (exp // 3) * 3
        s = 10.0 ** (-exp3)
        return y * s, f"L_z about COM (×10^{exp3} kg·m²/s)"
    return y, "L_z about COM (kg·m²/s)"


def scale_angular_momentum_origin_display(y_si: np.ndarray) -> tuple[np.ndarray, str]:
    """Total L_z about origin — same scaling rule as COM L_z."""
    y_plot, lab = scale_angular_momentum_display(y_si)
    if "about COM" in lab:
        return y_plot, lab.replace("L_z about COM", "Total L_z about origin")
    return y_plot, "Total L_z about origin (kg·m²/s)"


def scale_energy_display(y_si: np.ndarray) -> tuple[np.ndarray, str]:
    """Newtonian specific energy J/kg — use MJ/kg when |values| are huge."""
    y = np.asarray(y_si, dtype=np.float64)
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return y, "Newtonian specific energy (J/kg)"
    m = float(np.nanmax(np.abs(finite)))
    if m >= 1e6:
        return y / 1e6, "Newtonian specific energy (MJ/kg)"
    return y, "Newtonian specific energy (J/kg)"


def format_animation_time_caption(time_s: float, mode: str) -> str:
    """Short time string for animation frame title (SI time in seconds)."""
    t = float(time_s)
    if not np.isfinite(t):
        return "t = ?"
    if mode == "galaxy" and t >= 1000.0 * SECONDS_PER_JULIAN_YEAR:
        return f"t = {t / SECONDS_PER_JULIAN_YEAR:.3g} yr"
    if t >= 1e7:
        return f"t = {t:.4g} s"
    if t >= 1000:
        return f"t = {t:.4g} s"
    return f"t = {t:.6g} s"
