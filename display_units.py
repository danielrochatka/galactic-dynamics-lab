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
METERS_PER_AU = 1.495978707e11
METERS_PER_PC = 3.085677581e16
METERS_PER_KPC = 3.085677581e19
SECONDS_PER_JULIAN_YEAR = 365.25 * 86400.0
SECONDS_PER_DAY = 86400.0
SECONDS_PER_HOUR = 3600.0
SECONDS_PER_MINUTE = 60.0


@dataclass(frozen=True)
class DisplayUnitPrefs:
    """Parsed from run_info (cpp_sim postprocess); internal simulation stays SI."""

    distance: str  # normalized: auto | m | km | au | ly | pc | kpc
    time: str  # auto | s | min | hr | day | yr | kyr | myr
    velocity: str  # auto | m/s | km/s
    units_in_overlay: bool
    show_unit_reference: bool


def _norm_key(s: str) -> str:
    return str(s or "").strip().lower()


def display_prefs_from_run_info(run_info: dict | None) -> DisplayUnitPrefs:
    if not run_info:
        return DisplayUnitPrefs("auto", "auto", "auto", True, True)
    du = _norm_key(run_info.get("display_distance_unit", "auto"))
    tu = _norm_key(run_info.get("display_time_unit", "auto"))
    vu = _norm_key(run_info.get("display_velocity_unit", "auto"))
    if du not in ("auto", "m", "km", "au", "ly", "pc", "kpc"):
        du = "auto"
    if tu not in ("auto", "s", "min", "hr", "day", "yr", "kyr", "myr"):
        tu = "auto"
    if vu not in ("auto", "m/s", "km/s"):
        vu = "auto"

    def _boolish(v: object, default: bool) -> bool:
        if v is None:
            return default
        if isinstance(v, bool):
            return v
        if isinstance(v, (int, float)):
            return bool(v)
        s = str(v).strip().lower()
        if s in ("0", "false", "no", "off"):
            return False
        if s in ("1", "true", "yes", "on"):
            return True
        return default

    uio = _boolish(run_info.get("display_units_in_overlay"), True)
    sur = _boolish(run_info.get("display_show_unit_reference"), True)
    return DisplayUnitPrefs(du, tu, vu, uio, sur)


def distance_display_factor(unit: str) -> tuple[float, str]:
    """SI meters → display coordinate; label for axes."""
    u = _norm_key(unit)
    if u == "m":
        return 1.0, "m"
    if u == "km":
        return 1e-3, "km"
    if u == "au":
        return 1.0 / METERS_PER_AU, "AU"
    if u == "ly":
        return 1.0 / METERS_PER_LY, "ly"
    if u == "pc":
        return 1.0 / METERS_PER_PC, "pc"
    if u == "kpc":
        return 1.0 / METERS_PER_KPC, "kpc"
    raise ValueError(f"unsupported display distance unit: {unit!r}")


def time_display_factor(unit: str) -> tuple[float, str]:
    """SI seconds → display time; label for axes."""
    u = _norm_key(unit)
    if u == "s":
        return 1.0, "s"
    if u == "min":
        return 1.0 / SECONDS_PER_MINUTE, "min"
    if u == "hr":
        return 1.0 / SECONDS_PER_HOUR, "hr"
    if u == "day":
        return 1.0 / SECONDS_PER_DAY, "day"
    if u == "yr":
        return 1.0 / SECONDS_PER_JULIAN_YEAR, "yr"
    if u == "kyr":
        return 1.0 / (1000.0 * SECONDS_PER_JULIAN_YEAR), "kyr"
    if u == "myr":
        return 1.0 / (1e6 * SECONDS_PER_JULIAN_YEAR), "Myr"
    raise ValueError(f"unsupported display time unit: {unit!r}")


def velocity_display_factor(unit: str) -> tuple[float, str]:
    u = _norm_key(unit)
    if u == "m/s":
        return 1.0, "m/s"
    if u == "km/s":
        return 1e-3, "km/s"
    raise ValueError(f"unsupported display velocity unit: {unit!r}")


def spatial_display_explicit(distance_unit: str) -> SpatialDisplay:
    fac, lab = distance_display_factor(distance_unit)
    return SpatialDisplay(fac, lab)


def spatial_display_from_run_info(
    mode: str, half_axis_m: float, run_info: dict | None
) -> SpatialDisplay:
    prefs = display_prefs_from_run_info(run_info)
    if prefs.distance == "auto":
        return spatial_display_for_xy_plot(mode, half_axis_m)
    return spatial_display_explicit(prefs.distance)


def apply_display_prefs_to_series(
    base: SeriesDisplay,
    prefs: DisplayUnitPrefs,
) -> SeriesDisplay:
    """
    Override auto-chosen series display factors when run_info requests explicit units.
    """
    df, du = base.distance_factor, base.distance_unit
    if prefs.distance != "auto":
        df, du = distance_display_factor(prefs.distance)

    tf, tu = base.time_factor, base.time_unit
    if prefs.time != "auto":
        tf, tu = time_display_factor(prefs.time)

    sf, su = base.speed_factor, base.speed_unit
    if prefs.velocity != "auto":
        sf, su = velocity_display_factor(prefs.velocity)

    return SeriesDisplay(
        distance_factor=df,
        distance_unit=du,
        speed_factor=sf,
        speed_unit=su,
        time_factor=tf,
        time_unit=tu,
        suppress_tick_offset=base.suppress_tick_offset,
    )


def series_display_resolved_for_bulk(
    simulation_mode: str,
    two_body_secondary: bool,
    diagnostics: dict[str, np.ndarray],
    cutoff_radius: float,
    run_info: dict | None,
) -> SeriesDisplay:
    """Bulk radial diagnostics: auto series + optional run_info overrides."""
    max_r = max(
        float(np.max(diagnostics["median_r"])),
        float(np.max(diagnostics["mean_r"])),
        float(np.max(diagnostics["std_r"])),
        float(np.max(diagnostics["max_r"])),
        float(cutoff_radius),
    )
    max_t = float(np.max(diagnostics["time"]))
    if two_body_secondary and simulation_mode == "earth_moon_benchmark":
        base = series_display_for_two_body("earth_moon_benchmark", max_distance_m=max_r)
    elif simulation_mode == "galaxy":
        base = series_display_for_galaxy_diagnostics(max_r, max_t)
    else:
        base = series_display_generic_validation(max_r)

    prefs = display_prefs_from_run_info(run_info)
    return apply_display_prefs_to_series(base, prefs)


def series_display_resolved_for_two_body(
    mode: str,
    *,
    max_distance_m: float,
    run_info: dict | None,
) -> SeriesDisplay:
    base = series_display_for_two_body(mode, max_distance_m=max_distance_m)
    prefs = display_prefs_from_run_info(run_info)
    return apply_display_prefs_to_series(base, prefs)


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
    if ha >= 0.5 * METERS_PER_AU:
        return SpatialDisplay(1.0 / METERS_PER_AU, "AU")
    if ha >= 5.0e5:
        return SpatialDisplay(1e-3, "km")
    return SpatialDisplay(1.0, "m")


def _pick_astronomical_distance_scale(half_axis_m: float) -> SpatialDisplay:
    """Galaxy-scale runs: ly / kly / Mly from viewport size; AU for solar-system-ish; km/m for toy runs."""
    if half_axis_m >= 0.5 * 1e6 * METERS_PER_LY:
        return SpatialDisplay(1.0 / (1e6 * METERS_PER_LY), "Mly")
    if half_axis_m >= 0.5 * 1e3 * METERS_PER_LY:
        return SpatialDisplay(1.0 / (1e3 * METERS_PER_LY), "kly")
    if half_axis_m >= 0.5 * METERS_PER_LY:
        return SpatialDisplay(1.0 / METERS_PER_LY, "ly")
    if half_axis_m >= 0.5 * METERS_PER_AU:
        return SpatialDisplay(1.0 / METERS_PER_AU, "AU")
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


def _rotation_curve_display_auto(x_max_m: float) -> tuple[float, float, str, str]:
    """Auto policy when run_info display keys are all auto."""
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
    elif xm >= 0.5 * METERS_PER_AU:
        rf = 1.0 / METERS_PER_AU
        xl = "Distance from galactic center (AU)"
    elif xm >= 5.0e8:
        rf = 1e-3
        xl = "Distance from galactic center (km)"
    else:
        rf = 1.0
        xl = "Distance from galactic center (m)"

    if xm >= 1e9:
        vf = 1e-3
        yl = "Orbital speed (km/s)"
    else:
        vf = 1.0
        yl = "Orbital speed (m/s)"
    return rf, vf, xl, yl


def rotation_curve_display(
    x_max_m: float, run_info: dict | None = None
) -> tuple[float, float, str, str]:
    """
    Returns (r_factor, v_factor, x_label, y_label) for r (m) and v (m/s) inputs.
    Plotted values are r_plot = r_m * r_factor, v_plot = v_m_s * v_factor.
    """
    prefs = display_prefs_from_run_info(run_info)
    rf, vf, xl, yl = _rotation_curve_display_auto(x_max_m)
    if prefs.distance != "auto":
        rf, du = distance_display_factor(prefs.distance)
        xl = f"Distance from galactic center ({du})"
    if prefs.velocity != "auto":
        vf, vu = velocity_display_factor(prefs.velocity)
        yl = f"Orbital speed ({vu})"
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


def format_animation_time_caption(
    time_s: float,
    mode: str,
    *,
    run_info: dict | None = None,
    max_time_s: float | None = None,
) -> str:
    """
    Short time string for animation frame title (SI time in seconds in memory; caption uses display units).
    """
    prefs = display_prefs_from_run_info(run_info)
    if prefs.time != "auto":
        tf, tu = time_display_factor(prefs.time)
        tv = float(time_s) * tf
        if not np.isfinite(tv):
            return "t = ?"
        return f"t = {tv:.6g} {tu}"

    t = float(time_s)
    if not np.isfinite(t):
        return "t = ?"

    mt = float(max_time_s) if max_time_s is not None and np.isfinite(max_time_s) else t
    if mode == "galaxy" and mt >= 1.0e9 * SECONDS_PER_JULIAN_YEAR:
        return f"t = {t / (1e6 * SECONDS_PER_JULIAN_YEAR):.3g} Myr"
    if mode == "galaxy" and mt >= 1.0e3 * SECONDS_PER_JULIAN_YEAR:
        return f"t = {t / (1e3 * SECONDS_PER_JULIAN_YEAR):.3g} kyr"
    if mode == "galaxy" and t >= 1000.0 * SECONDS_PER_JULIAN_YEAR:
        return f"t = {t / SECONDS_PER_JULIAN_YEAR:.3g} yr"
    if mode == "galaxy" and mt >= 10 * SECONDS_PER_DAY and mt < 1000.0 * SECONDS_PER_JULIAN_YEAR:
        return f"t = {t / SECONDS_PER_DAY:.4g} day"
    if t >= 1e7:
        return f"t = {t:.4g} s"
    if t >= 1000:
        return f"t = {t:.4g} s"
    return f"t = {t:.6g} s"


def display_overlay_time_caption(
    time_s: float, run_info: dict | None, *, simulation_mode: str
) -> str:
    """Overlay line for current frame time (matches animation caption policy)."""
    mode = simulation_mode
    max_t: float | None = None
    if run_info:
        raw = run_info.get("total_simulated_time")
        if raw is not None:
            try:
                max_t = float(raw)
            except (TypeError, ValueError):
                max_t = None
    if max_t is None or not np.isfinite(max_t):
        max_t = float(time_s)
    return format_animation_time_caption(
        time_s, mode, run_info=run_info, max_time_s=max_t
    )


def display_unit_reference_lines(
    run_info: dict | None,
    *,
    spatial_unit: str,
) -> list[str]:
    """Small reference block: active display units (SI data scaled at draw time)."""
    prefs = display_prefs_from_run_info(run_info)
    lines: list[str] = []
    if not prefs.show_unit_reference:
        return lines
    du = prefs.distance if prefs.distance != "auto" else spatial_unit
    tu = prefs.time if prefs.time != "auto" else "auto"
    vu = prefs.velocity if prefs.velocity != "auto" else "auto"
    lines.append(f"distance display: {du}  |  time: {tu}  |  velocity: {vu}")
    lines.append("positions/CSV: SI (m, m/s, s)")
    return lines
