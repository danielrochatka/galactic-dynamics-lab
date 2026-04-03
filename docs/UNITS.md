# Units: SI internally, display units on plots

## Internal (simulation)

The C++ simulator (`cpp_sim`) integrates in **SI** (meters, seconds, kilograms, m/s). Snapshot CSV columns (`x`, `y`, `vx`, `vy`, `time`, …) and numeric physics fields in **`run_info.txt`** (e.g. `dt`, `galaxy_radius`, `bh_mass`, `softening`) are **SI** unless a key explicitly documents otherwise.

This does **not** change when you set display-unit options below.

## Display (postprocess only)

Plots, PNGs, MP4s, and diagnostics produced by **`plot_cpp_run.py`**, **`render.py`**, **`render_overlay.py`**, **`diagnostics.py`**, and **`plot_rotation_curve.py`** scale values **at draw time** so axes and captions show human-friendly numbers. Raw arrays on disk are unchanged.

## Config keys (written to `run_info.txt`)

Set these in your C++ run config (same pattern as `render_overlay_mode` or `diagnostic_cutoff_radius`):

| Key | Purpose |
|-----|---------|
| `display_distance_unit` | Axis scaling for positions / radii in plots |
| `display_time_unit` | Time axis and animation/overlay time captions |
| `display_velocity_unit` | Speed axes where applicable |
| `display_units_in_overlay` | If true (default), overlay lists which display units are active |
| `display_show_unit_reference` | If true (default), overlay adds a short SI vs display reminder |

### Allowed values

- **Distance:** `auto`, `m`, `km`, `AU`, `ly`, `pc`, `kpc`
- **Time:** `auto`, `s`, `min`, `hr`, `day`, `yr`, `kyr`, `Myr`
- **Velocity:** `auto`, `m/s`, `km/s`
- **Bools:** `true`/`false`, `1`/`0`, `yes`/`no` (as elsewhere in config)

Case is normalized (`AU` and `au` both accepted in config; run_info may store lowercase).

### Auto mode (default)

`auto` picks a sensible scale from the **viewport** (scatter/animation), **diagnostic extent**, or **run duration** heuristics:

- **Distance:** Earth–Moon style runs use **km**; compact tests use **m** or **km**; solar-system–ish half-axes use **AU**; galactic scales use **ly** / **kly** / **Mly**; explicit config overrides this.
- **Time:** Short runs use **s**; multi-day spans may use **day**; long galaxy runs use **yr**, **kyr**, or **Myr** from the simulated time span.
- **Velocity:** Two-body and galaxy series follow existing defaults (e.g. **km/s** for typical Earth–Moon and galaxy speeds); explicit `velocity` overrides.

## Related files

- **`display_units.py`** — Scaling factors and labels
- **`cpp_sim/config.hpp`**, **`cpp_sim/config.cpp`**, **`cpp_sim/output.cpp`** — Parse and echo keys into **`run_info.txt`**
