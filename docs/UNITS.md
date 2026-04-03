# Units policy: SI internals, configurable display units

## Core rule

- Simulation physics, integration state, and force calculations are in **SI**.
- Snapshot CSV numeric values stay in **SI** (`x,y` in m; `vx,vy` in m/s; time in s).
- Display units are applied only when making plots/renders/videos in postprocess.

## Config keys

Set in run config (`configs/*.cfg`):

- `display_distance_unit = auto | m | km | AU | ly | pc | kpc`
- `display_time_unit = auto | s | min | hr | day | yr | kyr | Myr`
- `display_velocity_unit = auto | m/s | km/s`
- `display_units_in_overlay = true|false`
- `display_show_unit_reference = true|false`

These are written by C++ into `run_info.txt` and read by Python postprocess (`plot_cpp_run.py`, diagnostics, rotation-curve plotting, render overlay).

## Auto mode policy

`auto` is intentionally simple and practical:

- Distance (from plot/render scale): `m` → `km` → `AU` → `ly` → `pc` → `kpc` as extent grows.
- Time (from run span): `s` → `min` → `hr` → `day` → `yr` → `kyr` → `Myr`.
- Velocity (from speed scale): `m/s` for small tests; `km/s` for larger/astro-style runs.

## Output behavior

- Axis labels include units (for example `x (AU)`, `Time (day)`, `Orbital speed (km/s)`).
- Overlays can include active display unit labels.
- Optional unit-reference text can be shown on PNG/MP4 frames.
- Raw CSV and SI run outputs are unchanged.
