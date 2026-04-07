# engine — C++ simulation runtime

`engine/` contains the **authoritative simulation runtime** for Galactic Dynamics Lab: integrator, mode dispatch, package dispatch, config resolution, and structured run outputs.

This is the only simulator runtime in the repository. Python is post-processing/tooling only.

## What lives in `engine/`

- `main.cpp`, `simulation.cpp`: runtime entrypoint and stepping loop.
- `config.*`, `scenario_defaults.cpp`: config model and scenario resolution.
- `physics/`: compiled runtime physics packages (for example, `Newtonian`, `TPFCore`).
- `tests/`: C++ unit/regression/integration test harnesses.

## Build

From repo root:

```bash
(cd engine && make)
```

Requires C++11 (`g++` or `clang++`).

Windows users should follow the WSL-focused setup guide: [`../docs/WINDOWS.md`](../docs/WINDOWS.md).

## Run

From repo root:

```bash
# default mode (galaxy)
(cd engine && ./galaxy_sim)

# explicit mode
(cd engine && ./galaxy_sim earth_moon_benchmark)

# mode + overrides
(cd engine && ./galaxy_sim galaxy --num_particles=2000 --time_step=1e13)
```

Notes:
- `simulation_mode` defaults to `galaxy` if omitted.
- CLI `--key=value` overrides config values.
- Unknown modes print the supported mode list.

### Compare mode (easy entry)

Use package compare keys in your run config to run side-by-side package comparisons from the same C++ runtime/IC seed:

```ini
physics_package = Newtonian
physics_package_compare = TPFCore
```

## Outputs

Runs write to root `outputs/<run_id>/` (when running from `engine/`, this resolves to `../outputs/<run_id>/`).

Common files:
- `run_info.txt` — resolved config + runtime branch labels.
- `resolved_scenario.txt` / `resolved_scenario.json` — effective scenario after resolution.
- `snapshot_*.csv` — particle snapshots.
- `render_manifest.json` (mode/config dependent) — render-time audit payload.

## Config precedence

Operational config surfaces:
- Built-in defaults: `engine/config.hpp` (`Config`).
- Package defaults: `engine/physics/<Package>/defaults.cfg`.
- Run config files: root `configs/` (for example `configs/example.cfg`).
- CLI overrides: `--key=value`.
- Scenario resolution by mode: `resolve_scenario` in `scenario_defaults.cpp`.

Effective precedence used at runtime:
1. Built-in defaults
2. Package defaults
3. Run config + CLI overrides
4. Scenario resolution (mode-aware effective IC/control values)

The resolved block is printed at startup and written to `run_info.txt`.

## Plot a run

The C++ binary does not render images/videos itself. Plot after the run from repo root:

```bash
python3 plot_cpp_run.py outputs/<run_id>
```

## Physics packages (high level)

A physics package is a compiled C++ module selected at runtime by name (`physics/registry.cpp`).

Current built-in packages:
- `Newtonian` (default orbital/gravity baseline)
- `TPFCore` (TPF-based runtime package)

To add a package: implement `physics/physics_package.hpp`, register a factory, and optionally provide `physics/<Name>/defaults.cfg`.

## TPFCore theory details

This README keeps TPFCore details intentionally short. For Xi/Theta/I structure, readout/VDSG paths, and package-specific diagnostics, use:

- [`physics/TPFCore/README.md`](physics/TPFCore/README.md)

## Testing from `engine/`

```bash
# all engine tests (doctest + integration scripts)
(cd engine && make test)

# doctest binary only
(cd engine && make test_unit)

# integration scripts only
(cd engine && make test_integration)
```

For repository-wide testing guidance and Python test commands, see [`../docs/TESTING.md`](../docs/TESTING.md).
