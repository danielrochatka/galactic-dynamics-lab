# Galactic Dynamics Lab

**Galactic Dynamics Lab** is an experimental simulation and analysis environment for galactic dynamics, N-body systems, and pluggable physics packages.

## What this repository is

- **Unified framework** — One codebase hosts the simulation engine, package registry, configs, diagnostics, and plotting scripts.
- **Pluggable physics** — The C++ engine supports multiple packages (for example `Newtonian` and `TPFCore`) through package dispatch.
- **Simulation + analysis** — `cpp_sim/` provides the main C++ simulator; Python scripts handle post-processing and plotting.

## Repo map

| Area | Role |
|------|------|
| **`cpp_sim/`** | C++ N-body engine, package registry, config parsing, outputs. |
| **`cpp_sim/physics/Newtonian/`** | Baseline Newtonian package. |
| **`cpp_sim/physics/TPFCore/`** | TPF package implementation and package-local documentation. |
| **`configs/`** | Run configs used by the simulator. |
| **`plot_cpp_run.py`** | Post-process C++ run directories into plots and animations. |
| **Python (`main.py`, …)** | Reference Python tools, helpers, and diagnostics. |

## Quick start

| Goal | Start here |
|------|------------|
| Build and run C++ simulator | [cpp_sim/README.md](cpp_sim/README.md) |
| TPF package docs | [cpp_sim/physics/TPFCore/README.md](cpp_sim/physics/TPFCore/README.md) |
| Example config keys | [configs/example.cfg](configs/example.cfg) |

For audit outputs, inspect `run_info.txt` and (when enabled) `render_manifest.json` in each run directory.
