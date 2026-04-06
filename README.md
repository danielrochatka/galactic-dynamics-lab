# Galactic Dynamics Lab

Galactic Dynamics Lab now has **one simulator runtime**: the native C++ engine in `engine/`.

## Repository layout

| Path | Purpose |
|---|---|
| `engine/` | Authoritative simulator runtime (build, stepping, package dispatch, manifests). |
| `configs/` | Root-level run configuration files (authoritative run-config surface). |
| `outputs/` | Root-level run outputs (authoritative output surface). |
| `docs/` | Project documentation. |
| `tools/` *(optional)* / root Python scripts | Plotting, diagnostics, analysis, and migration helpers. |
| `python_tests/` | Python tests for tooling/post-processing. |

## Runtime model

- There is no Python simulation runtime path in this repository anymore.
- Physics runtime packages remain compiled C++ packages under `engine/physics/`.
- Python is tooling-only: plotting, diagnostics, analysis, and tests.

## Quick start

- Build/run engine: [`engine/README.md`](engine/README.md)
- TPFCore package docs: [`engine/physics/TPFCore/README.md`](engine/physics/TPFCore/README.md)
- Config examples: [`configs/example.cfg`](configs/example.cfg)

For each run, inspect `run_info.txt` and (when enabled) `render_manifest.json` in the corresponding `outputs/<run_id>/` directory.
