# Galactic Dynamics Lab

**Galactic Dynamics Lab** is an experimental simulation and analysis environment for galactic dynamics, N-body systems, and pluggable physics models. This repository provides tooling to download, configure, run, and audit experiments—whether you need a Newtonian baseline, want to compare alternative force implementations, or explore included packages such as **TPFCore**.

## What this repository is

- **A unified lab / framework** — One codebase hosts the simulation engine, physics packages, configs, analysis scripts, and experiment workflows. The design favors **auditable runs**: resolved parameters in `run_info`, manifests, and diagnostics so results can be traced without guesswork.
- **Pluggable physics** — The C++ engine loads **Newtonian** gravity as a baseline and registers additional **physics packages** at runtime. **TPFCore** (Transformational Physics Framework–aligned field structure and experimental closures) is **one** included package—not the sole reason the repository exists. More packages can be added over time.
- **Simulation + analysis** — `cpp_sim/` is the primary **C++** N-body application. A **Python** stack (`main.py`, `config.py`, …) offers reference implementations, rendering, and validation modes. Post-processing (e.g. `plot_cpp_run.py`) turns run directories into plots and animations.

**Why one repository:** Engine improvements and physics modeling still evolve together. This project keeps **TPFCore inside the same repo** as the engine—no separate checkout or split-package install flow—so integrator, registry, and package code can be iterated jointly.

## cpp_sim in one line

`cpp_sim/` is the **simulation application**: build `galaxy_sim`, drive modes from root `configs/`, write `run_info.txt`, snapshots, and optional **`render_manifest.json`**. Details: **[cpp_sim/README.md](cpp_sim/README.md)**.

## Repo map

| Area | Role |
|------|------|
| **`cpp_sim/`** | C++ N-body engine, physics package registry, configs precedence, outputs. |
| **`cpp_sim/physics/Newtonian/`** | Baseline Newtonian package (default). |
| **`cpp_sim/physics/TPFCore/`** | TPFCore package: primitive Ξ/Θ/I, readout/VDSG paths, manuscript-aligned vs exploratory layers. **[cpp_sim/physics/TPFCore/README.md](cpp_sim/physics/TPFCore/README.md)** |
| **`configs/`** | Run configs (root only; see cpp_sim README). `example.cfg` is a reference list. |
| **`plot_cpp_run.py`** | Post-process C++ run directories → galaxy PNGs, animation, rotation curve. |
| **Python (`main.py`, …)** | Reference N-body + plots; optional path alongside the C++ audit trail. |

## How the pieces fit together

- **Simulation** — Particles are integrated with a chosen **physics package**; the engine exposes **resolved** parameters in `run_info` and manifests.
- **Experiments** — Long runs, sweeps, and diagnostics use configs and scripts; outputs separate **configured** options from **actual** dynamics branches where the engine records both (see **`cpp_sim/README.md`** and, for TPFCore, **`cpp_sim/physics/TPFCore/README.md`**).
- **TPF manuscript and package theory** — Detailed TPF scope, paper-vs-code mapping, and Ξ/Θ/I semantics live under **`cpp_sim/physics/TPFCore/`** (e.g. **`TPF_PAPER_V11_SCOPE.md`**). The root README does not duplicate that material; the code implements testable pieces and labeled stand-ins where theory is deferred.

## Quick start (pointers)

| Goal | Start here |
|------|------------|
| Build & run C++, configs, outputs, overlays | [cpp_sim/README.md](cpp_sim/README.md) |
| TPFCore package: Θ, I, λ, VDSG, readout vs dynamics | [cpp_sim/physics/TPFCore/README.md](cpp_sim/physics/TPFCore/README.md) |
| Paper vs code scope (v11), TPF-specific | [cpp_sim/physics/TPFCore/TPF_PAPER_V11_SCOPE.md](cpp_sim/physics/TPFCore/TPF_PAPER_V11_SCOPE.md) |
| Example keys (not exhaustive) | [configs/example.cfg](configs/example.cfg) |

For Python-only quick runs, see comments in `configs/example.cfg` and `main.py`. The **authoritative** C++ simulation workflow for audit trails is documented under **`cpp_sim/README.md`**.
