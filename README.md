# Galaxy — TPF simulation repository

This repository hosts **simulation and tooling** around **Triality Primordial Field (TPF)**–style galaxy dynamics: a 2D N-body engine, Python reference code, and post-processing, with a clear split between **manuscript-aligned structure**, **exploratory closures**, and **experiment workflows**.

## What this repo is

- **Theory / manuscript** — TPF structural claims and scope boundaries live in project docs (e.g. paper drafts and `cpp_sim/physics/TPFCore/TPF_PAPER_V11_SCOPE.md`). The code does not replace the paper; it implements testable pieces and labeled stand-ins where the theory is deferred.
- **Simulation** — `cpp_sim/` is the primary **C++** N-body application (Newtonian baseline and pluggable physics packages). A **Python** stack (`main.py`, `config.py`, …) provides a reference implementation, rendering, and validation modes.
- **Experiments** — Long runs, sweeps, and diagnostics are orchestrated via configs and scripts; outputs are designed to be **auditable** (run metadata, manifests, optional plot overlays).

## cpp_sim in one line

`cpp_sim/` is the **simulation application**: build `galaxy_sim`, drive modes from root `configs/`, write `run_info.txt`, snapshots, and optional **`render_manifest.json`** for render-time truth. Details: **[cpp_sim/README.md](cpp_sim/README.md)**.

## Repo map

| Area | Role |
|------|------|
| **`cpp_sim/`** | C++ simulator, physics packages, configs precedence, outputs. |
| **`cpp_sim/physics/TPFCore/`** | TPFCore package: primitive Ξ/Θ/I, readout/VDSG paths. **[cpp_sim/physics/TPFCore/README.md](cpp_sim/physics/TPFCore/README.md)** |
| **`configs/`** | Run configs (root only; see cpp_sim README). `example.cfg` is a reference list. |
| **`plot_cpp_run.py`** | Post-process C++ run directories → galaxy PNGs, animation, rotation curve. |
| **Python (`main.py`, …)** | Reference N-body + plots; not the only path to science results. |

## Paper vs simulation vs experiments

- **Paper** — Defines intended ontology, equations where closed, and what is *not* expanded (∆C, full nonlinear field solve, etc.).
- **Simulation** — Integrates particles with a chosen **physics package**; exposes **resolved** parameters in `run_info` / manifests.
- **Experiments** — Compare branches (e.g. Newtonian control vs TPFCore), tune closures, and record outcomes without conflating **configured** readout with **actual acceleration branch** (see manifests and TPFCore README).

## Quick start (pointers)

| Goal | Start here |
|------|------------|
| Build & run C++, configs, outputs, overlays | [cpp_sim/README.md](cpp_sim/README.md) |
| TPFCore: Θ, I, λ, VDSG, readout vs dynamics | [cpp_sim/physics/TPFCore/README.md](cpp_sim/physics/TPFCore/README.md) |
| Paper vs code scope (v11) | [cpp_sim/physics/TPFCore/TPF_PAPER_V11_SCOPE.md](cpp_sim/physics/TPFCore/TPF_PAPER_V11_SCOPE.md) |
| Example keys (not exhaustive) | [configs/example.cfg](configs/example.cfg) |

For Python-only quick runs, see comments in `configs/example.cfg` and `main.py`; the **authoritative** simulation workflow for this project’s audit trail is documented under **`cpp_sim/README.md`**.
