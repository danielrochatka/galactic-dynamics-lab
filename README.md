# Galactic Dynamics Lab

Galactic Dynamics Lab is a galactic dynamics simulation repository with a **single authoritative runtime**: the native C++ engine in [`engine/`](engine/).

## Quickstart (first run)

From the repository root:

```bash
# 1) build
(cd engine && make)

# 2) run (uses root configs; default mode is galaxy)
(cd engine && ./galaxy_sim)

# 3) plot a completed run
python3 plot_cpp_run.py outputs/<run_id>

# 4) test
./run_tests.sh
```

> Tip: each run writes `run_info.txt` (and usually `render_manifest.json`) under `outputs/<run_id>/`.

## Repository layout

| Path | Purpose |
|---|---|
| `engine/` | **Authoritative simulator runtime** (build, stepping, package dispatch, manifests). |
| `configs/` | **Authoritative run configuration surface** for operations. |
| `outputs/` | **Authoritative run output surface** for generated results. |
| `docs/` | Project docs (including test policy and test entrypoints). |
| `tools/` *(optional)* / root Python scripts | Tooling only: plotting, diagnostics, analysis, migration helpers. |
| `python_tests/` | Python tests for tooling and post-processing (not a simulation runtime). |

## Runtime truth (important)

- There is **one** simulator runtime in this repo: `engine/` (C++).
- Physics runtime packages are compiled C++ packages under `engine/physics/`.
- Python is **tooling only** (plotting, diagnostics, analysis, tests).
- There is no legacy Python simulation runtime path.

## Configs and outputs

- Put run configs in root `configs/` (for example, `configs/example.cfg`).
- Runs write outputs to root `outputs/<run_id>/`.
- If running from `engine/`, default output path is `../outputs/<run_id>/` (same root `outputs/` directory).

## Where to go next

- Engine build/run/config details: [`engine/README.md`](engine/README.md)
- TPFCore package theory + package-specific diagnostics: [`engine/physics/TPFCore/README.md`](engine/physics/TPFCore/README.md)
- Testing layers and commands: [`docs/TESTING.md`](docs/TESTING.md)
- Example run config: [`configs/example.cfg`](configs/example.cfg)
