# cpp_sim — Galaxy N-body application

C++ **simulation framework**: velocity Verlet, pluggable **physics packages**, layered configuration, and structured outputs. **Rendering** is not built into the binary; use **`plot_cpp_run.py`** from the repo root on a finished run directory.

This README describes the **application** (build, config, outputs, plotting). **TPFCore-specific theory** (Ξ, Θ, I, readout vs VDSG) lives in **[physics/TPFCore/README.md](physics/TPFCore/README.md)**.

---

## Build

```bash
cd cpp_sim
make
```

Requires C++11 (`g++` or `clang++`). No external libraries.

---

## Run

From **`cpp_sim/`**:

```bash
./galaxy_sim [simulation_mode]
```

`simulation_mode` defaults to **`galaxy`** if omitted. CLI overrides **`simulation_mode`** after loading config. Additional overrides use **`--key=value`** (see startup banner for resolved config).

**Modes** (high level):

| Mode | Purpose |
|------|---------|
| `galaxy` | Disk N-body (production morphology runs). |
| `two_body_orbit`, `symmetric_pair`, `small_n_conservation`, `timestep_convergence` | Validation / checks. |
| `tpf_single_source_inspect`, `tpf_symmetric_pair_inspect`, … | Package-specific **inspection** (require `physics_package = TPFCore`). |
| `tpf_two_body_sweep`, `tpf_weak_field_calibration`, `tpf_newtonian_force_compare`, `tpf_bound_orbit_sweep`, `tpf_diagnostic_consistency_audit` | TPFCore **sweeps / calibration / diagnostics** (see binary help and `SimulationMode` in `config.hpp`). |

Full mode list and errors for unknown modes are printed by the binary. For **TPFCore**, **dynamical** modes need **`tpfcore_enable_provisional_readout = true`**; inspection and sweep modes have their own requirements — see the TPFCore package README.

Outputs go under **`output_dir`** (default `outputs/<run_id>/`, `run_id` often `YYYYMMDD_HHMMSS`).

---

## Configuration

### Where files live

| Layer | Location |
|-------|----------|
| Built-in defaults | `config.hpp` (`Config`). |
| Package defaults | **`physics/<Package>/defaults.cfg`** only (e.g. `physics/TPFCore/defaults.cfg`). |
| Run config | **Repository root `configs/`** only — e.g. `../configs/my.local.cfg` when cwd is `cpp_sim/`. **Not** `cpp_sim/configs/`. |

**Precedence:** built-in → package defaults → root run config → CLI (`simulation_mode` and `--key=value`).

Startup prints **Run config selected**, **Loaded package defaults**, and a **Resolved config** summary. The same resolved block appears at the top of **`run_info.txt`**.

### Physics packages

Dispatch is by name at runtime (see **`physics/registry.cpp`**).

| Package | Role |
|---------|------|
| **Newtonian** | Default: BH at origin + optional star–star pairwise gravity with softening. |
| **TPFCore** | Primitive TPF field structure + optional provisional readout / VDSG dynamics path. **Details:** [physics/TPFCore/README.md](physics/TPFCore/README.md). |

Set in run config:

```ini
physics_package = Newtonian
# or
physics_package = TPFCore
```

Adding a new package: implement **`physics/physics_package.hpp`**, register in **`registry.cpp`**, optional **`physics/<Name>/defaults.cfg`**. See **`physics/Template/`** for a stub.

### Galaxy initialization (templates)

**Galaxy mode** uses a **named template** plus optional noise and structured seeds (`galaxy_init_*` keys in `config.hpp`). Valid **`galaxy_init_template`** names match **`galaxy_init.hpp`**: `symmetric_disk`, `symmetric_disk_noisy`, `clumpy_disk`, `weak_m2`, `weak_m3`, `weak_bar`, `preformed_spiral`. Resolved initialization is logged in **`run_info.txt`** (galaxy init audit block) and **`galaxy_init_diagnostics.txt`**.

**RNG:** `galaxy_init_seed` controls reproducibility.

This is **simulation-app** configuration, not TPF manuscript content; see parameter comments in **`config.hpp`**.

---

## Outputs (dynamical runs)

Typical **`outputs/<run_id>/`** contents:

| File | Content |
|------|---------|
| **`run_info.txt`** | Tab-separated resolved config, counters, **`active_dynamics_branch`**, **`active_metrics_branch`**, **`acceleration_code_path`**, TPF/IC keys when applicable. |
| **`render_manifest.json`** / **`render_manifest.txt`** | Full resolved audit for renders (galaxy mode when `save_run_info`): branches, coupling, cooling, IC parameters, aliases note. |
| **`snapshot_*.csv`** | State: `# step,…,time,…` then `i,x,y,vx,vy,mass`. |
| **`galaxy_init_diagnostics.txt`** | Initial radial / speed / L_z summaries (galaxy mode). |

**TPFCore** dynamical runs may additionally write readout/regime/trajectory diagnostics when enabled (see TPFCore README).

---

## Renders, overlays, and post-processing

The binary **does not** draw PNG/MP4. After a run:

```bash
# from repo root (use python3 or your venv’s python)
python3 plot_cpp_run.py cpp_sim/outputs/<run_id>
```

This produces e.g. **`galaxy_initial.png`**, **`galaxy_final.png`**, optional **`galaxy.mp4`**, **`rotation_curve.png`**, and may honor **`render_overlay_mode`** from **`run_info`** (or **`--render-overlay-mode`**): **`none`** | **`minimal`** | **`audit_full`**. Overlays read **`run_info`** and **`render_manifest.json`** so plots show **dynamics vs metrics** and coupling without opening source.

**`plot_animation_dynamic_zoom`** in run config controls animation viewport smoothing (see `plot_cpp_run.py` docstring).

---

## Application vs physics package

| Concern | Lives in |
|---------|----------|
| Integrator, snapshots, `run_info`, manifests, IC templates, registry | **This README / `config.hpp` / `main.cpp`** |
| Ξ, Θ, I, λ, readout closures, VDSG vs tensor path, paper alignment | **`physics/TPFCore/README.md`** |

---

## Python reference vs C++

The repo includes a **Python** simulator (`main.py`, …) for reference and validation. For comparable modes, results should agree within floating-point tolerance. **C++** is the path used for large galaxy runs and audit manifests.

---

## Compare C++ vs Python (sanity)

Same `simulation_mode` and aligned parameters: compare final **`snapshot_*.csv`** and **`run_info`** fields. Run C++ from `cpp_sim/`; run Python from repo root per `main.py` / `config.py`.
