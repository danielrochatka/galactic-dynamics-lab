# Galaxy N-body simulator (C++)

Phase 1 port of the Python galaxy simulator. Same Newtonian pairwise gravity and velocity Verlet integrator; no rendering (use Python for plots).

## Build

```bash
cd cpp_sim
make
```

Requires C++11 (e.g. `g++` or `clang++`). No external dependencies.

## Run

From the `cpp_sim/` directory:

```bash
./galaxy_sim [mode]
```

**Modes** (default: `galaxy` if no argument):

| Mode | Description |
|------|-------------|
| `galaxy` | Full disk (n_stars, inner/outer radius, enclosed-mass ICs). |
| `two_body_orbit` | One star orbiting fixed BH (validation). |
| `symmetric_pair` | Two mirrored stars (validation). |
| `small_n_conservation` | Small N-body, n=3..10 (validation). |
| `timestep_convergence` | Two-body at dt, dt/2, dt/4 (validation). |
| `tpf_single_source_inspect` | TPFCore: one source at origin, probe Theta/I along +x (requires `physics_package = TPFCore`). |
| `tpf_symmetric_pair_inspect` | TPFCore: symmetric pair at (±d,0), probe Theta/I along +x and +y (requires `physics_package = TPFCore`). |
| `tpf_single_source_optimize_c` | TPFCore: sweep c, fit against field-equation residual. Exploratory ansatz-tuning; fitted c is NOT a final constant. |

Example (Newtonian dynamics):

```bash
./galaxy_sim two_body_orbit
./galaxy_sim symmetric_pair
./galaxy_sim galaxy
```

Example (TPFCore inspection):

```bash
# In config: physics_package = TPFCore
./galaxy_sim tpf_single_source_inspect
./galaxy_sim tpf_symmetric_pair_inspect
```

Outputs go to `outputs/<run_id>/` (run_id = `YYYYMMDD_HHMMSS`).

## Output files

**Dynamical modes (Newtonian):**

- **`run_info.txt`** — dt, n_steps, softening, bh_mass, star_mass, enable_star_star_gravity, total_simulated_time, number_of_snapshots, n_stars, simulation_mode, physics_package.
- **`snapshot_00000.csv`**, **`snapshot_00010.csv`**, … — Per-snapshot state: header `# step,<step>,time,<time>`, then rows `i,x,y,vx,vy,mass`.
- **`validation_timestep_convergence.txt`** — Only for `timestep_convergence`: dt vs final_x, final_y, final_r, L_z, E_drift.

**TPFCore inspection modes** (`tpf_single_source_inspect`, `tpf_symmetric_pair_inspect`):

- **`theta_profile.csv`** — radius, x, y, xi_x, xi_y, theta_xx, theta_xy, theta_yy, theta_trace, invariant_I, residual_x, residual_y, residual_norm (symmetric pair adds `axis` column: x or y).
- **`invariant_profile.csv`** — radius, invariant_I, theta_trace, residual_norm (symmetric pair adds `axis` column).
- **`field_summary.txt`** — geometry, source positions, isotropic correction coefficient c, max residual magnitudes, symmetry checks (e.g. residual_y ≈ 0 on +x axis), provisional-ansatz flag.

**TPFCore c-sweep utility** (`tpf_single_source_optimize_c`): Exploratory ansatz-tuning. Numerically fits c against field-equation residual. **Fitted c is NOT a final paper-derived constant.**

- **`c_sweep.csv`** — c, max_residual_norm, mean_residual_norm, l2_residual_norm
- **`c_sweep_summary.txt`** — sweep range, steps, chosen objective, best c, best objective value

See `physics/TPFCore/README.md` for how to run the c sweep and what the objective metrics mean.

Snapshot CSVs can be loaded in Python for plotting/diagnostics.

## Physics packages

The simulator loads the physics model by **package name** from config. All packages are compiled into the binary; dispatch is by name at runtime.

**Available packages:** Newtonian, TPFCore. The old misleading TPF package (which collapsed to Newtonian-style scalar potential) has been **removed**.

### Structure

- **`physics/physics_package.hpp`** — Shared interface: package name, `compute_accelerations(...)`, optional `compute_potential_energy`, `init`, `init_from_config`, `validation_name`.
- **`physics/Newtonian/`** — Default package: Newtonian gravity (BH at origin + optional star–star with softening).
- **`physics/TPFCore/`** — Primitive TPF structure (Xi, Theta, I). Hessian-based provisional ansatz. Inspection-first; see below.
- **`physics/Template/`** — Stub package and README for adding a new package.
- **`physics/registry.cpp`** — Registry: maps package name → implementation. Add new packages here.

Example layout:

```
cpp_sim/physics/
  physics_package.hpp
  registry.cpp
  Newtonian/
  TPFCore/
  Template/
  MyCustomPhysics/   (your package)
```

### Selecting a package in config

In your config file (e.g. `configs/my.local.cfg`):

```
physics_package = Newtonian
```

- **Default**: If `physics_package` is omitted or empty, **Newtonian** is used.
- If the name is unknown, the program exits with a clear error and lists available packages (Newtonian, TPFCore).
- The chosen package name is written to **`run_info.txt`** as `physics_package\t<name>`.

### TPFCore (inspection-first)

**Two usage modes:**

1. **Inspection-only** (default): Set `physics_package = TPFCore`. Run inspection/utility modes:
   - `./galaxy_sim tpf_single_source_inspect` — one source at origin, probe along +x
   - `./galaxy_sim tpf_symmetric_pair_inspect` — sources at (±d,0), probe along +x and +y
   - `./galaxy_sim tpf_single_source_optimize_c` — sweep c, fit against residual (exploratory)
   - Dynamical modes (galaxy, two_body_orbit, etc.) will fail with a clear message.

2. **With provisional readout** (exploratory): Add `tpfcore_enable_provisional_readout = true` to config. Allows dynamical modes: `two_body_orbit`, `symmetric_pair`, `small_n_conservation`, `galaxy`. Motion is tensor-driven (Theta·r_hat), **not** Newtonian −grad(Φ). This is **NOT** the full derived TPF dynamics—see `physics/TPFCore/README.md`.

**Outputs:** Inspection: `theta_profile.csv`, `invariant_profile.csv`, `field_summary.txt`. C-sweep: `c_sweep.csv`, `c_sweep_summary.txt`. Dynamical: `run_info.txt`, `snapshot_*.csv`. When `tpfcore_dump_readout_debug=true`: `tpf_readout_debug.csv` with per-snapshot readout diagnostics (a_radial, a_inward, a_tangential, Theta, etc.) to diagnose runaway vs bound behavior. See `physics/TPFCore/README.md`.

### Adding a new package

1. Create a folder under `physics/`, e.g. **`physics/MyCustomPhysics/`**.
2. Implement the **physics package interface** (see `physics/physics_package.hpp` and `physics/Template/`):
   - **Required**: `name()` and `compute_accelerations(state, bh_mass, softening, star_star, ax, ay)`.
   - Optional: `compute_potential_energy`, `init()`, `init_from_config`, `validation_name()`.
3. **Register** the package in **`physics/registry.cpp`**: `#include "MyCustomPhysics/mycustom.hpp"`, add a static instance, and append it to the `s_packages` array.
4. Add **`physics/MyCustomPhysics/defaults.cfg`** for package-specific options (optional; see package config system).
5. Set `physics_package = MyCustomPhysics` in run config.

Packages are **compiled C++** implementing the shared interface; there are no Python or interpreted plugins. The integrator and simulation loop call only the interface, so behavior is unchanged except which implementation is selected.

## Ported behavior (Python reference)

- **Physics**: By default the **Newtonian** package: BH at origin + optional star–star pairwise, same softening. Other packages can be added under `physics/`.
- **Integrator**: Velocity Verlet; same step formula as Python.
- **Config**: Same meaning of dt, n_steps, snapshot_every, softening, bh_mass, star_mass, n_stars, enable_star_star_gravity.
- **Galaxy ICs**: Disk uniform in area, v_circ from enclosed mass (BH + stars inside r), tangential direction + velocity noise.
- **Validation ICs**: two_body (r0, v0 from speed_ratio×v_circ), symmetric_pair (±a, ±v), small_n (ring, seed 42).

## Configuration

The simulator uses a **layered config system**: built-in defaults → package defaults → user run config. No recompilation needed.

### Shared run config vs package-local defaults

| Layer | Location | Purpose |
|-------|----------|---------|
| **Built-in defaults** | `Config` in `config.hpp` | Fallback when no file supplies a value |
| **Package defaults** | `cpp_sim/physics/<PackageName>/defaults.cfg` | Package-specific options (e.g. TPFCore: `tpfcore_probe_radius_min`, etc.) |
| **User run config** | `configs/my.local.cfg` (or `configs/local/my.local.cfg`) | Your personal overrides; overrides everything |

### Config precedence

**Order (lowest to highest):** built-in defaults → package defaults → run config → **CLI override** (if applicable).

1. **Built-in defaults** — `Config` in `config.hpp`
2. **Package defaults** — `physics/<PackageName>/defaults.cfg` (chosen by `physics_package` probed from run config)
3. **Run config** — first existing of `configs/my.local.cfg`, `../configs/my.local.cfg`, etc. (depends on current working directory)
4. **CLI override** — `./galaxy_sim <mode>` overrides `simulation_mode` only

On startup the simulator prints **Run config selected: …**, **CLI override applied: …** (if any), and a **Resolved config** banner (RUN CONFIG, PACKAGE DEFAULTS, PHYSICS PACKAGE, SIMULATION MODE, OUTPUT DIR, and key values). The same resolved block is written at the top of **run_info.txt**. Use these to confirm which config and mode were actually used.

### Where config files live

- **Run config**: `configs/my.local.cfg`, `configs/local/my.local.cfg` (tried in order). Works from repo root or `cpp_sim/`.
- **Package defaults**: `cpp_sim/physics/Newtonian/defaults.cfg`, `cpp_sim/physics/TPFCore/defaults.cfg`, etc. Version-controlled with each package.

### Personal local config

1. Copy **`configs/example.cfg`** to **`configs/my.local.cfg`**
2. Edit `configs/my.local.cfg` — it is gitignored
3. Override any setting, including package-specific ones (e.g. `tpfcore_probe_radius_max`)

On startup you'll see:
```
Loaded package defaults: physics/TPFCore/defaults.cfg
Loaded run config: configs/my.local.cfg
```

### Adding a new package with defaults

1. Create `cpp_sim/physics/MyPackage/defaults.cfg` with your package's keys
2. Register the package in `physics/registry.cpp`
3. Users select it with `physics_package = MyPackage` in the run config
4. Your package defaults apply first; run config overrides them

- Command-line mode (e.g. `./galaxy_sim two_body_orbit`) overrides `simulation_mode` after the config is loaded.

## Not ported (use Python)

- Matplotlib rendering, diagnostic plots, MP4/GIF.
- Command-line config overrides (C++ uses defaults; edit `config.hpp` / `Config` to change).

## Compare C++ vs Python

Same mode and defaults should give the same results to floating-point tolerance. Example (two_body, 5000 steps, dt=0.01):

```bash
# C++
./galaxy_sim two_body_orbit
# Inspect outputs/YYYYMMDD_HHMMSS/snapshot_05000.csv

# Python (from repo root)
python main.py  # with simulation_mode = "two_body_orbit" in config
# Or: python -c "from config import SimulationConfig; from validation import run_two_body_orbit; ..."
```

Compare final x, y, r, L_z and run_info values.
