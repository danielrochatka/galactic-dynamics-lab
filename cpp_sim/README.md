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

Example:

```bash
./galaxy_sim two_body_orbit
./galaxy_sim symmetric_pair
./galaxy_sim small_n_conservation
./galaxy_sim timestep_convergence
./galaxy_sim galaxy
```

Outputs go to `outputs/<run_id>/` (run_id = `YYYYMMDD_HHMMSS`).

## Output files

- **`run_info.txt`** — dt, n_steps, softening, bh_mass, star_mass, enable_star_star_gravity, total_simulated_time, number_of_snapshots, n_stars, simulation_mode, physics_package.
- **`snapshot_00000.csv`**, **`snapshot_00010.csv`**, … — Per-snapshot state: header `# step,<step>,time,<time>`, then rows `i,x,y,vx,vy,mass`.
- **`validation_timestep_convergence.txt`** — Only for `timestep_convergence`: dt vs final_x, final_y, final_r, L_z, E_drift.

Snapshot CSVs can be loaded in Python for plotting/diagnostics.

## Physics packages

The simulator loads the physics model by **package name** from config. All packages are compiled into the binary; dispatch is by name at runtime.

### Structure

- **`physics/physics_package.hpp`** — Shared interface: package name, `compute_accelerations(...)`, optional `compute_potential_energy`, `init`, `init_from_config`, `validation_name`.
- **`physics/Newtonian/`** — Default package: Newtonian gravity (BH at origin + optional star–star with softening).
- **`physics/TPF/`** — TPF weak-field correspondence package (linearized sector only; see below).
- **`physics/Template/`** — Stub package and README for adding a new package.
- **`physics/registry.cpp`** — Registry: maps package name → implementation. Add new packages here.

Example layout:

```
cpp_sim/physics/
  physics_package.hpp
  registry.cpp
  Newtonian/
  TPF/
  Template/
  MyCustomPhysics/   (your package)
```

### Selecting a package in config

In your config file (e.g. `configs/my.local.cfg`):

```
physics_package = Newtonian
```

- **Default**: If `physics_package` is omitted or empty, **Newtonian** is used.
- If the name is unknown, the program exits with a clear error and lists available packages (Newtonian, TPF, etc.).
- The chosen package name is written to **`run_info.txt`** as `physics_package\t<name>`.

### TPF weak-field package

The **TPF** package implements the paper’s **weak-field correspondence sector only**. It is **not** a full variational TPF solver.

- **Selection**: Set `physics_package = TPF` in your run config.
- **Package defaults**: `physics/TPF/defaults.cfg` defines `tpf_alpha`, `tpf_match_newtonian_scale`, `tpf_softening`. Override from run config if desired.
- **Output**: When TPF is used, `run_info.txt` includes `tpf_alpha`, `tpf_match_newtonian_scale`, `tpf_softening`. Console prints `Physics: TPF weak-field correspondence package`.
- **Scope**: Linearized Poisson-type sector (nabla² phi = alpha×rho). No full nonlinear/dynamic TPF. See `physics/TPF/README.md`.

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
| **Package defaults** | `cpp_sim/physics/<PackageName>/defaults.cfg` | Package-specific options (e.g. TPF: `tpf_alpha`, `tpf_match_newtonian_scale`) |
| **User run config** | `configs/my.local.cfg` (or `configs/local/my.local.cfg`) | Your personal overrides; overrides everything |

### Load precedence

1. **Built-in defaults** (C++ `Config` struct)
2. **Package defaults** — when `physics_package = X` is selected, load `physics/X/defaults.cfg` if it exists
3. **User run config** — overrides all (including package defaults)

Package authors own their package defaults. Users override from the run config when needed.

### Where config files live

- **Run config**: `configs/my.local.cfg`, `configs/local/my.local.cfg` (tried in order). Works from repo root or `cpp_sim/`.
- **Package defaults**: `cpp_sim/physics/Newtonian/defaults.cfg`, `cpp_sim/physics/TPF/defaults.cfg`, etc. Version-controlled with each package.

### Personal local config

1. Copy **`configs/example.cfg`** to **`configs/my.local.cfg`**
2. Edit `configs/my.local.cfg` — it is gitignored
3. Override any setting, including package-specific ones (e.g. `tpf_alpha`)

On startup you'll see:
```
Loaded package defaults: physics/TPF/defaults.cfg
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
