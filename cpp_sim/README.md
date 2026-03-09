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

- **`physics/physics_package.hpp`** — Shared interface: package name, `compute_accelerations(...)`, optional `compute_potential_energy`, `init`, `validation_name`.
- **`physics/Newtonian/`** — Default package: Newtonian gravity (BH at origin + optional star–star with softening).
- **`physics/Template/`** — Stub package and README for adding a new package.
- **`physics/registry.cpp`** — Registry: maps package name → implementation (e.g. `"Newtonian"` → Newtonian package). Add new packages here.

Example layout:

```
cpp_sim/physics/
  physics_package.hpp
  registry.cpp
  Newtonian/
  Template/
  MyCustomPhysics/   (your package)
```

### Selecting a package in config

In your config file (e.g. `configs/my.local.cfg`):

```
physics_package = Newtonian
```

- **Default**: If `physics_package` is omitted or empty, **Newtonian** is used.
- If the name is unknown, the program exits with a clear error and lists that Newtonian is available (and that more can be added in `physics/registry.cpp`).
- The chosen package name is written to **`run_info.txt`** as `physics_package\t<name>`.

### Adding a new package

1. Create a folder under `physics/`, e.g. **`physics/MyCustomPhysics/`**.
2. Implement the **physics package interface** (see `physics/physics_package.hpp` and `physics/Template/`):
   - **Required**: `name()` and `compute_accelerations(state, bh_mass, softening, star_star, ax, ay)`.
   - Optional: `compute_potential_energy`, `init()`, `validation_name()`.
3. **Register** the package in **`physics/registry.cpp`**: `#include "MyCustomPhysics/mycustom.hpp"`, add a static instance, and append it to the `s_packages` array.
4. Set `physics_package = MyCustomPhysics` in config.

Packages are **compiled C++** implementing the shared interface; there are no Python or interpreted plugins. The integrator and simulation loop call only the interface, so behavior is unchanged except which implementation is selected.

## Ported behavior (Python reference)

- **Physics**: By default the **Newtonian** package: BH at origin + optional star–star pairwise, same softening. Other packages can be added under `physics/`.
- **Integrator**: Velocity Verlet; same step formula as Python.
- **Config**: Same meaning of dt, n_steps, snapshot_every, softening, bh_mass, star_mass, n_stars, enable_star_star_gravity.
- **Galaxy ICs**: Disk uniform in area, v_circ from enclosed mass (BH + stars inside r), tangential direction + velocity noise.
- **Validation ICs**: two_body (r0, v0 from speed_ratio×v_circ), symmetric_pair (±a, ±v), small_n (ring, seed 42).

## Configuration

- **Defaults** are in **`config.hpp`** (`Config`). You can override them without recompiling by using a local config file.
- **Config file**: If present, the binary loads **`configs/my.local.cfg`** (or `configs/local/my.local.cfg`). It tries, in order: `configs/my.local.cfg`, `../configs/my.local.cfg`, `configs/local/my.local.cfg`, `../configs/local/my.local.cfg` so it works whether you run from the repo root or from `cpp_sim/`. Keys match the option names (e.g. `n_steps`, `simulation_mode`, `dt`). On startup you’ll see `Config loaded from: ...` when a file was used.
- Copy **`configs/example.cfg`** to **`configs/my.local.cfg`** and edit; that file is gitignored. Command-line mode (e.g. `./galaxy_sim two_body_orbit`) still overrides `simulation_mode` after the file is loaded.

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
