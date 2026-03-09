# Galaxy N-body simulator

2D N-body galaxy simulation (Newtonian pairwise gravity, velocity Verlet). Python reference implementation with rendering and diagnostics; C++ engine for speed.

## Quick start

**Python** (from repo root):

```bash
pip install numpy matplotlib
python main.py
```

**C++** (faster, no built-in plots):

```bash
cd cpp_sim && make && ./galaxy_sim galaxy
```

Use `plot_cpp_run.py` to generate plots from C++ outputs (see below).

## Configuration

- **Defaults** — Python: `config.py` (`SimulationConfig`). C++: `cpp_sim/config.hpp` (`Config`).
- **Example config** — `configs/example.cfg` lists available options and sensible defaults. It is version-controlled as a reference.
- **Local overrides** — Keep personal or machine-specific settings out of version control:
  - Put local config files in **`configs/local/`** (entire directory is gitignored), or
  - Use any **`configs/*.local.cfg`** file (e.g. `configs/dev.local.cfg`).

Copy `configs/example.cfg` to `configs/my.local.cfg` or `configs/local/my.cfg` and edit as needed. **Python** `main.py` automatically loads (in order) `configs/my.local.cfg`, `configs/local/my.local.cfg`, and any `configs/local/*.cfg`; later files override earlier. Keys match option names (e.g. `n_steps`, `simulation_mode`). On startup, a line like `Config loaded from: ['configs/my.local.cfg']` confirms your file was used.

## Project layout

- **Python**: `main.py`, `config.py`, `simulation.py`, `physics.py`, `init_conditions.py`, `render.py`, `diagnostics.py`, `validation.py` — full run, rendering, validation modes.
- **C++**: `cpp_sim/` — simulation core only; see `cpp_sim/README.md` for build and run.
- **Post-processing**: `plot_cpp_run.py <run_dir>` — reads C++ snapshot CSVs and produces initial/final plots, optional animation and diagnostics.
- **Config**: `configs/example.cfg` (tracked), `configs/local/` and `configs/*.local.cfg` (gitignored).

## Outputs

- Python: `outputs/<run_id>/` (run_id = `YYYYMMDD_HHMMSS`).
- C++: `cpp_sim/outputs/<run_id>/` — `run_info.txt`, `snapshot_*.csv`; then run `python plot_cpp_run.py cpp_sim/outputs/<run_id>` to generate plots.
