# Physics package template

Copy this folder to create a new physics package (e.g. `MyCustomPhysics/`).

## You must implement

1. **name()** – Return the package identifier string (e.g. `"MyCustomPhysics"`). This is the name users set in config as `physics_package = MyCustomPhysics`.

2. **compute_accelerations(...)** – Compute accelerations for all particles. Signature:
   ```cpp
   void compute_accelerations(const State& state,
                              double bh_mass,
                              double softening,
                              bool star_star,
                              std::vector<double>& ax,
                              std::vector<double>& ay) const override;
   ```
   Fill `ax[i]` and `ay[i]` for each particle `i`. The integrator calls this twice per step (velocity Verlet). `state` has `.x`, `.y`, `.vx`, `.vy`, `.mass` and `.n()`.

## Optional overrides (generic hooks)

- **display_name()** – Human-facing package label.
- **runtime_metadata()** – Generic package runtime key/value metadata.
- **supports_utility_mode(mode)** – Capability check for utility-mode handling.
- **run_utility_mode(config, output_dir)** – Package utility-mode entry point; default returns false.
- **run_info_metadata(config)** – Additional run-info key/value metadata.
- **render_metadata(config)** – Additional render/report key/value metadata.
- **package_config_metadata()** – Package config schema/default metadata.
- **compute_potential_energy(state, bh_mass, softening)** – For diagnostics/validation. Default returns 0.
- **init()** – Called once before the run (e.g. load tables). Default no-op.
- **validation_name()** – For reporting. Default returns `name()`.

## Registration model (current)

Packages register via static self-registration in package translation units by calling
`register_physics_package_factory("PackageName", ...)`.

You should **not** edit `physics/registry.cpp` to add package entries; the registry is factory-based.

## Build inclusion note (current)

The engine build is still source-list based. Add your package `.cpp` files to `engine/Makefile` (`LIB_SRCS`) so they are compiled and linked.

## Optional package defaults

You may add `defaults.cfg` in your package folder with package-specific keys. Users can override from run config. See `physics/TPFCore/defaults.cfg` for an example.
