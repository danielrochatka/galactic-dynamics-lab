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

## Optional overrides

- **compute_potential_energy(state, bh_mass, softening)** – For diagnostics/validation. Default returns 0.
- **init()** – Called once before the run (e.g. load tables). Default no-op.
- **validation_name()** – For reporting. Default returns `name()`.

## Wiring your package into the binary

1. Add your header to `physics/registry.cpp`:  
   `#include "MyCustomPhysics/mycustom.hpp"`
2. Instantiate your package (e.g. a static instance) and add it to the `s_packages` array in `registry.cpp`.

Then set `physics_package = MyCustomPhysics` in your config.
