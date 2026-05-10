# Pushback Artifact: TPF package isolation reset (2026-05-10)

## Exact blocker
The current engine architecture hard-binds TPF behavior into the generic `Config` type and `SimulationMode` enum, and that generic state is consumed across core orchestration layers (`main`, `simulation`, `output`, `render_audit`, `force_compare`, `galaxy_init`, `resolved_scenario`, and generic unit tests). A safe one-shot move to package-owned ownership requires coordinated interface redesign touching parser, CLI mode routing, runtime hooks, and serialization contracts simultaneously.

Proceeding incrementally in this PR without that explicit architectural decision would create high risk of:
- config compatibility breakage for existing TPF config files,
- run_info/render metadata schema drift,
- behavior drift in TPF utility/benchmark modes,
- regressions in Newtonian/default paths due to shared generic structs.

## Exact files/functions preventing full isolation
### Generic config/parser ownership currently containing TPF specifics
- `engine/config.hpp`
  - `SimulationMode` contains TPF-only enum values.
  - `Config` owns large sets of `tpf_*`/Xi/Theta/VDSG fields.
  - `is_tpf_utility_mode(...)` and TPF-diagnostic function declarations.
- `engine/config.cpp`
  - `apply_config_kv(...)` hardcodes TPF keys/validation/defaulting.
  - `parse_simulation_mode(...)` maps TPF mode names.
  - `serialize_config_kv(...)` emits TPF keys.

### Generic orchestration code branching on TPF semantics
- `engine/main.cpp`
  - mode and package-specific runtime dispatch includes TPF-specific branches.
- `engine/simulation.cpp`
  - TPF/Xi-specific runtime policy gates (e.g. `xi_kernel_deformation_active`).
- `engine/output.cpp`
  - run_info/output metadata keys and logic tied to TPF config fields.
- `engine/render_audit.cpp`
  - render/audit emissions include TPF-specific handling.
- `engine/force_compare.cpp` + `engine/force_compare.hpp`
  - compare behavior includes TPF-targeted naming/logic.
- `engine/galaxy_init.cpp` + `engine/galaxy_init.hpp` + `engine/init_conditions.hpp`
  - comments/paths and initialization flow contain TPF-owned behavior and assumptions.
- `engine/resolved_scenario.cpp`
  - effective dynamics/mode resolution hardcodes TPF mode reinterpretation.
- `engine/accel_pipeline_stats.hpp`
  - runtime stats ownership currently treated as generic while carrying TPF semantics.

### Generic tests still asserting TPF semantics
- `engine/tests/unit/test_config.cpp` (heavy TPF key/mode parsing/serialization assertions).
- `engine/tests/unit/test_tpf_core_direct_tpf.cpp`
- `engine/tests/unit/test_tpf_core_xi_kernel_deformed.cpp`
- `engine/tests/unit/test_render_audit.cpp` / `test_compare_orchestration.cpp` (TPF-coupled expectations).
- integration scripts under `engine/tests/integration/` with TPF-specific expectations in generic test tree.

## Exact architectural decision needed
Adopt a package-owned extension interface that makes generic engine oblivious to TPF symbols:
1. **Config extension boundary:**
   - Generic `Config` keeps only generic fields plus `package_options` map (string->string).
   - Package API provides parse/validate/default/serialize hooks for package-owned keys.
2. **Mode ownership boundary:**
   - Generic `SimulationMode` reduced to generic modes only (or replaced by generic string token).
   - Package claims package modes through `owns_mode(mode_token)` and package-run entry points.
3. **Runtime policy boundary:**
   - Package hook contracts for integrator policy, startup cooling/step policy, diagnostics lifecycle, and compare behavior.
4. **Metadata/output boundary:**
   - Generic output/render writes package-supplied K/V metadata blobs, with no hardcoded TPF keys.
5. **Test ownership boundary:**
   - TPF behavior/config/mode/output tests moved under `engine/physics/TPFCore/tests/**`.

## Safest implementation sequence
1. Introduce generic package extension interfaces (config hooks, mode ownership, runtime/output hooks) with no behavior changes.
2. Migrate `Config` TPF keys into TPFCore-owned typed config built from generic `package_options`, preserving legacy key acceptance through TPFCore parser hook.
3. Migrate TPF mode parsing/dispatch from generic enum switches to package-owned mode resolution.
4. Move TPF runtime policy branches out of `simulation.cpp` into TPFCore runtime hooks.
5. Move TPF metadata/render/audit/compare output responsibilities to package hooks and package-owned artifacts.
6. Relocate TPF tests from generic `engine/tests/**` to `engine/physics/TPFCore/tests/**`; keep generic tests physics-package-agnostic.
7. Run `make test_unit` + TPF smoke/integration checks, then enforce grep gate with zero non-TPFCore matches.

## Why proceeding now would risk behavior/output/config breakage
Without first landing the interface boundary decisions above, directly deleting/moving references to satisfy grep would force either:
- temporary stubs that silently drop legacy config keys,
- altered mode resolution behavior for CLI/config modes,
- changed run_info/render schema keys relied on tooling/scripts,
- partially migrated runtime policy with mixed ownership.

That would violate the stated constraints (no scaffold-only, preserve existing behavior and config compatibility) and risks non-obvious simulation/output regressions.

## Current boundary scan evidence
Command run:

```bash
rg -n "TPF|tpf_|TPFCore|Xi|xi_|Theta|theta_|VDSG|direct_tpf|legacy_readout|geodesic_correspondence|v11_weak_field_truncation|C_mu_nu|C_μν" \
  engine \
  --glob '!engine/physics/TPFCore/**' \
  --glob '!build/**' \
  --glob '!outputs/**'
```

Result: extensive matches remain across generic engine sources and tests (1,000+ lines), confirming full isolation is not yet safely achievable as a small in-place edit.
