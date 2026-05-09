# TPF Package Boundary Audit (tpf-xi-theta-v1)

Date: 2026-05-08  
Branch: `tpf-xi-theta-v1`

## Scope and architecture rule

This audit checks for TPF-specific references **outside** `engine/physics/TPFCore/`.

Rule enforced:
- App/engine layers should know only generic package APIs and generic package selection.
- App/engine should not encode TPF theory objects (Xi/Theta), TPF route names, mode names, tensor names, or TPF-specific diagnostics pipelines.

## Inputs read before audit

- `engine/physics/TPFCore/TPF_FOUNDATIONS.md`
- `engine/physics/TPFCore/ARCHITECTURE_REFACTOR_PLAN.md`
- Package-loading/readme context:
  - `README.md`
  - `engine/README.md`
  - `engine/physics/Template/README.md`
  - `engine/physics/TPFCore/README.md`

## Search method (baseline, reproducible)

Searched outside the package boundary for:

- `TPF`, `tpf_`, `Xi`, `Ξ`, `Theta`, `Θ`, `VDSG`, `direct_tpf`, `xi_kernel_deformed`, `geodesic_correspondence`, `v11_weak_field_truncation`, `legacy_readout`, `C_mu_nu`, `C_μν`

Command used:

```bash
rg -n "TPF|tpf_|Xi|Ξ|Theta|Θ|VDSG|direct_tpf|xi_kernel_deformed|geodesic_correspondence|v11_weak_field_truncation|legacy_readout|C_mu_nu|C_μν" \
  --glob '!engine/physics/TPFCore/**' \
  --glob '!dev/**' \
  --glob '!outputs/**' \
  --glob '!.git/**' \
  --glob '!docs/audits/tpf_package_boundary_audit.md'
```

Baseline note: keep `dev/**` excluded for comparable leak tracking over time. Python/dev tooling is tracked separately below.

## Inventory summary

- **Raw matches outside boundary:** 2,342
- **Files with matches:** 73
- **Highest-risk files:** `engine/config.cpp` (247), `engine/output.cpp` (211), `engine/config.hpp` (173), `engine/render_audit.cpp` (118), `engine/main.cpp` (104), `engine/simulation.cpp` (9).

## Classification key

- **A** must move into TPFCore
- **B** must become generic package metadata/hook
- **C** acceptable package-owned reference
- **D** acceptable historical/audit/doc reference
- **E** stale or misleading reference to delete/update
- **F** test should move under TPFCore package tests

## Per-file inventory (boundary scan)

| path | raw match count | main matched terms | class | proposed action | risk | phase target |
|---|---:|---|:--:|---|---|---|
| `engine/config.cpp` | 247 | tpf_ | A | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | high | Phase 1-4 |
| `engine/tests/unit/test_config.cpp` | 245 | tpf_, geodesic_correspondence, TPF, direct_tpf | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `THEORY_IMPLEMENTATION_AUDIT.md` | 225 | Θ, direct_tpf, tpf_, TPF | D | Keep as historical/audit reference; optionally trim stale naming later. | low | Phase 6 |
| `engine/output.cpp` | 211 | tpf_, TPF, Xi, direct_tpf | A | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | high | Phase 1-4 |
| `engine/config.hpp` | 173 | tpf_, TPF, Xi, VDSG | A | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | high | Phase 1-4 |
| `PLAN_GEODESIC_CORRESPONDENCE.md` | 169 | tpf_, geodesic_correspondence, Ξ, TPF | D | Keep as historical/audit reference; optionally trim stale naming later. | low | Phase 6 |
| `engine/render_audit.cpp` | 118 | tpf_, TPF, Xi, legacy_readout | A | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | high | Phase 1-4 |
| `engine/main.cpp` | 104 | tpf_, TPF, direct_tpf, xi_kernel_deformed | A | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | high | Phase 1-4 |
| `engine/tests/unit/test_render_audit.cpp` | 100 | tpf_, TPF, Xi, direct_tpf | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `python_tests/unit/test_plot_cpp_compare.py` | 75 | tpf_, TPF, xi_kernel_deformed, direct_tpf | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `render_overlay.py` | 54 | tpf_, TPF, Xi, xi_kernel_deformed | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `docs/audits/tpfcore_legacy_runtime_cleanup_plan_2026-05-01.md` | 46 | tpf_, TPF, legacy_readout, xi_kernel_deformed | D | Keep as historical/audit reference; optionally trim stale naming later. | low | Phase 6 |
| `engine/tests/unit/test_resolved_scenario.cpp` | 40 | tpf_, TPF, direct_tpf, Theta | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `python_tests/unit/test_render_overlay.py` | 34 | tpf_, TPF, xi_kernel_deformed, legacy_readout | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/tests/integration/compare_declared_workflow_xi_kernel_deformed_smoke.sh` | 28 | xi_kernel_deformed, TPF, tpf_, direct_tpf | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `engine/tests/unit/test_tpf_xi_theta_v1.cpp` | 26 | tpf_, TPF, Xi, Theta | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `plot_cpp_compare.py` | 25 | tpf_, TPF, xi_kernel_deformed, direct_tpf | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `plot_cpp_run.py` | 21 | tpf_, legacy_readout, xi_kernel_deformed, direct_tpf | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `plot_tpf_4d_xi_motion_probe.py` | 20 | tpf_, Xi | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `docs/audits/tpfcore_calculation_audit_2026-04-29.md` | 19 | tpf_, Xi, TPF, xi_kernel_deformed | D | Keep as historical/audit reference; optionally trim stale naming later. | low | Phase 6 |
| `run_pure_vdsg_tests.sh` | 18 | TPF, tpf_, VDSG | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `plot_tpf_4d_static_residual.py` | 18 | tpf_, Xi | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/tests/integration/tpf_4d_static_residual_plot_harness_smoke.sh` | 18 | tpf_, TPF | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `plot_v11_earth_moon_line_benchmark.py` | 17 | tpf_, TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/force_compare.cpp` | 16 | TPF, tpf_ | B | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | medium | Phase 1-4 |
| `engine/tests/unit/test_simulation_integrator_gate.cpp` | 15 | tpf_, TPF, Xi | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `python_tests/unit/test_plot_cpp_run_sim_mode.py` | 14 | tpf_, legacy_readout, TPF, xi_kernel_deformed | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/galaxy_init.cpp` | 13 | tpf_, TPF | B | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | medium | Phase 1-4 |
| `validate/geodesic_correspondence_newtonian/compute_metrics.py` | 12 | geodesic_correspondence, tpf_, TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/tests/integration/compare_newtonian_equivalence_audit.sh` | 12 | tpf_, TPF, VDSG | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `diagnostics.py` | 11 | tpf_, TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `docs/audits/softening_audit_2026-04-29.md` | 11 | Xi, TPF, tpf_, xi_kernel_deformed | D | Keep as historical/audit reference; optionally trim stale naming later. | low | Phase 6 |
| `configs/example.cfg` | 10 | tpf_, TPF | C | Keep as package-scoped example config; align wording after hooks land. | low | Phase 6 |
| `tools_calibrate_tpfcore_readout_scale.py` | 9 | TPF, tpf_, legacy_readout | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/simulation.cpp` | 9 | tpf_, TPF | A | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | high | Phase 1-4 |
| `python_tests/unit/test_plot_tpf_4d_xi_motion_probe.py` | 8 | tpf_, TPF, Xi | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/resolved_scenario.cpp` | 8 | tpf_ | A | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | high | Phase 1-4 |
| `engine/Makefile` | 8 | TPF, tpf_ | B | Adjust build/package discovery toward package-local registration units without generic file edits. | medium | Phase 1 |
| `engine/tests/integration/tpf_4d_static_residual_plot_script_smoke.sh` | 8 | tpf_ | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `engine/force_compare.hpp` | 7 | tpf_, TPF, Theta, Xi | B | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | medium | Phase 1-4 |
| `engine/tests/integration/compare_declared_workflow_smoke.sh` | 7 | TPF, tpf_ | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `engine/tests/integration/tpf_4d_xi_motion_probe_plot_script_smoke.sh` | 7 | tpf_ | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `configs/tpf_xi_theta_v1_example.cfg` | 6 | tpf_, TPF, Theta, Xi | C | Keep as package-scoped example config; align wording after hooks land. | low | Phase 6 |
| `configs/tpf_legacy_readout_galaxy_example.cfg` | 6 | tpf_, legacy_readout, TPF, Theta | E | Update/archive stale preset naming and deprecated route semantics. | medium | Phase 6 |
| `configs/bh_orbit_validation_tpfcore_legacy_readout_vdsg_off.cfg` | 6 | legacy_readout, tpf_, TPF | E | Update/archive stale preset naming and deprecated route semantics. | medium | Phase 6 |
| `validate/geodesic_correspondence_newtonian/artifacts/metrics.json` | 6 | geodesic_correspondence, tpf_ | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/README.md` | 6 | TPF, Theta, VDSG, Xi | D | Keep as historical/audit reference; optionally trim stale naming later. | low | Phase 6 |
| `engine/tests/unit/test_tpf_core_xi_kernel_deformed.cpp` | 6 | tpf_, xi_kernel_deformed, TPF | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `engine/tests/unit/test_tpf_core_direct_tpf.cpp` | 6 | direct_tpf, tpf_, TPF | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `plot_rotation_curve.py` | 5 | TPF, tpf_ | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `validate/geodesic_correspondence_newtonian/configs/case3_smalln_starstar_tpf_gc.cfg` | 5 | geodesic_correspondence, tpf_, TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `validate/geodesic_correspondence_newtonian/configs/case1_bh_orbit_tpf_gc.cfg` | 5 | geodesic_correspondence, tpf_, TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `validate/geodesic_correspondence_newtonian/configs/case2_galaxy_nostarstar_tpf_gc.cfg` | 5 | geodesic_correspondence, tpf_, TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `hunt_coupling.py` | 4 | tpf_, TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `configs/tpf_direct_tpf_placeholder_example.cfg` | 4 | direct_tpf, tpf_, Theta, Xi | E | Update/archive stale preset naming and deprecated route semantics. | medium | Phase 6 |
| `configs/geodesic_correspondence_baseline.cfg` | 4 | geodesic_correspondence, tpf_, TPF | E | Update/archive stale preset naming and deprecated route semantics. | medium | Phase 6 |
| `docs/TESTING.md` | 3 | TPF | D | Keep as historical/audit reference; optionally trim stale naming later. | low | Phase 6 |
| `analyze_coherence.py` | 3 | TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `python_tests/integration/test_postprocess_truthfulness.py` | 3 | tpf_, TPF, xi_kernel_deformed | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/accel_pipeline_stats.hpp` | 3 | VDSG, TPF | B | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | medium | Phase 1-4 |
| `engine/tests/integration/preflight_compare_yes_bypass_smoke.sh` | 3 | TPF, tpf_ | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `python_tests/unit/test_plot_rotation_curve_run_info_precedence.py` | 2 | TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `engine/render_audit.hpp` | 2 | VDSG, TPF, legacy_readout, Theta | B | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | medium | Phase 1-4 |
| `engine/tests/run_integration.sh` | 2 | xi_kernel_deformed, tpf_ | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `engine/tests/integration/utility_run_info_provenance_and_particles_smoke.sh` | 2 | tpf_, TPF | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `engine/output.hpp` | 2 | tpf_, VDSG, TPF | B | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | medium | Phase 1-4 |
| `README.md` | 1 | TPF | D | Keep as historical/audit reference; optionally trim stale naming later. | low | Phase 6 |
| `python_tests/integration/test_compare_render_success_reporting.py` | 1 | TPF | B | Defer to tooling audit; treat as consumer of package metadata hooks and remove route hard-codes later. | medium | Phase 6 |
| `configs/smoke_test.cfg` | 1 | tpf_ | C | Keep as package-scoped example config; align wording after hooks land. | low | Phase 6 |
| `configs/bh_orbit_validation_newtonian_baseline.cfg` | 1 | direct_tpf, VDSG | C | Keep as package-scoped example config; align wording after hooks land. | low | Phase 6 |
| `engine/init_conditions.hpp` | 1 | TPF | B | Introduce generic package hook/metadata boundary; remove TPF-specific branching from engine layer incrementally. | medium | Phase 1-4 |
| `engine/tests/integration/newtonian_stability_matrix_audit.sh` | 1 | tpf_, TPF | F | Move/retarget to package test targets where TPF-specific; keep only generic contract tests in engine suite. | medium | Phase 5 |
| `engine/physics/Template/README.md` | 1 | TPF | B/E | Update template docs in Phase 1 to describe the current static self-registration model and the new generic package hooks; remove stale instructions about editing registry.cpp / s_packages. | medium | Phase 1 |

## Package registration / autoload audit (audit-only)

Files inspected:
- `engine/physics/physics_package.hpp`
- `engine/physics/registry.cpp`
- `engine/physics/Template/README.md`
- `engine/Makefile`
- package registration units under `engine/physics/*` (notably `Newtonian/newtonian.cpp`, `TPFCore/tpf_core_package.cpp`)

Findings:

1. **Does adding a package currently require editing generic engine files?**  
   **Yes.** At minimum, `engine/Makefile` must list package `.cpp` files in `LIB_SRCS`; otherwise they are not compiled/linked into the binary.

2. **Does adding a package require adding headers outside the package?**  
   **Not strictly in the current registrar pattern, but template guidance is stale.** Current runtime registry is factory-based via `register_physics_package_factory(...)` and package-local static registrars; however `engine/physics/Template/README.md` still instructs edits to `physics/registry.cpp` and references a non-current `s_packages` array model.

3. **Is registration manual, static self-registration, build-discovered, or dynamic?**  
   **Hybrid/manual today:**
   - Runtime registration is **static self-registration** (package translation unit calls registry function at static init).
   - Build inclusion is **manual** (explicit source list in Makefile).
   - Test discovery is partially build-discovered for package tests (`find physics ... tests/unit|tests/regression`).
   - No runtime plugin/dynamic loading currently exists.

4. **Smallest safe Phase 1 step toward “drop package folder in `engine/physics` and autoload it”?**  
   - Keep runtime behavior unchanged.
   - Add **generic package metadata/capability hook interface** in `PhysicsPackage` (no-op defaults).
   - Add a **generic package utility dispatch hook** and run-info/render metadata hook surfaces (no-op defaults).
   - Update template docs to match registrar reality.
   - Make one minimal build-system improvement for package onboarding (e.g., a package-local source include mechanism) without moving TPF behavior yet.

5. **What should be deferred to later?**  
   - Full build autoload/discovery for arbitrary package folders.
   - Runtime dynamic plugin loading.
   - Moving existing TPF behavior and labels out of engine runtime/output/render paths (Phase 2+).
   - Config schema isolation and migration logic (Phase 4).

## Separate Python/tooling note (deferred)

Baseline leak counting intentionally excludes `dev/**` for comparability.

However, Python/tooling appears to carry package architecture debt and should be audited in a dedicated pass (Phase 6), including likely candidates already surfaced by this scan:
- `render_overlay.py`, `plot_cpp_compare.py`, `plot_cpp_run.py`, `plot_rotation_curve.py`
- `diagnostics.py`, `hunt_coupling.py`, `tools_calibrate_tpfcore_readout_scale.py`
- `python_tests/unit/test_plot_cpp_compare.py`, `python_tests/unit/test_render_overlay.py`
- `validate/geodesic_correspondence_newtonian/*`

No Python/tooling behavior changes are made in this audit PR.


## Phase 1 status (scaffold added)

Implemented in Phase 1:
- Added generic `PhysicsPackage` hooks for display/runtime metadata, utility-mode capability/dispatch, run-info metadata, render metadata, and package config metadata (all default no-op/empty).
- Added minimal engine-level hook plumbing/tests to validate interface availability with Newtonian and TPFCore while preserving behavior.
- Updated the package template README to reflect static self-registration and current Makefile source-inclusion requirements.

Intentionally not moved in Phase 1:
- TPF utility dispatch remains in `engine/main.cpp` (no behavior migration yet).
- TPF-specific run-info/render field generation remains in existing engine paths.
- No config schema migration or TPF config field relocation.

Next phase target:
- Move TPF utility dispatch out of `main.cpp` into package hook dispatch.

## Actionable phased plan

### Phase 1 (next implementation PR only)

**Goal:** add generic package capability/metadata hooks with no behavior movement.

Scope (exact recommendation):
- Add package display/runtime metadata hook(s).
- Add package utility-mode capability query hook.
- Add package utility-mode dispatch hook.
- Add package output/run-info metadata hook.
- Add package render/diagnostic metadata hook.
- Add package config schema/defaults hook or package config adapter hook.
- If needed, add only trivial no-op adapters in TPFCore to satisfy interface compile/link requirements.

Non-goals in Phase 1:
- No TPF logic migration from engine runtime flow (except trivial no-op adapter plumbing).
- No output/render behavior rewiring.
- No config migration.
- No test relocation beyond minimal compile wiring if strictly required.

### Later phases

- **Phase 2:** move TPF utility dispatch out of `main.cpp`.
- **Phase 3:** move output/render metadata into package hooks.
- **Phase 4:** isolate package config.
- **Phase 5:** move TPF-specific tests into package tests.
- **Phase 6:** docs/config/Python cleanup.

## Acceptance criteria checklist

- [x] Per-file inventory table added.
- [x] Baseline `rg` command retained and reproducible.
- [x] Package-registration/autoload blockers identified explicitly.
- [x] Runtime code unchanged (audit-only).
- [x] No behavior changes.
