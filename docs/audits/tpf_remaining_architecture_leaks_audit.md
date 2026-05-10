# TPF Remaining Architecture Leaks Audit (Phase 5D-A)

Date: 2026-05-10  
Branch: `tpf-xi-theta-v1`  
Scope: audit-only (no runtime/output/config/test behavior changes)

## 1) Executive summary

Current status after Phases 5B/5C:
- Package hooks are in place and are already used for utility dispatch, post-run diagnostics, and part of run_info/render metadata.
- Remaining leaks are concentrated in five engine files (`main.cpp`, `simulation.cpp`, `output.cpp`, `render_audit.cpp`, `config.*`) plus hook-shape pressure in `physics_package.hpp`.
- The largest near-term risk is **mixed ownership blocks** where package metadata, schema compatibility, and global config provenance are interleaved.

Biggest remaining leaks:
1. `engine/output.cpp` still contains selected-package TPF run_info plus pipeline stats and compatibility/schema blocks mixed together.
2. `engine/simulation.cpp` still contains package-specific integrator policy gate (`xi_kernel_deformation_active`) in generic simulation flow.
3. `engine/main.cpp` still contains TPF-specific step-0 acceleration audit/decomposition and benchmark plot coupling.
4. `engine/render_audit.cpp` still contains TPF semantic interpretation and config-coupled fallback logic in generic manifest rendering.
5. `engine/config.hpp/.cpp` remains globally TPF-coupled (modes, fields, parse/serialize keys), which should not be altered in Phase 5D-B.

Recommended next implementation slice:
- **Phase 5D-B should only move A-class safe selected-package TPFCore run_info metadata in `engine/output.cpp` into package ownership** while preserving compatibility fallback behavior and existing output schema keys.

## 2) Per-file/function table

| file | function/area | current package-specific dependency | current owner | desired owner | classification | risk | recommended phase | notes |
|---|---|---|---|---|---|---|---|---|
| `engine/main.cpp` | `write_galaxy_step0_accel_audit(...)` | Direct `TPFCore` checks; `tpf_dynamics_mode`; Xi/Theta labeling; direct package fetch/init and decomposition recipe | generic engine | TPFCore diagnostics hook | B | High | 5F | Engine function embeds TPF route semantics and historical direct/vdsg decomposition text. |
| `engine/main.cpp` | Utility flow after package hook dispatch | Generic utility dispatch, but optional benchmark plotting hard-codes `tpf_4d_static_residual_benchmark` script/path/PNG names | mixed | generic + package-provided plotting metadata later | B | Medium | 5F (or 5E/6 tooling) | Dispatch is good; post-dispatch plotting path is still TPF-named compatibility behavior. |
| `engine/main.cpp` | direct package-name conditionals | `if (physics_package == "TPFCore")` remains in diagnostics path | generic engine | package hook or compatibility wrapper | B | Medium | 5F | Keep untouched in 5D-B to avoid cross-file migration risk. |
| `engine/simulation.cpp` | `xi_kernel_deformation_active` gate | TPF-specific config checks inside generic integrator loop; selects Euler vs Verlet | generic simulation | package runtime/integrator policy hook | B | High | 5G | This is package-specific runtime policy leak in generic code. Also appears currently unreachable for accepted TPFCore configs because non-off xi-kernel/nonzero coupling is rejected at package init. |
| `engine/simulation.cpp` | cooling usage (`cooling_active`, `apply_cooling_step`, `suppress_snapshot_for_cooling`) | Hook-based and generic | generic simulation | keep generic | D | Low | leave | Current pattern is aligned with package-oriented architecture. |
| `engine/output.cpp` | remaining selected-package TPFCore run_info block (`if physics_package == TPFCore`) | TPF dynamics keys, law/readout labels, xi-kernel fields, cooling and tpfcore notes in generic writer | generic output | TPFCore `run_info_metadata` (safe subset first) | A/B/C mixed | High | 5D-B (A only) | Split required by key group; avoid moving config/provenance coupled lines now. |
| `engine/output.cpp` | benchmark compatibility fallback (non-TPF package + TPF benchmark modes) | Schema/backward-compat run_info behavior in generic code | generic output | generic compatibility layer | D | Medium | leave | This fallback is deliberate compatibility debt from 5C; keep in generic layer for now. |
| `engine/output.cpp` | `tpf_last_*` pipeline stats (`AccelPipelineStats`) | Needs runtime step stats context not exposed via metadata hook | generic output + runtime struct | future package runtime stats hook/context | E | Medium | 5D-C | Requires hook/context shape change before move. |
| `engine/output.cpp` | config/provenance fields and schema compatibility keys | Interleaved with package labels and legacy aliases | generic output | mostly generic until config isolation | C/D | Medium | post-6 | Not safe for 5D-B; ties into global config model and compatibility commitments. |
| `engine/output.cpp` | false-guarded v11 blocks | stale/dead emission blocks | generic output | cleanup later | F | Low | later explicit cleanup | Keep untouched until explicit cleanup decision. |
| `engine/render_audit.cpp` | `compute_active_dynamics_branch`, `compute_active_metrics_branch`, `compute_acceleration_code_path` | TPF semantic label interpreter logic in generic manifest writer | generic render audit | likely package semantic metadata provider + generic formatter | B | High | 5E | Helper strings are semantic ownership leak, but migration needs careful manifest parity plan. |
| `engine/render_audit.cpp` | cooling active calc + alias/effective-mode helpers | TPF config-coupled logic (`legacy_readout`, geodesic/v11 alias behavior) | generic render audit | defer until config isolation | C | Medium | 6+ | Coupled to global config and historical mode alias logic. |
| `engine/render_audit.cpp` | non-TPF legacy-readout schema fallback | Explicit compatibility fallback branch | generic render audit | generic compatibility layer | D | Medium | leave | Keep until schema compatibility policy changes explicitly. |
| `engine/render_audit.cpp` | direct include `physics/TPFCore/derived_tpf_radial.hpp` | TPF type coupling in generic layer | generic render audit | package hook payload types / string metadata only | B | Medium | 5E | Indicates manifest semantic computation is not fully package-isolated. |
| `engine/render_audit.cpp` | false-guarded v11 blocks | stale/dead | generic | cleanup later | F | Low | later | Explicitly defer cleanup. |
| `engine/config.hpp` | `SimulationMode` includes many `tpf_*` utility/benchmark modes | global enum owns package-specific runtime modes | global config | config isolation layer + package mode descriptors later | C | High | 6 | Out of scope for 5D-B by design. |
| `engine/config.hpp` | many `Config` fields are TPFCore-only | global config schema tightly coupled to TPF | global config | package config metadata/isolation | C | High | 6 | Must wait for config isolation plan. |
| `engine/config.cpp` | `parse_mode` / `is_tpf_utility_mode` / `simulation_mode_requires_output_dir` | global parser/runtime routing includes TPF modes | global config | deferred config/package boundary redesign | C | High | 6 | Not safe for 5D-B. |
| `engine/config.cpp` | `apply_config_kv` TPF key parsing + branch guards | huge TPF key handling in generic parser | global config | config isolation + package parser metadata | C | High | 6 | Migration here would broaden scope and risk parity regressions. |
| `engine/config.cpp` | `serialize_config_kv` TPF key emission | global serialization schema includes TPF keys | global config | config isolation | C | High | 6 | Must remain stable through 5D sequence. |
| `engine/config.cpp` | `find_package_defaults_path` assumptions | package defaults path resolution remains generic string/path convention | generic config | can stay generic short-term | D | Low | leave | Later revisit for autoload/discovery phase. |
| `engine/physics/physics_package.hpp` | hook surface breadth (`run_info_metadata`, `render_metadata`, cooling, diagnostics, utility dispatch) | growing interface surface; metadata/runtime policy concerns spread across many methods | package interface | possible grouped provider/policy structs (future) | B | Medium | 5D-C/5E design note | Do not refactor now; note pressure for future struct grouping. |

## 3) Classification legend

- **A** = safe selected-package ownership move candidate  
- **B** = needs generic hook/context design before move  
- **C** = config-isolation-coupled; defer  
- **D** = generic or schema/provenance compatibility; leave for now  
- **E** = runtime/pipeline stats require extra context; defer  
- **F** = stale/dead/false-guarded; cleanup later only after explicit decision

## 4) Recommended next PR sequence

1. **Phase 5D-B**: move only A-class safe selected-package TPF dynamics run_info metadata in `engine/output.cpp` to `TPFCorePackage::run_info_metadata(...)`.
2. **Phase 5D-C**: design and add runtime/pipeline metadata context for `tpf_last_*` entries (without schema churn).
3. **Phase 5E**: render_audit semantic label ownership audit/design (`compute_active_*` ownership split and TPF semantic provider design).
4. **Phase 5F**: `main.cpp` diagnostic hook migration (`write_galaxy_step0_accel_audit` and residual optional plotting coupling review).
5. **Phase 5G**: simulation integrator policy hook/design (remove `xi_kernel_deformation_active` package logic from generic loop).
6. **Phase 6**: config isolation (mode/field parsing and serialization boundaries).
7. **Phase 7**: package autoload/discovery.
8. **Phase 8**: broader test relocation and Python/tooling cleanup.

## 5) Do-not-touch list for Phase 5D-B

Must not be touched in the next implementation PR:
- Config migration (`engine/config.hpp`, `engine/config.cpp` schema/parse/serialize changes)
- `simulation.cpp` integrator policy logic
- `render_audit.cpp` semantic labeling logic
- `main.cpp` diagnostics/plotting behavior
- package autoload/discovery work
- physics/math implementation behavior
- output schema key names/meanings/order changes

## 6) Search commands used

```bash
rg -n 'TPFCore|tpf_|Xi|xi_|Theta|theta_|VDSG|direct_tpf|legacy_readout|geodesic_correspondence|v11_weak_field_truncation|C_mu_nu|C_μν' engine/main.cpp engine/simulation.cpp engine/output.cpp engine/render_audit.cpp engine/config.hpp engine/config.cpp engine/physics/physics_package.hpp

rg -n 'physics_package ==|get_physics_package|run_info_metadata|render_metadata|write_post_run_diagnostics|cooling_active|apply_cooling_step|suppress_snapshot_for_cooling|semi_implicit_euler_step|velocity_verlet_step' engine/main.cpp engine/simulation.cpp engine/output.cpp engine/render_audit.cpp engine/physics/physics_package.hpp

rg -n 'tpf_last_|AccelPipelineStats|tpf_accel_pipeline|pipeline' engine/output.cpp engine/main.cpp engine/physics/TPFCore
```


## Phase 5D-B status (implemented)

Phase 5D-B completed for safe selected-package run_info supplement metadata ownership.

### Candidate inventory and disposition

| key/block | current condition | original output location/order | proposed owner | class | moved in 5D-B? |
|---|---|---|---|---|---|
| `tpf_dynamics_mode` | `physics_package==TPFCore` and non-benchmark branch | branch/physics supplement block, before route/readout status lines | `TPFCorePackage::run_info_supplement_metadata` | A | yes |
| `tpf_runtime_path_tier` | same, mode-conditional (`direct_tpf`/`legacy_readout`/`xi_kernel_deformed`) | immediately after `tpf_dynamics_mode` | package | A | yes |
| `tpf_runtime_deprecation_note` | `legacy_readout` | after `tpf_runtime_path_tier` in legacy path | package | A | yes |
| `tpf_core_law_mode` | `direct_tpf` or `xi_kernel_deformed` | route-status lines | package | A | yes |
| `tpf_core_dynamics_route` | `xi_kernel_deformed` | route-status lines | package | A | yes |
| `acceleration_formula`, `K_xi` | `xi_kernel_deformed` | route-status lines | package | A | yes |
| `tpf_4d_xi_motion_readout_scale` (supplement copy) | `xi_kernel_deformed` | xi route block | package | A | yes |
| `xi_kernel_mode`, `xi_kernel_label`, `xi_kernel_coupling`, `xi_kernel_factor_mode`, `xi_kernel_metric_min`, `xi_kernel_metric_max`, `xi_temporal_mode` | `xi_kernel_deformed` | xi route block | package | A | yes |
| `tpf_truncation_status`, `tpf_higher_order_status`, `tpf_extension_status`, `tpf_provisional_readout_status`, `tpf_readout_closure_knobs_status`, `tpf_stabilizer_status` | `direct_tpf` | direct route status lines | package | A | yes |
| `tpfcore_enable_provisional_readout` / `tpfcore_readout_mode` or `*_status` variants | legacy/non-legacy mode split | end of moved route/status subset | package | A | yes |
| `tpf_kappa`, `tpf_vdsg_*`, `tpf_poisson_*`, `tpf_cooling_fraction`, `tpf_global_accel_shunt_*`, `tpf_accel_pipeline_diagnostics_csv` | same enclosing TPF block | after moved subset | mixed config/pipeline/schema | C/E/D | deferred |
| `tpf_last_*` and `AccelPipelineStats` section | pipeline stats block | later dedicated section | runtime stats context | E | deferred |
| benchmark compatibility fallback | non-TPF + benchmark modes | earlier compatibility branch | generic compatibility | D | deferred |

## Phase 5D-C-A status (audit/design only)
- Phase 5D-C-A completed as audit/design only; no runtime or output behavior changed.
- Pipeline/runtime run_info metadata (`tpf_last_*` and related runtime counter context) was inventoried for future package ownership.
- Implementation is explicitly deferred to Phase 5D-C-B.
- See: `docs/audits/tpf_pipeline_run_info_metadata_audit.md`.
