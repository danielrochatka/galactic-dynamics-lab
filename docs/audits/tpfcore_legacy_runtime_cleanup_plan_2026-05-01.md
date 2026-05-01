# TPFCore legacy/runtime cleanup audit and staged removal plan (2026-05-01)

## Scope and active paths
Active/supported paths for current work are treated as:
1. Newtonian baseline.
2. TPFCore `xi_kernel_deformed` runtime path.
3. TPFCore `direct_tpf` paper-facing path (kept unless paper scope drops it).
4. Compare rendering/diagnostics needed to inspect those paths.

This audit classifies legacy/runtime-adjacent items as:
- **A Keep**: required for active workflow now.
- **B Quarantine**: still useful, but should not affect normal runs by default.
- **C Delete**: dead/legacy and not required by current active scope.

---

## Inventory and classification

### 1) `legacy_readout` runtime path and provisional gate
- **Files / symbols**
  - `engine/physics/TPFCore/tpf_core_package.cpp`
    - `TPFCorePackage::compute_accelerations` legacy branch behavior (legacy readout + VDSG + optional shunt), constructor default `tpf_dynamics_mode_="legacy_readout"`.
  - `engine/config.hpp`, `engine/config.cpp`
    - `Config::tpf_dynamics_mode` default is `legacy_readout`.
    - `tpfcore_enable_provisional_readout`, `tpfcore_readout_mode`, `tpfcore_readout_scale`, `tpfcore_theta_tt_scale`, `tpfcore_theta_tr_scale`.
  - `engine/physics/TPFCore/provisional_readout.cpp/.hpp`
    - `compute_provisional_readout_acceleration`, `compute_provisional_readout_with_diagnostics`, `write_readout_debug_csv`.
  - Tests heavily centered on legacy routing:
    - `engine/physics/TPFCore/tests/unit/test_routing_semantics.cpp`
    - `python_tests/unit/test_render_overlay.py`
    - `python_tests/unit/test_plot_cpp_run_sim_mode.py`
- **Classification**: **B Quarantine**.
- **Rationale**: still exercised and documented, but interferes by being default and by broad diagnostics/logging assumptions.
- **Removal plan**
  1. **Phase 1 (safe):** flip default runtime to `xi_kernel_deformed` in `Config` and `TPFCorePackage` constructor; require explicit opt-in for `legacy_readout`.
  2. Move legacy/provisional code into `engine/physics/TPFCore/legacy/` namespace segment and gate with `tpf_enable_legacy_runtime=false` default.
  3. Convert legacy-only tests into quarantined test target (e.g. `make test-legacy`) and remove from default CI.
  4. **Phase 2:** delete provisional readout runtime path once no active scenario references it.

### 2) `provisional_readout` closures and debug exports
- **Files / symbols**
  - `engine/physics/TPFCore/provisional_readout.cpp` modes:
    - `tensor_radial_projection`, `tensor_radial_projection_negated`, `experimental_radial_r_scaling`, derived radial mode dispatch.
  - `engine/config.hpp/.cpp`: `tpfcore_dump_readout_debug`, readout mode/scale keys.
  - Output hooks in `engine/output.cpp` and `engine/render_audit.cpp` emit readout-centric notes.
- **Classification**: **B Quarantine**.
- **Rationale**: useful for forensic benchmarking, but not part of preferred forward runtime.
- **Removal plan**
  - Introduce `tpf_readout_legacy_diagnostics_enabled=false` default to suppress CSV/debug emissions during normal runs.
  - Fence `write_readout_debug_csv` behind both legacy runtime selection and explicit flag.
  - Split readout diagnostics text into `legacy_runtime_notes` section so standard run_info stays xi/direct-focused.

### 3) Old VDSG additive branches
- **Files / symbols**
  - `engine/physics/TPFCore/extensions_vdsg.cpp/.hpp`: `accumulate_vdsg_velocity_modifier`, `vdsg_effective_coupling`.
  - Invocation points in legacy and direct flows (`tpf_core_package.cpp`) and extensive render labels (`render_audit.cpp`).
  - Tests:
    - `engine/physics/TPFCore/tests/integration/tpfcore_vdsg_additive_branch_regression.sh`
    - Many assertions in `test_routing_semantics.cpp`, `test_render_audit.cpp`, `test_render_overlay.py`.
- **Classification**: **A Keep** for `direct_tpf`; **B Quarantine** for `legacy_readout` coupling behavior.
- **Rationale**: user requested keeping direct path if paper-facing; direct path currently includes optional additive VDSG terms.
- **Removal plan**
  - Keep VDSG implementation but stop referencing it in xi path metadata beyond "not used".
  - Move legacy-specific VDSG messaging/tests to quarantine test suite.
  - If paper scope confirms direct path no longer needs additive VDSG, delete in a dedicated PR (code + tests together).

### 4) Stale inspection/debug routes
- **Files / symbols**
  - Simulation modes in `engine/config.hpp/.cpp`:
    - `tpf_single_source_inspect`, `tpf_symmetric_pair_inspect`, `tpf_source_field_benchmark`, `tpf_weak_field_calibration`, `tpf_newtonian_force_compare`, `tpf_diagnostic_consistency_audit`.
  - Diagnostic runner and outputs in `engine/force_compare.cpp` and `engine/output.cpp`.
  - Xi exterior inspection knobs used by optional experimental solver.
- **Classification**: **B Quarantine** (most); **A Keep** only where directly used by xi/direct/newtonian compare.
- **Rationale**: valuable dev tools, but should not appear as first-class runtime behavior for production progress.
- **Removal plan**
  - Tag these modes as `experimental/inspection` in parser help and render metadata.
  - Move related scripts/tests to non-default integration lane.
  - Delete `tpf_diagnostic_consistency_audit` mode if no owner/consumer in next milestone.

### 5) Unused/deprecated config keys and aliases
- **Files / symbols**
  - `engine/config.cpp`: accepts alias `tpf_gdd_coupling -> tpf_vdsg_coupling` and deprecated `weak_field_correspondence` mode alias.
  - `engine/tests/unit/test_config.cpp`: explicit tests for these aliases.
  - `engine/render_audit.cpp`: emits alias notes in JSON/TXT.
- **Classification**: **C Delete** (after one deprecation cycle) for `tpf_gdd_coupling` and weak-field alias acceptance.
- **Rationale**: pure compatibility baggage; keeps legacy naming alive in current run metadata.
- **Removal plan**
  1. Phase 1: warn-on-use, add strict mode failing on alias keys.
  2. Phase 2: remove parser acceptance blocks and alias text emission.
  3. Delete/update affected tests:
     - `engine/tests/unit/test_config.cpp` alias cases.
     - `engine/tests/unit/test_render_audit.cpp` alias-note expectations.

### 6) Temporary startup diagnostics/noise in run_info/render metadata
- **Files / symbols**
  - `engine/output.cpp`: large TPF branch/diagnostics supplements blocks, including legacy mode notes and repeated key dumps.
  - `engine/render_audit.cpp`: verbose legacy alias/dynamics annotations.
  - `engine/main.cpp`: reference decomposition columns and notes for direct/VDSG decomposition in snapshots.
- **Classification**: **B Quarantine** (verbosity); selective **C Delete** for clearly duplicated fields.
- **Rationale**: useful during transition, but currently noisy and broad enough to blur active-path intent.
- **Removal plan**
  - Introduce concise default metadata profile (active-path only), plus optional `tpf_verbose_audit_metadata=true`.
  - Remove duplicated/legacy-only lines from default `run_info.txt`.
  - Keep a machine-readable detailed artifact only when verbose flag is enabled.

### 7) Legacy-oriented configs/examples
- **Files**
  - `configs/tpf_legacy_readout_galaxy_example.cfg`
  - `configs/bh_orbit_validation_tpfcore_legacy_readout_vdsg_off.cfg`
  - `configs/tpf_direct_tpf_placeholder_example.cfg` (keep but rename once stable)
- **Classification**:
  - legacy examples: **B Quarantine** now, **C Delete** later.
  - direct placeholder: **A Keep** if paper-facing.
- **Removal plan**
  - Move legacy examples under `configs/legacy/` and prepend deprecation headers.
  - Add canonical `xi_kernel_deformed` + `direct_tpf` minimal configs as top-level examples.

---

## Concrete staged execution plan (no deletion in this pass)

### Stage 0: safety rails (next PR)
1. Change defaults to active-path-first:
   - `Config::tpf_dynamics_mode` default -> `xi_kernel_deformed`.
   - `TPFCorePackage` constructor default `tpf_dynamics_mode_` -> `xi_kernel_deformed`.
2. Add explicit opt-in key `tpf_enable_legacy_runtime=false` for `legacy_readout` branch.
3. Add startup warning whenever legacy runtime/alias keys are used.

### Stage 1: quarantine mechanics
1. Move legacy readout and its diagnostics/test fixtures into isolated tree (`TPFCore/legacy/`).
2. Split CI:
   - default: Newtonian + xi_kernel_deformed + direct_tpf + required compare diagnostics.
   - optional: legacy/quarantine suite.
3. Move legacy configs into `configs/legacy/`.

### Stage 2: alias and stale-mode removals
1. Remove parser support for:
   - `tpf_gdd_coupling`
   - `tpf_dynamics_mode=weak_field_correspondence` alias.
2. Remove stale inspection modes lacking active owners/consumers (`tpf_diagnostic_consistency_audit` first candidate).
3. Update `render_audit` and `output` metadata to drop legacy compatibility language.

### Stage 3: hard delete legacy runtime
1. Delete `provisional_readout.cpp/.hpp` and legacy branch wiring from `tpf_core_package.cpp`.
2. Delete now-obsolete config keys:
   - `tpfcore_enable_provisional_readout`
   - `tpfcore_readout_mode`
   - `tpfcore_readout_scale`
   - `tpfcore_theta_tt_scale`
   - `tpfcore_theta_tr_scale`
   - `tpfcore_dump_readout_debug`
3. Remove/update affected tests:
   - `engine/physics/TPFCore/tests/unit/test_routing_semantics.cpp` (legacy cases)
   - `python_tests/unit/test_render_overlay.py` legacy labels
   - `python_tests/unit/test_plot_cpp_run_sim_mode.py` legacy mode label checks
   - `engine/tests/unit/test_render_audit.cpp` legacy/alias notes
   - `engine/tests/unit/test_config.cpp` legacy alias + deprecated mode cases

---

## Risk controls before any delete PR
- Require one green run for:
  - Newtonian baseline unit/integration tests.
  - `engine/tests/unit/test_tpf_core_xi_kernel_deformed.cpp`.
  - `engine/tests/unit/test_tpf_core_direct_tpf.cpp` and `engine/physics/TPFCore/tests/unit/test_direct_tpf_baseline_structure.cpp`.
  - compare workflow smokes (`engine/tests/integration/compare_declared_workflow*.sh`).
- Freeze interface contract for `run_info.txt` fields consumed by current Python diagnostics prior to removing legacy entries.
