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

## Search method

Searched outside the package boundary for:

- `TPF`, `tpf_`, `Xi`, `Ξ`, `Theta`, `Θ`, `VDSG`, `direct_tpf`, `xi_kernel_deformed`, `geodesic_correspondence`, `v11_weak_field_truncation`, `legacy_readout`, `C_mu_nu`, `C_μν`

Command used:

```bash
rg -n "TPF|tpf_|Xi|Ξ|Theta|Θ|VDSG|direct_tpf|xi_kernel_deformed|geodesic_correspondence|v11_weak_field_truncation|legacy_readout|C_mu_nu|C_μν" \
  --glob '!engine/physics/TPFCore/**' --glob '!dev/**' --glob '!outputs/**' --glob '!.git/**'
```

## Inventory summary

- **Raw matches outside boundary:** 2,342
- **Files with matches:** 66
- Largest concentrations:
  1. `engine/config.cpp` (247)
  2. `engine/tests/unit/test_config.cpp` (245)
  3. `engine/output.cpp` (211)
  4. `engine/config.hpp` (173)
  5. `engine/render_audit.cpp` (118)
  6. `engine/main.cpp` (104)

## Classification key

- **A** must move into TPFCore
- **B** must become generic package metadata/hook
- **C** acceptable package-owned reference
- **D** acceptable historical/audit/doc reference
- **E** stale or misleading reference to delete/update
- **F** test should move under TPFCore package tests

## Findings by boundary area

### 1) Engine runtime/orchestration boundary leaks (highest severity)

#### 1.1 `engine/main.cpp` — **A/B/E**

Observed:
- Direct include of TPF package type and downcast usage.
- Hard-coded TPF utility simulation modes and route-specific strings.
- TPF-specific plot script invocation and fixed TPF artifact filenames.
- Inline route logic/labels for Xi/Theta and VDSG/reference decomposition.

Classification:
- **A**: TPF-specific runtime orchestration logic should be package-owned.
- **B**: replace engine branching with generic package capability hooks (`run_package_utility_mode`, `postprocess_artifacts`, `diagnostic_manifest`).
- **E**: comments/labels referencing removed/legacy route semantics are misleading in app layer.

Risk if moved:
- Utility mode CLI behavior and artifact naming can break unless hook contract is stabilized first.

#### 1.2 `engine/simulation.cpp` — **A/B**

Observed:
- Engine-level cooling and kernel gates keyed on `physics_package == "TPFCore"` and `tpf_*` fields.

Classification:
- **A** for TPF policy logic in generic stepping loop.
- **B** for converting to generic optional package runtime policy hook (e.g., package-provided startup damping policy).

Risk if moved:
- Step count/snapshot skipping behavior may shift if ordering differs from current logic.

### 2) Config surface contamination in engine core

#### 2.1 `engine/config.hpp` + `engine/config.cpp` + `engine/resolved_scenario.cpp` — **A/B**

Observed:
- Extensive TPF-only keys and mode enums are defined in global app config model.
- Parsing/validation logic enforces TPF route names and aliases at engine-global level.
- Scenario resolution contains TPF-only defaults/guards.

Classification:
- **A**: TPF theory and route knobs should leave generic app config structures.
- **B**: keep only generic package names + opaque package config map (or package schema adapters).

Risk if moved:
- Existing config files and CLI overrides depend on current key names; migration layer needed.

### 3) Output/audit/render coupling to TPF internals

#### 3.1 `engine/output.cpp`, `engine/output.hpp`, `engine/render_audit.cpp`, `engine/render_audit.hpp`, `engine/accel_pipeline_stats.hpp` — **A/B/E**

Observed:
- Engine output manifests contain TPF route names and theory labels.
- Render overlays infer branches from TPF keys (`legacy_readout`, `direct_tpf`, `xi_kernel_deformed`, VDSG semantics).
- TPF diagnostic filenames and branch IDs are hard-coded outside package.

Classification:
- **A** for package-specific branch semantics in engine output path.
- **B** for generic metadata envelope where package contributes self-described labels and artifacts.
- **E** for stale legacy labels where branch contract already narrowed to `tpf_xi_theta_v1`.

Risk if moved:
- Downstream Python tooling that parses current run_info/render fields may break; requires compatibility fields during transition.

### 4) Non-runtime docs/configs/tooling

#### 4.1 Root/engine readmes and historical plans — **D**

Files like `README.md`, `engine/README.md`, `PLAN_GEODESIC_CORRESPONDENCE.md`, `THEORY_IMPLEMENTATION_AUDIT.md`, and prior audit docs are acceptable as historical/project documentation.

#### 4.2 `configs/*.cfg` TPF examples and deprecated route presets — **C/E**

- **C**: explicit package example configs are acceptable if clearly package-scoped.
- **E**: stale presets referencing disallowed/removed modes (`geodesic_correspondence`, legacy route naming) should be archived or updated to current branch truth.

#### 4.3 Python plotting/tool scripts at repo root — **B/C/E**

- Scripts that introspect TPF run metadata from generic outputs are partly acceptable (**C**), but many currently encode route internals and branch labels (**B/A tendency**).
- Some labels mention removed modes as active runtime options (**E**).

### 5) Test placement and scope

#### 5.1 Engine unit/integration tests with TPF internals — **F/B**

Observed candidates:
- `engine/tests/unit/test_tpf_core_direct_tpf.cpp`
- `engine/tests/unit/test_tpf_core_xi_kernel_deformed.cpp`
- `engine/tests/unit/test_tpf_xi_theta_v1.cpp`
- Many TPF-specific assertions inside `test_config.cpp`, `test_render_audit.cpp`, and TPF-specific integration smoke scripts.

Classification:
- **F**: package-internal behavior tests should move under TPFCore package tests (or dedicated package test target).
- **B**: engine tests should validate only generic package contract and non-package-specific behavior.

Risk if moved:
- CI wiring may lose coverage unless new package test target is added before test migration.

## Worst boundary violations (priority)

1. **`engine/main.cpp` directly orchestrates TPF theory/utility modes and artifacts.** (A/B)  
2. **Global config model (`engine/config.*`) is saturated with TPF-specific keys and mode semantics.** (A/B)  
3. **Engine render/output layer (`engine/output.*`, `engine/render_audit.*`) encodes TPF route names/diagnostics.** (A/B/E)  
4. **Engine simulation loop has TPF package behavior branches (`simulation.cpp`).** (A/B)  

## Proposed Phase 1 order (no runtime behavior changes yet)

1. **Define generic package metadata/capability interfaces** (B):
   - package-provided diagnostics labels/artifacts
   - package utility-mode dispatch interface
   - package-specific config schema adapter
2. **Move TPF utility-mode execution from `engine/main.cpp` into TPFCore-owned entrypoints** behind generic hooks (A→B).
3. **Move TPF-specific output/render branch labeling to package metadata payloads** while preserving existing external fields for compatibility (A→B).
4. **Isolate TPF config keys behind package config namespace/adapter**, leaving app config with only generic keys (A→B).
5. **Relocate TPF package behavior tests from engine test suite into package test suite** (F).
6. **Clean stale route/config references** in example configs and labels (E).

## Risky move areas (behavior drift risk)

- Run-info/render-manifest field names consumed by Python tools.
- Utility-mode CLI expectations and output artifact filenames.
- Cooling/snapshot skip ordering in stepping loop.
- Compare workflows that currently infer branch identity from engine-side TPF heuristics.

## Per-class rollup (action-oriented)

- **A (must move into TPFCore):** engine main/simulation/config/output/render TPF internals.
- **B (become generic metadata/hook):** package capability dispatch, diagnostics labels, config schema integration, utility-mode launcher interface.
- **C (acceptable package-owned reference):** explicit package example configs and package selection mention in generic docs.
- **D (acceptable historical/audit/doc reference):** prior audits/plans and project-level documentation.
- **E (stale/misleading):** old route names/labels and deprecated correspondence preset references in active examples/tool labels.
- **F (tests to move):** TPF mode/route specific engine tests and integration scripts.
