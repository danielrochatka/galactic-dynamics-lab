# TPF Pipeline/Runtime run_info Metadata Audit (Phase 5D-C-A)

Date: 2026-05-10  
Branch: `tpf-xi-theta-v1`  
Scope: audit/design only (no runtime/output/schema/test behavior changes)

## Purpose
This document inventories remaining runtime/pipeline run_info metadata (`tpf_last_*`, pipeline stats, and adjacent runtime counters), classifies ownership/migration readiness, and proposes a generic runtime metadata hook/context design for Phase 5D-C-B.

## Commands executed
```bash
rg -n 'tpf_last_|AccelPipelineStats|tpf_accel_pipeline|pipeline|diagnostics_csv|compute_direct_tpf|xi_runtime|runtime_counters' engine/output.cpp engine/main.cpp engine/accel_pipeline_stats.hpp engine/physics/TPFCore engine/tests

rg -n 'run_info_metadata|run_info_supplement_metadata|render_metadata' engine

cd engine && make test_unit
```

## Files inspected
- `engine/output.cpp`
- `engine/main.cpp`
- `engine/accel_pipeline_stats.hpp`
- `engine/physics/TPFCore/tpf_core_package.hpp`
- `engine/physics/TPFCore/tpf_core_package.cpp`
- `engine/tests` matches for runtime/pipeline metadata assertions

## 1) Exact runtime/pipeline run_info inventory

| key name | current condition | current source value | current output location/order | current owner | desired owner | formatting requirement | dependency type | class | proposed future phase |
|---|---|---|---|---|---|---|---|---|---|
| `tpf_last_mean_baseline_accel_mag` | `physics_package==TPFCore && tpf_pipeline && tpf_pipeline->valid` | `AccelPipelineStats::mean_baseline_mag` | `run_info.txt`, inside `=== TPF acceleration pipeline (last integrator step) ===`, first key in section | `engine/output.cpp` + `AccelPipelineStats` | TPFCore package runtime metadata emission | numeric via `<<` stream default | runtime stats + package internals | E | 5D-C-B |
| `tpf_last_mean_vdsg_accel_mag` | same gate as above | `AccelPipelineStats::mean_vdsg_mag` | same section, second key | generic output | package-owned runtime metadata | numeric stream default | runtime stats + package internals | E | 5D-C-B |
| `tpf_last_vdsg_over_baseline_ratio` | same gate as above | `AccelPipelineStats::vdsg_over_baseline_ratio` | same section, third key | generic output | package-owned runtime metadata | numeric stream default | runtime stats + package internals | E | 5D-C-B |
| `tpf_last_shunt_events` | same gate as above | `AccelPipelineStats::shunt_events_last_step` | same section, fourth key | generic output | package-owned runtime metadata | integer stream default | runtime stats + package internals | E | 5D-C-B |
| `tpf_last_frac_capped` | same gate as above | `AccelPipelineStats::frac_capped_last_step` | same section, fifth key | generic output | package-owned runtime metadata | numeric stream default | runtime stats + package internals | E | 5D-C-B |
| `tpf_last_global_accel_shunt_enabled` | same gate as above | `(AccelPipelineStats::shunt_enabled ? 1 : 0)` | same section, sixth key | generic output | package-owned runtime metadata | bool encoded as `0/1` | runtime stats + package internals | E | 5D-C-B |
| `tpf_last_global_accel_shunt_fraction` | same gate as above | `AccelPipelineStats::shunt_fraction` | same section, seventh key | generic output | package-owned runtime metadata | numeric stream default | runtime stats + package internals | E | 5D-C-B |
| `xi_runtime_theta_evaluations` | `save_run_info && tpf_dynamics_mode=="xi_kernel_deformed"` | `XiRuntimeCounters::theta_evaluations` | appended to `run_info.txt` post-run in `write_post_run_diagnostics`, after diagnostics writes | TPFCore post-run append path | TPFCore package runtime metadata append path (can remain package-local) | integer stream default | runtime stats + package internals | A | keep (5D-C-B no move required) |
| `xi_runtime_invariant_I_evaluations` | same as above | `XiRuntimeCounters::invariant_I_evaluations` | same append block, second key | TPFCore | TPFCore | integer stream default | runtime stats + package internals | A | keep |
| `xi_last_call_pair_evaluations` | same as above | `XiRuntimeCounters::xi_last_call_pair_evaluations` | same append block, third key | TPFCore | TPFCore | integer stream default | runtime stats + package internals | A | keep |
| `xi_total_pair_evaluations` | same as above | `XiRuntimeCounters::xi_total_pair_evaluations` | same append block, fourth key | TPFCore | TPFCore | integer stream default | runtime stats + package internals | A | keep |

### Notes
- `AccelPipelineStats` currently contains extra fields not emitted into run_info (`vdsg_pairs_evaluated`, beta angle stats), so the present inventory is narrower than the struct surface and should remain so unless a new explicit schema decision is made.
- `tpf_accel_pipeline_diagnostics_csv` appears in run_info as config/provenance metadata, not runtime stat metadata; it is out of this migration target slice.

## 2) Classification summary
- **E (package-internal runtime stat; should become package-owned):** all `tpf_last_*` keys.
- **A (safe to move once hook/context exists or already package-owned):** `xi_runtime_*` and xi pair counters are already emitted from package-owned runtime path; no generic-output ownership leak remains there.
- **C/D defer:** config/provenance adjacent keys (e.g., `tpf_accel_pipeline_diagnostics_csv`) stay out of 5D-C-B.

## 3) Hook/context design options for Phase 5D-C-B

### Option 1 — Generic `RunInfoContext` + new runtime hook
**Shape:**
- Add generic context struct in engine layer (no TPF types):
  - common optional runtime fields (e.g., step counts, snapshot counts if needed)
  - opaque package runtime payload handle (type-erased token)
- Add hook:
  - `PhysicsPackage::run_info_runtime_metadata(const Config&, const RunInfoContext&)`

**Evaluation**
- Dependency direction: good (generic -> package interface; package decodes own payload).
- Complexity: medium (new hook + context wiring).
- Output parity risk: low/medium (if section order preserved).
- Keeps `AccelPipelineStats` out of generic interface: **yes**.
- Future packages: clean extensibility for other runtime-stat emitters.

### Option 2 — Package-owned runtime stats object with type-erased interface
**Shape:**
- Package exposes/accepts a type-erased runtime stats provider object; output calls package with that object.

**Evaluation**
- Dependency direction: acceptable but heavier abstraction footprint.
- Complexity: medium/high.
- Output parity risk: medium.
- Keeps `AccelPipelineStats` out of generic interface: **yes**.
- Future packages: flexible but likely over-engineered for current narrow need.

### Option 3 — Keep `output.cpp` fallback until broader config/runtime isolation
**Shape:**
- No new hook now; defer `tpf_last_*` movement.

**Evaluation**
- Dependency direction: unchanged (still leaking runtime type to generic output).
- Complexity: low now, higher debt later.
- Output parity risk: none now.
- Keeps `AccelPipelineStats` out of generic interface: **no**.
- Future packages: poor precedent.

### Recommended design (for 5D-C-B)
**Recommend Option 1**: add `RunInfoContext` + `run_info_runtime_metadata(...)`.

Why:
- Smallest generic addition that isolates runtime-stat transport without exposing `AccelPipelineStats` in generic code.
- Preserves long-term package-owned metadata direction.
- Supports future package runtime metadata without TPF-specific dependencies.

## 4) Exact Phase 5D-C-B implementation scope (narrow)
1. Add generic runtime metadata hook/context only (no config/schema changes).
2. Wire package-owned emission for `tpf_last_*` only.
3. Preserve exact run_info section boundaries, key order, and value formatting (`0/1` for bool key).
4. Keep `xi_runtime_*` append behavior unchanged unless explicitly consolidated with byte-for-byte parity.
5. Leave all config/provenance and compatibility/fallback keys untouched.
6. Do not modify `engine/main.cpp`, `engine/simulation.cpp`, `engine/config.*`, `engine/render_audit.cpp` behavior.

## 5) Parity guardrails for 5D-C-B
- Section header/footer text must remain identical:
  - `=== TPF acceleration pipeline (last integrator step) ===`
  - `=== End TPF acceleration pipeline ===`
- Key order must remain identical to current order in `output.cpp`.
- No key additions/removals/renames.
- No test semantic weakening.

## 6) Test/check result
- `cd engine && make test_unit` — pass in this phase after doc-only edits.

## Phase 5D-C-B status (implemented)
- Implemented `RunInfoContext` + `PhysicsPackage::run_info_runtime_metadata(const Config&, const RunInfoContext&)`.
- Generic `engine/output.cpp` now passes runtime context and emits package-provided runtime metadata entries.
- Moved exactly this section into TPFCore package runtime metadata ownership:
  - Header: `=== TPF acceleration pipeline (last integrator step) ===`
  - Keys in preserved order: `tpf_last_mean_baseline_accel_mag`, `tpf_last_mean_vdsg_accel_mag`, `tpf_last_vdsg_over_baseline_ratio`, `tpf_last_shunt_events`, `tpf_last_frac_capped`, `tpf_last_global_accel_shunt_enabled`, `tpf_last_global_accel_shunt_fraction`
  - Footer: `=== End TPF acceleration pipeline ===`
- Preserved emission condition parity: emitted only when runtime pipeline stats exist and are valid.
- Deferred unchanged: `xi_runtime_*`, `tpf_accel_pipeline_diagnostics_csv`, config/provenance fields.
