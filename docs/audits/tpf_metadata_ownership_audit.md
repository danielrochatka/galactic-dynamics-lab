# TPF Metadata Ownership Audit (Phase 5A/5B)

## Scope
Audit only. No runtime code or schema changes.

Inspected files:
- `engine/output.cpp`, `engine/output.hpp`
- `engine/render_audit.cpp`, `engine/render_audit.hpp`
- `engine/main.cpp`
- `engine/tests/unit/test_render_audit.cpp`
- `engine/tests/unit/test_resolved_scenario.cpp`
- `engine/tests/unit/test_physics_package_hooks.cpp`

## Commands used
```bash
rg -n "TPFCore|tpf_|Xi|xi_|Theta|theta_|VDSG|direct_tpf|legacy_readout|geodesic_correspondence|v11_weak_field_truncation|C_mu_nu|C_μν|run_info_metadata|render_metadata|runtime_metadata|package_config_metadata" \
  engine/output.cpp engine/output.hpp engine/render_audit.cpp engine/render_audit.hpp engine/main.cpp \
  engine/tests/unit/test_render_audit.cpp engine/tests/unit/test_resolved_scenario.cpp engine/tests/unit/test_physics_package_hooks.cpp

rg -n "f << \"tpf_|json_kv\(.*\"tpf_|tf << \"tpf_|TPFCore|legacy_readout|direct_tpf|v11_weak_field|geodesic_correspondence|Xi|Theta|VDSG|C_mu_nu|C_μν" \
  engine/output.cpp engine/render_audit.cpp engine/tests/unit/test_render_audit.cpp engine/tests/unit/test_resolved_scenario.cpp
```

## Classification legend
- **A**: safe package metadata emission to move in Phase 5B
- **B**: requires generic hook shape adjustment before move
- **C**: tied to global Config; wait for config isolation
- **D**: generic output/reporting; leave in engine
- **E**: test-only reference
- **F**: stale/dead/possibly removable later

## Classification table
| file | line/range | emitted key or behavior | current owner | proposed owner | class | move phase | parity risk | notes |
|---|---:|---|---|---|---|---|---|---|
| engine/output.cpp | 112-155,159-175 | `tpf_4d_*benchmark_*` run_info sections | generic output | TPFCore package metadata hook | A | 5B | High | Long contiguous key block/order-sensitive; many tests assert substrings. |
| engine/output.cpp | 177-278 | TPFCore run_info keys (`tpf_dynamics_mode*`, `tpf_kappa`, `tpf_vdsg_*`, readout/provisional labels) | generic output | TPFCore package metadata hook | A | 5B | High | Includes bool-as-0/1 formatting and mode-conditional keys. |
| engine/output.cpp | 367-375 | pipeline keys `tpf_last_*` | generic output + pipeline runtime struct | TPFCore runtime metadata hook | B | 5C | Medium | Existing hook signatures do not carry `AccelPipelineStats*`; shape adjustment needed. |
| engine/output.cpp | 379-445 | duplicate/extended TPFCore run_info block (`tpfcore_*`, repeat notes, diagnostics labels) | generic output | split: package metadata (A) + config/global (C) | C | post-5B | High | Mixture of package labels and global scenario/config provenance; should be split carefully. |
| engine/output.cpp | 86-107 + false-guarded v11 blocks | dead/stale v11 audit text | generic output | remove/defer | F | later | Low | Inactive `if (false)` code; not emitted currently. |
| engine/output.cpp | 54-66 | generic branch/metrics/path keys | generic output | keep generic | D | n/a | Low | Calls generic helper functions and basic run bookkeeping. |
| engine/render_audit.cpp | 241-269,321-339 | manifest TPF keys (`tpf_core_law_mode`, Xi kernel keys, `tpf_runtime_route_status`) in JSON/TXT | generic render | TPFCore render metadata hook | A | 5B | High | Must preserve JSON/TXT pair parity and insertion order. |
| engine/render_audit.cpp | 260-269,383-389 | legacy readout status keys and labels | generic render | TPFCore render metadata hook | A | 5B | Medium | Conditional branch by `legacy_readout`; needs identical status strings. |
| engine/render_audit.cpp | 184-188,370-374 | cooling active + `tpf_cooling_fraction` | generic render | likely config/global (not package-only) | C | later | Medium | Depends on global mode gating and historical dynamics alias logic. |
| engine/render_audit.cpp | 74-81 | alias/effective mode helpers (`geodesic_correspondence`, `v11_weak_field_truncation`) | generic render | keep generic until config isolation | C | later | Medium | Not direct metadata emission hook payload today. |
| engine/render_audit.cpp | 90-170 | `compute_active_*` helper strings | generic render helper | keep generic | D | n/a | Medium | Cross-package diagnostics labeling; not pure package-owned metadata entries yet. |
| engine/render_audit.cpp | false-guarded v11 emission blocks | inactive v11 manifest fields | generic render | remove/defer | F | later | Low | `if (false)` blocks not active. |
| engine/main.cpp | matches are package selection/runtime path, not metadata file emission | simulation orchestration | engine | keep | D | n/a | Low | Out of 5A scope for metadata movement. |
| engine/tests/unit/test_render_audit.cpp | multiple | tests assert TPF keys in render manifest outputs | tests | tests | E | follow 5B | High | Must update/retain expected keys exactly during migration. |
| engine/tests/unit/test_resolved_scenario.cpp | multiple | tests assert run_info TPF keys and notes | tests | tests | E | follow 5B | High | Strong parity guard for run_info lines and notes. |
| engine/tests/unit/test_physics_package_hooks.cpp | metadata hook default-empty checks | tests | tests | E | n/a | Low | Good baseline for default generic hooks. |

## Phase 5B proposed safe scope (A-only)
1. Move **render-manifest TPF-specific JSON/TXT key emission** from `render_audit.cpp` into `TPFCorePackage::render_metadata(...)` and iterate in generic renderer.
2. Move **run_info TPF-specific benchmark and dynamics metadata blocks** from `output.cpp` into `TPFCorePackage::run_info_metadata(...)` and iterate in generic writer.
3. Leave `tpf_last_*` pipeline keys for later (requires hook shape update).
4. Leave config-isolation-coupled alias/effective-mode logic for later.

## A-class parity requirements
- Preserve exact keys and literal values for:
  - `tpf_4d_static_residual_benchmark_*`
  - `tpf_4d_static_motion_readout_benchmark_*`
  - `tpf_4d_xi_motion_probe_benchmark_*`
  - `tpf_dynamics_mode*`, `tpf_runtime_path_tier`, `tpf_core_law_mode`, Xi-kernel keys, `tpf_kappa`, `tpf_vdsg_*`, `tpf_poisson_*`, `tpf_cooling_fraction`, shunt/diagnostic booleans.
- Formatting:
  - run_info bools remain `0/1` (not `true/false`).
  - JSON bools remain true booleans in manifest JSON helper path.
  - Numeric rendering must remain stream-default/scientific behavior as currently emitted per target.
- Ordering:
  - Preserve line order within each moved block; tests rely on substring presence and downstream tools may assume deterministic order.
- Output target parity:
  - `run_info.txt` entries remain tab-delimited key/value lines.
  - `render_manifest.json` and `render_manifest.txt` remain mutually consistent.

## Risks to call out before 5B
- Moving mixed blocks that combine package metadata and global config/provenance may introduce dropped or duplicated keys.
- JSON/TXT dual-output drift risk if one path migrates and the other does not.
- Boolean formatting mismatch risk (`0/1` vs `true/false`) between run_info and manifest outputs.
- Numeric formatting drift risk if package hook serialization differs from existing helper functions.
- Existing unit tests (`test_render_audit`, `test_resolved_scenario`) are sensitive to key names/phrasing.
- Downstream plotting/audit scripts may parse exact keys and stable ordering.

## 5A conclusion
Phase 5A is complete as an audit-only step. No runtime behavior should change in this PR.


## Phase 5B status (render_manifest A-class moved)

Completed in Phase 5B:
- Moved A-class render-manifest TPF metadata emission from `engine/render_audit.cpp` into `TPFCorePackage::render_metadata(...)`.
- Generic render-manifest emission now iterates package-provided metadata entries for both JSON and TXT outputs.
- Preserved key/value parity for moved keys, including legacy-readout status keys/labels.

Deferred to Phase 5C/5D:
- A-class `run_info` metadata migration in `engine/output.cpp`.
- Any config-isolation-dependent C-class items.

Remaining `render_audit` TPF references and class:
- C: cooling-active/cooling-fraction context and alias/effective-mode helpers.
- D: `compute_active_*` generic diagnostics strings.
- F: false-guarded v11 blocks.
