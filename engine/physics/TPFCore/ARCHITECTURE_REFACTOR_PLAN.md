# TPFCore architecture diagnosis and refactor plan (initial pass)

## Why this document exists
This plan is a corrective response to drift where TPFCore routing and mode flags ended up controlling *which formulas* run, not just *which parameters* are fed into canonical formulas.

This first pass is intentionally conservative: it reorganizes source routing and records canonicalization rules without changing declared physical models.

## A) Architecture diagnosis

### 1) Duplicate/near-duplicate math paths currently present

- **Source traversal logic duplicated in multiple files**:
  - `field_evaluation.cpp` had its own BH + star-star superposition loops.
  - `provisional_readout.cpp` had separate BH + star-star loops for tensor-radial readout.
  - `tpf_core_package.cpp` had another BH + star-star loop pair in VDSG modifier code.
  - These loops represent the same conceptual object: “iterate active gravitational sources for target i”.

- **Theta/I usage split across multiple branch-specific entry points**:
  - `compute_direct_tpf_accelerations` uses `evaluate_provisional_field_multi_source` + direct baseline.
  - readout closures perform per-mode recomputation or per-source evaluation.
  - weak-field correspondence has an explicitly separate implementation (kept separate intentionally; different model family).

### 2) Config and routing semantics doing too much

- `tpf_dynamics_mode` currently selects fundamentally different acceleration constructions:
  - `legacy_readout` (provisional readout + optional VDSG + optional shunt)
  - `direct_tpf` (principal tensor baseline + optional VDSG)
  - `v11_weak_field_truncation` (paper correspondence helper)
- This is valid for model-family selection, but historically drift happened because minor knobs were also treated like hidden route switches.

### 3) Branch-specific implementations that should collapse

- **Source iteration** should not be branch-owned. It is now canonicalized (initial pass) via a shared iterator helper.
- **Field object construction** is already partly canonical (`source_ansatz` + `field_evaluation`); follow-up should continue collapsing closure-local re-derivations where mathematically equivalent.

### 4) Places representing same conceptual object in multiple places

- “Active source set for particle i” used to be represented in three local loop implementations:
  - `field_evaluation.cpp`
  - `provisional_readout.cpp`
  - `tpf_core_package.cpp`
- This pass unifies that object under `source_iteration.hpp`.

## B) Target architecture proposal

The target structure should make the math layers explicit:

- **`source_*` layer** (canonical field/source primitives)
  - `source_ansatz.*`
  - `source_iteration.hpp` (new canonical source routing primitive)

- **`field_*` layer** (canonical mathematical objects)
  - `field_evaluation.*` for Xi/Theta/I packaging
  - direct baseline objects (`direct_tpf_baseline.*`) for principal tensor mapping

- **`readout_*` layer** (readout models only)
  - `provisional_readout.*` for explicit readout models
  - each readout mode should stay a named model, not a hidden branch

- **`extensions_*` layer** (optional explicit effects)
  - VDSG additive extension logic (currently in `tpf_core_package.cpp`, candidate extraction in follow-up)

- **`runtime/orchestration` layer**
  - `tpf_core_package.*` should orchestrate model selection and compose canonical pieces.

## C) Canonicalization rules

1. **One canonical function per mathematical object** (Xi/Theta/I/principal tensor/source iterator).
2. **Configs pass parameters, not stealth formula switches**, except explicit model-family selectors.
3. **Defaults are centralized once**; computational functions require explicit inputs.
4. **Routing composes canonical functions** instead of re-implementing formulas in branch code.
5. **New model families must be explicit and isolated** (e.g., correspondence helper, VDSG extension).
6. **If equivalence is uncertain, do not merge blindly**; label as variant or separate model until proven.

## D) Initial refactor pass completed in this change

- Introduced `source_iteration.hpp` as the canonical source enumerator for BH + star-star interactions.
- Rewired:
  - `field_evaluation.cpp`
  - tensor-radial path in `provisional_readout.cpp`
  - `accumulate_vdsg_velocity_modifier` in `tpf_core_package.cpp`
- Result: a single canonical source traversal primitive now feeds multiple math/readout/extension paths.

## E) Transparency report

### Moved
- Source enumeration logic moved out of local lambdas/loops into `source_iteration.hpp`.

### Merged
- Equivalent BH/star routing branches in field-evaluation, tensor-radial readout, and VDSG modifier now share one canonical iterator.

### Remaining duplicated or near-duplicated areas (explicitly left for follow-up)
- Derived radial closure and tensor-radial closure still compute different readout formulas and remain separate.
- Weak-field correspondence implementation remains separate from direct baseline by design (distinct model family).
- Some Theta-derived projections appear in multiple closure functions; these require careful equivalence proof before merge.

### Intentionally separate because math is genuinely distinct
- `direct_tpf` principal-part baseline vs `v11_weak_field_truncation` correspondence helper.
- Provisional readout models (`tensor_radial_projection`, `derived_tpf_radial_readout`, `experimental_radial_r_scaling`) as explicitly different readout models.

### Behavior-sensitive areas to revalidate after this pass
- VDSG additive modifier magnitude/sign for star-star interactions (routing unified, formulas intended unchanged).
- Tensor-radial mode acceleration parity against previous snapshots/regression baselines.
- Diagnostics depending on source ordering should be spot-checked (if any downstream code assumed ad hoc loop ordering).

## F) Step 3 readout-family separation (this pass)

### Short diagnosis of readout drift before this pass

- `provisional_readout.cpp` mixed three distinct readout model families in one implementation unit:
  - tensor-radial projection readout,
  - derived-radial (TR coherence / derived TPF radial) readout,
  - experimental radial r-scaling readout.
- Shared orchestration and per-family formulas lived side-by-side, which made it harder to audit what was:
  - canonical shared field packaging,
  - mode-family-specific math,
  - diagnostics-only bookkeeping.
- Helper boundaries were implicit (`static` functions in one file), so tracing “which model family is active” required navigating a monolithic file.

### Structural separation completed in this pass

- Separated readout families into explicit implementation units:
  - `readout_tensor_radial.cpp`
  - `readout_derived_radial.cpp`
  - `readout_experimental_radial_scaling.cpp`
- Added a small readout-family interface header:
  - `readout_model_families.hpp`
- Kept `provisional_readout.cpp` as readout orchestration + dispatch + debug CSV wiring.

### What remains canonical/shared (intentionally)

- Canonical source traversal remains shared through `source_iteration.hpp`.
- Canonical field packaging remains shared through `field_evaluation.*` (`as_canonical_field_objects`, `build_canonical_field_objects`, Theta invariants).
- Derived radial gravity profile construction remains shared in `derived_tpf_radial.*`.

### Explicitly left for next step

- `tpf_core_package.*` still contains broad orchestration and mode-selection/reporting responsibilities.
- Optional extension/model-selection plumbing (including VDSG coupling branches) is still primarily package-owned.
- Next step should continue shrinking `tpf_core_package.*` toward orchestration-only responsibilities and extract additional optional extension logic as safe.
