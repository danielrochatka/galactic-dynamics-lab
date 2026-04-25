# Static 4D residual contamination audit

- **Branch:** `spike/tpf-4d-field-core`
- **Mode audited:** `simulation_mode=tpf_4d_static_residual_benchmark`
- **Audit type:** contamination/sanity evidence record (documentation-only)

## Audited execution path

`main.cpp`
-> `TPFCorePackage::run_4d_static_residual_benchmark(...)`
-> `evaluate_static_configuration_residual_4d(...)`
-> `evaluate_static_sources_field_4d(...)`
-> CSV/summary/plot artifacts
-> return before integration

## Evidence notes

This record reflects code-path inspection and command-driven evidence collection. It does not assert theory validation and does not claim TPF is validated.

---

## A. Confirmed clean

1. No `compute_accelerations` path is used in this benchmark mode.
2. No Newtonian package dynamics are used in this benchmark path.
3. No `direct_tpf` acceleration routing is used in this benchmark path.
4. No kappa/readout scaling is used for particle accelerations in this benchmark path.
5. No VDSG/shunt/cooling path is used in this benchmark path.
6. No particle integration loop is entered for this benchmark mode.
7. Plotting is post-process only (artifact rendering from CSV inputs).
8. CSV field columns use `evaluate_static_sources_field_4d(...)`, not duplicated source formulas in benchmark CSV writer logic.
9. Benchmark path returns after static benchmark artifact writes; no particle integration handoff.

---

## B. Suspicious but harmless

1. `tpf_4d_static_residual_slice.csv` is a legacy alias of the central `xy` view-plane export.
2. `run_info` still lists inherited config knobs (readout/kappa/vdsg/shunt/cooling), but labels them as non-operative for this mode.

---

## C. Actual contamination / must fix

- None found in this audit.

---

## D. Unsupported / not found

1. No benchmark callsite into `compute_accelerations` found for `tpf_4d_static_residual_benchmark`.
2. No Newtonian package invocation found in the static residual evaluator path.

---

## Coupling/G status note (scope clarification)

1. `tpf_4d_static_residual_benchmark` does **not** use Newtonian `G` or `TPF_G_SI` in the static residual/tensor diagnostic path.
2. This benchmark therefore tests uncalibrated geometric/tensor consistency of the source ansatz on the sampled grid.
3. This must **not** be interpreted as a derivation or prediction of physical gravitational coupling `G`.
4. Any coupling-derivation claim would require a later units/source-ledger/observable-coupling benchmark path.

---

## Commands used (evidence collection)

- `rg -n "tpf_4d_static_residual_benchmark|run_4d_static_residual_benchmark|evaluate_static_configuration_residual_4d|static_residual"`
- `nl -ba engine/main.cpp | sed -n '540,640p'`
- `nl -ba engine/main.cpp | sed -n '880,1045p'`
- `nl -ba engine/physics/TPFCore/tpf_core_package.cpp | sed -n '1300,1515p'`
- `nl -ba engine/physics/TPFCore/tpf_4d_static_residual.cpp | sed -n '1,280p'`
- `nl -ba engine/physics/TPFCore/tpf_4d_static_field.cpp | sed -n '1,260p'`
- `nl -ba plot_tpf_4d_static_residual.py | sed -n '1,280p'`
- `nl -ba engine/output.cpp | sed -n '90,190p'`
- `nl -ba engine/resolved_scenario.cpp | sed -n '90,170p'`
- `nl -ba engine/render_audit.cpp | sed -n '70,170p'`
- `git status --short`

## Branch/process note

This is spike-only audit documentation and should not be treated as merge-readiness to `main`.
