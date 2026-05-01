# TPF 4D Field Core Plan (spike/tpf-4d-field-core)

## Current branch goal
- Stage 0 + Stage 1 + Stage 2 + Stage 3 + Stage 4 + Stage 5 (visualization/export only): add isolated 4D field/tensor math objects, a static 4D source-field evaluator, a static full spatial-support residual diagnostic, and a benchmark harness that writes run artifacts while keeping runtime behavior unchanged.

## Staged migration
- Stage 0 (hygiene): adjust actively misleading comments/docs to match branch status.
- Stage 1: add `Xi4`, `Theta4`, metric-sign helper, 4D trace/double-contraction, and 4D invariant using `LAMBDA_4D`, with unit tests.
- Stage 2: add static 4D source-field evaluation and static multi-source superposition helpers using full 4D tensor objects with spatial field support; no dynamics routing changes.
- Stage 3 (this change): add static full spatial-support configuration-equation residual diagnostics using the Stage 2 static evaluator; no dynamics/time-evolution/source-velocity/runtime-routing migration.
- Stage 4 (this change): add a `tpf_4d_static_residual_benchmark` harness that runs the Stage 3 static full spatial-support residual diagnostics and writes summary/slice artifacts; no dynamics/source-velocity/time-evolution/runtime-routing migration.
- Stage 5 (this change): add richer view-plane CSV exports, source marker CSV, and optional post-process plotting script wiring behind `--plot`; view-plane diagnostic renderings are produced from the full spatial-support static 4D residual evaluation (no residual/source-field math changes, no dynamics/source-worldline/time-evolution migration in this stage).
- Stage 4 observed benchmark behavior (current tests): monopole, equal-mass bonded pair, and unequal-mass bonded pair residual diagnostics are available and have shown controlled exterior residual behavior under current static benchmark runs. This is observed benchmark behavior, not proof/validation of full dynamics.
- Later stages (not in this task): source-velocity/time-evolution extensions, and integration into runtime routing.

## Not implemented in this branch stage
- No dynamics/time evolution changes.
- No source velocity handling changes.
- No source velocities in the Stage 4 benchmark path.
- No galaxy behavior changes.
- No dynamics/time evolution implementation in Stage 3 residual diagnostics.
- No migration of runtime solver/routing to the Stage 2 evaluator.
- No runtime migration to Xi4/Theta4 benchmark evaluator yet.
- No DeltaC closure implementation in this stage.
- No replacement of existing benchmark/runtime paths.
- Stage 4/5 benchmark harness provides static-field diagnostic evidence and residual benchmark evidence; dynamics/moving-sources/orbital/coupling/DeltaC benchmarks remain future work.
