# TPF 4D Field Core Plan (spike/tpf-4d-field-core)

## Current branch goal
- Stage 0 + Stage 1 + Stage 2 + Stage 3 + Stage 4 (benchmark/run-output harness): add isolated 4D field/tensor math objects, a static 4D source-field evaluator, a static full spatial-support residual diagnostic, and a benchmark harness that writes run artifacts while keeping runtime behavior unchanged.

## Staged migration
- Stage 0 (hygiene): adjust actively misleading comments/docs to match branch status.
- Stage 1: add `Xi4`, `Theta4`, metric-sign helper, 4D trace/double-contraction, and 4D invariant using `LAMBDA_4D`, with unit tests.
- Stage 2: add static 4D source-field evaluation and static multi-source superposition helpers using full 4D tensor objects with spatial field support; no dynamics routing changes.
- Stage 3 (this change): add static full spatial-support configuration-equation residual diagnostics using the Stage 2 static evaluator; no dynamics/time-evolution/source-velocity/runtime-routing migration.
- Stage 4 (this change): add a `tpf_4d_static_residual_benchmark` harness that runs the Stage 3 static full spatial-support residual diagnostics and writes summary/slice artifacts; no dynamics/source-velocity/time-evolution/runtime-routing migration.
- Later stages (not in this task): source-velocity/time-evolution extensions, and integration into runtime routing.

## Not implemented in this branch stage
- No dynamics/time evolution changes.
- No source velocity handling changes.
- No galaxy behavior changes.
- No dynamics/time evolution implementation in Stage 3 residual diagnostics.
- No migration of runtime solver/routing to the Stage 2 evaluator.
- No replacement of existing benchmark/runtime paths.
- Stage 4 benchmark harness is inspection-only and does not claim validation/proof.
