# TPF 4D Field Core Plan (spike/tpf-4d-field-core)

## Current branch goal
- Stage 0 + Stage 1 + Stage 2 (static evaluator only): add isolated 4D field/tensor math objects and a static 4D source-field evaluator while keeping runtime behavior unchanged.

## Staged migration
- Stage 0 (hygiene): adjust actively misleading comments/docs to match branch status.
- Stage 1: add `Xi4`, `Theta4`, metric-sign helper, 4D trace/double-contraction, and 4D invariant using `LAMBDA_4D`, with unit tests.
- Stage 2 (this change): add static 4D source-field evaluation and static multi-source superposition helpers using full 4D tensor objects with spatial field support; no dynamics routing changes.
- Later stages (not in this task): residual implementation, source-velocity/time-evolution extensions, and integration into runtime routing.

## Not implemented in this branch stage
- No dynamics/time evolution changes.
- No source velocity handling changes.
- No galaxy behavior changes.
- No residual implementation in Stage 2.
- No migration of runtime solver/routing to the Stage 2 evaluator.
- No replacement of existing benchmark/runtime paths.
