# Testing

This repository separates implementation checks from scientific interpretation. The automated suite verifies math consistency, config parsing, branch labeling, manifests, reproducibility, and small end-to-end runs.

## Layers

| Layer | Purpose | Typical contents |
|-------|---------|------------------|
| **Unit** | Pure or localized logic: parsing, helpers, deterministic math | `cpp_sim/tests/unit/`, `physics/*/tests/unit/`, `python_tests/unit/` |
| **Integration** | Tiny real runs, subprocess smoke, outputs on disk | `cpp_sim/tests/integration/*.sh`, `physics/*/tests/integration/*.sh`, `python_tests/integration/` |
| **Regression** | Frozen expectations to catch silent drift (stats, labels) | `physics/TPFCore/tests/regression/` |
| **Future: observation / validation** | Not part of this suite; would compare to external data or long campaigns | — |

## Layout

- **App / framework (C++):** `cpp_sim/tests/unit/`, `cpp_sim/tests/integration/`
- **Physics packages:** `cpp_sim/physics/<Name>/tests/{unit,integration,regression}/`
- **Python:** `python_tests/{unit,integration}/` (manifests, overlays, loaders — not `tests/` at repo root, to avoid confusion with C++)

New physics packages can add their own `tests/` trees; root commands remain the single entrypoint. `cpp_sim/Makefile` auto-discovers doctest sources from `cpp_sim/tests/{unit,regression}` plus `cpp_sim/physics/*/tests/{unit,regression}` via wildcards + `find`, and `cpp_sim/tests/run_integration.sh` auto-discovers shell scripts from both `cpp_sim/tests/integration` and `cpp_sim/physics/*/tests/integration`.


## Test ownership / location policy

- Package-owned behavior belongs under `cpp_sim/physics/<Package>/tests/` (unit/integration/regression).
- Root C++ tests in `cpp_sim/tests/` are reserved for app-level behavior: config parsing, package registry/dispatch, CLI/main behavior, cross-package workflows, and generic infrastructure.
- For TPFCore specifically, routing/closure labels, TPF-specific manifests/audit labels, and TPF-specific smoke/regression checks should live under `cpp_sim/physics/TPFCore/tests/`.

Classification rule: if removing one package makes a test meaningless, that test belongs in that package.

## Frameworks

- **C++:** [doctest](https://github.com/doctest/doctest) (single header: `cpp_sim/vendor/doctest.h`). Fast, header-only, low ceremony.
- **Python:** `unittest` (stdlib) — no extra pip dependency for CI.

## Run configs for CI

Integration scripts set **`GALAXY_RUN_CONFIG`** to `configs/smoke_test.cfg` so runs do not depend on a developer’s `configs/my.local.cfg`. The simulator checks this environment variable first in `find_run_config_path()`.

## Commands

| Goal | Command |
|------|---------|
| **Everything** | `./run_tests.sh` (from repo root) |
| **cpp_sim only (unit + regression + integration)** | `(cd cpp_sim && make test)` |
| **cpp_sim doctest binary only** | `(cd cpp_sim && make test_unit)` |
| **cpp_sim integration shell scripts only** | `(cd cpp_sim && make test_integration)` |
| **Newtonian package unit tests** | Built inside `galaxy_tests`; filter: `./cpp_sim/galaxy_tests -tc="*Newtonian*"` (doctest) |
| **TPFCore package tests** | Same binary; filter e.g. `-tc="*source_ansatz*"` or run full `./cpp_sim/galaxy_tests` |
| **Python only** | `PYTHONPATH=. python3 -m unittest discover -s python_tests/unit -p 'test_*.py' -v` then the same with `-s python_tests/integration` |

## Honesty

- Passing tests means the **code behaves as encoded** under the stated assumptions.
- Failing tests mean regression or inconsistency in the implemented code path.
- Treat test results as software verification only.
