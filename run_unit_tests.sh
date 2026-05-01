#!/usr/bin/env bash
# Run C++ unit/regression doctest target only (non-interactive).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "=== unit tests: engine doctest target only ==="
(cd "$ROOT/engine" && make test_unit)

echo "unit tests: ok"
