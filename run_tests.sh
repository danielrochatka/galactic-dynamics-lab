#!/usr/bin/env bash
# Run the full repository test suite (C++ unit/regression + integration + Python unittest).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "=== cpp_sim (make test) ==="
(cd "$ROOT/cpp_sim" && make test)

echo ""
echo "=== python_tests (unittest) ==="
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
python3 -m unittest discover -s python_tests/unit -p 'test_*.py' -v
python3 -m unittest discover -s python_tests/integration -p 'test_*.py' -v

echo ""
echo "All tests finished OK."
