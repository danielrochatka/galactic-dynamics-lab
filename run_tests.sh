#!/usr/bin/env bash
# Run full repository test harness via explicit unit + integration entrypoints.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

./run_unit_tests.sh
./run_integration_tests.sh
./run_python_tests.sh

echo ""
echo "All tests finished OK."
