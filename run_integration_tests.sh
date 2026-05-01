#!/usr/bin/env bash
# Run integration smoke scripts only (non-interactive).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

# Ensure non-interactive preflight handling for galaxy scenarios.
export GALAXY_PREFLIGHT_MODE="${GALAXY_PREFLIGHT_MODE:-advisory}"

echo "=== integration tests: engine smoke scripts only ==="
(cd "$ROOT/engine" && make test_integration)

echo "integration tests: ok"
