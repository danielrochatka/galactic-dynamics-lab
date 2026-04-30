#!/usr/bin/env bash
# Run Python unit + integration test discovery suites (non-interactive).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

echo "=== python tests: unit ==="
python3 -m unittest discover -s python_tests/unit -p 'test_*.py' -v

echo "=== python tests: integration ==="
python3 -m unittest discover -s python_tests/integration -p 'test_*.py' -v

echo "python tests: ok"
