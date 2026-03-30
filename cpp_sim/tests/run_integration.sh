#!/usr/bin/env bash
# Run all integration test scripts (app + physics packages). Sets CPP_SIM_ROOT for children.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export CPP_SIM_ROOT="$ROOT"
# Deterministic baseline; overrides developer my.local.cfg when set in main.cpp.
export GALAXY_RUN_CONFIG="$ROOT/../configs/smoke_test.cfg"
cd "$ROOT"

if [[ ! -x ./galaxy_sim ]]; then
  echo "error: build galaxy_sim first (cd cpp_sim && make)" >&2
  exit 1
fi

mapfile -t SCRIPTS < <(
  {
    find "$ROOT/tests/integration" -maxdepth 1 -type f -name '*.sh' ! -name '_env.sh' 2>/dev/null
    find "$ROOT/physics" -type f -path '*/tests/integration/*.sh' 2>/dev/null
  } | sort -u
)

if [[ ${#SCRIPTS[@]} -eq 0 ]]; then
  echo "run_integration.sh: no integration scripts found"
  exit 0
fi

echo "=== cpp_sim integration (${#SCRIPTS[@]} script(s)) ==="
for s in "${SCRIPTS[@]}"; do
  echo "--- $s ---"
  bash "$s"
done
echo "integration: ok"
