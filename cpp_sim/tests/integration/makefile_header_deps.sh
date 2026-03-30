#!/usr/bin/env bash
# Regression: Makefile must rebuild translation units when config.hpp changes (Config layout).
# Stale .o files across TUs caused undefined behavior (e.g. std::length_error / SIGSEGV in galaxy init).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Script lives in cpp_sim/tests/integration/ — Makefile is two levels up.
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
make clean >/dev/null
make >/dev/null
touch "$ROOT/config.hpp"
OUT=$(make -n 2>&1)
if ! echo "$OUT" | grep -q 'main[.]cpp'; then
  echo "makefile_header_deps: FAIL — expected main.cpp to rebuild after touching config.hpp" >&2
  echo "$OUT" >&2
  exit 1
fi
if ! echo "$OUT" | grep -q 'galaxy_init[.]cpp'; then
  echo "makefile_header_deps: FAIL — expected galaxy_init.cpp to rebuild after touching config.hpp" >&2
  echo "$OUT" >&2
  exit 1
fi
echo "makefile_header_deps: ok"
