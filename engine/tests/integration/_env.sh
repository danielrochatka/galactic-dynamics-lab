# Sourced by integration scripts; run_integration.sh sets ENGINE_ROOT.
: "${ENGINE_ROOT:?ENGINE_ROOT must be set (run via tests/run_integration.sh)}"
cd "$ENGINE_ROOT"
if [[ ! -x ./galaxy_sim ]]; then
  echo "error: galaxy_sim not built in $ENGINE_ROOT" >&2
  exit 1
fi
