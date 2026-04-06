# Sourced by integration scripts; run_integration.sh sets CPP_SIM_ROOT.
: "${CPP_SIM_ROOT:?CPP_SIM_ROOT must be set (run via tests/run_integration.sh)}"
cd "$CPP_SIM_ROOT"
if [[ ! -x ./galaxy_sim ]]; then
  echo "error: galaxy_sim not built in $CPP_SIM_ROOT" >&2
  exit 1
fi
