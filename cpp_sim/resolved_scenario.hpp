#ifndef GALAXY_RESOLVED_SCENARIO_HPP
#define GALAXY_RESOLVED_SCENARIO_HPP

#include "config.hpp"
#include "types.hpp"
#include <string>

namespace galaxy {

struct ResolvedScenario {
  Config config;
  State initial_state;
  std::string initializer_used;
  std::string mode_label;
  int effective_n_steps = 0;
  int effective_snapshot_every = 1;
  double effective_total_sim_time = 0.0;
  std::string timing_policy;
  std::string softening_policy;
};

ResolvedScenario resolve_scenario(const Config& input);

}  // namespace galaxy

#endif
