#pragma once

#include "config.hpp"
#include "resolved_scenario.hpp"

#include <string>
#include <vector>

namespace galaxy {

struct GalaxyPreflightSummary {
  double disk_mass = 0.0;
  double disk_to_bh_mass_ratio = 0.0;
  double total_sim_time = 0.0;
  double outer_orbital_period = 0.0;
  double simulated_outer_orbits = 0.0;
  double softening_to_outer_radius = 0.0;
  std::vector<std::string> warnings;
};

GalaxyPreflightSummary build_galaxy_preflight_summary(const Config& config,
                                                       const ResolvedScenario& resolved);

void append_galaxy_preflight_to_run_info(const std::string& output_dir,
                                         const GalaxyPreflightSummary& summary);

}
