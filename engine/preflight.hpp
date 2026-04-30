#pragma once

#include "config.hpp"
#include "resolved_scenario.hpp"

#include <string>
#include <vector>

namespace galaxy {

struct GalaxyPreflightSummary {
  int n_stars = 0;
  double star_mass = 0.0;
  double star_mass_solar = 0.0;
  double disk_mass = 0.0;
  double disk_mass_solar = 0.0;
  double bh_mass = 0.0;
  double bh_mass_solar = 0.0;
  double disk_to_bh_mass_ratio = 0.0;
  double inner_radius_m = 0.0;
  double inner_radius_kpc = 0.0;
  double outer_radius_m = 0.0;
  double outer_radius_kpc = 0.0;
  double galaxy_radius_m = 0.0;
  double galaxy_radius_kpc = 0.0;
  double effective_softening_m = 0.0;
  double effective_softening_kpc = 0.0;
  double softening_to_outer_radius = 0.0;
  double dt = 0.0;
  int n_steps = 0;
  double total_sim_time_s = 0.0;
  double total_sim_time_myr = 0.0;
  double outer_orbital_period_s = 0.0;
  double outer_orbital_period_myr = 0.0;
  double simulated_outer_orbits = 0.0;
  std::vector<std::string> warnings;
};

GalaxyPreflightSummary build_galaxy_preflight_summary(const Config& config,
                                                       const ResolvedScenario& resolved);
void print_galaxy_preflight_summary(const GalaxyPreflightSummary& summary);
void append_galaxy_preflight_to_run_info(const std::string& output_dir,
                                         const GalaxyPreflightSummary& summary);

}
