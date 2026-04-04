#ifndef GALAXY_SCENARIO_DEFAULTS_HPP
#define GALAXY_SCENARIO_DEFAULTS_HPP

#include "config.hpp"

namespace galaxy {

struct ModeScenarioDefaults {
  bool applies = false;
  double dt = 0.0;
  double softening = 0.0;
  double bh_mass = 0.0;
  bool enable_star_star_gravity = true;
  int n_steps = 0;
  int snapshot_every = 1;
  const char* timing_policy = "";
  const char* softening_policy = "";

  double validation_two_body_radius = 0.0;
  double validation_two_body_speed_ratio = 1.0;
  bool validation_symmetric_include_bh = true;
  double validation_symmetric_separation = 0.0;
  double validation_symmetric_speed = 0.0;
  int validation_small_n = 0;

  double validation_earth_mass = 0.0;
  double validation_moon_mass = 0.0;
  double validation_earth_moon_distance = 0.0;
  double validation_moon_tangential_speed = 0.0;
};

ModeScenarioDefaults scenario_defaults_for_mode(SimulationMode mode);

}  // namespace galaxy

#endif
