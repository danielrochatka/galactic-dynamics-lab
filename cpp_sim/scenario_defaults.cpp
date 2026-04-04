#include "scenario_defaults.hpp"

namespace galaxy {

ModeScenarioDefaults scenario_defaults_for_mode(SimulationMode mode) {
  ModeScenarioDefaults d;
  switch (mode) {
    case SimulationMode::two_body_orbit:
    case SimulationMode::earth_moon_benchmark:
      d.applies = true;
      d.bh_mass = 0.0;
      d.enable_star_star_gravity = true;
      d.n_steps = 5000;
      d.snapshot_every = 5;
      d.validation_earth_mass = kDefaultEarthMassKg;
      d.validation_moon_mass = kDefaultMoonMassKg;
      d.validation_earth_moon_distance = kDefaultEarthMoonDistanceM;
      d.validation_moon_tangential_speed = kDefaultMoonTangentialSpeedMps;
      break;
    case SimulationMode::bh_orbit_validation:
      d.applies = true;
      d.bh_mass = kDefaultBhMassKg;
      d.enable_star_star_gravity = false;
      d.n_steps = 5000;
      d.snapshot_every = 5;
      d.validation_two_body_radius = kDefaultValidationTwoBodyRadiusM;
      d.validation_two_body_speed_ratio = 1.0;
      break;
    case SimulationMode::symmetric_pair:
      d.applies = true;
      d.bh_mass = kDefaultBhMassKg;
      d.enable_star_star_gravity = true;
      d.n_steps = 5000;
      d.snapshot_every = 5;
      d.validation_symmetric_include_bh = true;
      d.validation_symmetric_separation = 7.48e10;
      d.validation_symmetric_speed = 3.0e4;
      break;
    case SimulationMode::small_n_conservation:
      d.applies = true;
      d.bh_mass = kDefaultBhMassKg;
      d.enable_star_star_gravity = true;
      d.n_steps = 5000;
      d.snapshot_every = 5;
      d.validation_small_n = 5;
      break;
    case SimulationMode::timestep_convergence:
      d.applies = true;
      d.bh_mass = kDefaultBhMassKg;
      d.enable_star_star_gravity = false;
      d.n_steps = 5000;
      d.snapshot_every = 5;
      d.validation_two_body_radius = kDefaultValidationTwoBodyRadiusM;
      d.validation_two_body_speed_ratio = 1.0;
      break;
    default:
      break;
  }
  return d;
}

}  // namespace galaxy
