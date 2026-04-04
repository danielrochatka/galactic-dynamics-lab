#include "scenario_defaults.hpp"

namespace galaxy {

ModeScenarioDefaults scenario_defaults_for_mode(SimulationMode mode) {
  ModeScenarioDefaults d;
  switch (mode) {
    case SimulationMode::galaxy:
      d.applies = true;
      d.dt = 0.01;
      d.softening = kDefaultSofteningM;
      d.bh_mass = kDefaultBhMassKg;
      d.enable_star_star_gravity = true;
      d.n_steps = 50000;
      d.snapshot_every = 50;
      d.timing_policy = "galaxy_scale_baseline";
      d.softening_policy = "galaxy_plummer_regularization";
      break;
    case SimulationMode::two_body_orbit:
    case SimulationMode::earth_moon_benchmark:
      d.applies = true;
      d.dt = 3600.0;
      d.softening = 0.0;
      d.bh_mass = 0.0;
      d.enable_star_star_gravity = true;
      d.n_steps = 1440;
      d.snapshot_every = 6;
      d.timing_policy = "earth_moon_hourly_step_60d_horizon";
      d.softening_policy = "exact_two_body_newtonian_no_softening";
      d.validation_earth_mass = kDefaultEarthMassKg;
      d.validation_moon_mass = kDefaultMoonMassKg;
      d.validation_earth_moon_distance = kDefaultEarthMoonDistanceM;
      d.validation_moon_tangential_speed = kDefaultMoonTangentialSpeedMps;
      break;
    case SimulationMode::bh_orbit_validation:
      d.applies = true;
      d.dt = 10000.0;
      d.softening = 0.0;
      d.bh_mass = kDefaultBhMassKg;
      d.enable_star_star_gravity = false;
      d.n_steps = 6000;
      d.snapshot_every = 10;
      d.timing_policy = "single_star_bh_orbit_hour_scale_step_multi_month_horizon";
      d.softening_policy = "single_test_particle_newtonian_no_softening";
      d.validation_two_body_radius = 1.0e13;
      d.validation_two_body_speed_ratio = 1.0;
      break;
    case SimulationMode::symmetric_pair:
      d.applies = true;
      d.dt = 3600.0;
      d.softening = 0.0;
      d.bh_mass = kDefaultBhMassKg;
      d.enable_star_star_gravity = true;
      d.n_steps = 3000;
      d.snapshot_every = 6;
      d.timing_policy = "au_scale_binary_hourly_step_multi_month_horizon";
      d.softening_policy = "symmetric_pair_newtonian_no_softening";
      d.validation_symmetric_include_bh = true;
      d.validation_symmetric_separation = 7.48e10;
      d.validation_symmetric_speed = 3.0e4;
      break;
    case SimulationMode::small_n_conservation:
      d.applies = true;
      d.dt = 1.0e-4;
      d.softening = 0.0;
      d.bh_mass = 1.0e18;
      d.enable_star_star_gravity = true;
      d.n_steps = 20000;
      d.snapshot_every = 20;
      d.timing_policy = "small_n_toy_scale_short_dt_for_orbital_resolution";
      d.softening_policy = "conservation_mode_defaults_to_no_softening";
      d.validation_small_n = 5;
      break;
    case SimulationMode::timestep_convergence:
      d.applies = true;
      d.dt = 10000.0;
      d.softening = 0.0;
      d.bh_mass = kDefaultBhMassKg;
      d.enable_star_star_gravity = false;
      d.n_steps = 6000;
      d.snapshot_every = 10;
      d.timing_policy = "bh_orbit_reference_for_dt_sweeps";
      d.softening_policy = "convergence_baseline_no_softening";
      d.validation_two_body_radius = 1.0e13;
      d.validation_two_body_speed_ratio = 1.0;
      break;
    default:
      break;
  }
  return d;
}

}  // namespace galaxy
