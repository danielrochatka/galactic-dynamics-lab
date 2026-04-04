#include "resolved_scenario.hpp"

#include "galaxy_init.hpp"
#include "init_conditions.hpp"
#include "scenario_defaults.hpp"

namespace galaxy {

namespace {

void apply_mode_defaults(Config& c, const Config& user_cfg) {
  const Config engine_defaults;
  const ModeScenarioDefaults d = scenario_defaults_for_mode(c.simulation_mode);
  if (!d.applies) return;
  if (user_cfg.bh_mass == engine_defaults.bh_mass) c.bh_mass = d.bh_mass;
  if (user_cfg.enable_star_star_gravity == engine_defaults.enable_star_star_gravity) c.enable_star_star_gravity = d.enable_star_star_gravity;
  if (user_cfg.validation_n_steps == engine_defaults.validation_n_steps) c.validation_n_steps = d.n_steps;
  if (user_cfg.validation_snapshot_every == engine_defaults.validation_snapshot_every) c.validation_snapshot_every = d.snapshot_every;
  if (d.validation_two_body_radius > 0.0 &&
      user_cfg.validation_two_body_radius == engine_defaults.validation_two_body_radius)
    c.validation_two_body_radius = d.validation_two_body_radius;
  if (user_cfg.validation_two_body_speed_ratio == engine_defaults.validation_two_body_speed_ratio)
    c.validation_two_body_speed_ratio = d.validation_two_body_speed_ratio;
  if (user_cfg.validation_symmetric_include_bh == engine_defaults.validation_symmetric_include_bh)
    c.validation_symmetric_include_bh = d.validation_symmetric_include_bh;
  if (d.validation_symmetric_separation > 0.0 &&
      user_cfg.validation_symmetric_separation == engine_defaults.validation_symmetric_separation)
    c.validation_symmetric_separation = d.validation_symmetric_separation;
  if (d.validation_symmetric_speed > 0.0 &&
      user_cfg.validation_symmetric_speed == engine_defaults.validation_symmetric_speed)
    c.validation_symmetric_speed = d.validation_symmetric_speed;
  if (d.validation_small_n > 0 && user_cfg.validation_small_n == engine_defaults.validation_small_n)
    c.validation_small_n = d.validation_small_n;
  if (d.validation_earth_mass > 0.0 && user_cfg.validation_earth_mass == engine_defaults.validation_earth_mass)
    c.validation_earth_mass = d.validation_earth_mass;
  if (d.validation_moon_mass > 0.0 && user_cfg.validation_moon_mass == engine_defaults.validation_moon_mass)
    c.validation_moon_mass = d.validation_moon_mass;
  if (d.validation_earth_moon_distance > 0.0 &&
      user_cfg.validation_earth_moon_distance == engine_defaults.validation_earth_moon_distance)
    c.validation_earth_moon_distance = d.validation_earth_moon_distance;
  if (d.validation_moon_tangential_speed > 0.0 &&
      user_cfg.validation_moon_tangential_speed == engine_defaults.validation_moon_tangential_speed)
    c.validation_moon_tangential_speed = d.validation_moon_tangential_speed;
}

}  // namespace

ResolvedScenario resolve_scenario(const Config& input) {
  ResolvedScenario r;
  r.config = input;
  apply_mode_defaults(r.config, input);
  r.mode_label = mode_to_string(r.config.simulation_mode);

  switch (r.config.simulation_mode) {
    case SimulationMode::galaxy:
      init_galaxy_disk(r.config, r.initial_state);
      sync_config_galaxy_init_from_last_audit(r.config);
      r.initializer_used = "init_galaxy_disk";
      r.effective_n_steps = r.config.n_steps;
      r.effective_snapshot_every = r.config.snapshot_every;
      break;
    case SimulationMode::two_body_orbit:
    case SimulationMode::earth_moon_benchmark:
      init_two_body(r.config, r.initial_state);
      r.initializer_used = "init_two_body";
      r.effective_n_steps = r.config.validation_n_steps;
      r.effective_snapshot_every = r.config.validation_snapshot_every;
      break;
    case SimulationMode::bh_orbit_validation:
      init_two_body_star_around_bh(r.config, r.initial_state);
      r.initializer_used = "init_two_body_star_around_bh";
      r.effective_n_steps = r.config.validation_n_steps;
      r.effective_snapshot_every = r.config.validation_snapshot_every;
      break;
    case SimulationMode::symmetric_pair:
      init_symmetric_pair(r.config, r.initial_state);
      if (!r.config.validation_symmetric_include_bh) r.config.bh_mass = 0.0;
      r.initializer_used = "init_symmetric_pair";
      r.effective_n_steps = r.config.validation_n_steps;
      r.effective_snapshot_every = r.config.validation_snapshot_every;
      break;
    case SimulationMode::small_n_conservation:
      init_small_n(r.config, r.initial_state);
      r.config.n_stars = r.initial_state.n();
      r.initializer_used = "init_small_n";
      r.effective_n_steps = r.config.validation_n_steps;
      r.effective_snapshot_every = r.config.validation_snapshot_every;
      break;
    default:
      r.initializer_used = "not_applicable";
      r.effective_n_steps = r.config.n_steps;
      r.effective_snapshot_every = r.config.snapshot_every;
      break;
  }

  return r;
}

}  // namespace galaxy
