#include "doctest.h"

#include "config.hpp"
#include "preflight.hpp"
#include "resolved_scenario.hpp"

#include <algorithm>
#include <cmath>

TEST_CASE("galaxy preflight calculations") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::galaxy;
  c.n_stars = 1300;
  c.star_mass = 3.0e34;
  c.bh_mass = 8.55e36;
  c.inner_radius = 1.0e19;
  c.outer_radius = 4.6e20;
  c.galaxy_radius = 4.6e20;
  c.dt = 1.0e13;
  c.softening = 1.254e19;

  galaxy::ResolvedScenario r;
  r.config = c;
  r.config.softening = 1.254e19;
  r.effective_n_steps = 1500;

  const auto p = galaxy::build_galaxy_preflight_summary(c, r);
  CHECK(p.disk_mass == doctest::Approx(3.9e37));
  CHECK(p.total_sim_time_s == doctest::Approx(1.5e16));
  CHECK(p.softening_to_outer_radius == doctest::Approx(1.254e19 / 4.6e20));
  CHECK(p.outer_orbital_period_s > 0.0);
  CHECK(p.total_sim_time_myr == doctest::Approx(p.total_sim_time_s / (31557600.0 * 1.0e6)));
}

TEST_CASE("galaxy preflight warning generation") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::galaxy;
  c.n_stars = 100;
  c.star_mass = 2.0;
  c.bh_mass = 1.0;
  c.dt = 0.01;
  c.n_steps = 1;
  c.outer_radius = 1.0;
  c.softening = 0.2;

  galaxy::ResolvedScenario r;
  r.config = c;
  r.effective_n_steps = c.n_steps;

  const auto p = galaxy::build_galaxy_preflight_summary(c, r);
  CHECK(std::find(p.warnings.begin(), p.warnings.end(), "simulated_outer_orbits < 0.25") != p.warnings.end());
  CHECK(std::find(p.warnings.begin(), p.warnings.end(), "softening_to_outer_radius > 0.05") != p.warnings.end());
  CHECK(std::find(p.warnings.begin(), p.warnings.end(), "star_mass > bh_mass and n_stars > 10") != p.warnings.end());
}
