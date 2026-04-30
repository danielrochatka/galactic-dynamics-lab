#include "doctest.h"

#include "config.hpp"
#include "preflight.hpp"
#include "resolved_scenario.hpp"

#include <algorithm>

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
  CHECK(p.disk_mass == doctest::Approx(200.0));
  CHECK(p.disk_to_bh_mass_ratio == doctest::Approx(200.0));
  CHECK(p.softening_to_outer_radius == doctest::Approx(0.2));

  CHECK(std::find(p.warnings.begin(), p.warnings.end(), "simulated_outer_orbits < 0.25") != p.warnings.end());
  CHECK(std::find(p.warnings.begin(), p.warnings.end(), "softening_to_outer_radius > 0.05") != p.warnings.end());
  CHECK(std::find(p.warnings.begin(), p.warnings.end(), "star_mass > bh_mass and n_stars > 10") != p.warnings.end());
  CHECK(std::find(p.warnings.begin(), p.warnings.end(), "disk_to_bh_mass_ratio > 1000") == p.warnings.end());
}
