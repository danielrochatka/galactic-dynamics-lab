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


namespace {
bool has_warning(const galaxy::GalaxyPreflightSummary& p, const char* warning) {
  return std::find(p.warnings.begin(), p.warnings.end(), warning) != p.warnings.end();
}
}

TEST_CASE("collisionless auto softening compact low-N avoids softening preflight warning") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::galaxy;
  c.softening_mode = "auto";
  c.softening_auto_profile = "collisionless";
  c.auto_softening_dimension = 2;
  c.auto_softening_factor = 1.8;
  c.n_stars = 200;
  c.inner_radius = 2.5e17;
  c.outer_radius = 4.6e18;
  c.galaxy_radius = 4.6e18;

  const galaxy::ResolvedScenario r = galaxy::resolve_scenario(c);
  const auto p = galaxy::build_galaxy_preflight_summary(c, r);

  CHECK(r.config.softening / c.outer_radius <= doctest::Approx(0.05));
  CHECK_FALSE(has_warning(p, "softening_to_outer_radius > 0.05"));
}

TEST_CASE("collisionless auto softening typical galaxy avoids softening preflight warning") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::galaxy;
  c.softening_mode = "auto";
  c.softening_auto_profile = "collisionless";
  c.auto_softening_dimension = 2;
  c.auto_softening_factor = 1.8;
  c.n_stars = 2800;
  c.inner_radius = 2.5e19;
  c.outer_radius = 4.6e20;
  c.galaxy_radius = 4.6e20;

  const galaxy::ResolvedScenario r = galaxy::resolve_scenario(c);
  const auto p = galaxy::build_galaxy_preflight_summary(c, r);

  CHECK(r.config.softening > 0.0);
  CHECK(r.config.softening / c.outer_radius <= doctest::Approx(0.05));
  CHECK_FALSE(has_warning(p, "softening_to_outer_radius > 0.05"));
}

TEST_CASE("collisionless auto explicit max controls preflight softening warning behavior") {
  SUBCASE("safe explicit auto_softening_max") {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::galaxy;
    c.softening_mode = "auto";
    c.softening_auto_profile = "collisionless";
    c.auto_softening_dimension = 2;
    c.auto_softening_factor = 1.8;
    c.n_stars = 200;
    c.inner_radius = 2.5e17;
    c.outer_radius = 4.6e18;
    c.galaxy_radius = 4.6e18;
    c.auto_softening_max = 0.03 * c.outer_radius;

    const galaxy::ResolvedScenario r = galaxy::resolve_scenario(c);
    const auto p = galaxy::build_galaxy_preflight_summary(c, r);

    CHECK(r.config.softening <= doctest::Approx(c.auto_softening_max));
    CHECK(r.config.softening / c.outer_radius <= doctest::Approx(0.05));
    CHECK_FALSE(has_warning(p, "softening_to_outer_radius > 0.05"));
  }

  SUBCASE("unsafe explicit auto_softening_max can trigger warning") {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::galaxy;
    c.softening_mode = "auto";
    c.softening_auto_profile = "collisionless";
    c.auto_softening_dimension = 2;
    c.auto_softening_factor = 1.8;
    c.n_stars = 10;
    c.inner_radius = 2.5e17;
    c.outer_radius = 4.6e18;
    c.galaxy_radius = 4.6e18;
    c.auto_softening_max = 0.10 * c.outer_radius;

    const galaxy::ResolvedScenario r = galaxy::resolve_scenario(c);
    const auto p = galaxy::build_galaxy_preflight_summary(c, r);

    CHECK(r.config.softening == doctest::Approx(c.auto_softening_max));
    CHECK(r.config.softening / c.outer_radius > 0.05);
    CHECK(has_warning(p, "softening_to_outer_radius > 0.05"));
  }
}

TEST_CASE("manual softening above threshold still triggers preflight warning") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::galaxy;
  c.softening_mode = "manual";
  c.n_stars = 2800;
  c.inner_radius = 2.5e19;
  c.outer_radius = 4.6e20;
  c.galaxy_radius = 4.6e20;
  c.softening = 0.10 * c.outer_radius;
  c.explicit_overrides.softening = true;

  const galaxy::ResolvedScenario r = galaxy::resolve_scenario(c);
  const auto p = galaxy::build_galaxy_preflight_summary(c, r);

  CHECK(r.config.softening == doctest::Approx(c.softening));
  CHECK(has_warning(p, "softening_to_outer_radius > 0.05"));
}
