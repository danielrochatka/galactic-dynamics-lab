#include "config.hpp"
#include "doctest.h"
#include "softening_policy.hpp"

#include <cmath>

TEST_CASE("plummer_softening_scale matches historical pow formulation") {
  const double r_sq = 9.0;
  const double eps = 2.0;
  const double eps2 = eps * eps;
  const double softened = r_sq + eps2;
  const double expected = std::pow(r_sq / softened, 1.5);

  const double got = galaxy::plummer_softening_scale(r_sq, eps);
  CHECK(got == doctest::Approx(expected));
}

TEST_CASE("plummer_softening_scale handles non-positive softening and radius-squared guard") {
  CHECK(galaxy::plummer_softening_scale(4.0, 0.0) == doctest::Approx(1.0));
  CHECK(galaxy::plummer_softening_scale(4.0, -1.0) == doctest::Approx(1.0));

  CHECK(galaxy::plummer_softening_scale(0.0, 1.0) == doctest::Approx(0.0));
  CHECK(galaxy::plummer_softening_scale(-4.0, 1.0) == doctest::Approx(0.0));
}

TEST_CASE("collisionless auto softening defaults to 2D when dimension is unset") {
  galaxy::Config cfg;
  cfg.softening_mode = "auto";
  cfg.softening_auto_profile = "collisionless";
  cfg.simulation_mode = galaxy::SimulationMode::galaxy;
  cfg.auto_softening_dimension = 0;
  cfg.galaxy_radius = 100.0;
  cfg.n_stars = 1000;

  const galaxy::ResolvedSoftening r = galaxy::resolve_softening(cfg, galaxy::State{});
  CHECK(r.dimension == 2);
}

TEST_CASE("collisionless auto softening explicit dimension overrides default") {
  galaxy::Config cfg;
  cfg.softening_mode = "auto";
  cfg.softening_auto_profile = "collisionless";
  cfg.auto_softening_dimension = 4;
  cfg.galaxy_radius = 100.0;
  cfg.n_stars = 1000;

  const galaxy::ResolvedSoftening r = galaxy::resolve_softening(cfg, galaxy::State{});
  CHECK(r.dimension == 4);
}

TEST_CASE("collisionless auto softening applies contextual cap for compact low-N disk") {
  galaxy::Config cfg;
  cfg.softening_mode = "auto";
  cfg.softening_auto_profile = "collisionless";
  cfg.inner_radius = 0.0;
  cfg.outer_radius = 10.0;
  cfg.galaxy_radius = 10.0;
  cfg.n_stars = 1;

  const galaxy::ResolvedSoftening r = galaxy::resolve_softening(cfg, galaxy::State{});
  CHECK(r.effective_softening <= doctest::Approx(0.05 * r.radius_outer_used));
  CHECK(r.max_capped);
  CHECK(r.max_cap_source == "contextual");
}

TEST_CASE("auto_softening_max smaller than structural result caps downward") {
  galaxy::Config cfg;
  cfg.softening_mode = "auto";
  cfg.softening_auto_profile = "collisionless";
  cfg.simulation_mode = galaxy::SimulationMode::galaxy;
  cfg.outer_radius = 10.0;
  cfg.galaxy_radius = 10.0;
  cfg.inner_radius = 1.0;
  cfg.n_stars = 1;
  cfg.auto_softening_max = 0.02;

  const galaxy::ResolvedSoftening r = galaxy::resolve_softening(cfg, galaxy::State{});
  CHECK(r.effective_softening == doctest::Approx(cfg.auto_softening_max));
  CHECK(r.max_cap_source == "contextual+explicit");
}

TEST_CASE("normal collisionless auto softening remains unchanged below cap") {
  galaxy::Config cfg;
  cfg.softening_mode = "auto";
  cfg.softening_auto_profile = "collisionless";
  cfg.simulation_mode = galaxy::SimulationMode::two_body_orbit;
  cfg.galaxy_radius = 100.0;
  cfg.outer_radius = 100.0;
  cfg.inner_radius = 1.0;
  cfg.n_stars = 1000000;

  const galaxy::ResolvedSoftening r = galaxy::resolve_softening(cfg, galaxy::State{});
  CHECK(r.effective_softening < 0.05 * r.radius_outer_used);
  CHECK_FALSE(r.max_capped);
  CHECK(r.max_cap_source.empty());
}

TEST_CASE("auto_softening_min floors value and remains visible in metadata") {
  galaxy::Config cfg;
  cfg.softening_mode = "auto";
  cfg.softening_auto_profile = "collisionless";
  cfg.galaxy_radius = 100.0;
  cfg.outer_radius = 100.0;
  cfg.inner_radius = 1.0;
  cfg.n_stars = 1000000;
  cfg.auto_softening_min = 10.0;

  const galaxy::ResolvedSoftening r = galaxy::resolve_softening(cfg, galaxy::State{});
  CHECK(r.effective_softening == doctest::Approx(10.0));
  CHECK(r.min_floored);
  CHECK(r.source.find("min_floor=explicit") != std::string::npos);
}

TEST_CASE("collisionless auto softening respects inner structural cap for typical galaxy") {
  galaxy::Config cfg;
  cfg.softening_mode = "auto";
  cfg.softening_auto_profile = "collisionless";
  cfg.simulation_mode = galaxy::SimulationMode::galaxy;
  cfg.auto_softening_dimension = 2;
  cfg.auto_softening_factor = 1.8;
  cfg.inner_radius = 2.5e19;
  cfg.outer_radius = 4.6e20;
  cfg.galaxy_radius = cfg.outer_radius;
  cfg.n_stars = 1100;

  const galaxy::ResolvedSoftening r = galaxy::resolve_softening(cfg, galaxy::State{});
  CHECK(r.effective_softening <= doctest::Approx(0.05 * cfg.inner_radius));
}

TEST_CASE("stellar_physical and nuclear_cluster auto profiles are not inner-radius capped") {
  galaxy::Config stellar;
  stellar.softening_mode = "auto";
  stellar.softening_auto_profile = "stellar_physical";
  stellar.simulation_mode = galaxy::SimulationMode::galaxy;
  stellar.auto_softening_factor = 10.0;
  stellar.inner_radius = 1.0;
  const galaxy::ResolvedSoftening s = galaxy::resolve_softening(stellar, galaxy::State{});
  CHECK(s.effective_softening == doctest::Approx(10.0 * 6.957e8));

  galaxy::Config cluster;
  cluster.softening_mode = "auto";
  cluster.softening_auto_profile = "nuclear_cluster";
  cluster.simulation_mode = galaxy::SimulationMode::galaxy;
  cluster.auto_softening_factor = 3.0;
  cluster.inner_radius = 1.0;
  const galaxy::ResolvedSoftening n = galaxy::resolve_softening(cluster, galaxy::State{});
  CHECK(n.effective_softening == doctest::Approx(3.0 * 6.957e8));
}

TEST_CASE("auto_softening_max larger than structural result does not raise softening") {
  galaxy::Config cfg;
  cfg.softening_mode = "auto";
  cfg.softening_auto_profile = "collisionless";
  cfg.simulation_mode = galaxy::SimulationMode::galaxy;
  cfg.auto_softening_dimension = 2;
  cfg.auto_softening_factor = 1.8;
  cfg.inner_radius = 2.5e19;
  cfg.outer_radius = 4.6e20;
  cfg.galaxy_radius = cfg.outer_radius;
  cfg.n_stars = 1100;
  cfg.auto_softening_max = 0.10 * cfg.outer_radius;

  const galaxy::ResolvedSoftening r = galaxy::resolve_softening(cfg, galaxy::State{});
  CHECK(r.effective_softening <= doctest::Approx(0.05 * cfg.inner_radius));
  CHECK(r.max_cap_source == "contextual");
}
