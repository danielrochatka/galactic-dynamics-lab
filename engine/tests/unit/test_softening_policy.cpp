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
