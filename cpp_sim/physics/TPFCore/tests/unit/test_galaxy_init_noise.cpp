#include "config.hpp"
#include "doctest.h"
#include "galaxy_init.hpp"
#include "types.hpp"

#include <cmath>

TEST_CASE("same seed reproduces galaxy state") {
  galaxy::Config c;
  c.galaxy_init_template = "symmetric_disk";
  c.galaxy_init_seed = 99991u;
  c.n_stars = 64;
  c.galaxy_radius = 40.0;
  c.galaxy_init_position_noise = 0.0;
  c.galaxy_init_velocity_angle_noise = 0.0;
  c.galaxy_init_velocity_magnitude_noise = 0.0;

  galaxy::State a, b;
  galaxy::initialize_galaxy_disk(c, a, nullptr);
  galaxy::initialize_galaxy_disk(c, b, nullptr);
  REQUIRE(a.n() == b.n());
  for (int i = 0; i < a.n(); ++i) {
    CHECK(a.x[i] == doctest::Approx(b.x[i]));
    CHECK(a.y[i] == doctest::Approx(b.y[i]));
  }
}

TEST_CASE("preformed_spiral: higher master_chaos scales eff noise") {
  galaxy::Config low;
  low.galaxy_init_template = "preformed_spiral";
  low.galaxy_init_seed = 42u;
  low.n_stars = 32;
  low.galaxy_radius = 40.0;
  low.galaxy_init_position_noise = 0.02;
  low.galaxy_init_master_chaos = 1.0;

  galaxy::Config high = low;
  high.galaxy_init_master_chaos = 4.0;

  galaxy::State s1, s2;
  galaxy::GalaxyInitAudit a1, a2;
  galaxy::initialize_galaxy_disk(low, s1, &a1);
  galaxy::initialize_galaxy_disk(high, s2, &a2);

  CHECK(a1.eff_position_noise == doctest::Approx(0.02));
  CHECK(a2.eff_position_noise == doctest::Approx(0.08));
  CHECK(a2.eff_position_noise == doctest::Approx(4.0 * a1.eff_position_noise));
}
