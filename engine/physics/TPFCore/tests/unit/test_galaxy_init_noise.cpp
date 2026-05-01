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

TEST_CASE("BH-only initializer ignores star mass when star-star gravity disabled") {
  galaxy::Config low_mass;
  low_mass.galaxy_init_template = "preformed_spiral";
  low_mass.physics_package = "Newtonian";
  low_mass.enable_star_star_gravity = false;
  low_mass.galaxy_init_seed = 12345u;
  low_mass.n_stars = 48;
  low_mass.galaxy_radius = 100.0;
  low_mass.bh_mass = 1.0e30;
  low_mass.star_mass = 1.0e12;
  low_mass.galaxy_init_position_noise = 0.0;
  low_mass.galaxy_init_velocity_angle_noise = 0.0;
  low_mass.galaxy_init_velocity_magnitude_noise = 0.0;
  low_mass.velocity_noise = 0.0;

  galaxy::Config high_mass = low_mass;
  high_mass.star_mass = 1.0e24;

  galaxy::State low_state, high_state;
  galaxy::GalaxyInitAudit low_audit, high_audit;
  galaxy::initialize_galaxy_disk(low_mass, low_state, &low_audit);
  galaxy::initialize_galaxy_disk(high_mass, high_state, &high_audit);

  REQUIRE(low_state.n() == high_state.n());
  CHECK(low_audit.velocity_mass_model == "bh_only");
  CHECK(high_audit.velocity_mass_model == "bh_only");
  CHECK(low_audit.velocity_uses_star_mass == false);
  CHECK(high_audit.velocity_uses_star_mass == false);

  for (int i = 0; i < low_state.n(); ++i) {
    CHECK(low_state.x[i] == doctest::Approx(high_state.x[i]));
    CHECK(low_state.y[i] == doctest::Approx(high_state.y[i]));
    CHECK(low_state.vx[i] == doctest::Approx(high_state.vx[i]));
    CHECK(low_state.vy[i] == doctest::Approx(high_state.vy[i]));
  }
}

TEST_CASE("initializer includes enclosed stellar mass when star-star gravity enabled") {
  galaxy::Config low_mass;
  low_mass.galaxy_init_template = "preformed_spiral";
  low_mass.physics_package = "Newtonian";
  low_mass.enable_star_star_gravity = true;
  low_mass.galaxy_init_seed = 67890u;
  low_mass.n_stars = 48;
  low_mass.galaxy_radius = 100.0;
  low_mass.bh_mass = 1.0e25;
  low_mass.star_mass = 1.0e22;
  low_mass.galaxy_init_position_noise = 0.0;
  low_mass.galaxy_init_velocity_angle_noise = 0.0;
  low_mass.galaxy_init_velocity_magnitude_noise = 0.0;
  low_mass.velocity_noise = 0.0;

  galaxy::Config high_mass = low_mass;
  high_mass.star_mass = 1.0e24;

  galaxy::State low_state, high_state;
  galaxy::GalaxyInitAudit low_audit, high_audit;
  galaxy::initialize_galaxy_disk(low_mass, low_state, &low_audit);
  galaxy::initialize_galaxy_disk(high_mass, high_state, &high_audit);

  REQUIRE(low_state.n() == high_state.n());
  CHECK(low_audit.velocity_mass_model == "bh_plus_enclosed_stars");
  CHECK(high_audit.velocity_mass_model == "bh_plus_enclosed_stars");
  CHECK(low_audit.velocity_uses_star_mass == true);
  CHECK(high_audit.velocity_uses_star_mass == true);

  bool found_velocity_difference = false;
  for (int i = 0; i < low_state.n(); ++i) {
    if (std::abs(low_state.vx[i] - high_state.vx[i]) > 1e-12 ||
        std::abs(low_state.vy[i] - high_state.vy[i]) > 1e-12) {
      found_velocity_difference = true;
      break;
    }
  }
  CHECK(found_velocity_difference);
}

TEST_CASE("derived TPF initializer audit truthfully reports star-mass usage vs toggle state") {
  galaxy::Config c;
  c.galaxy_init_template = "preformed_spiral";
  c.physics_package = "TPFCore";
  c.tpfcore_enable_provisional_readout = true;
  c.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c.galaxy_init_seed = 2026u;
  c.n_stars = 24;
  c.galaxy_radius = 100.0;
  c.galaxy_init_position_noise = 0.0;
  c.galaxy_init_velocity_angle_noise = 0.0;
  c.galaxy_init_velocity_magnitude_noise = 0.0;
  c.velocity_noise = 0.0;

  galaxy::State s_false, s_true;
  galaxy::GalaxyInitAudit a_false, a_true;

  c.enable_star_star_gravity = false;
  galaxy::initialize_galaxy_disk(c, s_false, &a_false);
  CHECK(a_false.velocity_uses_star_mass == true);
  CHECK(a_false.velocity_respects_star_star_flag == false);
  CHECK(a_false.velocity_mass_model == "derived_tpf_profile_includes_stars_toggle_not_applied");

  c.enable_star_star_gravity = true;
  galaxy::initialize_galaxy_disk(c, s_true, &a_true);
  CHECK(a_true.velocity_uses_star_mass == true);
  CHECK(a_true.velocity_respects_star_star_flag == true);
  CHECK(a_true.velocity_mass_model == "derived_tpf_profile_includes_stars");
}
