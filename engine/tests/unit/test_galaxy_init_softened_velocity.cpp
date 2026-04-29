#include "doctest.h"

#include "config.hpp"
#include "galaxy_init.hpp"
#include "physics/Newtonian/newtonian.hpp"
#include "types.hpp"

#include <algorithm>
#include <cmath>

using galaxy::Config;
using galaxy::State;

namespace {
constexpr double G_SI = 6.6743e-11;

Config base_config() {
  Config c;
  c.n_stars = 1;
  c.galaxy_radius = 1.0e3;
  c.bh_mass = 2.0e30;
  c.star_mass = 3.0e25;
  c.initial_velocity_scale = 1.0;
  c.enable_star_star_gravity = false;
  c.velocity_noise = 0.0;
  c.galaxy_init_position_noise = 0.0;
  c.galaxy_init_velocity_angle_noise = 0.0;
  c.galaxy_init_velocity_magnitude_noise = 0.0;
  c.galaxy_init_clumpiness = 0.0;
  c.galaxy_init_m2_amplitude = 0.0;
  c.galaxy_init_m3_amplitude = 0.0;
  c.galaxy_init_bar_amplitude = 0.0;
  c.galaxy_init_spiral_amplitude = 0.0;
  c.galaxy_init_seed = 42;
  return c;
}

double speed0(const State& s) { return std::hypot(s.vx[0], s.vy[0]); }

double radius0(const State& s) { return std::hypot(s.x[0], s.y[0]); }

}  // namespace

TEST_CASE("galaxy init softened enclosed-mass speed reduces to legacy at eps=0") {
  Config c = base_config();
  c.softening = 0.0;
  State s;
  galaxy::initialize_galaxy_disk(c, s, nullptr);
  const double r = radius0(s);
  const double v_expected = std::sqrt(G_SI * c.bh_mass / r);
  CHECK(speed0(s) == doctest::Approx(v_expected).epsilon(1e-12));
}

TEST_CASE("galaxy init speed lowers monotonically with larger softening") {
  Config c = base_config();
  State s0, s1, s2;

  c.softening = 0.0;
  galaxy::initialize_galaxy_disk(c, s0, nullptr);
  c.softening = 10.0;
  galaxy::initialize_galaxy_disk(c, s1, nullptr);
  c.softening = 100.0;
  galaxy::initialize_galaxy_disk(c, s2, nullptr);

  const double v0 = speed0(s0), v1 = speed0(s1), v2 = speed0(s2);
  CHECK(v1 < v0);
  CHECK(v2 < v1);
}

TEST_CASE("galaxy init single-BH circular speed matches runtime softened radial acceleration") {
  Config c = base_config();
  c.softening = 50.0;
  c.enable_star_star_gravity = false;

  State s;
  galaxy::initialize_galaxy_disk(c, s, nullptr);

  std::vector<double> ax, ay;
  galaxy::NewtonianPackage pkg;
  pkg.compute_accelerations(s, c.bh_mass, c.softening, false, ax, ay);

  const double r = radius0(s);
  const double v = speed0(s);
  const double ar = -((s.x[0] * ax[0] + s.y[0] * ay[0]) / std::max(r, 1e-30));
  CHECK((v * v / r) == doctest::Approx(ar).epsilon(1e-12));
}

TEST_CASE("pairwise_radial_equilibrium BH-only matches runtime softened BH circular speed") {
  Config c = base_config();
  c.softening = 50.0;
  c.enable_star_star_gravity = false;
  c.galaxy_init_velocity_mode = "pairwise_radial_equilibrium";
  State s;
  galaxy::GalaxyInitAudit audit;
  galaxy::initialize_galaxy_disk(c, s, &audit);
  CHECK(audit.pairwise_radial_equilibrium_used == true);
  CHECK(audit.pairwise_radial_equilibrium_fallback_count == 0);
  std::vector<double> ax, ay;
  galaxy::NewtonianPackage pkg;
  pkg.compute_accelerations(s, c.bh_mass, c.softening, false, ax, ay);
  const double r = radius0(s);
  const double inward = -((s.x[0] * ax[0] + s.y[0] * ay[0]) / std::max(r, 1e-30));
  CHECK((speed0(s) * speed0(s) / r) == doctest::Approx(inward).epsilon(1e-12));
}

TEST_CASE("galaxy init preserves PR109 star-star toggle semantics") {
  Config c = base_config();
  c.n_stars = 5;
  c.softening = 20.0;

  State s_false_a, s_false_b;
  c.enable_star_star_gravity = false;
  c.star_mass = 1.0e20;
  galaxy::initialize_galaxy_disk(c, s_false_a, nullptr);
  c.star_mass = 1.0e35;
  galaxy::initialize_galaxy_disk(c, s_false_b, nullptr);
  for (int i = 0; i < c.n_stars; ++i) CHECK(std::hypot(s_false_a.vx[i], s_false_a.vy[i]) == doctest::Approx(std::hypot(s_false_b.vx[i], s_false_b.vy[i])).epsilon(1e-12));

  State s_true_a, s_true_b;
  c.enable_star_star_gravity = true;
  c.star_mass = 1.0e20;
  galaxy::initialize_galaxy_disk(c, s_true_a, nullptr);
  c.star_mass = 1.0e35;
  galaxy::initialize_galaxy_disk(c, s_true_b, nullptr);

  bool found_diff = false;
  for (int i = 0; i < c.n_stars; ++i) {
    if (std::abs(std::hypot(s_true_a.vx[i], s_true_a.vy[i]) - std::hypot(s_true_b.vx[i], s_true_b.vy[i])) > 1e-6) {
      found_diff = true;
      break;
    }
  }
  CHECK(found_diff);
}

TEST_CASE("pairwise_radial_equilibrium star-star toggle semantics") {
  Config c = base_config();
  c.n_stars = 8;
  c.softening = 20.0;
  c.galaxy_init_velocity_mode = "pairwise_radial_equilibrium";
  State s_false_a, s_false_b;
  c.enable_star_star_gravity = false;
  c.star_mass = 1.0e20;
  galaxy::initialize_galaxy_disk(c, s_false_a, nullptr);
  c.star_mass = 1.0e35;
  galaxy::initialize_galaxy_disk(c, s_false_b, nullptr);
  for (int i = 0; i < c.n_stars; ++i) CHECK(std::hypot(s_false_a.vx[i], s_false_a.vy[i]) == doctest::Approx(std::hypot(s_false_b.vx[i], s_false_b.vy[i])).epsilon(1e-12));
  State s_true_a, s_true_b;
  c.enable_star_star_gravity = true;
  c.star_mass = 1.0e20;
  galaxy::initialize_galaxy_disk(c, s_true_a, nullptr);
  c.star_mass = 1.0e35;
  galaxy::initialize_galaxy_disk(c, s_true_b, nullptr);
  bool found_diff = false;
  for (int i = 0; i < c.n_stars; ++i) if (std::abs(std::hypot(s_true_a.vx[i], s_true_a.vy[i]) - std::hypot(s_true_b.vx[i], s_true_b.vy[i])) > 1e-6) { found_diff = true; break; }
  CHECK(found_diff);
}

TEST_CASE("galaxy init softened speed finite at tiny radius") {
  Config c = base_config();
  c.softening = 100.0;
  c.galaxy_radius = 1e-20;
  State s;
  galaxy::initialize_galaxy_disk(c, s, nullptr);
  CHECK(std::isfinite(speed0(s)));
}

TEST_CASE("pairwise_radial_equilibrium enforces v^2/r against runtime radial acceleration for multi-particle disk") {
  Config c = base_config();
  c.n_stars = 16;
  c.enable_star_star_gravity = true;
  c.softening = 10.0;
  c.galaxy_init_velocity_mode = "pairwise_radial_equilibrium";
  c.velocity_noise = 0.0;
  State s;
  galaxy::GalaxyInitAudit audit;
  galaxy::initialize_galaxy_disk(c, s, &audit);
  std::vector<double> ax, ay;
  galaxy::NewtonianPackage pkg;
  pkg.compute_accelerations(s, c.bh_mass, c.softening, true, ax, ay);
  for (int i = 0; i < s.n(); ++i) {
    const double r = std::hypot(s.x[i], s.y[i]);
    const double r_safe = std::max(r, 1e-30);
    const double inward = -((s.x[i] * ax[i] + s.y[i] * ay[i]) / r_safe);
    const double v_t = (-s.y[i] * s.vx[i] + s.x[i] * s.vy[i]) / r_safe;
    CHECK((v_t * v_t / r_safe) == doctest::Approx(inward).epsilon(1e-7));
  }
  CHECK(audit.pairwise_radial_equilibrium_fallback_count == 0);
}
