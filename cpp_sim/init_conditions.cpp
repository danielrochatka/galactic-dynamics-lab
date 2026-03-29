#include "init_conditions.hpp"
#include "galaxy_init.hpp"
#include <cmath>
#include <random>
#include <vector>

namespace galaxy {

namespace {
constexpr double PI = 3.14159265358979323846;
}

void init_galaxy_disk(const Config& config, State& state) {
  initialize_galaxy_disk(config, state, nullptr);
}

void init_two_body(const Config& config, State& state) {
  (void)config;
  /* Hardcoded Earth–Moon benchmark (SI). State has no fixed-body flags; both bodies evolve.
   * Use bh_mass = 0 and enable_star_star_gravity = true so only pairwise masses apply.
   * Note: NewtonianPackage uses a G=1 convention; these are literal SI storage values. */
  constexpr double kEarthMass = 5.972e24;
  constexpr double kMoonMass = 7.348e22;
  constexpr double kMoonX = 3.844e8;
  constexpr double kMoonVy = 1022.0;

  state.resize(0);
  state.resize(2);

  state.mass[0] = kEarthMass;
  state.x[0] = 0.0;
  state.y[0] = 0.0;
  state.vx[0] = 0.0;
  state.vy[0] = 0.0;

  state.mass[1] = kMoonMass;
  state.x[1] = kMoonX;
  state.y[1] = 0.0;
  state.vx[1] = 0.0;
  state.vy[1] = kMoonVy;
}

void init_two_body_star_around_bh(const Config& config, State& state) {
  state.resize(1);
  double r0 = config.validation_two_body_radius;
  double v_circ = std::sqrt(config.bh_mass / r0);
  double v0 = config.initial_velocity_scale * config.validation_two_body_speed_ratio * v_circ;

  state.x[0] = r0;
  state.y[0] = 0.0;
  state.vx[0] = 0.0;
  state.vy[0] = v0;
  state.mass[0] = config.star_mass;
}

void init_symmetric_pair(const Config& config, State& state) {
  state.resize(2);
  double a = config.validation_symmetric_separation;
  double v = config.initial_velocity_scale * config.validation_symmetric_speed;
  double m = config.star_mass;

  state.x[0] = -a; state.y[0] = 0.0;
  state.vx[0] = 0.0; state.vy[0] = v;
  state.mass[0] = m;

  state.x[1] = a; state.y[1] = 0.0;
  state.vx[1] = 0.0; state.vy[1] = -v;
  state.mass[1] = m;
}

void init_small_n(const Config& config, State& state) {
  int n = config.validation_small_n;
  if (n < 3) n = 3;
  if (n > 10) n = 10;
  state.resize(n);

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> u(-2, 2);

  const double bh_mass = config.bh_mass;
  const double star_mass = config.star_mass;

  std::vector<double> r(n), theta(n);
  for (int i = 0; i < n; ++i) {
    theta[i] = (2 * PI * i) / n + 0.1;
    r[i] = 20.0 + u(rng);
  }

  const double v_scale = config.initial_velocity_scale;
  for (int i = 0; i < n; ++i) {
    state.x[i] = r[i] * std::cos(theta[i]);
    state.y[i] = r[i] * std::sin(theta[i]);
    double v_circ = std::sqrt(bh_mass / r[i]);
    state.vx[i] = v_scale * (-v_circ * std::sin(theta[i]));
    state.vy[i] = v_scale * (v_circ * std::cos(theta[i]));
    state.mass[i] = star_mass;
  }
}

}  // namespace galaxy
