#include "init_conditions.hpp"
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>

namespace galaxy {

namespace {

const double PI = 3.14159265358979323846;

}  // namespace

void init_galaxy_disk(const Config& config, State& state, unsigned seed) {
  const int n = config.n_stars;
  state.resize(n);

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u01(0, 1);
  std::uniform_real_distribution<double> u02pi(0, 2 * PI);
  std::normal_distribution<double> normal(0, 1);

  const double R = config.galaxy_radius;
  const double r_min = 0.05 * R;
  const double bh_mass = config.bh_mass;
  const double star_mass = config.star_mass;
  const double noise = config.velocity_noise;

  std::vector<double> radii(n), theta(n);
  for (int i = 0; i < n; ++i) {
    double u = u01(rng);
    /* Uniform in area: r^2 uniform between r_min^2 and R^2 */
    double r_sq = r_min * r_min + u * (R * R - r_min * r_min);
    radii[i] = std::sqrt(r_sq);
    theta[i] = u02pi(rng);
  }

  // Sort by radius to compute n_inside (enclosed stellar count)
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&radii](int i, int j) { return radii[i] < radii[j]; });

  std::vector<int> n_inside(n);
  for (int k = 0; k < n; ++k)
    n_inside[order[k]] = k;

  for (int i = 0; i < n; ++i) {
    double r = radii[i];
    double th = theta[i];
    state.x[i] = r * std::cos(th);
    state.y[i] = r * std::sin(th);

    double enclosed_mass = bh_mass + n_inside[i] * star_mass;
    double a_newtonian = (6.6743e-11 * enclosed_mass) / (r * r);
    const double a0 = 1.2e-10;
    double a_tpf = (a_newtonian < a0) ? std::sqrt(a_newtonian * a0) : a_newtonian;
    double v_tpf = std::sqrt(a_tpf * r);

    double vx = -std::sin(th) * v_tpf;
    double vy = std::cos(th) * v_tpf;
    if (noise > 0) {
      double scale = noise * v_tpf;
      vx += scale * normal(rng);
      vy += scale * normal(rng);
    }
    state.vx[i] = config.initial_velocity_scale * vx;
    state.vy[i] = config.initial_velocity_scale * vy;
    state.mass[i] = star_mass;
  }
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
