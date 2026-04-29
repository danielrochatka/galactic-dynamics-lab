#include "config.hpp"
#include "doctest.h"
#include "physics/Newtonian/newtonian.hpp"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "types.hpp"

#include <cmath>
#include <vector>

namespace {
constexpr double G_SI = 6.6743e-11;

galaxy::Config base_xi_cfg() {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_kernel_coupling = 0.0;
  c.tpf_4d_xi_motion_readout_scale = 1.0;
  c.tpf_global_accel_shunt_enable = false;
  return c;
}

galaxy::State three_body_state(double m1, double m2) {
  galaxy::State s;
  s.resize(3);
  s.x = {2.0, -1.0, 3.0};
  s.y = {0.0, 0.0, 0.0};
  s.vx = {0.0, 0.0, 0.0};
  s.vy = {0.0, 0.0, 0.0};
  s.mass = {1.0, m1, m2};
  return s;
}

}

TEST_CASE("xi off mode with unit K_xi matches softened Newtonian shape/sign") {
  auto c = base_xi_cfg();
  c.tpf_4d_xi_motion_readout_scale = 1.0;
  galaxy::TPFCorePackage p;
  p.init_from_config(c);

  galaxy::State s;
  s.resize(1);
  s.x[0] = 4.0; s.y[0] = 3.0; s.vx[0] = 0.0; s.vy[0] = 0.0; s.mass[0] = 1.0;
  std::vector<double> ax, ay;
  const double bh_mass = 10.0;
  const double eps = 0.5;
  p.compute_accelerations(s, bh_mass, eps, false, ax, ay);

  const double r2 = s.x[0]*s.x[0] + s.y[0]*s.y[0] + eps*eps;
  const double inv_r3 = 1.0 / (r2 * std::sqrt(r2));
  CHECK(ax[0] == doctest::Approx(-1.0 * bh_mass * s.x[0] * inv_r3));
  CHECK(ay[0] == doctest::Approx(-1.0 * bh_mass * s.y[0] * inv_r3));
  CHECK(ax[0] * s.x[0] + ay[0] * s.y[0] < 0.0);
}

TEST_CASE("xi kernel star_star=false ignores stellar source masses") {
  auto c = base_xi_cfg();
  galaxy::TPFCorePackage p;
  p.init_from_config(c);

  std::vector<double> ax_a, ay_a, ax_b, ay_b;
  auto s_a = three_body_state(2.0, 5.0);
  auto s_b = three_body_state(2000.0, 5000.0);
  p.compute_accelerations(s_a, 7.0, 0.1, false, ax_a, ay_a);
  p.compute_accelerations(s_b, 7.0, 0.1, false, ax_b, ay_b);

  CHECK(ax_a[0] == doctest::Approx(ax_b[0]));
  CHECK(ay_a[0] == doctest::Approx(ay_b[0]));
}

TEST_CASE("xi kernel star_star=true responds to stellar source masses") {
  auto c = base_xi_cfg();
  galaxy::TPFCorePackage p;
  p.init_from_config(c);

  std::vector<double> ax_a, ay_a, ax_b, ay_b;
  auto s_a = three_body_state(2.0, 5.0);
  auto s_b = three_body_state(2000.0, 5000.0);
  p.compute_accelerations(s_a, 0.0, 0.1, true, ax_a, ay_a);
  p.compute_accelerations(s_b, 0.0, 0.1, true, ax_b, ay_b);

  CHECK(std::fabs(ax_a[0] - ax_b[0]) > 1e-12);
}

TEST_CASE("xi off with K_xi=G matches Newtonian for BH-only and star_star paths") {
  auto c = base_xi_cfg();
  c.tpf_4d_xi_motion_readout_scale = G_SI;
  galaxy::TPFCorePackage tpf;
  tpf.init_from_config(c);
  galaxy::NewtonianPackage newton;

  galaxy::State s;
  s.resize(3);
  s.x = {3.0, -1.5, 0.25};
  s.y = {-0.5, 2.0, -1.75};
  s.vx = {0.0, 0.0, 0.0};
  s.vy = {0.0, 0.0, 0.0};
  s.mass = {2.0, 5.0, 11.0};

  const double bh = 13.0;
  for (double eps : {0.0, 0.2, 1.3}) {
    std::vector<double> ax_t, ay_t, ax_n, ay_n;
    tpf.compute_accelerations(s, bh, eps, false, ax_t, ay_t);
    newton.compute_accelerations(s, bh, eps, false, ax_n, ay_n);
    for (std::size_t i = 0; i < ax_t.size(); ++i) {
      CHECK(ax_t[i] == doctest::Approx(ax_n[i]).epsilon(1e-12));
      CHECK(ay_t[i] == doctest::Approx(ay_n[i]).epsilon(1e-12));
    }

    tpf.compute_accelerations(s, bh, eps, true, ax_t, ay_t);
    newton.compute_accelerations(s, bh, eps, true, ax_n, ay_n);
    for (std::size_t i = 0; i < ax_t.size(); ++i) {
      CHECK(ax_t[i] == doctest::Approx(ax_n[i]).epsilon(1e-12));
      CHECK(ay_t[i] == doctest::Approx(ay_n[i]).epsilon(1e-12));
    }
  }
}

TEST_CASE("xi kernel uses tpfcore_source_softening override and fallback") {
  galaxy::State s;
  s.resize(1);
  s.x[0] = 1.0; s.y[0] = 0.0; s.vx[0] = 0.0; s.vy[0] = 0.0; s.mass[0] = 1.0;
  const double bh = 2.0;

  auto c_override = base_xi_cfg();
  c_override.tpf_4d_xi_motion_readout_scale = 1.0;
  c_override.tpfcore_source_softening = 0.75;
  galaxy::TPFCorePackage p_override;
  p_override.init_from_config(c_override);

  auto c_fallback = base_xi_cfg();
  c_fallback.tpf_4d_xi_motion_readout_scale = 1.0;
  c_fallback.tpfcore_source_softening = 0.0;
  galaxy::TPFCorePackage p_fallback;
  p_fallback.init_from_config(c_fallback);

  std::vector<double> ax, ay;

  // Override: should ignore call softening and use 0.75.
  p_override.compute_accelerations(s, bh, 0.1, false, ax, ay);
  const double r2_override = 1.0 + 0.75 * 0.75;
  const double expected_override = -bh * (1.0 / (r2_override * std::sqrt(r2_override)));
  CHECK(ax[0] == doctest::Approx(expected_override).epsilon(1e-12));

  // Fallback: source_softening=0 should use the call softening.
  p_fallback.compute_accelerations(s, bh, 0.1, false, ax, ay);
  const double r2_fallback = 1.0 + 0.1 * 0.1;
  const double expected_fallback = -bh * (1.0 / (r2_fallback * std::sqrt(r2_fallback)));
  CHECK(ax[0] == doctest::Approx(expected_fallback).epsilon(1e-12));
}
