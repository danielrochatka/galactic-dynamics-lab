#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "types.hpp"

#include <cmath>
#include <vector>

namespace {

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

TEST_CASE("xi off mode with K_xi=G matches Newtonian BH direction/sign with shared softening") {
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
