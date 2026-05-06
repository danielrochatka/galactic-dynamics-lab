#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "types.hpp"

#include <cmath>
#include <vector>

namespace {
galaxy::State one_body(double x, double y, double vx, double vy, double m) {
  galaxy::State s;
  s.resize(1);
  s.x[0] = x; s.y[0] = y;
  s.vx[0] = vx; s.vy[0] = vy;
  s.mass[0] = m;
  return s;
}
}

TEST_CASE("tpf_xi_theta_v1 acceleration is -K_xi * Xi_total from canonical sample") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_motion_readout_scale = 2.5e-12;
  galaxy::TPFCorePackage p;
  p.init_from_config(c);
  auto s = one_body(5.0, 0.0, 0.0, 0.0, 1.0);
  std::vector<double> ax, ay;
  p.compute_accelerations(s, 10.0, 0.1, false, ax, ay);
  const auto sample = p.last_xi_theta_v1_sample();
  CHECK(ax[0] == doctest::Approx(-c.tpf_4d_xi_motion_readout_scale * sample.xi_x));
  CHECK(ay[0] == doctest::Approx(-c.tpf_4d_xi_motion_readout_scale * sample.xi_y));
}

TEST_CASE("tpf_xi_theta_v1 cancellation can yield Xi_total ~ 0 with nonzero unsymmetrized Theta") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";
  c.tpf_4d_xi_kernel_mode = "off";
  c.enable_star_star_gravity = true;
  galaxy::TPFCorePackage p;
  p.init_from_config(c);
  galaxy::State s;
  s.resize(3);
  // probe
  s.x[0] = 0.0; s.y[0] = 0.0; s.vx[0] = s.vy[0] = 0.0; s.mass[0] = 1.0;
  // symmetric sources about origin
  s.x[1] = 1.0; s.y[1] = 0.0; s.vx[1] = s.vy[1] = 0.0; s.mass[1] = 2.0;
  s.x[2] = -1.0; s.y[2] = 0.0; s.vx[2] = s.vy[2] = 0.0; s.mass[2] = 2.0;
  std::vector<double> ax, ay;
  p.compute_accelerations(s, 0.0, 1e-6, true, ax, ay);
  const auto sample = p.last_xi_theta_v1_sample();
  CHECK(std::abs(sample.xi_x) < 1e-10);
  CHECK(std::abs(sample.xi_y) < 1e-10);
  const double theta_norm =
      std::abs(sample.theta[0][0]) + std::abs(sample.theta[0][1]) + std::abs(sample.theta[1][0]) +
      std::abs(sample.theta[1][1]) + std::abs(sample.theta[2][2]);
  CHECK(theta_norm > 0.0);
}

