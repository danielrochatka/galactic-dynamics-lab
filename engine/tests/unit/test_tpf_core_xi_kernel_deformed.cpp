#include "doctest.h"

#include "config.hpp"
#include "physics/TPFCore/tpf_core_package.hpp"

#include <cmath>
#include <vector>

namespace {

galaxy::State tiny_state() {
  galaxy::State s;
  s.resize(2);
  s.x[0] = 2.0; s.y[0] = 0.0; s.vx[0] = 1.0e7; s.vy[0] = 0.0; s.mass[0] = 5.0;
  s.x[1] = -2.0; s.y[1] = 0.0; s.vx[1] = -0.5e7; s.vy[1] = 0.2e7; s.mass[1] = 7.0;
  return s;
}

void run_xi_mode(const galaxy::Config& cfg,
                 const galaxy::State& s,
                 double bh_mass,
                 bool star_star,
                 std::vector<double>& ax,
                 std::vector<double>& ay) {
  galaxy::TPFCorePackage pkg;
  pkg.init_from_config(cfg);
  pkg.compute_accelerations(s, bh_mass, 0.0, star_star, ax, ay);
}

void check_bh_inward_cardinal_signs(const galaxy::Config& c) {
  constexpr double bh_mass = 10.0;
  constexpr double r = 2.0;
  const struct {
    double x;
    double y;
    double ax_sign;
    double ay_sign;
  } cases[] = {
      {+r, 0.0, -1.0, 0.0},
      {-r, 0.0, +1.0, 0.0},
      {0.0, +r, 0.0, -1.0},
      {0.0, -r, 0.0, +1.0},
  };

  for (const auto& tc : cases) {
    galaxy::State s;
    s.resize(1);
    s.x[0] = tc.x;
    s.y[0] = tc.y;
    s.vx[0] = 0.0;
    s.vy[0] = 0.0;
    s.mass[0] = 1.0;

    std::vector<double> ax, ay;
    run_xi_mode(c, s, bh_mass, false, ax, ay);
    REQUIRE(ax.size() == 1);
    REQUIRE(ay.size() == 1);
    if (tc.ax_sign < 0.0) CHECK(ax[0] < 0.0);
    if (tc.ax_sign > 0.0) CHECK(ax[0] > 0.0);
    if (tc.ax_sign == 0.0) CHECK(ax[0] == doctest::Approx(0.0));
    if (tc.ay_sign < 0.0) CHECK(ay[0] < 0.0);
    if (tc.ay_sign > 0.0) CHECK(ay[0] > 0.0);
    if (tc.ay_sign == 0.0) CHECK(ay[0] == doctest::Approx(0.0));
  }
}

}  // namespace

TEST_CASE("xi_kernel_deformed runtime path executes for tiny galaxy state") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_motion_readout_scale = 2.0e-12;

  std::vector<double> ax, ay;
  run_xi_mode(c, tiny_state(), 10.0, true, ax, ay);
  CHECK(ax.size() == 2);
  CHECK(ay.size() == 2);
  CHECK(std::isfinite(ax[0]));
  CHECK(std::isfinite(ay[0]));
}

TEST_CASE("xi_kernel_deformed off mode matches baseline BH Xi acceleration shape") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_motion_readout_scale = 3.0e-12;

  galaxy::State s;
  s.resize(1);
  s.x[0] = 2.0; s.y[0] = 0.0; s.vx[0] = 0.0; s.vy[0] = 0.0; s.mass[0] = 1.0;

  std::vector<double> ax, ay;
  run_xi_mode(c, s, 10.0, false, ax, ay);

  const double xi_x_base = 10.0 * 2.0 / std::pow(2.0, 3.0);
  CHECK(ax[0] == doctest::Approx(-c.tpf_4d_xi_motion_readout_scale * xi_x_base));
  CHECK(ay[0] == doctest::Approx(0.0));
}

TEST_CASE("xi_kernel_deformed BH-only inward cardinal directions hold in off mode") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_motion_readout_scale = 6.67430e-11;
  check_bh_inward_cardinal_signs(c);
}

TEST_CASE("xi_kernel_deformed scalar_beta uses Xi_eff and differs from Xi_base with nonzero relative speed") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "scalar_beta";
  c.tpf_4d_xi_kernel_coupling = 1.0;
  c.tpf_4d_xi_kernel_beta_power = 1.0;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_motion_readout_scale = 1.0e-12;

  galaxy::State s;
  s.resize(1);
  s.x[0] = 2.0; s.y[0] = 0.0; s.vx[0] = 1.0e7; s.vy[0] = 0.0; s.mass[0] = 1.0;

  std::vector<double> ax, ay;
  run_xi_mode(c, s, 10.0, false, ax, ay);

  const double xi_base = 10.0 * 2.0 / std::pow(2.0, 3.0);
  const double beta = std::fabs(s.vx[0]) / 299792458.0;
  const double factor = c.tpf_4d_xi_kernel_coupling * beta;
  const double xi_eff = xi_base * (1.0 + factor);
  CHECK(factor > 0.0);
  CHECK(xi_eff != doctest::Approx(xi_base));
  CHECK(ax[0] == doctest::Approx(-c.tpf_4d_xi_motion_readout_scale * xi_eff));
  CHECK(ay[0] == doctest::Approx(0.0));
}

TEST_CASE("xi_kernel_deformed metric_velocity applies non-identity metric scale and Xi_eff acceleration") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "metric_velocity";
  c.tpf_4d_xi_kernel_coupling = 1.0;
  c.tpf_4d_xi_kernel_beta_power = 1.0;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_motion_readout_scale = 1.0e-12;

  galaxy::State s;
  s.resize(1);
  s.x[0] = 2.0; s.y[0] = 1.0; s.vx[0] = 1.0e7; s.vy[0] = 0.0; s.mass[0] = 1.0;

  std::vector<double> ax_off, ay_off, ax_metric, ay_metric;

  galaxy::Config c_off = c;
  c_off.tpf_4d_xi_kernel_mode = "off";
  run_xi_mode(c_off, s, 10.0, false, ax_off, ay_off);
  run_xi_mode(c, s, 10.0, false, ax_metric, ay_metric);

  const double beta = std::fabs(s.vx[0]) / 299792458.0;
  const double metric_scale = std::max(c.tpf_4d_xi_kernel_metric_min,
                                       std::min(c.tpf_4d_xi_kernel_metric_max, 1.0 + c.tpf_4d_xi_kernel_coupling * beta));
  CHECK(metric_scale != doctest::Approx(1.0));
  CHECK(std::fabs(ax_metric[0] - ax_off[0]) > 1.0e-18);
  CHECK(std::isfinite(ax_metric[0]));
  CHECK(std::isfinite(ay_metric[0]));
}

TEST_CASE("xi_kernel_deformed BH-only inward cardinal directions hold in metric_velocity mode") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "metric_velocity";
  c.tpf_4d_xi_kernel_coupling = 10.0;
  c.tpf_4d_xi_kernel_beta_power = 1.0;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_motion_readout_scale = 6.67430e-11;
  check_bh_inward_cardinal_signs(c);
}

TEST_CASE("xi_kernel_deformed metric_transverse_wake uses transverse wake strength with radial metric axis") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "metric_transverse_wake";
  c.tpf_4d_xi_kernel_coupling = 10.0;
  c.tpf_4d_xi_kernel_beta_power = 1.0;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_motion_readout_scale = 1.0e-12;

  galaxy::State s;
  s.resize(1);
  s.x[0] = 2.0;
  s.y[0] = 0.0;
  s.vx[0] = 0.0;
  s.vy[0] = 1.0e7;
  s.mass[0] = 1.0;

  std::vector<double> ax_off, ay_off, ax_wake, ay_wake;
  galaxy::Config c_off = c;
  c_off.tpf_4d_xi_kernel_mode = "off";
  run_xi_mode(c_off, s, 10.0, false, ax_off, ay_off);
  run_xi_mode(c, s, 10.0, false, ax_wake, ay_wake);

  CHECK(std::fabs(ax_wake[0] - ax_off[0]) > 1.0e-18);
  CHECK(ay_wake[0] == doctest::Approx(0.0));
}

TEST_CASE("xi_kernel_deformed BH-only inward cardinal directions hold in metric_transverse_wake mode") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "metric_transverse_wake";
  c.tpf_4d_xi_kernel_coupling = 10.0;
  c.tpf_4d_xi_kernel_beta_power = 1.0;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_motion_readout_scale = 6.67430e-11;
  check_bh_inward_cardinal_signs(c);
}

TEST_CASE("xi_kernel_deformed route does not use additive VDSG or direct/provisional readout paths") {
  galaxy::Config c0;
  c0.tpf_dynamics_mode = "xi_kernel_deformed";
  c0.tpf_4d_xi_kernel_mode = "off";
  c0.tpf_4d_xi_motion_readout_scale = 1.0e-12;
  c0.tpf_vdsg_coupling = 0.0;
  c0.tpfcore_enable_provisional_readout = false;
  c0.tpfcore_readout_mode = "tensor_radial_projection";

  galaxy::Config c1 = c0;
  c1.tpf_vdsg_coupling = 1.0e12;
  c1.tpfcore_enable_provisional_readout = true;
  c1.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c1.tpfcore_readout_scale = 1234.0;

  std::vector<double> ax0, ay0, ax1, ay1;
  const galaxy::State s = tiny_state();
  run_xi_mode(c0, s, 10.0, true, ax0, ay0);
  run_xi_mode(c1, s, 10.0, true, ax1, ay1);

  REQUIRE(ax0.size() == ax1.size());
  for (std::size_t i = 0; i < ax0.size(); ++i) {
    CHECK(ax1[i] == doctest::Approx(ax0[i]));
    CHECK(ay1[i] == doctest::Approx(ay0[i]));
  }
}
