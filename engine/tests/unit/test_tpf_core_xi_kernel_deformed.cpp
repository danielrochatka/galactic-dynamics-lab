#include "doctest.h"

#include "config.hpp"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "physics/TPFCore/runtime_package_helpers.hpp"

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

TEST_CASE("xi_kernel_deformed metric_transverse_wake positive coupling maps to compressive metric_scale") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "metric_transverse_wake";
  c.tpf_4d_xi_kernel_coupling = 10.0;
  c.tpf_4d_xi_kernel_beta_power = 1.0;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";

  const double beta_pass = 1.0e7 / 299792458.0;
  constexpr double wake_threshold = 0.10;
  constexpr double wake_width = 0.05;
  const double radial_ratio = 0.20;
  const double wake_gate = 0.5 * (1.0 + std::tanh((radial_ratio - wake_threshold) / wake_width));
  const double beta_effective = wake_gate * beta_pass;
  const double factor_raw = c.tpf_4d_xi_kernel_coupling * beta_effective;
  const double metric_scale = std::max(c.tpf_4d_xi_kernel_metric_min,
                                       std::min(c.tpf_4d_xi_kernel_metric_max, 1.0 / (1.0 + factor_raw)));
  CHECK(beta_effective > 0.0);
  CHECK(metric_scale < 1.0);
}

TEST_CASE("xi_kernel_deformed metric_transverse_wake positive coupling strengthens Xi_eff and keeps inward acceleration") {
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
  s.vx[0] = 2.0e6;
  s.vy[0] = 1.0e7;
  s.mass[0] = 1.0;

  std::vector<double> ax_off, ay_off, ax_wake, ay_wake;
  galaxy::Config c_off = c;
  c_off.tpf_4d_xi_kernel_mode = "off";
  run_xi_mode(c_off, s, 10.0, false, ax_off, ay_off);
  run_xi_mode(c, s, 10.0, false, ax_wake, ay_wake);

  const double xi_base_mag = std::fabs(ax_off[0]) / c.tpf_4d_xi_motion_readout_scale;
  const double xi_wake_mag = std::fabs(ax_wake[0]) / c.tpf_4d_xi_motion_readout_scale;
  CHECK(xi_wake_mag > xi_base_mag);
  CHECK(ax_wake[0] < 0.0);
  CHECK(ay_wake[0] == doctest::Approx(0.0));
}

TEST_CASE("xi_kernel_deformed metric_transverse_wake stronger positive coupling increases Xi_eff until clamp") {
  galaxy::Config c_lo;
  c_lo.tpf_dynamics_mode = "xi_kernel_deformed";
  c_lo.tpf_4d_xi_kernel_mode = "metric_transverse_wake";
  c_lo.tpf_4d_xi_kernel_coupling = 10.0;
  c_lo.tpf_4d_xi_kernel_beta_power = 1.0;
  c_lo.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c_lo.tpf_4d_xi_motion_readout_scale = 1.0e-12;

  galaxy::Config c_hi = c_lo;
  c_hi.tpf_4d_xi_kernel_coupling = 100.0;

  galaxy::Config c_clamp = c_lo;
  c_clamp.tpf_4d_xi_kernel_coupling = 1000.0;

  galaxy::Config c_more_clamp = c_lo;
  c_more_clamp.tpf_4d_xi_kernel_coupling = 5000.0;

  galaxy::State s;
  s.resize(1);
  s.x[0] = 2.0;
  s.y[0] = 0.0;
  s.vx[0] = 2.0e6;
  s.vy[0] = 1.0e7;
  s.mass[0] = 1.0;

  std::vector<double> ax_lo, ay_lo, ax_hi, ay_hi, ax_clamp, ay_clamp, ax_more_clamp, ay_more_clamp;
  run_xi_mode(c_lo, s, 10.0, false, ax_lo, ay_lo);
  run_xi_mode(c_hi, s, 10.0, false, ax_hi, ay_hi);
  run_xi_mode(c_clamp, s, 10.0, false, ax_clamp, ay_clamp);
  run_xi_mode(c_more_clamp, s, 10.0, false, ax_more_clamp, ay_more_clamp);

  const double xi_lo = std::fabs(ax_lo[0]) / c_lo.tpf_4d_xi_motion_readout_scale;
  const double xi_hi = std::fabs(ax_hi[0]) / c_hi.tpf_4d_xi_motion_readout_scale;
  const double xi_clamp = std::fabs(ax_clamp[0]) / c_clamp.tpf_4d_xi_motion_readout_scale;
  const double xi_more_clamp = std::fabs(ax_more_clamp[0]) / c_more_clamp.tpf_4d_xi_motion_readout_scale;

  CHECK(xi_hi > xi_lo);
  CHECK(xi_clamp > xi_hi);
  CHECK(xi_more_clamp == doctest::Approx(xi_clamp));
}

TEST_CASE("xi_kernel_deformed metric_transverse_wake negative coupling weakens Xi_eff") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "metric_transverse_wake";
  c.tpf_4d_xi_kernel_coupling = -0.5;
  c.tpf_4d_xi_kernel_beta_power = 1.0;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_motion_readout_scale = 1.0e-12;

  galaxy::State s;
  s.resize(1);
  s.x[0] = 2.0;
  s.y[0] = 0.0;
  s.vx[0] = 2.0e6;
  s.vy[0] = 1.0e7;
  s.mass[0] = 1.0;

  std::vector<double> ax_off, ay_off, ax_wake, ay_wake;
  galaxy::Config c_off = c;
  c_off.tpf_4d_xi_kernel_mode = "off";
  run_xi_mode(c_off, s, 10.0, false, ax_off, ay_off);
  run_xi_mode(c, s, 10.0, false, ax_wake, ay_wake);

  const double xi_base_mag = std::fabs(ax_off[0]) / c.tpf_4d_xi_motion_readout_scale;
  const double xi_wake_mag = std::fabs(ax_wake[0]) / c.tpf_4d_xi_motion_readout_scale;
  CHECK(xi_wake_mag < xi_base_mag);
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

TEST_CASE("xi_kernel_deformed off with K_xi=G matches Newtonian BH law for eps=0 and eps>0") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_motion_readout_scale = galaxy::tpfcore::TPF_G_SI;

  galaxy::State s;
  s.resize(1);
  s.x[0] = 3.0; s.y[0] = 4.0; s.vx[0] = 0.0; s.vy[0] = 0.0; s.mass[0] = 1.0;
  const double bh_mass = 9.0;
  const double r2 = 25.0;

  std::vector<double> ax0, ay0;
  galaxy::TPFCorePackage p0;
  p0.init_from_config(c);
  p0.compute_accelerations(s, bh_mass, 0.0, false, ax0, ay0);
  CHECK(ax0[0] == doctest::Approx(-galaxy::tpfcore::TPF_G_SI * bh_mass * s.x[0] / std::pow(r2, 1.5)));
  CHECK(ay0[0] == doctest::Approx(-galaxy::tpfcore::TPF_G_SI * bh_mass * s.y[0] / std::pow(r2, 1.5)));

  std::vector<double> axs, ays;
  p0.compute_accelerations(s, bh_mass, 2.0, false, axs, ays);
  const double soft_denom = std::pow(r2 + 4.0, 1.5);
  CHECK(axs[0] == doctest::Approx(-galaxy::tpfcore::TPF_G_SI * bh_mass * s.x[0] / soft_denom));
  CHECK(ays[0] == doctest::Approx(-galaxy::tpfcore::TPF_G_SI * bh_mass * s.y[0] / soft_denom));
}

TEST_CASE("xi_kernel_deformed metric_velocity uses r_eff2 basis for softening scale") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "metric_velocity";
  c.tpf_4d_xi_kernel_coupling = 4.0;
  c.tpf_4d_xi_kernel_beta_power = 1.0;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_motion_readout_scale = 1.0;
  c.tpf_4d_xi_source_speed_x = 0.0;
  c.tpf_4d_xi_source_speed_y = 0.0;
  c.tpf_4d_xi_source_speed_z = 0.0;

  galaxy::State s;
  s.resize(1);
  s.x[0] = 2.0; s.y[0] = 1.0; s.vx[0] = 1.0e7; s.vy[0] = 2.0e6; s.mass[0] = 1.0;
  const double bh_mass = 10.0;
  const double eps = 0.7;

  std::vector<double> ax, ay;
  galaxy::TPFCorePackage p;
  p.init_from_config(c);
  p.compute_accelerations(s, bh_mass, eps, false, ax, ay);

  const double dx = s.x[0], dy = s.y[0];
  const double vrel = std::sqrt(s.vx[0] * s.vx[0] + s.vy[0] * s.vy[0]);
  const galaxy::tpfcore::XiWakeKinematics wake =
      galaxy::tpfcore::compute_xi_wake_kinematics(dx, dy, 0.0, s.vx[0], s.vy[0], 0.0, 299792458.0, false);
  const double beta = wake.beta_effective;
  const double metric_scale = std::max(c.tpf_4d_xi_kernel_metric_min,
                                       std::min(c.tpf_4d_xi_kernel_metric_max, 1.0 + c.tpf_4d_xi_kernel_coupling * beta));
  const double nx = s.vx[0] / vrel;
  const double ny = s.vy[0] / vrel;
  const double alpha = metric_scale - 1.0;
  const double nd = nx * dx + ny * dy;
  const double gx = dx + alpha * nx * nd;
  const double gy = dy + alpha * ny * nd;
  const double r_eff2 = dx * gx + dy * gy;
  const double denom = std::pow(r_eff2 + eps * eps, 1.5);
  CHECK(ax[0] == doctest::Approx(-bh_mass * gx / denom));
  CHECK(ay[0] == doctest::Approx(-bh_mass * gy / denom));
}

TEST_CASE("xi_kernel_deformed source_softening override wins and star_star toggle behaves") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_motion_readout_scale = 1.0;
  c.tpfcore_source_softening = 2.0;

  galaxy::State s;
  s.resize(2);
  s.x[0] = 1.0; s.y[0] = 0.0; s.vx[0] = 0.0; s.vy[0] = 0.0; s.mass[0] = 1.0;
  s.x[1] = 2.0; s.y[1] = 0.0; s.vx[1] = 0.0; s.vy[1] = 0.0; s.mass[1] = 3.0;

  galaxy::TPFCorePackage p;
  p.init_from_config(c);
  std::vector<double> ax_false, ay_false, ax_true, ay_true;
  p.compute_accelerations(s, /*bh_mass=*/5.0, /*softening=*/0.1, /*star_star=*/false, ax_false, ay_false);
  p.compute_accelerations(s, /*bh_mass=*/5.0, /*softening=*/0.1, /*star_star=*/true, ax_true, ay_true);

  // BH term with source softening override eps=2.0 (not global softening=0.1).
  CHECK(ax_false[0] == doctest::Approx(-5.0 / std::pow(1.0 + 4.0, 1.5)));
  CHECK(ay_false[0] == doctest::Approx(0.0));
  // Enabling star-star adds inward pull from star[1] at dx=-1, same eps override.
  const double star_term = 3.0 / std::pow(1.0 + 4.0, 1.5);
  CHECK(ax_true[0] == doctest::Approx(ax_false[0] + star_term));
}

TEST_CASE("xi_kernel_deformed refactor preserves BH-only and star_star=true accelerations") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "metric_transverse_wake";
  c.tpf_4d_xi_kernel_coupling = 7.5;
  c.tpf_4d_xi_kernel_beta_power = 1.25;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_motion_readout_scale = 2.5e-12;

  galaxy::State s;
  s.resize(3);
  s.x[0] = 2.0;  s.y[0] = 1.0;  s.vx[0] = 1.2e6;  s.vy[0] = -2.5e6; s.mass[0] = 2.0;
  s.x[1] = -1.0; s.y[1] = 2.5;  s.vx[1] = -0.8e6; s.vy[1] = 0.9e6;  s.mass[1] = 3.5;
  s.x[2] = 0.5;  s.y[2] = -1.5; s.vx[2] = 0.7e6;  s.vy[2] = 1.1e6;  s.mass[2] = 1.5;

  auto compute_reference = [&](bool star_star, std::vector<double>& ax_out, std::vector<double>& ay_out) {
    const double c_light = 299792458.0;
    ax_out.assign(3, 0.0);
    ay_out.assign(3, 0.0);
    for (int i = 0; i < 3; ++i) {
      auto accumulate = [&](double sx, double sy, double svx, double svy, double mass) {
        const double dx = s.x[i] - sx;
        const double dy = s.y[i] - sy;
        const double vx_rel = s.vx[i] - svx;
        const double vy_rel = s.vy[i] - svy;
        const galaxy::tpfcore::XiWakeKinematics wake =
            galaxy::tpfcore::compute_xi_wake_kinematics(dx, dy, 0.0, vx_rel, vy_rel, 0.0, c_light, true);
        const double factor = c.tpf_4d_xi_kernel_coupling * std::pow(wake.beta_effective, c.tpf_4d_xi_kernel_beta_power);
        const double metric_scale =
            std::max(c.tpf_4d_xi_kernel_metric_min, std::min(c.tpf_4d_xi_kernel_metric_max, 1.0 / (1.0 + factor)));
        const double r = std::sqrt(dx * dx + dy * dy);
        const double nx = dx / r;
        const double ny = dy / r;
        const double alpha = metric_scale - 1.0;
        const double nd = nx * dx + ny * dy;
        const double gx = dx + alpha * nx * nd;
        const double gy = dy + alpha * ny * nd;
        const double r_eff2 = dx * gx + dy * gy;
        const double inv_r_eff3 = 1.0 / (r_eff2 * std::sqrt(r_eff2));
        ax_out[i] += -c.tpf_4d_xi_motion_readout_scale * mass * gx * inv_r_eff3;
        ay_out[i] += -c.tpf_4d_xi_motion_readout_scale * mass * gy * inv_r_eff3;
      };
      accumulate(0.0, 0.0, 0.0, 0.0, 6.0);
      if (star_star) {
        for (int j = 0; j < 3; ++j) {
          if (j == i) continue;
          accumulate(s.x[j], s.y[j], s.vx[j], s.vy[j], s.mass[j]);
        }
      }
    }
  };

  galaxy::TPFCorePackage p;
  p.init_from_config(c);

  std::vector<double> ax_bh, ay_bh, ax_ss, ay_ss;
  p.compute_accelerations(s, /*bh_mass=*/6.0, /*softening=*/0.0, /*star_star=*/false, ax_bh, ay_bh);
  p.compute_accelerations(s, /*bh_mass=*/6.0, /*softening=*/0.0, /*star_star=*/true, ax_ss, ay_ss);

  std::vector<double> ax_bh_ref, ay_bh_ref, ax_ss_ref, ay_ss_ref;
  compute_reference(false, ax_bh_ref, ay_bh_ref);
  compute_reference(true, ax_ss_ref, ay_ss_ref);

  for (int i = 0; i < 3; ++i) {
    CHECK(ax_bh[i] == doctest::Approx(ax_bh_ref[i]));
    CHECK(ay_bh[i] == doctest::Approx(ay_bh_ref[i]));
    CHECK(ax_ss[i] == doctest::Approx(ax_ss_ref[i]));
    CHECK(ay_ss[i] == doctest::Approx(ay_ss_ref[i]));
  }
}

TEST_CASE("xi_kernel_deformed pair counters remain accurate without diagnostics hot-loop increments") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_xi_kernel_dump_field_diagnostics = false;

  galaxy::State s;
  s.resize(3);
  s.x[0] = 1.0; s.y[0] = 0.0; s.vx[0] = 0.0; s.vy[0] = 0.0; s.mass[0] = 1.0;
  s.x[1] = 0.0; s.y[1] = 1.0; s.vx[1] = 0.0; s.vy[1] = 0.0; s.mass[1] = 2.0;
  s.x[2] = -1.0; s.y[2] = 0.0; s.vx[2] = 0.0; s.vy[2] = 0.0; s.mass[2] = 3.0;

  galaxy::TPFCorePackage p;
  p.init_from_config(c);
  std::vector<double> ax, ay;

  p.compute_accelerations(s, /*bh_mass=*/5.0, /*softening=*/0.0, /*star_star=*/false, ax, ay);
  auto c1 = p.xi_runtime_counters();
  CHECK(c1.xi_last_call_pair_evaluations == 3);
  CHECK(c1.xi_total_pair_evaluations == 3);

  p.compute_accelerations(s, /*bh_mass=*/5.0, /*softening=*/0.0, /*star_star=*/true, ax, ay);
  auto c2 = p.xi_runtime_counters();
  CHECK(c2.xi_last_call_pair_evaluations == 9);
  CHECK(c2.xi_total_pair_evaluations == 12);

  p.compute_accelerations(s, /*bh_mass=*/0.0, /*softening=*/0.0, /*star_star=*/true, ax, ay);
  auto c3 = p.xi_runtime_counters();
  CHECK(c3.xi_last_call_pair_evaluations == 6);
  CHECK(c3.xi_total_pair_evaluations == 18);
}
