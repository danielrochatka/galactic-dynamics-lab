#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/field_evaluation.hpp"
#include "physics/TPFCore/source_ansatz.hpp"
#include "physics/TPFCore/tpf_core_package.hpp"
#include <cmath>
#include <vector>

namespace {

galaxy::State sample_state() {
  galaxy::State s;
  s.resize(2);
  s.x[0] = 2.0e9;
  s.y[0] = 1.5e9;
  s.vx[0] = 12.0;
  s.vy[0] = -7.0;
  s.mass[0] = 5.0e20;

  s.x[1] = -1.4e9;
  s.y[1] = 1.1e9;
  s.vx[1] = -9.0;
  s.vy[1] = 3.0;
  s.mass[1] = 8.0e20;
  return s;
}

galaxy::TPFCorePackage make_package(const std::string& mode,
                                    double alpha_si,
                                    double kappa,
                                    double vdsg_coupling,
                                    double source_softening = 1.1e7) {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = mode;
  c.tpf_weak_field_correspondence_alpha_si = alpha_si;
  c.tpf_kappa = kappa;
  c.tpf_vdsg_coupling = vdsg_coupling;
  c.tpfcore_source_softening = source_softening;
  c.tpfcore_enable_provisional_readout = false;
  c.tpfcore_readout_mode = "tensor_radial_projection";
  c.tpfcore_readout_scale = 1.0;
  c.tpfcore_theta_tt_scale = 1.0;
  c.tpfcore_theta_tr_scale = 1.0;
  c.tpf_global_accel_shunt_enable = false;
  c.tpf_cooling_fraction = 0.0;

  galaxy::TPFCorePackage pkg;
  pkg.init_from_config(c);
  return pkg;
}

void run_accel(const galaxy::TPFCorePackage& pkg,
               const galaxy::State& s,
               double bh_mass,
               bool star_star,
               std::vector<double>& ax,
               std::vector<double>& ay) {
  pkg.compute_accelerations(s, bh_mass, 0.0, star_star, ax, ay);
}

}  // namespace

TEST_CASE("direct_tpf uses tensor principal-part formula and matches Theta/I/Cij projection") {
  const galaxy::State s = sample_state();
  const double bh_mass = 2.0e30;
  const double kappa = 3.7e4;
  const double eps = 1.1e7;

  galaxy::TPFCorePackage pkg = make_package("direct_tpf", -6.67430e-11, kappa, 0.0, eps);

  std::vector<double> ax, ay;
  run_accel(pkg, s, bh_mass, true, ax, ay);

  REQUIRE(ax.size() == static_cast<size_t>(s.n()));
  REQUIRE(ay.size() == static_cast<size_t>(s.n()));

  for (int i = 0; i < s.n(); ++i) {
    const galaxy::tpfcore::FieldAtPoint field =
        galaxy::tpfcore::evaluate_provisional_field_multi_source(s, i, bh_mass, true, eps);
    const galaxy::tpfcore::Theta3D& theta = field.theta;
    const double theta_trace = theta.trace();
    const double I = galaxy::tpfcore::compute_invariant_I(theta);

    const double c_xx = kappa * (theta.xx * theta.xx + theta.xy * theta.xy + theta.xz * theta.xz -
                                 galaxy::tpfcore::LAMBDA_4D * theta_trace * theta.xx - 0.5 * I);
    const double c_xy = kappa * (theta.xx * theta.xy + theta.xy * theta.yy + theta.xz * theta.yz -
                                 galaxy::tpfcore::LAMBDA_4D * theta_trace * theta.xy);
    const double c_yy = kappa * (theta.xy * theta.xy + theta.yy * theta.yy + theta.yz * theta.yz -
                                 galaxy::tpfcore::LAMBDA_4D * theta_trace * theta.yy - 0.5 * I);

    const double xi_norm = std::sqrt(field.xi.x * field.xi.x + field.xi.y * field.xi.y);
    REQUIRE(xi_norm > 1e-300);
    const double ux = field.xi.x / xi_norm;
    const double uy = field.xi.y / xi_norm;
    const double expected_ax = -(c_xx * ux + c_xy * uy);
    const double expected_ay = -(c_xy * ux + c_yy * uy);

    CHECK(ax[i] == doctest::Approx(expected_ax));
    CHECK(ay[i] == doctest::Approx(expected_ay));
  }
}

TEST_CASE("direct_tpf Earth-Moon step-0 is attractive along x-axis") {
  galaxy::State s;
  s.resize(2);
  s.x[0] = 0.0;             // Earth
  s.y[0] = 0.0;
  s.vx[0] = 0.0;
  s.vy[0] = 0.0;
  s.mass[0] = 5.97219e24;
  s.x[1] = 3.844e8;         // Moon on +x
  s.y[1] = 0.0;
  s.vx[1] = 0.0;
  s.vy[1] = 0.0;
  s.mass[1] = 7.34767309e22;

  galaxy::TPFCorePackage pkg = make_package("direct_tpf", -6.67430e-11, 2.0e4, 0.0, 1.0e7);
  std::vector<double> ax, ay;
  run_accel(pkg, s, 0.0, true, ax, ay);

  REQUIRE(ax.size() == 2);
  CHECK(ax[0] > 0.0);
  CHECK(ax[1] < 0.0);
}

TEST_CASE("direct_tpf gives nonzero origin response for off-origin source") {
  galaxy::State s;
  s.resize(2);
  s.x[0] = 0.0;
  s.y[0] = 0.0;
  s.vx[0] = 0.0;
  s.vy[0] = 0.0;
  s.mass[0] = 1.0e20;
  s.x[1] = 2.0e9;
  s.y[1] = 0.0;
  s.vx[1] = 0.0;
  s.vy[1] = 0.0;
  s.mass[1] = 8.0e20;

  galaxy::TPFCorePackage pkg = make_package("direct_tpf", -6.67430e-11, 2.0e4, 0.0, 1.0e7);
  std::vector<double> ax, ay;
  run_accel(pkg, s, 0.0, true, ax, ay);

  CHECK(std::abs(ax[0]) > 1e-30);
}

TEST_CASE("direct_tpf symmetric equal masses are equal-and-opposite") {
  galaxy::State s;
  s.resize(2);
  s.x[0] = -1.0e9;
  s.y[0] = 0.0;
  s.vx[0] = 0.0;
  s.vy[0] = 0.0;
  s.mass[0] = 5.0e20;
  s.x[1] = 1.0e9;
  s.y[1] = 0.0;
  s.vx[1] = 0.0;
  s.vy[1] = 0.0;
  s.mass[1] = 5.0e20;

  galaxy::TPFCorePackage pkg = make_package("direct_tpf", -6.67430e-11, 2.0e4, 0.0, 1.0e7);
  std::vector<double> ax, ay;
  run_accel(pkg, s, 0.0, true, ax, ay);

  CHECK(ax[0] > 0.0);
  CHECK(ax[1] < 0.0);
  CHECK(ax[0] == doctest::Approx(-ax[1]).epsilon(1e-12));
  CHECK(ay[0] == doctest::Approx(-ay[1]).epsilon(1e-12));
}

TEST_CASE("direct_tpf translation invariance for two-body relative acceleration") {
  galaxy::State s0;
  s0.resize(2);
  s0.x[0] = -4.0e8;
  s0.y[0] = 1.5e8;
  s0.vx[0] = 0.0;
  s0.vy[0] = 0.0;
  s0.mass[0] = 6.0e20;
  s0.x[1] = 7.0e8;
  s0.y[1] = -2.2e8;
  s0.vx[1] = 0.0;
  s0.vy[1] = 0.0;
  s0.mass[1] = 4.0e20;

  galaxy::State s1 = s0;
  const double tx = 2.0e9;
  const double ty = -3.0e9;
  for (int i = 0; i < s1.n(); ++i) {
    s1.x[i] += tx;
    s1.y[i] += ty;
  }

  galaxy::TPFCorePackage pkg = make_package("direct_tpf", -6.67430e-11, 2.0e4, 0.0, 1.0e7);
  std::vector<double> ax0, ay0, ax1, ay1;
  run_accel(pkg, s0, 0.0, true, ax0, ay0);
  run_accel(pkg, s1, 0.0, true, ax1, ay1);

  const double rel0x = ax0[1] - ax0[0];
  const double rel0y = ay0[1] - ay0[0];
  const double rel1x = ax1[1] - ax1[0];
  const double rel1y = ay1[1] - ay1[0];
  CHECK(rel1x == doctest::Approx(rel0x).epsilon(1e-12));
  CHECK(rel1y == doctest::Approx(rel0y).epsilon(1e-12));
}

TEST_CASE("direct_tpf is independent of correspondence alpha while v11_weak_field_truncation depends on alpha") {
  const galaxy::State s = sample_state();
  const double bh_mass = 2.0e30;

  std::vector<double> ax_d1, ay_d1, ax_d2, ay_d2;
  galaxy::TPFCorePackage direct_a = make_package("direct_tpf", -6.67430e-11, 2.0e4, 0.0);
  galaxy::TPFCorePackage direct_b = make_package("direct_tpf", -1.33486e-10, 2.0e4, 0.0);
  run_accel(direct_a, s, bh_mass, true, ax_d1, ay_d1);
  run_accel(direct_b, s, bh_mass, true, ax_d2, ay_d2);

  for (int i = 0; i < s.n(); ++i) {
    CHECK(ax_d1[i] == doctest::Approx(ax_d2[i]));
    CHECK(ay_d1[i] == doctest::Approx(ay_d2[i]));
  }

  std::vector<double> ax_v1, ay_v1, ax_v2, ay_v2;
  galaxy::TPFCorePackage v11_a = make_package("v11_weak_field_truncation", -6.67430e-11, 2.0e4, 0.0);
  galaxy::TPFCorePackage v11_b = make_package("v11_weak_field_truncation", -1.33486e-10, 2.0e4, 0.0);
  run_accel(v11_a, s, bh_mass, true, ax_v1, ay_v1);
  run_accel(v11_b, s, bh_mass, true, ax_v2, ay_v2);

  bool any_changed = false;
  for (int i = 0; i < s.n(); ++i) {
    if (std::abs(ax_v1[i] - ax_v2[i]) > 1e-30 || std::abs(ay_v1[i] - ay_v2[i]) > 1e-30) {
      any_changed = true;
      break;
    }
  }
  CHECK(any_changed);
}

TEST_CASE("direct_tpf VDSG contribution is additive and converges continuously to zero with coupling") {
  const galaxy::State s = sample_state();
  const double bh_mass = 2.0e30;
  const double kappa = 2.0e4;

  std::vector<double> ax_base, ay_base, ax_small, ay_small, ax_large, ay_large;
  galaxy::TPFCorePackage p_base = make_package("direct_tpf", -6.67430e-11, kappa, 0.0);
  galaxy::TPFCorePackage p_small = make_package("direct_tpf", -6.67430e-11, kappa, 1.0e-20);
  galaxy::TPFCorePackage p_large = make_package("direct_tpf", -6.67430e-11, kappa, 2.0e-20);

  run_accel(p_base, s, bh_mass, true, ax_base, ay_base);
  run_accel(p_small, s, bh_mass, true, ax_small, ay_small);
  run_accel(p_large, s, bh_mass, true, ax_large, ay_large);

  for (int i = 0; i < s.n(); ++i) {
    const double dax_small = ax_small[i] - ax_base[i];
    const double day_small = ay_small[i] - ay_base[i];
    const double dax_large = ax_large[i] - ax_base[i];
    const double day_large = ay_large[i] - ay_base[i];

    CHECK(std::isfinite(ax_base[i]));
    CHECK(std::isfinite(ay_base[i]));
    CHECK(std::abs(dax_small) <= std::abs(dax_large) + 1e-30);
    CHECK(std::abs(day_small) <= std::abs(day_large) + 1e-30);

    CHECK(dax_large == doctest::Approx(2.0 * dax_small).epsilon(1e-6));
    CHECK(day_large == doctest::Approx(2.0 * day_small).epsilon(1e-6));
  }
}

TEST_CASE("direct_tpf responds linearly to config tpf_kappa mapping") {
  const galaxy::State s = sample_state();
  const double bh_mass = 2.0e30;

  std::vector<double> ax_a, ay_a, ax_b, ay_b;
  galaxy::TPFCorePackage p_a = make_package("direct_tpf", -6.67430e-11, 1.3e4, 0.0);
  galaxy::TPFCorePackage p_b = make_package("direct_tpf", -6.67430e-11, 2.6e4, 0.0);
  run_accel(p_a, s, bh_mass, true, ax_a, ay_a);
  run_accel(p_b, s, bh_mass, true, ax_b, ay_b);

  for (int i = 0; i < s.n(); ++i) {
    CHECK(ax_b[i] == doctest::Approx(2.0 * ax_a[i]).epsilon(1e-12));
    CHECK(ay_b[i] == doctest::Approx(2.0 * ay_a[i]).epsilon(1e-12));
  }
}
