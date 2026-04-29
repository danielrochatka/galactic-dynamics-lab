#include "doctest.h"
#include "physics/Newtonian/newtonian.hpp"
#include "types.hpp"

#include <cmath>

namespace {
constexpr double G_SI = 6.6743e-11;
}

TEST_CASE("Newtonian BH acceleration: radial inward on +x axis") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(1);
  st.x[0] = 10.0;
  st.y[0] = 0.0;
  st.vx[0] = st.vy[0] = 0.0;
  st.mass[0] = 0.1;
  double bh = 1000.0;
  double eps = 1.0;
  std::vector<double> ax, ay;
  n.compute_accelerations(st, bh, eps, false, ax, ay);
  REQUIRE(ax.size() == 1);
  CHECK(ax[0] < 0.0);
  CHECK(ay[0] == doctest::Approx(0.0));
}

TEST_CASE("Newtonian BH circular-orbit SI sanity: a = v^2/r = G*M/r^2") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(1);
  const double r = 1.0e10;
  const double bh_mass = 5.0e35;
  const double v_circ = std::sqrt(G_SI * bh_mass / r);
  st.x[0] = r;
  st.y[0] = 0.0;
  st.vx[0] = 0.0;
  st.vy[0] = v_circ;
  st.mass[0] = 1.0;
  std::vector<double> ax, ay;
  n.compute_accelerations(st, bh_mass, 0.0, false, ax, ay);
  const double a_mag = std::sqrt(ax[0] * ax[0] + ay[0] * ay[0]);
  const double expected = (v_circ * v_circ) / r;
  CHECK(a_mag == doctest::Approx(expected).epsilon(1e-12));
  CHECK(ax[0] < 0.0);
  CHECK(std::abs(ay[0]) < 1e-20);
}

TEST_CASE("Newtonian force magnitude uses SI G for star-star interaction") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(2);
  const double m0 = 3.0e30;
  const double m1 = 8.0e29;
  const double d = 2.0e9;
  st.x[0] = -0.5 * d; st.y[0] = 0.0; st.mass[0] = m0;
  st.x[1] =  0.5 * d; st.y[1] = 0.0; st.mass[1] = m1;
  st.vx[0] = st.vy[0] = st.vx[1] = st.vy[1] = 0.0;
  std::vector<double> ax, ay;
  n.compute_accelerations(st, 0.0, 0.0, true, ax, ay);
  const double expected_a0 = G_SI * m1 / (d * d);
  const double expected_a1 = G_SI * m0 / (d * d);
  CHECK(ax[0] == doctest::Approx(expected_a0).epsilon(1e-12));
  CHECK(ax[1] == doctest::Approx(-expected_a1).epsilon(1e-12));
  CHECK(ay[0] == doctest::Approx(0.0));
  CHECK(ay[1] == doctest::Approx(0.0));
}

TEST_CASE("Newtonian potential energy uses SI G for BH and pairwise terms") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(2);
  const double bh_mass = 1.0e36;
  const double m0 = 2.0e30;
  const double m1 = 3.0e30;
  st.x[0] = 4.0e9; st.y[0] = 0.0; st.mass[0] = m0;
  st.x[1] = 0.0;   st.y[1] = 3.0e9; st.mass[1] = m1;
  st.vx[0] = st.vy[0] = st.vx[1] = st.vy[1] = 0.0;
  const double pe = n.compute_potential_energy(st, bh_mass, 0.0, true);
  const double r0 = 4.0e9;
  const double r1 = 3.0e9;
  const double r01 = 5.0e9;
  const double expected = -G_SI * ((bh_mass * m0) / r0 + (bh_mass * m1) / r1 + (m0 * m1) / r01);
  CHECK(pe == doctest::Approx(expected).epsilon(1e-12));
}

TEST_CASE("Newtonian symmetry: opposite points have opposite ax") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(2);
  st.x[0] = 5.0;
  st.y[0] = 3.0;
  st.x[1] = -5.0;
  st.y[1] = -3.0;
  st.mass[0] = st.mass[1] = 0.05;
  std::vector<double> ax, ay;
  n.compute_accelerations(st, 100.0, 1.0, false, ax, ay);
  REQUIRE(ax.size() == 2);
  CHECK(ax[0] == doctest::Approx(-ax[1]));
  CHECK(ay[0] == doctest::Approx(-ay[1]));
}

TEST_CASE("Newtonian softened BH acceleration matches analytic formula exactly") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(1);
  st.x[0] = 3.0;
  st.y[0] = -4.0;
  st.vx[0] = st.vy[0] = 0.0;
  st.mass[0] = 1.0;
  const double bh_mass = 9.0e10;
  const double eps = 0.7;

  std::vector<double> ax, ay;
  n.compute_accelerations(st, bh_mass, eps, false, ax, ay);

  const double r2 = st.x[0] * st.x[0] + st.y[0] * st.y[0] + eps * eps;
  const double inv_r3 = 1.0 / (r2 * std::sqrt(r2));
  CHECK(ax[0] == doctest::Approx(-G_SI * bh_mass * st.x[0] * inv_r3).epsilon(1e-14));
  CHECK(ay[0] == doctest::Approx(-G_SI * bh_mass * st.y[0] * inv_r3).epsilon(1e-14));
}

TEST_CASE("Newtonian potential gradient matches acceleration (BH-only finite-difference)") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(1);
  st.x[0] = 2.3;
  st.y[0] = -1.7;
  st.vx[0] = st.vy[0] = 0.0;
  st.mass[0] = 5.0;
  const double bh_mass = 2.0e6;
  const double eps = 0.25;
  const double h = 1.0e-6;

  std::vector<double> ax, ay;
  n.compute_accelerations(st, bh_mass, eps, false, ax, ay);

  auto pe_at = [&](double x, double y) {
    galaxy::State t = st;
    t.x[0] = x;
    t.y[0] = y;
    return n.compute_potential_energy(t, bh_mass, eps, false);
  };
  const double dUdx = (pe_at(st.x[0] + h, st.y[0]) - pe_at(st.x[0] - h, st.y[0])) / (2.0 * h);
  const double dUdy = (pe_at(st.x[0], st.y[0] + h) - pe_at(st.x[0], st.y[0] - h)) / (2.0 * h);

  CHECK(ax[0] == doctest::Approx(-dUdx / st.mass[0]).epsilon(1e-7));
  CHECK(ay[0] == doctest::Approx(-dUdy / st.mass[0]).epsilon(1e-7));
}

TEST_CASE("Newtonian potential gradient matches acceleration (star-star finite-difference)") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(2);
  st.x[0] = -0.8; st.y[0] = 0.4; st.mass[0] = 2.0;
  st.x[1] = 1.2; st.y[1] = -0.1; st.mass[1] = 7.0;
  st.vx[0] = st.vy[0] = st.vx[1] = st.vy[1] = 0.0;
  const double eps = 0.3;
  const double h = 1.0e-6;

  std::vector<double> ax, ay;
  n.compute_accelerations(st, 0.0, eps, true, ax, ay);

  auto pe_at = [&](double x0, double y0) {
    galaxy::State t = st;
    t.x[0] = x0;
    t.y[0] = y0;
    return n.compute_potential_energy(t, 0.0, eps, true);
  };
  const double dUdx0 = (pe_at(st.x[0] + h, st.y[0]) - pe_at(st.x[0] - h, st.y[0])) / (2.0 * h);
  const double dUdy0 = (pe_at(st.x[0], st.y[0] + h) - pe_at(st.x[0], st.y[0] - h)) / (2.0 * h);

  CHECK(ax[0] == doctest::Approx(-dUdx0 / st.mass[0]).epsilon(1e-7));
  CHECK(ay[0] == doctest::Approx(-dUdy0 / st.mass[0]).epsilon(1e-7));
}

TEST_CASE("Newtonian star-star pair force conserves total momentum (unequal masses, eps>0)") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(2);
  st.x[0] = -1.3; st.y[0] = 0.2; st.mass[0] = 3.0;
  st.x[1] = 2.1; st.y[1] = -0.4; st.mass[1] = 11.0;
  st.vx[0] = st.vy[0] = st.vx[1] = st.vy[1] = 0.0;

  std::vector<double> ax, ay;
  n.compute_accelerations(st, 0.0, 0.6, true, ax, ay);
  const double px_dot = st.mass[0] * ax[0] + st.mass[1] * ax[1];
  const double py_dot = st.mass[0] * ay[0] + st.mass[1] * ay[1];
  CHECK(std::fabs(px_dot) < 1e-20);
  CHECK(std::fabs(py_dot) < 1e-20);
}

TEST_CASE("Newtonian softening limits are finite and monotone") {
  galaxy::NewtonianPackage n;
  galaxy::State st;
  st.resize(1);
  st.x[0] = 1.5;
  st.y[0] = 0.0;
  st.mass[0] = 1.0;
  st.vx[0] = st.vy[0] = 0.0;
  const double bh = 5.0e8;

  std::vector<double> ax0, ay0, ax1, ay1, ax2, ay2;
  n.compute_accelerations(st, bh, 0.0, false, ax0, ay0);
  n.compute_accelerations(st, bh, 0.5, false, ax1, ay1);
  n.compute_accelerations(st, bh, 5.0, false, ax2, ay2);

  const double mag0 = std::hypot(ax0[0], ay0[0]);
  const double mag1 = std::hypot(ax1[0], ay1[0]);
  const double mag2 = std::hypot(ax2[0], ay2[0]);
  CHECK(mag0 > mag1);
  CHECK(mag1 > mag2);
  CHECK(std::isfinite(mag2));

  // eps=0 case should recover inverse-square exactly for this axis-aligned setup.
  const double expected = G_SI * bh / (st.x[0] * st.x[0]);
  CHECK(std::fabs(ax0[0] + expected) < 1e-12 * expected);

  // r=0 with eps>0 should be finite and exactly zero vector by symmetry.
  st.x[0] = 0.0;
  st.y[0] = 0.0;
  n.compute_accelerations(st, bh, 2.0, false, ax2, ay2);
  CHECK(std::isfinite(ax2[0]));
  CHECK(std::isfinite(ay2[0]));
  CHECK(ax2[0] == doctest::Approx(0.0));
  CHECK(ay2[0] == doctest::Approx(0.0));
}

TEST_CASE("Newtonian BH-only softening factor audit finalizes and reports distances") {
  galaxy::Config c;
  c.softening_audit_enable = true;
  galaxy::NewtonianPackage n;
  n.init_from_config(c);
  galaxy::State st;
  st.resize(2);
  st.x = {3.0, -4.0};
  st.y = {4.0, 0.0};
  st.mass = {1.0, 2.0};
  st.vx = {0.0, 0.0};
  st.vy = {0.0, 0.0};
  const double eps = 0.5;
  std::vector<double> ax, ay;
  n.compute_accelerations(st, 10.0, eps, false, ax, ay);
  n.compute_accelerations(st, 10.0, eps, false, ax, ay);
  const auto& a = n.softening_audit_stats();
  CHECK(a.run_total_pair_count >= 4);
  CHECK(a.run_total_bh_pair_count > 0);
  CHECK(a.call_count == 2);
  CHECK(a.run_total_violation_count == 0);
  CHECK(a.run_min_r > 0.0);
  CHECK(a.last_call_median_r > 0.0);
  CHECK(a.eps_used == doctest::Approx(eps));
}

TEST_CASE("Newtonian star-star softening factor audit counts star pairs with no violations") {
  galaxy::Config c;
  c.softening_audit_enable = true;
  galaxy::NewtonianPackage n;
  n.init_from_config(c);
  galaxy::State st;
  st.resize(2);
  st.x = {-1.0, 2.0};
  st.y = {0.2, -0.3};
  st.mass = {3.0, 11.0};
  st.vx = {0.0, 0.0};
  st.vy = {0.0, 0.0};
  std::vector<double> ax, ay;
  n.compute_accelerations(st, 0.0, 0.6, true, ax, ay);
  n.compute_accelerations(st, 0.0, 0.6, true, ax, ay);
  const auto& a = n.softening_audit_stats();
  CHECK(a.run_total_star_pair_count >= 4);
  CHECK(a.call_count == 2);
  CHECK(a.run_total_violation_count == 0);
  CHECK(a.run_total_force_vector_pair_hardening_count == 0);
  CHECK(a.run_total_force_vector_pair_direction_flip_count == 0);
  CHECK(a.run_total_force_vector_pair_nan_inf_count == 0);
  CHECK(a.run_total_force_vector_pair_inward_bh_violation_count == 0);
  CHECK(std::isfinite(a.net_ratio_median));
  CHECK(std::isfinite(a.net_ratio_p95));
  CHECK(std::isfinite(a.net_ratio_max));
}

TEST_CASE("Newtonian BH exact sign/magnitude softened vs unsoftened") {
  galaxy::NewtonianPackage n;
  galaxy::State st; st.resize(1);
  const double r = 5.0, bh = 2.0e9;
  st.x[0] = r; st.y[0] = 0.0; st.mass[0] = 1.0; st.vx[0] = st.vy[0] = 0.0;
  std::vector<double> ax0, ay0, ax1, ay1;
  n.compute_accelerations(st, bh, 0.0, false, ax0, ay0);
  n.compute_accelerations(st, bh, 1.2, false, ax1, ay1);
  CHECK(ax0[0] == doctest::Approx(-G_SI * bh / (r * r)).epsilon(1e-12));
  CHECK(ax1[0] == doctest::Approx(-G_SI * bh * r / std::pow(r * r + 1.44, 1.5)).epsilon(1e-12));
  CHECK(ax1[0] < 0.0); CHECK(std::fabs(ay1[0]) < 1e-20);
  CHECK(std::fabs(ax1[0]) <= std::fabs(ax0[0]));
}
