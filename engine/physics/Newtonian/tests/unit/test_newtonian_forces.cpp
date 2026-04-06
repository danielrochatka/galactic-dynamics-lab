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
