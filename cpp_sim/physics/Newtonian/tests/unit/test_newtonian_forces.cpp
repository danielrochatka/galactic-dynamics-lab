#include "doctest.h"
#include "physics/Newtonian/newtonian.hpp"
#include "types.hpp"

#include <cmath>

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
