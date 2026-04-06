#include "doctest.h"
#include "physics/TPFCore/field_evaluation.hpp"
#include "physics/TPFCore/source_ansatz.hpp"
#include "types.hpp"

#include <cmath>

using galaxy::State;
using galaxy::tpfcore::compute_invariant_I;
using galaxy::tpfcore::evaluate_provisional_field_multi_source;
using galaxy::tpfcore::evaluate_provisional_field_single_source;
using galaxy::tpfcore::LAMBDA_4D;
using galaxy::tpfcore::provisional_point_source_residual;

TEST_CASE("compute_invariant_I diagonal Theta") {
  galaxy::tpfcore::Theta3D t{};
  t.xx = t.yy = t.zz = 1.0;
  double I = compute_invariant_I(t);
  CHECK(I == doctest::Approx(3.0 - LAMBDA_4D * 9.0));
}

TEST_CASE("provisional_point_source_field is finite away from origin") {
  auto pf = galaxy::tpfcore::provisional_point_source_field(0.0, 0.0, 100.0, 3.0, 4.0, 0.5);
  CHECK(std::isfinite(pf.xi.x));
  CHECK(std::isfinite(pf.xi.y));
  CHECK(std::isfinite(pf.theta.xx));
}

TEST_CASE("residual norm small for moderate softening") {
  auto r = provisional_point_source_residual(0.0, 0.0, 1.0, 2.0, 0.0, 0.3);
  CHECK(std::isfinite(r.norm()));
  CHECK(r.norm() < 1.0);
}

TEST_CASE("single-source vs BH-only multi-source field match") {
  double eps = 0.4;
  double M = 888.0;
  State st;
  st.resize(1);
  st.x[0] = 12.0;
  st.y[0] = -5.0;
  st.mass[0] = 1.0;
  auto single = evaluate_provisional_field_single_source(0.0, 0.0, M, st.x[0], st.y[0], eps);
  auto multi = evaluate_provisional_field_multi_source(st, 0, M, false, eps);
  CHECK(single.invariant_I == doctest::Approx(multi.invariant_I));
  CHECK(single.xi.x == doctest::Approx(multi.xi.x));
  CHECK(single.xi.y == doctest::Approx(multi.xi.y));
}
