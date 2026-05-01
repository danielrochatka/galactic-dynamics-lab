#include "doctest.h"
#include "physics/TPFCore/tpf_4d_field.hpp"
#include "physics/TPFCore/tpf_4d_static_field.hpp"

#include <cmath>
#include <vector>

using galaxy::tpfcore::evaluate_static_source_field_4d;
using galaxy::tpfcore::evaluate_static_sources_field_4d;
using galaxy::tpfcore::Field4DAtPoint;
using galaxy::tpfcore::invariant_I_4d;
using galaxy::tpfcore::StaticProbePoint4D;
using galaxy::tpfcore::StaticSourcePoint4D;
using galaxy::tpfcore::trace_contraction_4d;

TEST_CASE("4D static monopole xi components on x-axis") {
  const StaticSourcePoint4D source{5.0, 0.0, 0.0, 0.0};
  const StaticProbePoint4D probe{2.0, 0.0, 0.0};
  const Field4DAtPoint out = evaluate_static_source_field_4d(source, probe, 0.0);

  const double expected_x = source.mass / (2.0 * 2.0);
  CHECK(out.xi.x == doctest::Approx(expected_x));
  CHECK(out.xi.x > 0.0);
  CHECK(out.xi.y == doctest::Approx(0.0));
  CHECK(out.xi.z == doctest::Approx(0.0));
  CHECK(out.xi.t == doctest::Approx(0.0));

  CHECK(out.theta.tt == doctest::Approx(0.0));
  CHECK(out.theta.tx == doctest::Approx(0.0));
  CHECK(out.theta.ty == doctest::Approx(0.0));
  CHECK(out.theta.tz == doctest::Approx(0.0));
  CHECK(out.theta.xt == doctest::Approx(0.0));
  CHECK(out.theta.yt == doctest::Approx(0.0));
  CHECK(out.theta.zt == doctest::Approx(0.0));
}

TEST_CASE("4D static monopole gradient on x-axis") {
  const double m = 6.0;
  const double r = 3.0;
  const Field4DAtPoint out = evaluate_static_source_field_4d(
      StaticSourcePoint4D{m, 0.0, 0.0, 0.0}, StaticProbePoint4D{r, 0.0, 0.0}, 0.0);

  const double mr3 = m / (r * r * r);
  CHECK(out.theta.xx == doctest::Approx(-2.0 * mr3));
  CHECK(out.theta.yy == doctest::Approx(mr3));
  CHECK(out.theta.zz == doctest::Approx(mr3));

  CHECK(out.theta.xy == doctest::Approx(0.0));
  CHECK(out.theta.xz == doctest::Approx(0.0));
  CHECK(out.theta.yx == doctest::Approx(0.0));
  CHECK(out.theta.yz == doctest::Approx(0.0));
  CHECK(out.theta.zx == doctest::Approx(0.0));
  CHECK(out.theta.zy == doctest::Approx(0.0));

  CHECK(out.theta_trace_4d == doctest::Approx(0.0));
}

TEST_CASE("4D static evaluator writes ordered spatial tensor components explicitly") {
  const Field4DAtPoint out = evaluate_static_source_field_4d(
      StaticSourcePoint4D{2.0, 0.5, -0.25, 1.0}, StaticProbePoint4D{1.75, 2.5, -0.75}, 0.0);

  CHECK(out.theta.xy == doctest::Approx(out.theta.yx));
  CHECK(out.theta.xz == doctest::Approx(out.theta.zx));
  CHECK(out.theta.yz == doctest::Approx(out.theta.zy));
}

TEST_CASE("4D static evaluator rotational consistency on principal axes") {
  const double m = 4.0;
  const double r = 2.0;

  const Field4DAtPoint x_axis = evaluate_static_source_field_4d(
      StaticSourcePoint4D{m, 0.0, 0.0, 0.0}, StaticProbePoint4D{r, 0.0, 0.0}, 0.0);
  const Field4DAtPoint y_axis = evaluate_static_source_field_4d(
      StaticSourcePoint4D{m, 0.0, 0.0, 0.0}, StaticProbePoint4D{0.0, r, 0.0}, 0.0);
  const Field4DAtPoint z_axis = evaluate_static_source_field_4d(
      StaticSourcePoint4D{m, 0.0, 0.0, 0.0}, StaticProbePoint4D{0.0, 0.0, r}, 0.0);

  const double expected_mag = m / (r * r);
  CHECK(x_axis.xi.x == doctest::Approx(expected_mag));
  CHECK(y_axis.xi.y == doctest::Approx(expected_mag));
  CHECK(z_axis.xi.z == doctest::Approx(expected_mag));

  CHECK(std::sqrt(x_axis.xi.x * x_axis.xi.x + x_axis.xi.y * x_axis.xi.y + x_axis.xi.z * x_axis.xi.z) == doctest::Approx(expected_mag));
  CHECK(std::sqrt(y_axis.xi.x * y_axis.xi.x + y_axis.xi.y * y_axis.xi.y + y_axis.xi.z * y_axis.xi.z) == doctest::Approx(expected_mag));
  CHECK(std::sqrt(z_axis.xi.x * z_axis.xi.x + z_axis.xi.y * z_axis.xi.y + z_axis.xi.z * z_axis.xi.z) == doctest::Approx(expected_mag));
}

TEST_CASE("4D static multi-source superposition recomputes trace and invariant") {
  const StaticProbePoint4D probe{0.25, -0.5, 0.75};
  const double eps = 0.1;
  const StaticSourcePoint4D s1{3.0, -1.0, 0.0, 0.5};
  const StaticSourcePoint4D s2{2.0, 0.5, 1.0, -0.25};

  const Field4DAtPoint one = evaluate_static_source_field_4d(s1, probe, eps);
  const Field4DAtPoint two = evaluate_static_source_field_4d(s2, probe, eps);
  const Field4DAtPoint sum = evaluate_static_sources_field_4d(std::vector<StaticSourcePoint4D>{s1, s2}, probe, eps);

  CHECK(sum.xi.t == doctest::Approx(one.xi.t + two.xi.t));
  CHECK(sum.xi.x == doctest::Approx(one.xi.x + two.xi.x));
  CHECK(sum.xi.y == doctest::Approx(one.xi.y + two.xi.y));
  CHECK(sum.xi.z == doctest::Approx(one.xi.z + two.xi.z));

  CHECK(sum.theta.tt == doctest::Approx(one.theta.tt + two.theta.tt));
  CHECK(sum.theta.tx == doctest::Approx(one.theta.tx + two.theta.tx));
  CHECK(sum.theta.ty == doctest::Approx(one.theta.ty + two.theta.ty));
  CHECK(sum.theta.tz == doctest::Approx(one.theta.tz + two.theta.tz));
  CHECK(sum.theta.xt == doctest::Approx(one.theta.xt + two.theta.xt));
  CHECK(sum.theta.xx == doctest::Approx(one.theta.xx + two.theta.xx));
  CHECK(sum.theta.xy == doctest::Approx(one.theta.xy + two.theta.xy));
  CHECK(sum.theta.xz == doctest::Approx(one.theta.xz + two.theta.xz));
  CHECK(sum.theta.yt == doctest::Approx(one.theta.yt + two.theta.yt));
  CHECK(sum.theta.yx == doctest::Approx(one.theta.yx + two.theta.yx));
  CHECK(sum.theta.yy == doctest::Approx(one.theta.yy + two.theta.yy));
  CHECK(sum.theta.yz == doctest::Approx(one.theta.yz + two.theta.yz));
  CHECK(sum.theta.zt == doctest::Approx(one.theta.zt + two.theta.zt));
  CHECK(sum.theta.zx == doctest::Approx(one.theta.zx + two.theta.zx));
  CHECK(sum.theta.zy == doctest::Approx(one.theta.zy + two.theta.zy));
  CHECK(sum.theta.zz == doctest::Approx(one.theta.zz + two.theta.zz));

  CHECK(sum.theta_trace_4d == doctest::Approx(trace_contraction_4d(sum.theta)));
  CHECK(sum.invariant_I_4d == doctest::Approx(invariant_I_4d(sum.theta)));
}

TEST_CASE("4D static evaluator softening keeps center finite") {
  const Field4DAtPoint center = evaluate_static_source_field_4d(
      StaticSourcePoint4D{5.0, 1.0, -2.0, 3.0}, StaticProbePoint4D{1.0, -2.0, 3.0}, 0.5);

  CHECK(std::isfinite(center.xi.x));
  CHECK(std::isfinite(center.xi.y));
  CHECK(std::isfinite(center.xi.z));

  CHECK(std::isfinite(center.theta.xx));
  CHECK(std::isfinite(center.theta.yy));
  CHECK(std::isfinite(center.theta.zz));
  CHECK(std::isfinite(center.theta_trace_4d));
  CHECK(std::isfinite(center.invariant_I_4d));
}
