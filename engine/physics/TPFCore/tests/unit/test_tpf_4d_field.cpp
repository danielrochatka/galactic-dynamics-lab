#include "doctest.h"
#include "physics/TPFCore/source_ansatz.hpp"
#include "physics/TPFCore/tpf_4d_field.hpp"

#include <cstddef>

using galaxy::tpfcore::invariant_I_4d;
using galaxy::tpfcore::LAMBDA_4D;
using galaxy::tpfcore::metric_signature_sign;
using galaxy::tpfcore::theta_double_contraction_4d;
using galaxy::tpfcore::Theta4;
using galaxy::tpfcore::trace_contraction_4d;
using galaxy::tpfcore::Xi4;

TEST_CASE("4D metric signature signs") {
  CHECK(metric_signature_sign(0) == doctest::Approx(-1.0));
  CHECK(metric_signature_sign(1) == doctest::Approx(1.0));
  CHECK(metric_signature_sign(2) == doctest::Approx(1.0));
  CHECK(metric_signature_sign(3) == doctest::Approx(1.0));
}

TEST_CASE("4D trace contraction uses explicit signature") {
  const Theta4 theta{2.0, 0.0, 0.0, 0.0,
                     5.0, 0.0, 0.0,
                     7.0, 0.0,
                     11.0};
  CHECK(trace_contraction_4d(theta) == doctest::Approx(-2.0 + 5.0 + 7.0 + 11.0));
}

TEST_CASE("4D invariant uses lambda_4d") {
  const Theta4 theta{2.0, 3.0, 5.0, 7.0,
                     11.0, 13.0, 17.0,
                     19.0, 23.0,
                     29.0};
  const double trace = -2.0 + 11.0 + 19.0 + 29.0;
  const double signs[4] = {-1.0, 1.0, 1.0, 1.0};
  double contraction = 0.0;
  for (int mu = 0; mu < 4; ++mu) {
    for (int nu = 0; nu < 4; ++nu) {
      const double v = theta.component(static_cast<std::size_t>(mu), static_cast<std::size_t>(nu));
      contraction += signs[mu] * signs[nu] * v * v;
    }
  }
  CHECK(theta_double_contraction_4d(theta) == doctest::Approx(contraction));
  CHECK(invariant_I_4d(theta) == doctest::Approx(contraction - LAMBDA_4D * trace * trace));
}

TEST_CASE("4D zero/static component sanity") {
  Xi4 xi{0.0, 0.0, 0.0, 0.0};
  CHECK(xi.component(0) == doctest::Approx(0.0));
  CHECK(xi.component(1) == doctest::Approx(0.0));
  CHECK(xi.component(2) == doctest::Approx(0.0));
  CHECK(xi.component(3) == doctest::Approx(0.0));

  Theta4 theta{};
  CHECK(trace_contraction_4d(theta) == doctest::Approx(0.0));
  CHECK(theta_double_contraction_4d(theta) == doctest::Approx(0.0));
  CHECK(invariant_I_4d(theta) == doctest::Approx(0.0));
}

TEST_CASE("4D simple diagonal tensor sanity") {
  const Theta4 theta{1.5, 0.0, 0.0, 0.0,
                     2.0, 0.0, 0.0,
                     3.0, 0.0,
                     4.0};
  const double trace = -1.5 + 2.0 + 3.0 + 4.0;
  const double contraction = 1.5 * 1.5 + 2.0 * 2.0 + 3.0 * 3.0 + 4.0 * 4.0;
  CHECK(trace_contraction_4d(theta) == doctest::Approx(trace));
  CHECK(theta_double_contraction_4d(theta) == doctest::Approx(contraction));
  CHECK(invariant_I_4d(theta) == doctest::Approx(contraction - LAMBDA_4D * trace * trace));
}
