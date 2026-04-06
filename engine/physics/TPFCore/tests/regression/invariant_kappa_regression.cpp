/**
 * Regression: manuscript I vs legacy Frobenius scalar; κ-ledger profile smoke.
 * Mirrors legacy test_tpf_regression.cpp expectations.
 */
#include "doctest.h"
#include "physics/TPFCore/derived_tpf_radial.hpp"
#include "physics/TPFCore/source_ansatz.hpp"
#include "types.hpp"

using galaxy::State;
using galaxy::tpfcore::build_tpf_gravity_profile;
using galaxy::tpfcore::compute_invariant_I;
using galaxy::tpfcore::derived_invariant_I_contracted;
using galaxy::tpfcore::LAMBDA_4D;

TEST_CASE("regression: compute_invariant_I vs derived_invariant_I_contracted") {
  galaxy::tpfcore::Theta3D t{};
  t.xx = t.yy = t.zz = 1.0;
  double I_ms = compute_invariant_I(t);
  double I_frob = derived_invariant_I_contracted(t);
  CHECK(I_ms == doctest::Approx(3.0 - LAMBDA_4D * 9.0));
  CHECK(I_frob == doctest::Approx(3.0));
}

TEST_CASE("regression: BH-only profile bin layout") {
  State st;
  st.resize(0);
  double bh_kg = 1.98847e30;
  auto prof = build_tpf_gravity_profile(st, bh_kg, 1e20, 10, 1e12, 1.0);
  CHECK(prof.bins == 10);
  CHECK(prof.M_eff_enc.size() == 10u);
  CHECK(std::isfinite(prof.M_eff_enc.back()));
}
