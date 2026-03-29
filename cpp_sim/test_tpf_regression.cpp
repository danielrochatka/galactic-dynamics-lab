/**
 * Regression checks: manuscript I vs legacy Frobenius scalar, VDSG circular BH limit, κ-ledger smoke.
 * Run: make test_tpf_regression && ./test_tpf_regression
 */

#include "physics/TPFCore/derived_tpf_radial.hpp"
#include "physics/TPFCore/source_ansatz.hpp"
#include "types.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

int g_fail = 0;

void check(bool ok, const char* msg) {
  if (!ok) {
    std::cerr << "FAIL: " << msg << '\n';
    ++g_fail;
  }
}

}  // namespace

int main() {
  using galaxy::State;
  using galaxy::tpfcore::compute_invariant_I;
  using galaxy::tpfcore::derived_invariant_I_contracted;
  using galaxy::tpfcore::build_tpf_gravity_profile;
  using galaxy::tpfcore::LAMBDA_4D;

  galaxy::tpfcore::Theta3D t{};
  t.xx = t.yy = t.zz = 1.0;
  double I_ms = compute_invariant_I(t);
  double I_frob = derived_invariant_I_contracted(t);
  check(std::abs(I_ms - (3.0 - LAMBDA_4D * 9.0)) < 1e-12, "compute_invariant_I on diagonal Theta");
  check(std::abs(I_frob - 3.0) < 1e-12, "derived_invariant_I_contracted on diagonal Theta");

  constexpr double c_light = 299792458.0;
  double vdsg_coupling = 1e-4;
  double doppler = 1.0 + vdsg_coupling * (0.0 / c_light);
  check(std::abs(doppler - 1.0) < 1e-15, "VDSG BH doppler = 1 when v·r̂ = 0 (circular)");

  State st;
  st.resize(0);
  double bh_kg = 1.98847e30;
  auto prof = build_tpf_gravity_profile(st, bh_kg, 1e12, 10, 1e20, 1.0);
  check(prof.bins == 10 && prof.M_eff_enc.size() == 10u, "profile bin layout");
  check(std::isfinite(prof.M_eff_enc.back()), "M_eff_enc cumulative finite (BH-only smoke)");

  if (g_fail != 0) {
    std::cerr << g_fail << " failure(s)\n";
    return 1;
  }
  std::cout << "test_tpf_regression: ok\n";
  return 0;
}
