#include "doctest.h"

#include "../../../../types.hpp"
#include "../../extensions_vdsg.hpp"
#include "../../tpf_core_package.hpp"
#include "../../../../config.hpp"

namespace {

galaxy::State one_target(double x, double y, double vx, double vy, double m) {
  galaxy::State s;
  s.resize(1);
  s.x[0] = x; s.y[0] = y; s.vx[0] = vx; s.vy[0] = vy; s.mass[0] = m;
  return s;
}

void compute_delta(const galaxy::State& s, const std::string& mode, double coupling,
                   double m0, bool weak_gate, double weak_a0, double soft,
                   std::vector<double>& ax, std::vector<double>& ay) {
  galaxy::tpfcore::accumulate_vdsg_velocity_modifier(
      s, 1.98847e30, soft, false, coupling, 0.0, mode,
      m0, 1.0, 0.25, weak_gate, weak_a0, 1.0, 0.25, ax, ay, nullptr);
}

}

TEST_CASE("invalid tpf_vdsg_mode fails loudly in TPFCore init") {
  galaxy::Config c;
  c.tpf_vdsg_mode = "bad_mode";
  galaxy::TPFCorePackage pkg;
  CHECK_THROWS_WITH(pkg.init_from_config(c), doctest::Contains("invalid tpf_vdsg_mode"));
}

TEST_CASE("radial mode tangential motion yields zero vdsg correction") {
  auto s = one_target(1.0e7, 0.0, 0.0, 1.0e3, 1.0e30);
  std::vector<double> ax, ay;
  compute_delta(s, "radial_doppler_rational", 1.0e3, 1.0e20, false, 1.0e-10, 1.0e6, ax, ay);
  CHECK(std::abs(ax[0]) < 1e-30);
  CHECK(std::abs(ay[0]) < 1e-30);
}

TEST_CASE("radial receding strengthens inward pull; approaching weakens") {
  std::vector<double> ax_r, ay_r, ax_a, ay_a;
  compute_delta(one_target(1.0e7, 0.0, 1.0e3, 0.0, 1.0e30), "radial_doppler_exp", 5.0e3, 1.0e20, false, 1.0e-10, 0.0, ax_r, ay_r);
  compute_delta(one_target(1.0e7, 0.0, -1.0e3, 0.0, 1.0e30), "radial_doppler_exp", 5.0e3, 1.0e20, false, 1.0e-10, 0.0, ax_a, ay_a);
  CHECK(ax_r[0] < 0.0);
  CHECK(ax_a[0] > 0.0);
}

TEST_CASE("radial reduced-mass gate suppresses low-mass target") {
  std::vector<double> ax_big, ay_big, ax_small, ay_small;
  compute_delta(one_target(1.0e7, 0.0, 1.0e3, 0.0, 1.0e30), "radial_doppler_bounded", 1.0e4, 1.0e24, false, 1.0e-10, 0.0, ax_big, ay_big);
  compute_delta(one_target(1.0e7, 0.0, 1.0e3, 0.0, 1.0e3), "radial_doppler_bounded", 1.0e4, 1.0e24, false, 1.0e-10, 0.0, ax_small, ay_small);
  CHECK(std::abs(ax_small[0]) < std::abs(ax_big[0]) * 1e-6);
}

TEST_CASE("weak-field gate suppresses compact high-acceleration pairs") {
  auto s = one_target(1.0, 0.0, 1.0e3, 0.0, 1.0e30);
  std::vector<double> ax_off, ay_off, ax_on, ay_on;
  compute_delta(s, "radial_doppler_rational", 1.0e4, 1.0e20, false, 1.0e-10, 0.0, ax_off, ay_off);
  compute_delta(s, "radial_doppler_rational", 1.0e4, 1.0e20, true, 1.0e-10, 0.0, ax_on, ay_on);
  CHECK(std::abs(ax_on[0]) < std::abs(ax_off[0]) * 1e-6);
}
