#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"

#include <cmath>
#include <vector>

namespace {

galaxy::TPFCorePackage make_pkg(const std::string& source, double kappa, double vdsg) {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_direct_field_source = source;
  c.tpf_kappa = kappa;
  c.tpf_vdsg_coupling = vdsg;
  c.tpf_cooling_fraction = 0.0;
  c.tpf_global_accel_shunt_enable = false;
  c.tpfcore_enable_provisional_readout = false;
  c.tpf_xi_constraint_grid_n = 65;
  c.tpf_xi_constraint_half_extent = 10.0;
  c.tpf_xi_constraint_inner_radius = 0.5;
  c.tpf_xi_constraint_max_iters = 80;
  c.tpf_xi_constraint_tol = 1e-6;

  galaxy::TPFCorePackage pkg;
  pkg.init_from_config(c);
  return pkg;
}

galaxy::State one_body(double x, double y, double vx = 0.0, double vy = 0.0) {
  galaxy::State s;
  s.resize(1);
  s.x[0] = x;
  s.y[0] = y;
  s.vx[0] = vx;
  s.vy[0] = vy;
  s.mass[0] = 1.0e20;
  return s;
}

}  // namespace

TEST_CASE("direct_tpf default field source matches explicit provisional_ansatz") {
  galaxy::Config c_default;
  c_default.physics_package = "TPFCore";
  c_default.tpf_dynamics_mode = "direct_tpf";
  c_default.tpf_kappa = 3.0e4;
  c_default.tpf_vdsg_coupling = 0.0;
  c_default.tpf_cooling_fraction = 0.0;
  c_default.tpf_global_accel_shunt_enable = false;
  c_default.tpfcore_enable_provisional_readout = false;

  galaxy::Config c_explicit = c_default;
  c_explicit.tpf_direct_field_source = "provisional_ansatz";

  galaxy::TPFCorePackage p0;
  galaxy::TPFCorePackage p1;
  p0.init_from_config(c_default);
  p1.init_from_config(c_explicit);

  const galaxy::State s = one_body(2.0, 1.0);
  std::vector<double> ax0, ay0, ax1, ay1;
  p0.compute_accelerations(s, /*bh_mass=*/5.0, /*softening=*/0.1, /*star_star=*/false, ax0, ay0);
  p1.compute_accelerations(s, /*bh_mass=*/5.0, /*softening=*/0.1, /*star_star=*/false, ax1, ay1);

  REQUIRE(ax0.size() == 1);
  REQUIRE(ax1.size() == 1);
  CHECK(ax0[0] == doctest::Approx(ax1[0]));
  CHECK(ay0[0] == doctest::Approx(ay1[0]));
}

TEST_CASE("direct_tpf xi_constraint runtime source rejects star_star") {
  const galaxy::State s = one_body(2.0, 1.0);
  galaxy::TPFCorePackage pkg = make_pkg("xi_constraint_exterior_single_source", 2.0e4, 0.0);
  std::vector<double> ax, ay;
  CHECK_THROWS_WITH(pkg.compute_accelerations(s, 5.0, 0.1, true, ax, ay),
                    doctest::Contains("experimental planar single-source-only runtime Xi solve"));
}

TEST_CASE("direct_tpf xi_constraint runtime source rejects outside domain and inner excision") {
  galaxy::TPFCorePackage pkg = make_pkg("xi_constraint_exterior_single_source", 2.0e4, 0.0);
  std::vector<double> ax, ay;

  const galaxy::State outside = one_body(11.0, 0.0);
  CHECK_THROWS_WITH(pkg.compute_accelerations(outside, 5.0, 0.1, false, ax, ay),
                    doctest::Contains("outside solver grid domain"));

  const galaxy::State inner = one_body(0.1, 0.0);
  CHECK_THROWS_WITH(pkg.compute_accelerations(inner, 5.0, 0.1, false, ax, ay),
                    doctest::Contains("inside inner excision radius"));
}

TEST_CASE("direct_tpf xi_constraint runtime source returns finite accelerations for supported single-source case") {
  galaxy::TPFCorePackage pkg = make_pkg("xi_constraint_exterior_single_source", 2.0e4, 0.0);
  const galaxy::State s = one_body(2.0, 1.0);
  std::vector<double> ax, ay;
  pkg.compute_accelerations(s, /*bh_mass=*/5.0, /*softening=*/0.1, /*star_star=*/false, ax, ay);

  REQUIRE(ax.size() == 1);
  REQUIRE(ay.size() == 1);
  CHECK(std::isfinite(ax[0]));
  CHECK(std::isfinite(ay[0]));
}

TEST_CASE("direct_tpf xi_constraint runtime source keeps VDSG additive after baseline") {
  const galaxy::State s = one_body(3.0, 0.0, 0.0, 1200.0);

  galaxy::TPFCorePackage base = make_pkg("xi_constraint_exterior_single_source", /*kappa=*/0.0, /*vdsg=*/0.0);
  galaxy::TPFCorePackage with_vdsg = make_pkg("xi_constraint_exterior_single_source", /*kappa=*/0.0, /*vdsg=*/1.0e5);

  std::vector<double> ax0, ay0, ax1, ay1;
  base.compute_accelerations(s, /*bh_mass=*/5.0, /*softening=*/0.1, /*star_star=*/false, ax0, ay0);
  with_vdsg.compute_accelerations(s, /*bh_mass=*/5.0, /*softening=*/0.1, /*star_star=*/false, ax1, ay1);

  REQUIRE(ax0.size() == 1);
  REQUIRE(ax1.size() == 1);
  CHECK(ax0[0] == doctest::Approx(0.0));
  CHECK(ay0[0] == doctest::Approx(0.0));
  CHECK(std::fabs(ax1[0]) + std::fabs(ay1[0]) > 0.0);
}
