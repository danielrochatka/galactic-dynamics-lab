#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/derived_tpf_radial.hpp"
#include "physics/TPFCore/provisional_readout.hpp"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "types.hpp"

#include <cmath>
#include <vector>

namespace {

galaxy::State one_body_state(double x, double y, double vx, double vy, double m) {
  galaxy::State s;
  s.resize(1);
  s.x[0] = x;
  s.y[0] = y;
  s.vx[0] = vx;
  s.vy[0] = vy;
  s.mass[0] = m;
  return s;
}

}  // namespace

TEST_CASE("VDSG coupling bypasses readout mode for accelerations") {
  galaxy::Config c_a;
  c_a.tpfcore_enable_provisional_readout = true;
  c_a.tpfcore_readout_mode = "tensor_radial_projection";
  c_a.tpf_vdsg_coupling = 1.0;  // activate VDSG branch
  c_a.dt = 0.01;
  c_a.star_mass = 2.0;

  galaxy::Config c_b = c_a;
  c_b.tpfcore_readout_mode = "unknown_mode_that_readout_would_ignore";

  galaxy::TPFCorePackage p_a, p_b;
  p_a.init_from_config(c_a);
  p_b.init_from_config(c_b);

  galaxy::State s = one_body_state(10.0, 0.0, 0.0, 0.0, 2.0);
  std::vector<double> ax_a, ay_a, ax_b, ay_b;
  p_a.compute_accelerations(s, /*bh_mass=*/100.0, /*softening=*/1.0, /*star_star=*/false, ax_a, ay_a);
  p_b.compute_accelerations(s, /*bh_mass=*/100.0, /*softening=*/1.0, /*star_star=*/false, ax_b, ay_b);

  REQUIRE(ax_a.size() == 1);
  REQUIRE(ax_b.size() == 1);
  CHECK(ax_a[0] == doctest::Approx(ax_b[0]));
  CHECK(ay_a[0] == doctest::Approx(ay_b[0]));
}

TEST_CASE("tr_coherence_readout and derived_tpf_radial_readout share acceleration closure") {
  galaxy::Config c_tr;
  c_tr.tpfcore_enable_provisional_readout = true;
  c_tr.tpfcore_readout_mode = "tr_coherence_readout";
  c_tr.tpf_vdsg_coupling = 0.0;  // force readout branch
  c_tr.tpf_kappa = 123.0;
  c_tr.tpf_poisson_bins = 64;
  c_tr.galaxy_radius = 50.0;

  galaxy::Config c_dr = c_tr;
  c_dr.tpfcore_readout_mode = "derived_tpf_radial_readout";

  galaxy::TPFCorePackage p_tr, p_dr;
  p_tr.init_from_config(c_tr);
  p_dr.init_from_config(c_dr);

  galaxy::State s = one_body_state(7.0, 3.0, 1.0, -0.2, 5.0);
  std::vector<double> ax_tr, ay_tr, ax_dr, ay_dr;
  p_tr.compute_accelerations(s, /*bh_mass=*/200.0, /*softening=*/0.5, /*star_star=*/false, ax_tr, ay_tr);
  p_dr.compute_accelerations(s, /*bh_mass=*/200.0, /*softening=*/0.5, /*star_star=*/false, ax_dr, ay_dr);

  REQUIRE(ax_tr.size() == 1);
  REQUIRE(ax_dr.size() == 1);
  CHECK(ax_tr[0] == doctest::Approx(ax_dr[0]));
  CHECK(ay_tr[0] == doctest::Approx(ay_dr[0]));
}

TEST_CASE("derived-radial readout acceleration equals radial_acceleration_scalar_derived projection") {
  galaxy::State s = one_body_state(12.0, 5.0, 0.0, 0.0, 2.0);
  const double bh_mass = 50.0;
  const double softening = 0.7;
  galaxy::tpfcore::DerivedTpfPoissonConfig cfg;
  cfg.kappa = 77.0;
  cfg.bins = 32;
  cfg.max_radius = 40.0;
  cfg.galaxy_radius = 40.0;
  auto profile = galaxy::tpfcore::build_tpf_gravity_profile(s, bh_mass, cfg, softening);

  double ax = 0.0, ay = 0.0;
  galaxy::tpfcore::compute_provisional_readout_acceleration(
      s,
      /*i=*/0,
      bh_mass,
      /*star_star=*/false,
      softening,
      /*source_softening=*/softening,
      "derived_tpf_radial_readout",
      /*readout_scale=*/1.0,
      /*theta_tt_scale=*/1.0,
      /*theta_tr_scale=*/1.0,
      ax,
      ay,
      &cfg,
      &profile);

  const double x = s.x[0];
  const double y = s.y[0];
  const double r_soft = std::sqrt(x * x + y * y + softening * softening);
  const double r_cyl = std::hypot(x, y);
  const double a_s = galaxy::tpfcore::radial_acceleration_scalar_derived(s, bh_mass, profile, r_cyl, softening);
  const double expected_ax = a_s * (x / r_soft);
  const double expected_ay = a_s * (y / r_soft);

  CHECK(ax == doctest::Approx(expected_ax));
  CHECK(ay == doctest::Approx(expected_ay));
}
