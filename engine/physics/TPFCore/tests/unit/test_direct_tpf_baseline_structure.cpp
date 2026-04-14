#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/direct_tpf_baseline.hpp"
#include "physics/TPFCore/field_evaluation.hpp"
#include "physics/TPFCore/source_ansatz.hpp"
#include "physics/TPFCore/tpf_core_package.hpp"

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

TEST_CASE("direct_tpf baseline helper derives invariant I from Theta prior to Xi readout") {
  galaxy::State s = one_body_state(4.0, 3.0, 0.0, 0.0, 1.0);
  const galaxy::tpfcore::FieldAtPoint field =
      galaxy::tpfcore::evaluate_provisional_field_multi_source(s, 0, /*bh_mass=*/5.0, /*star_star=*/false, /*eps=*/0.1);

  const auto artifacts = galaxy::tpfcore::compute_direct_tpf_baseline_artifacts(
      field, /*kappa=*/2.0, galaxy::tpfcore::LAMBDA_4D);
  const auto readout = galaxy::tpfcore::compute_xi_directed_tensor_readout(artifacts.xi, artifacts.principal_cij);
  const auto baseline = galaxy::tpfcore::compute_direct_tpf_baseline_acceleration(
      field, /*kappa=*/2.0, galaxy::tpfcore::LAMBDA_4D);

  const auto derived_from_theta = galaxy::tpfcore::compute_invariant_I(artifacts.theta);
  CHECK(artifacts.invariant_I.value == doctest::Approx(derived_from_theta.value));
  CHECK(artifacts.invariant_I.value == doctest::Approx(field.invariant_I));
  CHECK(baseline.ax == doctest::Approx(readout.ax));
  CHECK(baseline.ay == doctest::Approx(readout.ay));
}

TEST_CASE("direct_tpf baseline helper keeps DeltaC explicit as zero placeholder") {
  const auto delta_c = galaxy::tpfcore::compute_deltaC_placeholder_zero();
  CHECK(delta_c.xx == doctest::Approx(0.0));
  CHECK(delta_c.xy == doctest::Approx(0.0));
  CHECK(delta_c.yy == doctest::Approx(0.0));
  CHECK(delta_c.zz == doctest::Approx(0.0));
}

TEST_CASE("direct_tpf baseline acceleration remains independent of correspondence alpha config") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_weak_field_correspondence_alpha_si = -2.0;
  c.tpfcore_enable_provisional_readout = false;
  c.tpf_vdsg_coupling = 0.0;
  c.tpf_cooling_fraction = 0.0;
  c.tpf_global_accel_shunt_enable = false;

  galaxy::Config c2 = c;
  c2.tpf_weak_field_correspondence_alpha_si = -9.0;

  galaxy::TPFCorePackage p1;
  galaxy::TPFCorePackage p2;
  p1.init_from_config(c);
  p2.init_from_config(c2);

  galaxy::State s = one_body_state(3.0, -5.0, 0.0, 0.0, 1.0);
  std::vector<double> ax1, ay1, ax2, ay2;
  p1.compute_accelerations(s, /*bh_mass=*/7.0, /*softening=*/0.0, /*star_star=*/false, ax1, ay1);
  p2.compute_accelerations(s, /*bh_mass=*/7.0, /*softening=*/0.0, /*star_star=*/false, ax2, ay2);

  REQUIRE(ax1.size() == 1);
  REQUIRE(ax2.size() == 1);
  CHECK(ax1[0] == doctest::Approx(ax2[0]));
  CHECK(ay1[0] == doctest::Approx(ay2[0]));
}

TEST_CASE("direct_tpf one-source softened central response remains attractive") {
  galaxy::State s = one_body_state(0.1, 0.0, 0.0, 0.0, 1.0);
  const double bh_mass = 5.0;
  const double softening = 0.1;
  const auto field = galaxy::tpfcore::evaluate_provisional_field_multi_source(
      s, 0, bh_mass, /*star_star=*/false, softening);
  const auto baseline = galaxy::tpfcore::compute_direct_tpf_baseline_acceleration(
      field, /*kappa=*/1.0, galaxy::tpfcore::LAMBDA_4D);
  const double radial_dot = baseline.ax * s.x[0] + baseline.ay * s.y[0];
  CHECK(radial_dot < 0.0);
}

TEST_CASE("direct_tpf applies VDSG after Eq10 baseline acceleration") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpfcore_enable_provisional_readout = false;
  c.tpf_kappa = 0.0;
  c.tpf_vdsg_coupling = 1.0e5;
  c.tpf_cooling_fraction = 0.0;
  c.tpf_global_accel_shunt_enable = false;

  galaxy::TPFCorePackage p;
  p.init_from_config(c);

  galaxy::State s = one_body_state(8.0, 0.0, 0.0, 6000.0, 2.0);
  const auto field = galaxy::tpfcore::evaluate_provisional_field_multi_source(s, 0, 20.0, false, 0.2);
  const auto baseline = galaxy::tpfcore::compute_direct_tpf_baseline_acceleration(
      field, /*kappa=*/0.0, galaxy::tpfcore::LAMBDA_4D);

  std::vector<double> ax, ay;
  p.compute_accelerations(s, /*bh_mass=*/20.0, /*softening=*/0.2, /*star_star=*/false, ax, ay);

  REQUIRE(ax.size() == 1);
  CHECK(baseline.ax == doctest::Approx(0.0));
  CHECK(baseline.ay == doctest::Approx(0.0));
  CHECK(std::fabs(ax[0]) + std::fabs(ay[0]) > 0.0);
}
