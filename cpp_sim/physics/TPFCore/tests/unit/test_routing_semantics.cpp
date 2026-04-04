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

/** Clean audit default: global |a| shunt off (independent of λ). */
void apply_clean_tpf_defaults(galaxy::Config& c) {
  c.tpf_global_accel_shunt_enable = false;
}

}  // namespace

TEST_CASE("Task A: coupling 0 vs tiny nonzero — shunt off, zero shunt events, near-identical ax when v=0") {
  galaxy::Config c0;
  c0.tpfcore_enable_provisional_readout = true;
  c0.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c0.tpf_vdsg_coupling = 0.0;
  c0.dt = 0.01;
  c0.star_mass = 2.0;
  c0.tpf_kappa = 1.0e10;
  c0.tpf_poisson_bins = 32;
  c0.galaxy_radius = 100.0;
  apply_clean_tpf_defaults(c0);

  galaxy::Config c1 = c0;
  c1.tpf_vdsg_coupling = 1e-8;

  galaxy::TPFCorePackage p0, p1;
  p0.init_from_config(c0);
  p1.init_from_config(c1);

  galaxy::State s = one_body_state(10.0, 0.0, 0.0, 0.0, 2.0);
  std::vector<double> ax0, ay0, ax1, ay1;

  galaxy::tpf_test_reset_global_accel_shunt_events();
  p0.compute_accelerations(s, 100.0, 1.0, false, ax0, ay0);
  const unsigned n0 = galaxy::tpf_test_global_accel_shunt_events();

  galaxy::tpf_test_reset_global_accel_shunt_events();
  p1.compute_accelerations(s, 100.0, 1.0, false, ax1, ay1);
  const unsigned n1 = galaxy::tpf_test_global_accel_shunt_events();

  CHECK(n0 == 0u);
  CHECK(n1 == 0u);
  CHECK(ax1[0] == doctest::Approx(ax0[0]).epsilon(1e-9));
  CHECK(ay1[0] == doctest::Approx(ay0[0]).epsilon(1e-9));
}

TEST_CASE("Task C: λ=0 vs λ=1e-8 with shunt disabled — differences at float noise (VDSG term negligible at v=0)") {
  galaxy::Config c0;
  apply_clean_tpf_defaults(c0);
  c0.tpfcore_enable_provisional_readout = true;
  c0.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c0.tpf_vdsg_coupling = 0.0;
  c0.dt = 0.01;
  c0.star_mass = 2.0;
  c0.tpf_kappa = 1.0e10;
  c0.tpf_poisson_bins = 32;
  c0.galaxy_radius = 100.0;

  galaxy::Config c1 = c0;
  c1.tpf_vdsg_coupling = 1e-8;

  galaxy::TPFCorePackage p0, p1;
  p0.init_from_config(c0);
  p1.init_from_config(c1);

  galaxy::State s = one_body_state(10.0, 0.0, 0.0, 0.0, 2.0);
  std::vector<double> ax0, ay0, ax1, ay1;
  p0.compute_accelerations(s, 100.0, 1.0, false, ax0, ay0);
  p1.compute_accelerations(s, 100.0, 1.0, false, ax1, ay1);

  CHECK(ax1[0] == doctest::Approx(ax0[0]).epsilon(1e-9));
  CHECK(ay1[0] == doctest::Approx(ay0[0]).epsilon(1e-9));
  CHECK(p0.last_accel_pipeline_stats().shunt_events_last_step == 0u);
  CHECK(p1.last_accel_pipeline_stats().shunt_events_last_step == 0u);
}

TEST_CASE("Task C: λ=0 vs λ=1e-8 with shunt enabled — same policy; shunt events match at v=0") {
  galaxy::Config c0;
  c0.tpfcore_enable_provisional_readout = true;
  c0.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c0.tpf_vdsg_coupling = 0.0;
  c0.dt = 0.01;
  c0.star_mass = 2.0;
  c0.tpf_kappa = 1.0e10;
  c0.tpf_poisson_bins = 32;
  c0.galaxy_radius = 100.0;
  c0.tpf_global_accel_shunt_enable = true;
  c0.tpf_global_accel_shunt_fraction = 0.001;

  galaxy::Config c1 = c0;
  c1.tpf_vdsg_coupling = 1e-8;

  galaxy::TPFCorePackage p0, p1;
  p0.init_from_config(c0);
  p1.init_from_config(c1);

  galaxy::State s = one_body_state(10.0, 0.0, 0.0, 0.0, 2.0);
  std::vector<double> ax0, ay0, ax1, ay1;

  galaxy::tpf_test_reset_global_accel_shunt_events();
  p0.compute_accelerations(s, 100.0, 1.0, false, ax0, ay0);
  const unsigned n0 = galaxy::tpf_test_global_accel_shunt_events();

  galaxy::tpf_test_reset_global_accel_shunt_events();
  p1.compute_accelerations(s, 100.0, 1.0, false, ax1, ay1);
  const unsigned n1 = galaxy::tpf_test_global_accel_shunt_events();

  CHECK(n0 == n1);
  CHECK(ax1[0] == doctest::Approx(ax0[0]).epsilon(1e-9));
  CHECK(ay1[0] == doctest::Approx(ay0[0]).epsilon(1e-9));
  CHECK(p0.last_accel_pipeline_stats().shunt_enabled == true);
  CHECK(p1.last_accel_pipeline_stats().shunt_enabled == true);
}

TEST_CASE("Task B: tangential motion — VDSG modifier scales with |v|/c (nonzero for circular-like flow)") {
  galaxy::Config c0;
  apply_clean_tpf_defaults(c0);
  c0.tpfcore_enable_provisional_readout = true;
  c0.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c0.tpf_vdsg_coupling = 0.0;
  c0.dt = 0.01;
  c0.star_mass = 2.0;
  c0.tpf_kappa = 1.0e10;
  c0.tpf_poisson_bins = 32;
  c0.galaxy_radius = 100.0;

  galaxy::Config c1 = c0;
  /* Large enough λ with |v|/c so Δax is above double noise (old (v·r̂)/c gave 0 here). */
  c1.tpf_vdsg_coupling = 0.02;

  galaxy::TPFCorePackage p0, p1;
  p0.init_from_config(c0);
  p1.init_from_config(c1);

  /* On +x axis from BH; velocity purely tangential (+y) => old (v·r̂) would give zero excess. */
  galaxy::State s = one_body_state(10.0, 0.0, 0.0, 8000.0, 2.0);
  std::vector<double> ax0, ay0, ax1, ay1;
  p0.compute_accelerations(s, 100.0, 1.0, false, ax0, ay0);
  p1.compute_accelerations(s, 100.0, 1.0, false, ax1, ay1);

  REQUIRE(ax0.size() == 1);
  CHECK(std::fabs(ax1[0] - ax0[0]) > 1e-20);
}

TEST_CASE("VDSG coupling is additive: tiny coupling gives accelerations near baseline-only") {
  galaxy::Config c0;
  apply_clean_tpf_defaults(c0);
  c0.tpfcore_enable_provisional_readout = true;
  c0.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c0.tpf_vdsg_coupling = 0.0;
  c0.dt = 0.01;
  c0.star_mass = 2.0;
  c0.tpf_kappa = 1.0e10;
  c0.tpf_poisson_bins = 32;
  c0.galaxy_radius = 100.0;

  galaxy::Config c1 = c0;
  c1.tpf_vdsg_coupling = 1e-45;

  galaxy::TPFCorePackage p0, p1;
  p0.init_from_config(c0);
  p1.init_from_config(c1);

  galaxy::State s = one_body_state(10.0, 0.0, 0.0, 0.0, 2.0);
  std::vector<double> ax0, ay0, ax1, ay1;
  p0.compute_accelerations(s, /*bh_mass=*/100.0, /*softening=*/1.0, /*star_star=*/false, ax0, ay0);
  p1.compute_accelerations(s, /*bh_mass=*/100.0, /*softening=*/1.0, /*star_star=*/false, ax1, ay1);

  REQUIRE(ax0.size() == 1);
  REQUIRE(ax1.size() == 1);
  CHECK(ax1[0] == doctest::Approx(ax0[0]).epsilon(1e-9));
  CHECK(ay1[0] == doctest::Approx(ay0[0]).epsilon(1e-9));
}

TEST_CASE("readout mode affects accelerations when VDSG coupling is shared (additive)") {
  galaxy::Config c_tensor;
  c_tensor.tpfcore_enable_provisional_readout = true;
  c_tensor.tpfcore_readout_mode = "tensor_radial_projection";
  c_tensor.tpf_vdsg_coupling = 1e-20;
  c_tensor.dt = 0.01;
  c_tensor.star_mass = 2.0;
  apply_clean_tpf_defaults(c_tensor);

  galaxy::Config c_derived = c_tensor;
  c_derived.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c_derived.tpf_kappa = 1.0e10;
  c_derived.tpf_poisson_bins = 32;
  c_derived.galaxy_radius = 100.0;

  galaxy::TPFCorePackage p_t, p_d;
  p_t.init_from_config(c_tensor);
  p_d.init_from_config(c_derived);

  /* Generic off-axis point so tensor vs derived closures do not coincide. */
  galaxy::State s = one_body_state(12.0, 5.0, 1.0, -0.2, 5.0);
  std::vector<double> ax_t, ay_t, ax_d, ay_d;
  p_t.compute_accelerations(s, 200.0, 0.5, false, ax_t, ay_t);
  p_d.compute_accelerations(s, 200.0, 0.5, false, ax_d, ay_d);

  REQUIRE(ax_t.size() == 1);
  REQUIRE(ax_d.size() == 1);
  const double diff = std::hypot(ax_t[0] - ax_d[0], ay_t[0] - ay_d[0]);
  CHECK(diff > 1e-20);
}

TEST_CASE("tr_coherence_readout and derived_tpf_radial_readout share acceleration closure") {
  galaxy::Config c_tr;
  c_tr.tpfcore_enable_provisional_readout = true;
  c_tr.tpfcore_readout_mode = "tr_coherence_readout";
  c_tr.tpf_vdsg_coupling = 0.0;  // force readout branch
  c_tr.tpf_kappa = 123.0;
  c_tr.tpf_poisson_bins = 64;
  c_tr.galaxy_radius = 50.0;
  apply_clean_tpf_defaults(c_tr);

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

TEST_CASE("v11_weak_field_truncation dynamics: one-source acceleration matches alpha*M/r^2 radial law") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "v11_weak_field_truncation";
  c.tpf_weak_field_correspondence_alpha_si = -2.0;
  c.tpfcore_enable_provisional_readout = false;
  c.tpf_vdsg_coupling = 0.0;
  c.tpf_cooling_fraction = 0.0;
  c.tpf_global_accel_shunt_enable = false;

  galaxy::TPFCorePackage p;
  p.init_from_config(c);

  galaxy::State s = one_body_state(3.0, 4.0, 0.0, 0.0, 1.0);
  std::vector<double> ax, ay;
  const double bh_mass = 5.0;
  p.compute_accelerations(s, bh_mass, /*softening=*/0.0, /*star_star=*/false, ax, ay);

  REQUIRE(ax.size() == 1);
  const double r = 5.0;
  const double coeff = c.tpf_weak_field_correspondence_alpha_si * bh_mass / (r * r * r);
  CHECK(ax[0] == doctest::Approx(coeff * 3.0));
  CHECK(ay[0] == doctest::Approx(coeff * 4.0));
}

TEST_CASE("v11_weak_field_truncation dynamics: pair superposition reproduces Newtonian-like symmetry for alpha=-G") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "v11_weak_field_truncation";
  c.tpf_weak_field_correspondence_alpha_si = -galaxy::tpfcore::TPF_G_SI;
  c.tpfcore_enable_provisional_readout = false;
  c.tpf_vdsg_coupling = 0.0;
  c.tpf_cooling_fraction = 0.0;
  c.tpf_global_accel_shunt_enable = false;

  galaxy::TPFCorePackage p;
  p.init_from_config(c);

  galaxy::State s;
  s.resize(2);
  s.x[0] = -1.0; s.y[0] = 0.0; s.mass[0] = 10.0; s.vx[0] = s.vy[0] = 0.0;
  s.x[1] =  1.0; s.y[1] = 0.0; s.mass[1] = 10.0; s.vx[1] = s.vy[1] = 0.0;
  std::vector<double> ax, ay;
  p.compute_accelerations(s, /*bh_mass=*/0.0, /*softening=*/0.0, /*star_star=*/true, ax, ay);
  REQUIRE(ax.size() == 2);
  CHECK(ax[0] == doctest::Approx(-ax[1]).epsilon(1e-12));
  CHECK(ay[0] == doctest::Approx(0.0));
  CHECK(ay[1] == doctest::Approx(0.0));
}

TEST_CASE("v11_weak_field_truncation dynamics rejects exploratory/provisional/stabilizer knobs") {
  galaxy::State s = one_body_state(10.0, 0.0, 0.0, 0.0, 1.0);
  std::vector<double> ax, ay;

  auto make_base = []() {
    galaxy::Config c;
    c.tpf_dynamics_mode = "v11_weak_field_truncation";
    c.tpfcore_enable_provisional_readout = false;
    c.tpf_vdsg_coupling = 0.0;
    c.tpf_cooling_fraction = 0.0;
    c.tpf_global_accel_shunt_enable = false;
    return c;
  };

  {
    auto c = make_base();
    c.tpf_vdsg_coupling = 1e-9;
    galaxy::TPFCorePackage p; p.init_from_config(c);
    CHECK_THROWS(p.compute_accelerations(s, 1.0, 0.0, false, ax, ay));
  }
  {
    auto c = make_base();
    c.tpfcore_enable_provisional_readout = true;
    galaxy::TPFCorePackage p; p.init_from_config(c);
    CHECK_THROWS(p.compute_accelerations(s, 1.0, 0.0, false, ax, ay));
  }
  {
    auto c = make_base();
    c.tpf_global_accel_shunt_enable = true;
    galaxy::TPFCorePackage p; p.init_from_config(c);
    CHECK_THROWS(p.compute_accelerations(s, 1.0, 0.0, false, ax, ay));
  }
  {
    auto c = make_base();
    c.tpf_cooling_fraction = 0.1;
    galaxy::TPFCorePackage p; p.init_from_config(c);
    CHECK_THROWS(p.compute_accelerations(s, 1.0, 0.0, false, ax, ay));
  }
}
