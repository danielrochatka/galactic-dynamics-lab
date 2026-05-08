#include "config.hpp"
#include "doctest.h"
#include "force_compare.hpp"
#include <algorithm>
#include <fstream>

using galaxy::apply_config_kv;
using galaxy::Config;

TEST_CASE("config defaults") {
  Config c;
  CHECK(c.simulation_mode == galaxy::SimulationMode::galaxy);
  CHECK(c.physics_package == "Newtonian");
  CHECK(c.physics_package_compare == "");
  CHECK(c.tpf_vdsg_coupling == doctest::Approx(0.0));
  CHECK(c.tpf_global_accel_shunt_enable == false);
  CHECK(c.tpf_global_accel_shunt_fraction == doctest::Approx(0.001));
  CHECK(c.tpf_cooling_fraction == doctest::Approx(0.0));
  CHECK(c.tpf_accel_pipeline_diagnostics_csv == true);
  CHECK(c.tpf_source_benchmark_shape == "monopole");
  CHECK(c.tpf_source_benchmark_total_mass == doctest::Approx(galaxy::kDefaultBhMassKg));
  CHECK(c.tpf_source_benchmark_mass1 == doctest::Approx(0.0));
  CHECK(c.tpf_source_benchmark_mass2 == doctest::Approx(0.0));
  CHECK(c.tpf_source_benchmark_separation == doctest::Approx(1.0e11));
  CHECK(c.tpf_source_benchmark_orientation_deg == doctest::Approx(0.0));
  CHECK(c.tpf_source_probe_grid_half_extent == doctest::Approx(2.0e11));
  CHECK(c.tpf_source_probe_grid_n == 121);
  CHECK(c.tpf_source_residual_exclusion_radius == doctest::Approx(0.25));
  CHECK(c.tpf_4d_residual_grid_n == 33);
  CHECK(c.tpf_4d_residual_grid_half_extent == doctest::Approx(20.0));
  CHECK(c.tpf_4d_residual_source_exclusion_radius == doctest::Approx(2.0));
  CHECK(c.tpf_4d_residual_field_softening == doctest::Approx(0.1));
  CHECK(c.tpf_4d_residual_bin_count == 32);
  CHECK(c.tpf_4d_residual_bin_radius_max == doctest::Approx(0.0));
  CHECK(c.tpf_4d_motion_probe_grid_n == 33);
  CHECK(c.tpf_4d_motion_probe_grid_half_extent == doctest::Approx(20.0));
  CHECK(c.tpf_4d_motion_source_exclusion_radius == doctest::Approx(3.0));
  CHECK(c.tpf_4d_motion_field_softening == doctest::Approx(0.1));
  CHECK(c.tpf_4d_motion_kappa == doctest::Approx(1.0));
  CHECK(c.tpf_4d_motion_readout_scale == doctest::Approx(1.0));
  CHECK(c.tpf_4d_motion_bin_count == 32);
  CHECK(c.tpf_4d_xi_motion_dt == doctest::Approx(0.001));
  CHECK(c.tpf_4d_xi_motion_steps == 2000);
  CHECK(c.tpf_4d_xi_motion_readout_scale == doctest::Approx(1.0e-12));
  CHECK(c.tpf_4d_xi_motion_field_softening == doctest::Approx(0.1));
  CHECK(c.tpf_4d_xi_motion_source_exclusion_radius == doctest::Approx(0.5));
  CHECK(c.tpf_4d_xi_motion_probe_layout == "ring");
  CHECK(c.tpf_4d_xi_motion_probe_count == 24);
  CHECK(c.tpf_4d_xi_motion_probe_radius == doctest::Approx(10.0));
  CHECK(c.tpf_4d_xi_motion_probe_speed == doctest::Approx(0.0));
  CHECK(c.tpf_4d_xi_motion_integrator == "velocity_verlet");
  CHECK(c.tpf_4d_xi_motion_dump_every == 1);
  CHECK(c.tpf_4d_xi_kernel_mode == "off");
  CHECK(c.tpf_4d_xi_kernel_coupling == doctest::Approx(0.0));
  CHECK(c.tpf_4d_xi_kernel_beta_power == doctest::Approx(1.0));
  CHECK(c.tpf_4d_xi_kernel_factor_mode == "beta_power");
  CHECK(c.tpf_4d_xi_kernel_metric_min == doctest::Approx(0.1));
  CHECK(c.tpf_4d_xi_kernel_metric_max == doctest::Approx(10.0));
  CHECK(c.tpf_4d_xi_temporal_mode == "off");
  CHECK(c.tpf_4d_xi_temporal_coupling == doctest::Approx(0.0));
  CHECK(c.tpf_4d_xi_source_speed_x == doctest::Approx(0.0));
  CHECK(c.tpf_4d_xi_source_speed_y == doctest::Approx(0.0));
  CHECK(c.tpf_4d_xi_source_speed_z == doctest::Approx(0.0));
  CHECK(c.tpf_xi_constraint_exterior_inspect == false);
  CHECK(c.tpf_xi_constraint_grid_n == 65);
  CHECK(c.tpf_xi_constraint_half_extent == doctest::Approx(10.0));
  CHECK(c.tpf_xi_constraint_inner_radius == doctest::Approx(0.5));
  CHECK(c.tpf_xi_constraint_max_iters == 250);
  CHECK(c.tpf_xi_constraint_tol == doctest::Approx(1e-8));
  CHECK(c.render_overlay_mode == "none");
  CHECK(c.display_distance_unit == "auto");
  CHECK(c.display_time_unit == "auto");
  CHECK(c.display_velocity_unit == "auto");
  CHECK(c.display_units_in_overlay == true);
  CHECK(c.display_show_unit_reference == true);
  CHECK(c.softening == doctest::Approx(0.0));
  CHECK(c.galaxy_init_velocity_mode == "enclosed_mass");
}

TEST_CASE("physics_package_compare parsing") {
  Config c;
  CHECK(apply_config_kv("physics_package_compare", "Newtonian", c));
  CHECK(c.physics_package_compare == "Newtonian");
}

TEST_CASE("legacy alias tpf_gdd_coupling is rejected") {
  Config c;
  CHECK(!apply_config_kv("tpf_gdd_coupling", "3.125e-10", c));
}

TEST_CASE("render_overlay_mode parsing") {
  Config c;
  CHECK(apply_config_kv("render_overlay_mode", "minimal", c));
  CHECK(c.render_overlay_mode == "minimal");
  CHECK(apply_config_kv("render_overlay_mode", "audit_full", c));
  CHECK(c.render_overlay_mode == "audit_full");
}

TEST_CASE("plot diagnostics cutoff parsing") {
  Config c;
  CHECK(apply_config_kv("diagnostic_cutoff_radius", "123.5", c));
  CHECK(c.diagnostic_cutoff_radius == doctest::Approx(123.5));
}

TEST_CASE("display unit config parsing") {
  Config c;
  CHECK(apply_config_kv("display_distance_unit", "AU", c));
  CHECK(c.display_distance_unit == "AU");
  CHECK(apply_config_kv("display_time_unit", "day", c));
  CHECK(c.display_time_unit == "day");
  CHECK(apply_config_kv("display_velocity_unit", "km/s", c));
  CHECK(c.display_velocity_unit == "km/s");
  CHECK(apply_config_kv("display_units_in_overlay", "false", c));
  CHECK(c.display_units_in_overlay == false);
  CHECK(apply_config_kv("display_show_unit_reference", "0", c));
  CHECK(c.display_show_unit_reference == false);
}

TEST_CASE("display unit config rejects invalid values") {
  Config c;
  CHECK_THROWS(apply_config_kv("display_distance_unit", "mile", c));
  CHECK_THROWS(apply_config_kv("display_time_unit", "week", c));
  CHECK_THROWS(apply_config_kv("display_velocity_unit", "mph", c));
}

TEST_CASE("galaxy init keys") {
  Config c;
  CHECK(apply_config_kv("galaxy_init_template", "preformed_spiral", c));
  CHECK(c.galaxy_init_template == "preformed_spiral");
  CHECK(apply_config_kv("galaxy_init_seed", "424242", c));
  CHECK(c.galaxy_init_seed == 424242u);
  CHECK(apply_config_kv("galaxy_init_master_chaos", "2.5", c));
  CHECK(c.galaxy_init_master_chaos == doctest::Approx(2.5));
  CHECK(apply_config_kv("galaxy_init_velocity_mode", "enclosed_mass", c));
  CHECK(c.galaxy_init_velocity_mode == "enclosed_mass");
  CHECK(apply_config_kv("galaxy_init_velocity_mode", "pairwise_radial_equilibrium", c));
  CHECK(c.galaxy_init_velocity_mode == "pairwise_radial_equilibrium");
  CHECK_THROWS(apply_config_kv("galaxy_init_velocity_mode", "bad_mode", c));
}

TEST_CASE("explicit tpf_vdsg_coupling still works") {
  Config c;
  CHECK(apply_config_kv("tpf_vdsg_coupling", "1e-99", c));
  CHECK(c.tpf_vdsg_coupling == doctest::Approx(1e-99));
}

TEST_CASE("tpf_kappa key parses unchanged and closure alias remains supported") {
  Config c;
  CHECK(apply_config_kv("tpf_kappa", "7.5e11", c));
  CHECK(c.tpf_kappa == doctest::Approx(7.5e11));
  CHECK(apply_config_kv("tpfcore_closure_kappa", "4.5e12", c));
  CHECK(c.tpf_kappa == doctest::Approx(4.5e12));
}

TEST_CASE("tpfcore_closure_kappa key maps to same tpf_kappa storage") {
  Config c;
  CHECK(apply_config_kv("tpfcore_closure_kappa", "4.5e12", c));
  CHECK(c.tpf_kappa == doctest::Approx(4.5e12));
}


TEST_CASE("Config defaults TPFCore dynamics to tpf_xi_theta_v1") {
  galaxy::Config c;
  CHECK(c.tpf_dynamics_mode == "tpf_xi_theta_v1");
}

TEST_CASE("geodesic_correspondence_baseline preset is rejected on tpf_xi_theta_v1 branch") {
  Config c;
  CHECK_THROWS(galaxy::load_config_file("../configs/geodesic_correspondence_baseline.cfg", c));
}
TEST_CASE("tpf_dynamics_mode accepts only tpf_xi_theta_v1 on this branch") {
  Config c;
  CHECK(apply_config_kv("tpf_dynamics_mode", "tpf_xi_theta_v1", c));
  CHECK(c.tpf_dynamics_mode == "tpf_xi_theta_v1");
  CHECK_THROWS(apply_config_kv("tpf_dynamics_mode", "weak_field_correspondence", c));
  CHECK_THROWS(apply_config_kv("tpf_dynamics_mode", "direct_tpf", c));
  CHECK_THROWS(apply_config_kv("tpf_dynamics_mode", "xi_kernel_deformed", c));
  CHECK(apply_config_kv("tpf_weak_field_correspondence_alpha_si", "-6.0e-11", c));
  CHECK(c.tpf_weak_field_correspondence_alpha_si == doctest::Approx(-6.0e-11));
}


TEST_CASE("tpfcore_enable_provisional_readout parses true for non-galaxy diagnostic compatibility") {
  Config c;
  CHECK(apply_config_kv("tpfcore_enable_provisional_readout", "true", c));
  CHECK(c.tpfcore_enable_provisional_readout == true);
}
TEST_CASE("tpf_xi_constraint_exterior inspection config keys parse") {
  Config c;
  CHECK(apply_config_kv("tpf_xi_constraint_exterior_inspect", "true", c));
  CHECK(c.tpf_xi_constraint_exterior_inspect == true);
  CHECK(apply_config_kv("tpf_xi_constraint_grid_n", "49", c));
  CHECK(c.tpf_xi_constraint_grid_n == 49);
  CHECK(apply_config_kv("tpf_xi_constraint_half_extent", "12.5", c));
  CHECK(c.tpf_xi_constraint_half_extent == doctest::Approx(12.5));
  CHECK(apply_config_kv("tpf_xi_constraint_inner_radius", "0.75", c));
  CHECK(c.tpf_xi_constraint_inner_radius == doctest::Approx(0.75));
  CHECK(apply_config_kv("tpf_xi_constraint_max_iters", "300", c));
  CHECK(c.tpf_xi_constraint_max_iters == 300);
  CHECK(apply_config_kv("tpf_xi_constraint_tol", "2e-7", c));
  CHECK(c.tpf_xi_constraint_tol == doctest::Approx(2e-7));
}

TEST_CASE("tpf_source_field_benchmark config keys parse") {
  Config c;
  CHECK(apply_config_kv("tpf_source_benchmark_shape", "bonded_pair", c));
  CHECK(c.tpf_source_benchmark_shape == "bonded_pair");
  CHECK(apply_config_kv("tpf_source_benchmark_total_mass", "9.99e21", c));
  CHECK(c.tpf_source_benchmark_total_mass == doctest::Approx(9.99e21));
  CHECK(apply_config_kv("tpf_source_benchmark_mass1", "7.5e11", c));
  CHECK(c.tpf_source_benchmark_mass1 == doctest::Approx(7.5e11));
  CHECK(apply_config_kv("tpf_source_benchmark_mass2", "2.5e11", c));
  CHECK(c.tpf_source_benchmark_mass2 == doctest::Approx(2.5e11));
  CHECK(apply_config_kv("tpf_source_benchmark_separation", "42.0", c));
  CHECK(c.tpf_source_benchmark_separation == doctest::Approx(42.0));
  CHECK(apply_config_kv("tpf_source_benchmark_orientation_deg", "27.5", c));
  CHECK(c.tpf_source_benchmark_orientation_deg == doctest::Approx(27.5));
  CHECK(apply_config_kv("tpf_source_probe_grid_half_extent", "1234.5", c));
  CHECK(c.tpf_source_probe_grid_half_extent == doctest::Approx(1234.5));
  CHECK(apply_config_kv("tpf_source_probe_grid_n", "75", c));
  CHECK(c.tpf_source_probe_grid_n == 75);
  CHECK(apply_config_kv("tpf_source_residual_exclusion_radius", "0.42", c));
  CHECK(c.tpf_source_residual_exclusion_radius == doctest::Approx(0.42));
  CHECK_THROWS(apply_config_kv("tpf_source_benchmark_shape", "triangle", c));
}

TEST_CASE("tpf_source_field_benchmark keys serialize") {
  Config c;
  c.tpf_source_benchmark_mass1 = 3.0;
  c.tpf_source_benchmark_mass2 = 4.0;
  const auto kv = galaxy::serialize_config_kv(c);
  const auto has_key = [&kv](const char* key) {
    return std::any_of(
        kv.begin(), kv.end(),
        [key](const std::pair<std::string, std::string>& entry) { return entry.first == key; });
  };
  CHECK(has_key("tpf_source_benchmark_mass1"));
  CHECK(has_key("tpf_source_benchmark_mass2"));
  CHECK(has_key("tpf_source_residual_exclusion_radius"));
}

TEST_CASE("tpf_4d_static_residual_benchmark config keys parse and serialize") {
  Config c;
  CHECK(apply_config_kv("tpf_4d_residual_grid_n", "37", c));
  CHECK(c.tpf_4d_residual_grid_n == 37);
  CHECK(apply_config_kv("tpf_4d_residual_grid_half_extent", "12.5", c));
  CHECK(c.tpf_4d_residual_grid_half_extent == doctest::Approx(12.5));
  CHECK(apply_config_kv("tpf_4d_residual_source_exclusion_radius", "1.25", c));
  CHECK(c.tpf_4d_residual_source_exclusion_radius == doctest::Approx(1.25));
  CHECK(apply_config_kv("tpf_4d_residual_field_softening", "0.05", c));
  CHECK(c.tpf_4d_residual_field_softening == doctest::Approx(0.05));
  CHECK(apply_config_kv("tpf_4d_residual_bin_count", "48", c));
  CHECK(c.tpf_4d_residual_bin_count == 48);
  CHECK(apply_config_kv("tpf_4d_residual_bin_radius_max", "9.5", c));
  CHECK(c.tpf_4d_residual_bin_radius_max == doctest::Approx(9.5));

  const auto kv = galaxy::serialize_config_kv(c);
  const auto has_key = [&kv](const char* key) {
    return std::any_of(
        kv.begin(), kv.end(),
        [key](const std::pair<std::string, std::string>& entry) { return entry.first == key; });
  };
  CHECK(has_key("tpf_4d_residual_grid_n"));
  CHECK(has_key("tpf_4d_residual_grid_half_extent"));
  CHECK(has_key("tpf_4d_residual_source_exclusion_radius"));
  CHECK(has_key("tpf_4d_residual_field_softening"));
  CHECK(has_key("tpf_4d_residual_bin_count"));
  CHECK(has_key("tpf_4d_residual_bin_radius_max"));
}

TEST_CASE("tpf_4d_static_motion_readout_benchmark config keys parse and serialize") {
  Config c;
  CHECK(apply_config_kv("tpf_4d_motion_probe_grid_n", "35", c));
  CHECK(c.tpf_4d_motion_probe_grid_n == 35);
  CHECK(apply_config_kv("tpf_4d_motion_probe_grid_half_extent", "11.5", c));
  CHECK(c.tpf_4d_motion_probe_grid_half_extent == doctest::Approx(11.5));
  CHECK(apply_config_kv("tpf_4d_motion_source_exclusion_radius", "2.75", c));
  CHECK(c.tpf_4d_motion_source_exclusion_radius == doctest::Approx(2.75));
  CHECK(apply_config_kv("tpf_4d_motion_field_softening", "0.15", c));
  CHECK(c.tpf_4d_motion_field_softening == doctest::Approx(0.15));
  CHECK(apply_config_kv("tpf_4d_motion_kappa", "4.5", c));
  CHECK(c.tpf_4d_motion_kappa == doctest::Approx(4.5));
  CHECK(apply_config_kv("tpf_4d_motion_readout_scale", "0.25", c));
  CHECK(c.tpf_4d_motion_readout_scale == doctest::Approx(0.25));
  CHECK(apply_config_kv("tpf_4d_motion_bin_count", "12", c));
  CHECK(c.tpf_4d_motion_bin_count == 12);

  const auto kv = galaxy::serialize_config_kv(c);
  const auto has_key = [&kv](const char* key) {
    return std::any_of(
        kv.begin(), kv.end(),
        [key](const std::pair<std::string, std::string>& entry) { return entry.first == key; });
  };
  CHECK(has_key("tpf_4d_motion_probe_grid_n"));
  CHECK(has_key("tpf_4d_motion_probe_grid_half_extent"));
  CHECK(has_key("tpf_4d_motion_source_exclusion_radius"));
  CHECK(has_key("tpf_4d_motion_field_softening"));
  CHECK(has_key("tpf_4d_motion_kappa"));
  CHECK(has_key("tpf_4d_motion_readout_scale"));
  CHECK(has_key("tpf_4d_motion_bin_count"));
}

TEST_CASE("tpf_4d_xi_motion_probe_benchmark config keys parse and serialize") {
  Config c;
  CHECK(apply_config_kv("tpf_4d_xi_motion_dt", "0.0025", c));
  CHECK(c.tpf_4d_xi_motion_dt == doctest::Approx(0.0025));
  CHECK(apply_config_kv("tpf_4d_xi_motion_steps", "123", c));
  CHECK(c.tpf_4d_xi_motion_steps == 123);
  CHECK(apply_config_kv("tpf_4d_xi_motion_readout_scale", "9e-13", c));
  CHECK(c.tpf_4d_xi_motion_readout_scale == doctest::Approx(9e-13));
  CHECK(apply_config_kv("tpf_4d_xi_motion_field_softening", "0.2", c));
  CHECK(c.tpf_4d_xi_motion_field_softening == doctest::Approx(0.2));
  CHECK(apply_config_kv("tpf_4d_xi_motion_source_exclusion_radius", "0.75", c));
  CHECK(c.tpf_4d_xi_motion_source_exclusion_radius == doctest::Approx(0.75));
  CHECK(apply_config_kv("tpf_4d_xi_motion_probe_layout", "axis", c));
  CHECK(c.tpf_4d_xi_motion_probe_layout == "axis");
  CHECK(apply_config_kv("tpf_4d_xi_motion_probe_count", "64", c));
  CHECK(c.tpf_4d_xi_motion_probe_count == 64);
  CHECK(apply_config_kv("tpf_4d_xi_motion_probe_radius", "15.0", c));
  CHECK(c.tpf_4d_xi_motion_probe_radius == doctest::Approx(15.0));
  CHECK(apply_config_kv("tpf_4d_xi_motion_probe_speed", "2.0", c));
  CHECK(c.tpf_4d_xi_motion_probe_speed == doctest::Approx(2.0));
  CHECK(apply_config_kv("tpf_4d_xi_motion_integrator", "semi_implicit_euler", c));
  CHECK(c.tpf_4d_xi_motion_integrator == "semi_implicit_euler");
  CHECK(apply_config_kv("tpf_4d_xi_motion_dump_every", "3", c));
  CHECK(c.tpf_4d_xi_motion_dump_every == 3);
  CHECK(apply_config_kv("tpf_4d_xi_kernel_mode", "metric_velocity", c));
  CHECK(c.tpf_4d_xi_kernel_mode == "metric_velocity");
  CHECK(apply_config_kv("tpf_4d_xi_kernel_coupling", "11.0", c));
  CHECK(c.tpf_4d_xi_kernel_coupling == doctest::Approx(11.0));
  CHECK(apply_config_kv("tpf_4d_xi_kernel_beta_power", "1.5", c));
  CHECK(c.tpf_4d_xi_kernel_beta_power == doctest::Approx(1.5));
  CHECK(apply_config_kv("tpf_4d_xi_kernel_factor_mode", "gamma_minus_one", c));
  CHECK(c.tpf_4d_xi_kernel_factor_mode == "gamma_minus_one");
  CHECK(apply_config_kv("tpf_4d_xi_kernel_metric_min", "0.2", c));
  CHECK(c.tpf_4d_xi_kernel_metric_min == doctest::Approx(0.2));
  CHECK(apply_config_kv("tpf_4d_xi_kernel_metric_max", "12.0", c));
  CHECK(c.tpf_4d_xi_kernel_metric_max == doctest::Approx(12.0));
  CHECK(apply_config_kv("tpf_4d_xi_temporal_mode", "norm_scaled", c));
  CHECK(c.tpf_4d_xi_temporal_mode == "norm_scaled");
  CHECK(apply_config_kv("tpf_4d_xi_temporal_coupling", "0.25", c));
  CHECK(c.tpf_4d_xi_temporal_coupling == doctest::Approx(0.25));
  CHECK(apply_config_kv("tpf_4d_xi_source_speed_x", "1.0", c));
  CHECK(c.tpf_4d_xi_source_speed_x == doctest::Approx(1.0));
  CHECK(apply_config_kv("tpf_4d_xi_source_speed_y", "-2.0", c));
  CHECK(c.tpf_4d_xi_source_speed_y == doctest::Approx(-2.0));
  CHECK(apply_config_kv("tpf_4d_xi_source_speed_z", "3.5", c));
  CHECK(c.tpf_4d_xi_source_speed_z == doctest::Approx(3.5));

  const auto kv = galaxy::serialize_config_kv(c);
  const auto has_key = [&kv](const char* key) {
    return std::any_of(
        kv.begin(), kv.end(),
        [key](const std::pair<std::string, std::string>& entry) { return entry.first == key; });
  };
  CHECK(has_key("tpf_4d_xi_motion_dt"));
  CHECK(has_key("tpf_4d_xi_motion_steps"));
  CHECK(has_key("tpf_4d_xi_motion_readout_scale"));
  CHECK(has_key("tpf_4d_xi_motion_field_softening"));
  CHECK(has_key("tpf_4d_xi_motion_source_exclusion_radius"));
  CHECK(has_key("tpf_4d_xi_motion_probe_layout"));
  CHECK(has_key("tpf_4d_xi_motion_probe_count"));
  CHECK(has_key("tpf_4d_xi_motion_probe_radius"));
  CHECK(has_key("tpf_4d_xi_motion_probe_speed"));
  CHECK(has_key("tpf_4d_xi_motion_integrator"));
  CHECK(has_key("tpf_4d_xi_motion_dump_every"));
  CHECK(has_key("tpf_4d_xi_kernel_mode"));
  CHECK(has_key("tpf_4d_xi_kernel_coupling"));
  CHECK(has_key("tpf_4d_xi_kernel_beta_power"));
  CHECK(has_key("tpf_4d_xi_kernel_factor_mode"));
  CHECK(has_key("tpf_4d_xi_kernel_metric_min"));
  CHECK(has_key("tpf_4d_xi_kernel_metric_max"));
  CHECK(has_key("tpf_4d_xi_temporal_mode"));
  CHECK(has_key("tpf_4d_xi_temporal_coupling"));
  CHECK(has_key("tpf_4d_xi_source_speed_x"));
  CHECK(has_key("tpf_4d_xi_source_speed_y"));
  CHECK(has_key("tpf_4d_xi_source_speed_z"));
}


TEST_CASE("utility simulation modes remain parseable") {
  const char* modes[] = {
      "tpf_single_source_inspect", "tpf_symmetric_pair_inspect", "tpf_source_field_benchmark",
      "tpf_4d_static_residual_benchmark", "tpf_4d_static_motion_readout_benchmark",
      "tpf_4d_xi_motion_probe_benchmark", "tpf_weak_field_calibration",
      "tpf_newtonian_force_compare", "tpf_diagnostic_consistency_audit"};
  for (const char* m : modes) {
    galaxy::Config c;
    CHECK_NOTHROW(galaxy::apply_config_kv("simulation_mode", m, c));
  }
}

TEST_CASE("unsupported utility simulation modes are rejected") {
  galaxy::Config c;
  CHECK_THROWS(galaxy::apply_config_kv("simulation_mode", "tpf_two_body_sweep", c));
  CHECK_THROWS(galaxy::apply_config_kv("simulation_mode", "tpf_bound_orbit_sweep", c));
}

TEST_CASE("tpf utility mode dispatch classification") {
  CHECK(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::tpf_single_source_inspect));
  CHECK(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::tpf_symmetric_pair_inspect));
  CHECK(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::tpf_source_field_benchmark));
  CHECK(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::tpf_4d_static_residual_benchmark));
  CHECK(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::tpf_4d_static_motion_readout_benchmark));
  CHECK(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::tpf_4d_xi_motion_probe_benchmark));
  CHECK(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::tpf_weak_field_calibration));
  CHECK(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::tpf_newtonian_force_compare));
  CHECK(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::tpf_diagnostic_consistency_audit));
  CHECK_FALSE(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::galaxy));
  CHECK_FALSE(galaxy::is_tpf_utility_mode(galaxy::SimulationMode::earth_moon_benchmark));
}

TEST_CASE("diagnostic_consistency_audit reports unavailable on v1 branch") {
  galaxy::Config c;
  CHECK_FALSE(galaxy::run_tpf_diagnostic_consistency_audit(c, "/tmp"));
}
