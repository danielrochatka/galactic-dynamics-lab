#include "config.hpp"
#include "doctest.h"

using galaxy::apply_config_kv;
using galaxy::Config;

TEST_CASE("config defaults") {
  Config c;
  CHECK(c.simulation_mode == galaxy::SimulationMode::galaxy);
  CHECK(c.physics_package == "Newtonian");
  CHECK(c.physics_package_compare == "");
  CHECK(c.compare_parallel == false);
  CHECK(c.tpf_vdsg_coupling == doctest::Approx(1.0e-20));
  CHECK(c.tpf_global_accel_shunt_enable == false);
  CHECK(c.tpf_global_accel_shunt_fraction == doctest::Approx(0.001));
  CHECK(c.tpf_accel_pipeline_diagnostics_csv == true);
  CHECK(c.render_overlay_mode == "none");
  CHECK(c.display_distance_unit == "auto");
  CHECK(c.display_time_unit == "auto");
  CHECK(c.display_velocity_unit == "auto");
  CHECK(c.display_units_in_overlay == true);
  CHECK(c.display_show_unit_reference == true);
  CHECK(c.softening == doctest::Approx(0.0));
}

TEST_CASE("physics_package_compare parsing") {
  Config c;
  CHECK(apply_config_kv("physics_package_compare", "Newtonian", c));
  CHECK(c.physics_package_compare == "Newtonian");
}

TEST_CASE("compare_parallel parsing") {
  Config c;
  CHECK(apply_config_kv("compare_parallel", "true", c));
  CHECK(c.compare_parallel == true);
}

TEST_CASE("legacy alias tpf_gdd_coupling maps to tpf_vdsg_coupling") {
  Config c;
  CHECK(apply_config_kv("tpf_gdd_coupling", "3.125e-10", c));
  CHECK(c.tpf_vdsg_coupling == doctest::Approx(3.125e-10));
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
}

TEST_CASE("explicit tpf_vdsg_coupling still works") {
  Config c;
  CHECK(apply_config_kv("tpf_vdsg_coupling", "1e-99", c));
  CHECK(c.tpf_vdsg_coupling == doctest::Approx(1e-99));
}

TEST_CASE("tpfcore_closure_kappa key maps to closure ledger coefficient") {
  Config c;
  CHECK(apply_config_kv("tpfcore_closure_kappa", "4.5e12", c));
  CHECK(c.tpf_kappa == doctest::Approx(4.5e12));
}

TEST_CASE("tpf_dynamics_mode accepts v11_weak_field_truncation and legacy weak_field_correspondence alias") {
  Config c;
  CHECK(apply_config_kv("tpf_dynamics_mode", "v11_weak_field_truncation", c));
  CHECK(c.tpf_dynamics_mode == "v11_weak_field_truncation");
  CHECK(apply_config_kv("tpf_dynamics_mode", "weak_field_correspondence", c));
  CHECK(c.tpf_dynamics_mode == "v11_weak_field_truncation");
  CHECK(apply_config_kv("tpf_dynamics_mode", "direct_tpf", c));
  CHECK(c.tpf_dynamics_mode == "direct_tpf");
  CHECK(apply_config_kv("tpf_weak_field_correspondence_alpha_si", "-6.0e-11", c));
  CHECK(c.tpf_weak_field_correspondence_alpha_si == doctest::Approx(-6.0e-11));
}

TEST_CASE("tpf_analysis_mode and simulation_mode tpf_v11_weak_field_correspondence") {
  Config c;
  CHECK(apply_config_kv("tpf_analysis_mode", "v11_weak_field_correspondence", c));
  CHECK(c.tpf_analysis_mode == "v11_weak_field_correspondence");
  CHECK(apply_config_kv("simulation_mode", "tpf_v11_weak_field_correspondence", c));
  CHECK(c.simulation_mode == galaxy::SimulationMode::tpf_v11_weak_field_correspondence);
}

TEST_CASE("v11_weak_field_correspondence_benchmark and Earth-Moon SI keys") {
  Config c;
  CHECK(c.v11_weak_field_correspondence_benchmark == "axis_monopole");
  CHECK(apply_config_kv("v11_weak_field_correspondence_benchmark", "earth_moon_line_of_centers", c));
  CHECK(c.v11_weak_field_correspondence_benchmark == "earth_moon_line_of_centers");
  CHECK(apply_config_kv("v11_em_mean_distance_m", "3.844e8", c));
  CHECK(c.v11_em_mean_distance_m == doctest::Approx(3.844e8));
  CHECK(apply_config_kv("v11_em_calib_surface_g_m_s2", "9.81", c));
  CHECK(c.v11_em_calib_surface_g_m_s2 == doctest::Approx(9.81));
}
