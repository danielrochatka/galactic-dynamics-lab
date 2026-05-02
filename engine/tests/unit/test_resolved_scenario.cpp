#include "config.hpp"
#include "doctest.h"
#include "output.hpp"
#include "resolved_scenario.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace {

std::string slurp(const std::string& path) {
  std::ifstream f(path);
  std::ostringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

}  // namespace

TEST_CASE("resolve earth_moon_benchmark centralizes effective conditions") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::earth_moon_benchmark;

  galaxy::ResolvedScenario r = galaxy::resolve_scenario(c);

  CHECK(r.initializer_used == "init_two_body");
  CHECK(r.config.bh_mass == doctest::Approx(0.0));
  CHECK(r.config.enable_star_star_gravity == true);
  CHECK(r.config.dt == doctest::Approx(3600.0));
  CHECK(r.config.softening == doctest::Approx(0.0));
  CHECK(r.effective_n_steps == 1440);
  CHECK(r.effective_snapshot_every == 6);
  CHECK(r.effective_total_sim_time == doctest::Approx(3600.0 * 1440.0));
  REQUIRE(r.initial_state.n() == 2);
  CHECK(r.initial_state.mass[0] == doctest::Approx(galaxy::kDefaultEarthMassKg));
  CHECK(r.initial_state.mass[1] == doctest::Approx(galaxy::kDefaultMoonMassKg));
}

TEST_CASE("resolve bh_orbit_validation and symmetric_pair mode defaults") {
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::bh_orbit_validation;
    c.enable_star_star_gravity = true;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.initializer_used == "init_two_body_star_around_bh");
    CHECK(r.config.enable_star_star_gravity == false);
    CHECK(r.config.dt == doctest::Approx(10000.0));
    CHECK(r.config.softening == doctest::Approx(0.0));
    CHECK(r.config.validation_two_body_radius == doctest::Approx(1.0e13));
    CHECK(r.effective_n_steps == 6000);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::symmetric_pair;
    c.validation_symmetric_include_bh = false;
    c.explicit_overrides.validation_symmetric_include_bh = true;
    c.bh_mass = 42.0;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.initializer_used == "init_symmetric_pair");
    CHECK(r.config.bh_mass == doctest::Approx(0.0));
    CHECK(r.config.dt == doctest::Approx(3600.0));
    CHECK(r.config.softening == doctest::Approx(0.0));
    CHECK(r.effective_snapshot_every == 6);
  }
}

TEST_CASE("resolve small_n_conservation and timestep_convergence use explicit mode timing/softening") {
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::small_n_conservation;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.config.dt == doctest::Approx(1.0e-4));
    CHECK(r.config.softening == doctest::Approx(0.0));
    CHECK(r.config.bh_mass == doctest::Approx(1.0e18));
    CHECK(r.effective_n_steps == 20000);
    CHECK(r.effective_snapshot_every == 20);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::timestep_convergence;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.config.dt == doctest::Approx(10000.0));
    CHECK(r.config.softening == doctest::Approx(0.0));
    CHECK(r.config.validation_two_body_radius == doctest::Approx(1.0e13));
    CHECK(r.effective_n_steps == 6000);
    CHECK(r.effective_snapshot_every == 10);
  }
}

TEST_CASE("user overrides survive scenario resolution") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::earth_moon_benchmark;
  c.validation_n_steps = 111;
  c.validation_snapshot_every = 7;
  c.bh_mass = 123.0;
  c.dt = 12.5;
  c.softening = 2.0;
  c.explicit_overrides.validation_n_steps = true;
  c.explicit_overrides.validation_snapshot_every = true;
  c.explicit_overrides.bh_mass = true;
  c.explicit_overrides.dt = true;
  c.explicit_overrides.softening = true;

  auto r = galaxy::resolve_scenario(c);
  CHECK(r.config.bh_mass == doctest::Approx(123.0));
  CHECK(r.config.dt == doctest::Approx(12.5));
  CHECK(r.config.softening == doctest::Approx(2.0));
  CHECK(r.effective_n_steps == 111);
  CHECK(r.effective_snapshot_every == 7);
}

TEST_CASE("snapshot_target resolves snapshot cadence automatically") {
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::galaxy;
    c.n_steps = 200000;
    c.snapshot_every = 77;
    c.snapshot_target = 200;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.effective_snapshot_every == 1000);
    CHECK(r.snapshot_target_active == 1);
    CHECK(r.configured_snapshot_target == 200);
    CHECK(r.resolved_snapshot_every == 1000);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::galaxy;
    c.n_steps = 50000;
    c.snapshot_target = 200;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.effective_snapshot_every == 250);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::galaxy;
    c.n_steps = 100;
    c.snapshot_target = 200;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.effective_snapshot_every == 1);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::galaxy;
    c.n_steps = 10;
    c.snapshot_target = 6;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.effective_snapshot_every == 2);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::galaxy;
    c.n_steps = 4321;
    c.snapshot_every = 37;
    c.snapshot_target = 0;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.effective_snapshot_every == 37);
    CHECK(r.snapshot_target_active == 0);
    CHECK(r.resolved_snapshot_every == 37);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::galaxy;
    c.n_steps = 100;
    c.snapshot_target = 1000;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.effective_snapshot_every == 1);
  }
}

TEST_CASE("galaxy explicit zero softening and source softening are preserved") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::galaxy;
  REQUIRE(galaxy::apply_config_kv("softening", "0", c));
  REQUIRE(galaxy::apply_config_kv("tpfcore_source_softening", "0", c));

  const auto r = galaxy::resolve_scenario(c);
  CHECK(r.config.softening == doctest::Approx(0.0));
  CHECK(r.config.tpfcore_source_softening == doctest::Approx(0.0));
}

TEST_CASE("resolved scenario artifact writer outputs expected keys") {
  galaxy::Config c;
  c.output_dir = "../outputs/test_resolved_scenario_artifact";
  c.simulation_mode = galaxy::SimulationMode::earth_moon_benchmark;
  c.physics_package = "Newtonian";
  galaxy::ResolvedScenario r = galaxy::resolve_scenario(c);

#ifdef _WIN32
  _mkdir("outputs");
  _mkdir(c.output_dir.c_str());
#else
  const int mk_ok = std::system((std::string("mkdir -p ") + c.output_dir).c_str());
  (void)mk_ok;
#endif

  galaxy::write_resolved_scenario_artifacts(c.output_dir, r);
  const std::string txt = slurp(c.output_dir + "/resolved_scenario.txt");
  const std::string js = slurp(c.output_dir + "/resolved_scenario.json");

  CHECK(txt.find("simulation_mode\tearth_moon_benchmark") != std::string::npos);
  CHECK(txt.find("effective_bh_mass\t0") != std::string::npos);
  CHECK(txt.find("effective_dt\t3600") != std::string::npos);
  CHECK(txt.find("effective_softening\t0") != std::string::npos);
  CHECK(txt.find("effective_total_sim_time\t") != std::string::npos);
  CHECK(txt.find("timing_policy\tearth_moon_hourly_step_60d_horizon") != std::string::npos);
  CHECK(js.find("\"initializer_used\": \"init_two_body\"") != std::string::npos);
  CHECK(js.find("\"n_steps\": 1440") != std::string::npos);
  CHECK(js.find("\"total_sim_time\": 5184000") != std::string::npos);
  CHECK(js.find("\"softening\": 0") != std::string::npos);
  CHECK(js.find("\"timing\": \"earth_moon_hourly_step_60d_horizon\"") != std::string::npos);
}

TEST_CASE("run_info audit includes configured and effective sections and resolved consistency") {
  {
    galaxy::Config configured;
    configured.output_dir = "../outputs/test_run_info_earth_moon";
    configured.simulation_mode = galaxy::SimulationMode::earth_moon_benchmark;
    configured.snapshot_every = 50;
    configured.validation_snapshot_every = 10;
    configured.validation_n_steps = 111;
    configured.explicit_overrides.validation_snapshot_every = true;
    configured.explicit_overrides.validation_n_steps = true;
    configured.physics_package = "Newtonian";
    auto resolved = galaxy::resolve_scenario(configured);

    const int mk_ok = std::system((std::string("mkdir -p ") + configured.output_dir).c_str());
    (void)mk_ok;
    galaxy::write_run_info(configured.output_dir, resolved.config, resolved.effective_n_steps, 3, resolved.initial_state.n(),
                           "configs/my.local.cfg", "physics/Newtonian/defaults.cfg", &configured, &resolved);
    galaxy::write_resolved_scenario_artifacts(configured.output_dir, resolved);
    const std::string run_info = slurp(configured.output_dir + "/run_info.txt");
    const std::string resolved_txt = slurp(configured.output_dir + "/resolved_scenario.txt");

    CHECK(run_info.find("=== Configured values (post-layering, pre-resolution) ===") != std::string::npos);
    CHECK(run_info.find("=== Effective resolved runtime ===") != std::string::npos);
    CHECK(run_info.find("configured_snapshot_every\t50") != std::string::npos);
    CHECK(run_info.find("configured_validation_snapshot_every\t10") != std::string::npos);
    CHECK(run_info.find("configured_validation_n_steps\t111") != std::string::npos);
    CHECK(run_info.find("effective_snapshot_every\t10") != std::string::npos);
    CHECK(run_info.find("effective_n_steps\t111") != std::string::npos);
    CHECK(run_info.find("effective_initializer_used\tinit_two_body") != std::string::npos);
    CHECK(run_info.find("effective_particle_count\t2") != std::string::npos);
    CHECK(run_info.find("configured_display_distance_unit\tauto") != std::string::npos);
    CHECK(run_info.find("configured_tpfcore_readout_mode\ttensor_radial_projection") != std::string::npos);
    CHECK(resolved_txt.find("effective_snapshot_every\t10") != std::string::npos);
    CHECK(resolved_txt.find("effective_n_steps\t111") != std::string::npos);
  }

  {
    galaxy::Config configured;
    configured.output_dir = "../outputs/test_run_info_bh_orbit";
    configured.simulation_mode = galaxy::SimulationMode::bh_orbit_validation;
    configured.validation_n_steps = 321;
    configured.validation_snapshot_every = 9;
    configured.explicit_overrides.validation_n_steps = true;
    configured.explicit_overrides.validation_snapshot_every = true;
    auto resolved = galaxy::resolve_scenario(configured);
    const int mk_ok = std::system((std::string("mkdir -p ") + configured.output_dir).c_str());
    (void)mk_ok;
    galaxy::write_run_info(configured.output_dir, resolved.config, resolved.effective_n_steps, 5, resolved.initial_state.n(),
                           "configs/my.local.cfg", "physics/Newtonian/defaults.cfg", &configured, &resolved);
    const std::string run_info = slurp(configured.output_dir + "/run_info.txt");
    CHECK(run_info.find("configured_simulation_mode\tbh_orbit_validation") != std::string::npos);
    CHECK(run_info.find("configured_validation_n_steps\t321") != std::string::npos);
    CHECK(run_info.find("configured_validation_snapshot_every\t9") != std::string::npos);
    CHECK(run_info.find("effective_snapshot_every\t9") != std::string::npos);
    CHECK(run_info.find("effective_initializer_used\tinit_two_body_star_around_bh") != std::string::npos);
  }

  {
    galaxy::Config configured;
    configured.output_dir = "../outputs/test_run_info_tpf_4d_static_residual";
    configured.simulation_mode = galaxy::SimulationMode::tpf_4d_static_residual_benchmark;
    configured.physics_package = "TPFCore";
    configured.tpf_dynamics_mode = "direct_tpf";
    configured.tpfcore_enable_provisional_readout = true;
    configured.tpfcore_readout_mode = "derived_tpf_radial_readout";
    auto resolved = galaxy::resolve_scenario(configured);
    const int mk_ok = std::system((std::string("mkdir -p ") + configured.output_dir).c_str());
    (void)mk_ok;
    galaxy::write_run_info(configured.output_dir, resolved.config, resolved.effective_n_steps, 0, 0,
                           "configs/my.local.cfg", "physics/TPFCore/defaults.cfg", &configured, &resolved);
    const std::string run_info = slurp(configured.output_dir + "/run_info.txt");
    CHECK(run_info.find("TPF_4D_static_residual_benchmark (diagnostic-only; no particle integration; no acceleration "
                        "path)") != std::string::npos);
    CHECK(run_info.find("static_4D_field_residual: Xi4/Theta4 full spatial-support diagnostic") != std::string::npos);
    CHECK(run_info.find("none (tpf_4d_static_residual_benchmark uses evaluate_static_configuration_residual_4d; no "
                        "compute_accelerations call)") != std::string::npos);
    CHECK(run_info.find("=== TPF 4D static residual benchmark (diagnostic only; no integrator dynamics) ===") !=
          std::string::npos);
    CHECK(run_info.find("tpf_4d_static_residual_benchmark_artifacts\t"
                        "tpf_4d_static_residual_summary.txt;tpf_4d_static_residual_slice.csv;"
                        "tpf_4d_static_residual_slice_xy.csv;tpf_4d_static_residual_slice_xz.csv;"
                        "tpf_4d_static_residual_slice_yz.csv;tpf_4d_static_residual_sources.csv;"
                        "tpf_4d_static_residual_bins_nearest_source.csv;tpf_4d_static_residual_bins_origin.csv") !=
          std::string::npos);
    CHECK(run_info.find("tpf_4d_static_residual_benchmark_visualization_note\t"
                        "view_plane_renderings_derived_from_full_spatial_support_static_4d_residual_evaluation") !=
          std::string::npos);
    CHECK(run_info.find("tpf_4d_static_residual_benchmark_bins_note\t"
                        "binned_residual_files_provide_quantitative_static_field_evidence_for_the_4d_residual_"
                        "benchmark_and_support_regression_comparison_across_runs") != std::string::npos);
    CHECK(run_info.find("Eq. 10") == std::string::npos);
    CHECK(run_info.find("effective_tpf_dynamics_mode\tnone_static_residual_diagnostic_only") != std::string::npos);
    CHECK(run_info.find("effective_tpf_dynamics_mode\tdirect_tpf") == std::string::npos);
    CHECK(run_info.find("configured_tpf_dynamics_mode\tdirect_tpf") != std::string::npos);
    CHECK(run_info.find("tpf_dynamics_mode_effective_for_this_run\tnone_static_residual_diagnostic_only") !=
          std::string::npos);
    CHECK(run_info.find("tpf_4d_static_residual_benchmark_scope\t"
                        "static_residual_benchmark_exercises_static_Xi4_ordered_Theta4_full_spatial_support_x_y_z_"
                        "stencil") != std::string::npos);
    CHECK(run_info.find("tpf_4d_static_residual_benchmark_unexercised_scope\t"
                        "dynamics_moving_sources_physical_coupling_orbital_behavior_and_DeltaC_closure_are_not_"
                        "exercised_by_this_benchmark") != std::string::npos);
    CHECK(run_info.find("inspection only") == std::string::npos);
    CHECK(run_info.find("does not add physics validation") == std::string::npos);
    CHECK(run_info.find("tpf_dynamics_mode_effective_for_this_run\tdirect_tpf") == std::string::npos);
  }

  {
    const std::string plot_script = slurp("../plot_tpf_4d_static_residual.py");
    CHECK(plot_script.find("not full spatial-support physics domain") == std::string::npos);
    CHECK(plot_script.find("rendered from full spatial-support static 4D residual evaluation") != std::string::npos);
  }

  {
    const std::string tpfcore_readme = slurp("../engine/physics/TPFCore/README.md");
    const std::string tpfcore_plan = slurp("../engine/physics/TPFCore/TPF_4D_FIELD_CORE_PLAN.md");
    CHECK(tpfcore_readme.find("does not add physics validation") == std::string::npos);
    CHECK(tpfcore_readme.find("outside this benchmark") != std::string::npos);
    CHECK(tpfcore_plan.find("inspection-only") == std::string::npos);
    CHECK(tpfcore_plan.find("residual benchmark evidence") != std::string::npos);
  }

  {
    galaxy::Config configured;
    configured.output_dir = "../outputs/test_run_info_tpf_4d_static_residual_defaults_provenance_note";
    configured.simulation_mode = galaxy::SimulationMode::tpf_4d_static_residual_benchmark;
    configured.physics_package = "TPFCore";
    auto resolved = galaxy::resolve_scenario(configured);
    const int mk_ok = std::system((std::string("mkdir -p ") + configured.output_dir).c_str());
    (void)mk_ok;
    galaxy::write_run_info(configured.output_dir, resolved.config, resolved.effective_n_steps, 0, 0,
                           "configs/my.local.cfg", "physics/Newtonian/defaults.cfg", &configured, &resolved);
    const std::string run_info = slurp(configured.output_dir + "/run_info.txt");
    CHECK(run_info.find("package_defaults\tphysics/Newtonian/defaults.cfg") != std::string::npos);
    CHECK(run_info.find("tpf_4d_static_residual_benchmark_package_defaults_note\t"
                        "residual_evaluator_path_is_tpfcore_static_4d_evaluate_static_configuration_residual_4d; "
                        "layered_package_defaults_provenance_line_is_inherited_loader_selection_and_not_the_active_"
                        "residual_evaluator_path") != std::string::npos);
  }

  {
    galaxy::Config configured;
    configured.output_dir = "../outputs/test_run_info_tpf_direct_runtime_mode";
    configured.simulation_mode = galaxy::SimulationMode::tpf_source_field_benchmark;
    configured.physics_package = "TPFCore";
    configured.tpf_dynamics_mode = "direct_tpf";
    auto resolved = galaxy::resolve_scenario(configured);
    const int mk_ok = std::system((std::string("mkdir -p ") + configured.output_dir).c_str());
    (void)mk_ok;
    galaxy::write_run_info(configured.output_dir, resolved.config, resolved.effective_n_steps, 0, 0,
                           "configs/my.local.cfg", "physics/TPFCore/defaults.cfg", &configured, &resolved);
    const std::string run_info = slurp(configured.output_dir + "/run_info.txt");
    CHECK(run_info.find("effective_tpf_dynamics_mode\tdirect_tpf") != std::string::npos);
  }
}
