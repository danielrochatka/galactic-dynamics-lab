#include "config.hpp"
#include "doctest.h"
#include "render_audit.hpp"

using galaxy::Config;

TEST_CASE("compute_active_dynamics_branch: TPF legacy_readout labels VDSG in branch name") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "legacy_readout";
  c.tpfcore_enable_provisional_readout = true;
  c.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c.tpf_vdsg_coupling = 0.0;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_runtime_path_tier=deprecated_legacy; tpf_dynamics_mode=legacy_readout; provisional readout; mode=derived_tpf_radial_readout; vdsg_coupling=0.000000e+00");
  c.tpf_vdsg_coupling = 1e-20;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_runtime_path_tier=deprecated_legacy; tpf_dynamics_mode=legacy_readout; provisional readout; mode=derived_tpf_radial_readout; vdsg_coupling=1.000000e-20");
}

TEST_CASE("compute_active_dynamics_branch: TPF direct") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_vdsg_coupling = 0.0;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_runtime_path_tier=paper_facing; tpf_dynamics_mode=direct_tpf; Theta/I/kappa; DeltaC omitted; Xi-directed readout; vdsg_coupling=0.000000e+00");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf metrics; Theta/I/kappa; DeltaC omitted; vdsg_coupling=0.000000e+00");
  c.tpf_vdsg_coupling = 1e-20;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_runtime_path_tier=paper_facing; tpf_dynamics_mode=direct_tpf; Theta/I/kappa; DeltaC omitted; Xi-directed readout; vdsg_coupling=1.000000e-20");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf metrics; Theta/I/kappa; DeltaC omitted; vdsg_coupling=1.000000e-20");
}

TEST_CASE("compute_active_dynamics_branch: correspondence dynamics labels") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "geodesic_correspondence";
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_runtime_path_tier=paper_facing; tpf_dynamics_mode=geodesic_correspondence; rho_Xi -> Poisson analytic -> geodesic");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "geodesic correspondence metrics; rho_Xi -> Poisson analytic -> geodesic");
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "geodesic_correspondence: rho_Xi -> Poisson analytic -> geodesic");

  c.tpf_dynamics_mode = "v11_weak_field_truncation";
  CHECK(galaxy::compute_active_dynamics_branch(c).find("deprecated_legacy") != std::string::npos);
  CHECK(galaxy::compute_active_metrics_branch(c).find("legacy free-parameter") != std::string::npos);
}

TEST_CASE("compute_active_metrics_branch: metrics vs dynamics when VDSG on") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "legacy_readout";
  c.tpfcore_enable_provisional_readout = true;
  c.tpfcore_readout_mode = "derived_tpf_radial_readout";
  c.tpf_vdsg_coupling = 1e-15;
  CHECK(galaxy::compute_active_metrics_branch(c) == "tpfcore_readout:derived_tpf_radial_readout");
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "TPFCorePackage::compute_provisional_readout_acceleration + derived_tpf_radial_profile"
        " + accumulate_vdsg_velocity_modifier (global |a| shunt OFF)");
}

TEST_CASE("compute_acceleration_code_path: same pipeline string when vdsg coupling is zero") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "legacy_readout";
  c.tpfcore_enable_provisional_readout = true;
  c.tpfcore_readout_mode = "tensor_radial_projection";
  c.tpf_vdsg_coupling = 0.0;
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "TPFCorePackage::compute_provisional_readout_acceleration (tensor_radial_projection)"
        " + accumulate_vdsg_velocity_modifier (global |a| shunt OFF)");
}

TEST_CASE("compute_acceleration_code_path: includes shunt when explicitly enabled") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "legacy_readout";
  c.tpfcore_enable_provisional_readout = true;
  c.tpfcore_readout_mode = "tensor_radial_projection";
  c.tpf_global_accel_shunt_enable = true;
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "TPFCorePackage::compute_provisional_readout_acceleration (tensor_radial_projection)"
        " + accumulate_vdsg_velocity_modifier + apply_global_accel_magnitude_shunt (when tpf_global_accel_shunt_enable)");
}

TEST_CASE("compute_active_dynamics_branch: v11 weak-field correspondence audit mode") {
  Config c;
  c.simulation_mode = galaxy::SimulationMode::tpf_v11_weak_field_correspondence;
  c.physics_package = "TPFCore";
  CHECK(galaxy::compute_active_dynamics_branch(c).find("TPF_v11_weak_field_correspondence_audit") == 0);
  CHECK(galaxy::compute_active_metrics_branch(c).find("v11_paper_tensors") == 0);
  CHECK(galaxy::compute_acceleration_code_path(c).find("audit-only") != std::string::npos);
}

TEST_CASE("compute_active_* branch labels: tpf_4d_static_residual_benchmark is diagnostic-only") {
  Config c;
  c.simulation_mode = galaxy::SimulationMode::tpf_4d_static_residual_benchmark;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";

  const std::string dynamics = galaxy::compute_active_dynamics_branch(c);
  const std::string metrics = galaxy::compute_active_metrics_branch(c);
  const std::string accel = galaxy::compute_acceleration_code_path(c);

  CHECK(dynamics.find("diagnostic-only") != std::string::npos);
  CHECK(dynamics.find("no particle integration") != std::string::npos);
  CHECK(dynamics.find("no acceleration path") != std::string::npos);
  CHECK(dynamics.find("direct_tpf") == std::string::npos);
  CHECK(dynamics.find("compute_direct_tpf_accelerations") == std::string::npos);

  CHECK(metrics.find("static_4D_field_residual") != std::string::npos);
  CHECK(metrics.find("Xi4") != std::string::npos);
  CHECK(metrics.find("Theta4") != std::string::npos);

  CHECK(accel.find("none") != std::string::npos);
  CHECK(accel.find("evaluate_static_configuration_residual_4d") != std::string::npos);
  CHECK(accel.find("no compute_accelerations call") != std::string::npos);
  CHECK(accel.find("compute_direct_tpf_accelerations") == std::string::npos);
}

TEST_CASE("compute_active_* branch labels: tpf_4d_static_motion_readout_benchmark is benchmark-only") {
  Config c;
  c.simulation_mode = galaxy::SimulationMode::tpf_4d_static_motion_readout_benchmark;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";

  const std::string dynamics = galaxy::compute_active_dynamics_branch(c);
  const std::string metrics = galaxy::compute_active_metrics_branch(c);
  const std::string accel = galaxy::compute_acceleration_code_path(c);

  CHECK(dynamics.find("TPF_4D_static_motion_readout_benchmark") != std::string::npos);
  CHECK(metrics.find("GravityStaticMotionReadout_v1") != std::string::npos);
  CHECK(accel.find("evaluate_static_sources_field_4d") != std::string::npos);
  CHECK(accel.find("no compute_accelerations call") != std::string::npos);
  CHECK(accel.find("compute_direct_tpf_accelerations") == std::string::npos);
}

TEST_CASE("compute_active_* branch labels: tpf_4d_xi_motion_probe_benchmark is isolated Xi-direct benchmark") {
  Config c;
  c.simulation_mode = galaxy::SimulationMode::tpf_4d_xi_motion_probe_benchmark;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";

  const std::string dynamics = galaxy::compute_active_dynamics_branch(c);
  const std::string metrics = galaxy::compute_active_metrics_branch(c);
  const std::string accel = galaxy::compute_acceleration_code_path(c);

  CHECK(dynamics.find("TPF_4D_xi_motion_probe_benchmark") != std::string::npos);
  CHECK(metrics.find("GravityXiMotionReadout_v1") != std::string::npos);
  CHECK(accel.find("evaluate_static_sources_field_4d") != std::string::npos);
  CHECK(accel.find("GravityXiMotionReadout_v1") != std::string::npos);
  CHECK(accel.find("no compute_accelerations call") != std::string::npos);
  CHECK(accel.find("compute_direct_tpf_accelerations") == std::string::npos);
}

TEST_CASE("compute_acceleration_code_path: direct_tpf tensor principal-part route") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_vdsg_coupling = 0.0;
  CHECK(galaxy::compute_acceleration_code_path(c).find("field_evaluation -> Theta3D -> principal_Cij -> Xi_directed_tensor_readout") !=
        std::string::npos);
  CHECK(galaxy::compute_acceleration_code_path(c).find("compute_v11_weak_field_truncation_accelerations") ==
        std::string::npos);
  c.tpf_vdsg_coupling = 1e-10;
  CHECK(galaxy::compute_acceleration_code_path(c).find("optional additive VDSG extension") != std::string::npos);
}
