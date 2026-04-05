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
  CHECK(galaxy::compute_active_dynamics_branch(c) == "TPF_PROVISIONAL_legacy_readout:derived_tpf_radial_readout");
  c.tpf_vdsg_coupling = 1e-20;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "TPF_PROVISIONAL_legacy_readout_plus_EXPLORATORY_VDSG:derived_tpf_radial_readout");
}

TEST_CASE("compute_active_dynamics_branch: TPF direct") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_vdsg_coupling = 0.0;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "TPF_direct_tpf_tensor_principal_part_DeltaC_omitted_"
        "VDSG_off_provisional_readout_off_shunt_off_cooling_off");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf_metrics_tensor_principal_part_DeltaC_omitted_"
        "VDSG_off_provisional_readout_off_shunt_off_cooling_off");
  c.tpf_vdsg_coupling = 1e-20;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "TPF_direct_tpf_tensor_principal_part_DeltaC_omitted_"
        "VDSG_on_provisional_readout_off_shunt_off_cooling_off");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf_metrics_tensor_principal_part_DeltaC_omitted_"
        "VDSG_on_provisional_readout_off_shunt_off_cooling_off");
}

TEST_CASE("compute_active_dynamics_branch: v11 weak-field truncation dynamics") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "v11_weak_field_truncation";
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "TPF_v11_weak_field_truncation_dynamics_limited_static_quasistatic");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "v11_weak_field_truncation_dynamics_metrics_limited_scope");
  CHECK(galaxy::compute_acceleration_code_path(c).find("Eq.42-44") != std::string::npos);
}

TEST_CASE("compute_active_dynamics_branch: legacy weak_field_correspondence string aliases to truncation labels") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "weak_field_correspondence";
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "TPF_v11_weak_field_truncation_dynamics_limited_static_quasistatic");
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
        " + accumulate_vdsg_velocity_modifier (global |a| shunt OFF — clean readout+VDSG path without velocity cap)");
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
        " + accumulate_vdsg_velocity_modifier (global |a| shunt OFF — clean readout+VDSG path without velocity cap)");
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

TEST_CASE("compute_acceleration_code_path: direct_tpf canonical low-order truncation route") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_vdsg_coupling = 0.0;
  CHECK(galaxy::compute_acceleration_code_path(c).find("field_evaluation -> Theta3D -> principal_Cij -> tensor_projection") !=
        std::string::npos);
  CHECK(galaxy::compute_acceleration_code_path(c).find("compute_v11_weak_field_truncation_accelerations") ==
        std::string::npos);
  c.tpf_vdsg_coupling = 1e-10;
  CHECK(galaxy::compute_acceleration_code_path(c).find("optional additive VDSG extension") != std::string::npos);
}
