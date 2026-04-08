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
        "tpf_dynamics_mode=legacy_readout; provisional readout; mode=derived_tpf_radial_readout; vdsg_coupling=0.000000e+00");
  c.tpf_vdsg_coupling = 1e-20;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_dynamics_mode=legacy_readout; provisional readout; mode=derived_tpf_radial_readout; vdsg_coupling=1.000000e-20");
}

TEST_CASE("compute_active_dynamics_branch: TPF direct") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_vdsg_coupling = 0.0;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_dynamics_mode=direct_tpf; source=provisional_ansatz; Theta/I/kappa; DeltaC omitted; Xi-directed readout; vdsg_coupling=0.000000e+00");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf metrics; source=provisional_ansatz; DeltaC omitted; vdsg_coupling=0.000000e+00");
  c.tpf_vdsg_coupling = 1e-20;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_dynamics_mode=direct_tpf; source=provisional_ansatz; Theta/I/kappa; DeltaC omitted; Xi-directed readout; vdsg_coupling=1.000000e-20");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf metrics; source=provisional_ansatz; DeltaC omitted; vdsg_coupling=1.000000e-20");
}

TEST_CASE("compute_active_dynamics_branch: v11 weak-field truncation dynamics") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "v11_weak_field_truncation";
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_dynamics_mode=v11_weak_field_truncation; correspondence implementation; alpha_si");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "v11 correspondence metrics; alpha_si");
  CHECK(galaxy::compute_acceleration_code_path(c).find("Eq.42-44") != std::string::npos);
}

TEST_CASE("compute_active_dynamics_branch: deprecated weak_field_correspondence alias resolves to correspondence helper labels") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "weak_field_correspondence";
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_dynamics_mode=v11_weak_field_truncation; correspondence implementation; alpha_si");
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

TEST_CASE("compute_active_dynamics_branch: direct_tpf experimental xi_constraint runtime source is labeled") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_direct_field_source = "xi_constraint_exterior_single_source";
  const std::string dyn = galaxy::compute_active_dynamics_branch(c);
  CHECK(dyn.find("experimental planar single-source runtime Xi solve") != std::string::npos);
  CHECK(dyn.find("not full Eq. (10)") != std::string::npos);
}
