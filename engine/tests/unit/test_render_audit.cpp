#include "config.hpp"
#include "doctest.h"
#include "render_audit.hpp"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using galaxy::Config;

namespace {

std::string slurp_file(const std::string& path) {
  std::ifstream in(path.c_str());
  std::ostringstream buf;
  buf << in.rdbuf();
  return buf.str();
}

}  // namespace

TEST_CASE("compute_active_dynamics_branch: Newtonian") {
  Config c;
  c.physics_package = "Newtonian";
  CHECK(galaxy::compute_active_dynamics_branch(c) == "Newtonian_pairwise_G_SI");
}

TEST_CASE("compute_active_dynamics_branch: unknown package is tagged non-TPFCore") {
  Config c;
  c.physics_package = "CustomLabModel";
  CHECK(galaxy::compute_active_dynamics_branch(c) == "CustomLabModel (non-TPFCore)");
  CHECK(galaxy::compute_active_metrics_branch(c) == "unknown");
  CHECK(galaxy::compute_acceleration_code_path(c) == "unknown_package");
}

TEST_CASE("compute_active_dynamics_branch: direct_tpf reports VDSG extension status") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_vdsg_coupling = 0.0;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_runtime_path_tier=paper_facing; tpf_dynamics_mode=direct_tpf; Theta/I/kappa; DeltaC omitted; Xi-directed readout; vdsg_coupling=0.000000e+00");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf metrics; Theta/I/kappa; DeltaC omitted; vdsg_coupling=0.000000e+00");
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "TPFCorePackage::compute_direct_tpf_accelerations (principal-part implementation: field_evaluation -> Theta3D -> principal_Cij -> "
        "Xi_directed_tensor_readout; Theta/I/kappa baseline; DeltaC omitted in current implementation scope; readout/shunt/cooling "
        "rejected) + accumulate_vdsg_velocity_modifier (continuous zero contribution at tpf_vdsg_coupling == 0)");

  c.tpf_vdsg_coupling = 1e-12;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_runtime_path_tier=paper_facing; tpf_dynamics_mode=direct_tpf; Theta/I/kappa; DeltaC omitted; Xi-directed readout; vdsg_coupling=1.000000e-12");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf metrics; Theta/I/kappa; DeltaC omitted; vdsg_coupling=1.000000e-12");
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "TPFCorePackage::compute_direct_tpf_accelerations (principal-part implementation: field_evaluation -> Theta3D -> principal_Cij -> "
        "Xi_directed_tensor_readout; Theta/I/kappa baseline; DeltaC omitted in current implementation scope; readout/shunt/cooling "
        "rejected) + accumulate_vdsg_velocity_modifier (optional additive VDSG extension)");
}

TEST_CASE("compute_active_dynamics_branch: xi_kernel_deformed reports Xi runtime route semantics") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_motion_readout_scale = 1.0e-12;
  c.tpf_4d_xi_kernel_mode = "metric_velocity";
  c.tpf_4d_xi_kernel_coupling = 2.0;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_temporal_mode = "off";
  const std::string dyn = galaxy::compute_active_dynamics_branch(c);
  CHECK(dyn.find("tpf_dynamics_mode=xi_kernel_deformed") != std::string::npos);
  CHECK(dyn.find("a=-K_xi*Xi_eff_spatial") != std::string::npos);
  CHECK(dyn.find("additive_vdsg=off") != std::string::npos);
  CHECK(dyn.find("principal_c=off") != std::string::npos);
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "xi_kernel_deformed metrics; Xi_eff readout a=-K_xi*Xi_eff_spatial");
  const std::string code_path = galaxy::compute_acceleration_code_path(c);
  CHECK(code_path.find("compute_xi_kernel_deformed_accelerations") != std::string::npos);
  CHECK(code_path.find("no additive VDSG helper") != std::string::npos);
}

TEST_CASE("write_render_manifest TXT tpf_extension_mode semantics align with JSON semantics") {
  char dir_template[] = "/tmp/render_audit_txt_ext_mode_XXXXXX";
  char* out_dir_c = mkdtemp(dir_template);
  REQUIRE(out_dir_c != nullptr);
  const std::string out_dir(out_dir_c);

  Config c;
  c.physics_package = "TPFCore";
  c.simulation_mode = galaxy::SimulationMode::galaxy;
  c.tpf_dynamics_mode = "direct_tpf";
  c.tpf_vdsg_coupling = 0.0;
  galaxy::write_render_manifest(out_dir, c, 1, 1, 8, nullptr);
  std::string txt = slurp_file(out_dir + "/render_manifest.txt");
  CHECK(txt.find("tpf_extension_mode\tnone\n") != std::string::npos);

  c.tpf_vdsg_coupling = 1e-12;
  galaxy::write_render_manifest(out_dir, c, 1, 1, 8, nullptr);
  txt = slurp_file(out_dir + "/render_manifest.txt");
  CHECK(txt.find("tpf_extension_mode\tvdsg\n") != std::string::npos);

  c.tpf_dynamics_mode = "v11_weak_field_truncation";
  c.tpf_vdsg_coupling = 1e-12;
  galaxy::write_render_manifest(out_dir, c, 1, 1, 8, nullptr);
  txt = slurp_file(out_dir + "/render_manifest.txt");
  CHECK(txt.find("tpf_extension_mode\tnone_vdsg_rejected\n") != std::string::npos);

  std::remove((out_dir + "/render_manifest.txt").c_str());
  std::remove((out_dir + "/render_manifest.json").c_str());
  rmdir(out_dir.c_str());
}

TEST_CASE("write_render_manifest: xi_kernel_deformed reports xi law metadata (not legacy provisional)") {
  char dir_template[] = "/tmp/render_audit_xi_kernel_XXXXXX";
  char* out_dir_c = mkdtemp(dir_template);
  REQUIRE(out_dir_c != nullptr);
  const std::string out_dir(out_dir_c);

  Config c;
  c.physics_package = "TPFCore";
  c.simulation_mode = galaxy::SimulationMode::galaxy;
  c.tpf_dynamics_mode = "xi_kernel_deformed";
  c.tpf_4d_xi_kernel_mode = "scalar_beta";
  c.tpf_4d_xi_kernel_coupling = 3.0e-3;
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_temporal_mode = "off";
  galaxy::write_render_manifest(out_dir, c, 1, 1, 8, nullptr);

  const std::string txt = slurp_file(out_dir + "/render_manifest.txt");
  const std::string json = slurp_file(out_dir + "/render_manifest.json");
  CHECK(txt.find("tpf_core_law_mode\txi_kernel_deformed") != std::string::npos);
  CHECK(txt.find("K_xi\ttpf_4d_xi_motion_readout_scale") != std::string::npos);
  CHECK(txt.find("acceleration_formula\ta=-K_xi*Xi_eff_spatial") != std::string::npos);
  CHECK(txt.find("xi_kernel_factor_mode\tbeta_power") != std::string::npos);
  CHECK(txt.find("xi_temporal_mode\toff") != std::string::npos);
  CHECK(txt.find("xi_kernel_metric_min\t") != std::string::npos);
  CHECK(txt.find("xi_kernel_metric_max\t") != std::string::npos);
  CHECK(txt.find("tpf_extension_mode\t") == std::string::npos);
  CHECK(txt.find("tpf_stabilizer_mode\t") == std::string::npos);
  CHECK(txt.find("legacy_readout_provisional") == std::string::npos);
  CHECK(json.find("\"tpf_core_law_mode\": \"xi_kernel_deformed\"") != std::string::npos);
  CHECK(json.find("\"K_xi\": \"tpf_4d_xi_motion_readout_scale\"") != std::string::npos);
  CHECK(json.find("\"acceleration_formula\": \"a=-K_xi*Xi_eff_spatial\"") != std::string::npos);
  CHECK(json.find("\"xi_kernel_factor_mode\": \"beta_power\"") != std::string::npos);
  CHECK(json.find("\"xi_temporal_mode\": \"off\"") != std::string::npos);
  CHECK(json.find("\"tpf_extension_mode\"") == std::string::npos);
  CHECK(json.find("\"tpf_stabilizer_mode\"") == std::string::npos);

  std::remove((out_dir + "/render_manifest.txt").c_str());
  std::remove((out_dir + "/render_manifest.json").c_str());
  rmdir(out_dir.c_str());
}

TEST_CASE("write_render_manifest: tpf_4d_static_residual_benchmark uses diagnostic-only audit labels") {
  char dir_template[] = "/tmp/render_audit_static4d_XXXXXX";
  char* out_dir_c = mkdtemp(dir_template);
  REQUIRE(out_dir_c != nullptr);
  const std::string out_dir(out_dir_c);

  Config c;
  c.physics_package = "TPFCore";
  c.simulation_mode = galaxy::SimulationMode::tpf_4d_static_residual_benchmark;
  c.tpf_dynamics_mode = "direct_tpf";
  galaxy::write_render_manifest(out_dir, c, 0, 0, 0, nullptr);

  const std::string txt = slurp_file(out_dir + "/render_manifest.txt");
  const std::string json = slurp_file(out_dir + "/render_manifest.json");
  CHECK(txt.find("active_dynamics_branch\tTPF_4D_static_residual_benchmark (diagnostic-only; no particle integration; "
                 "no acceleration path)") != std::string::npos);
  CHECK(txt.find("active_metrics_branch\tstatic_4D_field_residual: Xi4/Theta4 full spatial-support diagnostic") !=
        std::string::npos);
  CHECK(txt.find("acceleration_code_path\tnone (tpf_4d_static_residual_benchmark uses "
                 "evaluate_static_configuration_residual_4d; no compute_accelerations call)") != std::string::npos);
  CHECK(json.find("\"active_dynamics_branch\": \"TPF_4D_static_residual_benchmark (diagnostic-only; no particle "
                  "integration; no acceleration path)\"") != std::string::npos);
  CHECK(json.find("\"active_metrics_branch\": \"static_4D_field_residual: Xi4/Theta4 full spatial-support diagnostic\"") !=
        std::string::npos);
  CHECK(json.find("\"acceleration_code_path\": \"none (tpf_4d_static_residual_benchmark uses "
                  "evaluate_static_configuration_residual_4d; no compute_accelerations call)\"") != std::string::npos);
  CHECK(json.find("compute_direct_tpf_accelerations") == std::string::npos);

  std::remove((out_dir + "/render_manifest.txt").c_str());
  std::remove((out_dir + "/render_manifest.json").c_str());
  rmdir(out_dir.c_str());
}

TEST_CASE("write_render_manifest: tpf_4d_static_motion_readout_benchmark uses benchmark-only audit labels") {
  char dir_template[] = "/tmp/render_audit_static4d_motion_XXXXXX";
  char* out_dir_c = mkdtemp(dir_template);
  REQUIRE(out_dir_c != nullptr);
  const std::string out_dir(out_dir_c);

  Config c;
  c.physics_package = "TPFCore";
  c.simulation_mode = galaxy::SimulationMode::tpf_4d_static_motion_readout_benchmark;
  c.tpf_dynamics_mode = "direct_tpf";
  galaxy::write_render_manifest(out_dir, c, 0, 0, 0, nullptr);

  const std::string txt = slurp_file(out_dir + "/render_manifest.txt");
  const std::string json = slurp_file(out_dir + "/render_manifest.json");
  CHECK(txt.find("active_dynamics_branch\tTPF_4D_static_motion_readout_benchmark") != std::string::npos);
  CHECK(txt.find("active_metrics_branch\tstatic_4D_motion_readout: GravityStaticMotionReadout_v1 over probe grid") !=
        std::string::npos);
  CHECK(txt.find("acceleration_code_path\tnone (tpf_4d_static_motion_readout_benchmark uses "
                 "evaluate_static_sources_field_4d + GravityStaticMotionReadout_v1; no compute_accelerations call)") !=
        std::string::npos);
  CHECK(json.find("\"active_dynamics_branch\": \"TPF_4D_static_motion_readout_benchmark") != std::string::npos);
  CHECK(json.find("\"active_metrics_branch\": \"static_4D_motion_readout: GravityStaticMotionReadout_v1 over probe "
                  "grid\"") != std::string::npos);
  CHECK(json.find("\"acceleration_code_path\": \"none (tpf_4d_static_motion_readout_benchmark uses "
                  "evaluate_static_sources_field_4d + GravityStaticMotionReadout_v1; no compute_accelerations call)\"") !=
        std::string::npos);

  std::remove((out_dir + "/render_manifest.txt").c_str());
  std::remove((out_dir + "/render_manifest.json").c_str());
  rmdir(out_dir.c_str());
}

TEST_CASE("write_render_manifest: tpf_4d_xi_motion_probe_benchmark uses Xi-direct benchmark-only audit labels") {
  char dir_template[] = "/tmp/render_audit_static4d_ximotion_XXXXXX";
  char* out_dir_c = mkdtemp(dir_template);
  REQUIRE(out_dir_c != nullptr);
  const std::string out_dir(out_dir_c);

  Config c;
  c.physics_package = "TPFCore";
  c.simulation_mode = galaxy::SimulationMode::tpf_4d_xi_motion_probe_benchmark;
  c.tpf_dynamics_mode = "direct_tpf";
  galaxy::write_render_manifest(out_dir, c, 0, 0, 0, nullptr);

  const std::string txt = slurp_file(out_dir + "/render_manifest.txt");
  const std::string json = slurp_file(out_dir + "/render_manifest.json");
  CHECK(txt.find("active_dynamics_branch\tTPF_4D_xi_motion_probe_benchmark") != std::string::npos);
  CHECK(txt.find("active_metrics_branch\tdynamic_4D_xi_motion_readout: GravityXiMotionReadout_v1 over moving probes") !=
        std::string::npos);
  CHECK(txt.find("acceleration_code_path\tnone (tpf_4d_xi_motion_probe_benchmark uses "
                 "evaluate_static_sources_field_4d + GravityXiMotionReadout_v1; no compute_accelerations call)") !=
        std::string::npos);
  CHECK(json.find("\"active_dynamics_branch\": \"TPF_4D_xi_motion_probe_benchmark") != std::string::npos);
  CHECK(json.find("\"active_metrics_branch\": \"dynamic_4D_xi_motion_readout: GravityXiMotionReadout_v1 over moving "
                  "probes\"") != std::string::npos);
  CHECK(json.find("\"acceleration_code_path\": \"none (tpf_4d_xi_motion_probe_benchmark uses "
                  "evaluate_static_sources_field_4d + GravityXiMotionReadout_v1; no compute_accelerations call)\"") !=
        std::string::npos);

  std::remove((out_dir + "/render_manifest.txt").c_str());
  std::remove((out_dir + "/render_manifest.json").c_str());
  rmdir(out_dir.c_str());
}
