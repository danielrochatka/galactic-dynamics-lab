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
        "tpf_dynamics_mode=direct_tpf; Theta/I/kappa; DeltaC omitted; Xi-directed readout; vdsg_coupling=0.000000e+00");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf metrics; Theta/I/kappa; DeltaC omitted; vdsg_coupling=0.000000e+00");
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "TPFCorePackage::compute_direct_tpf_accelerations (principal-part implementation: field_evaluation -> Theta3D -> principal_Cij -> "
        "Xi_directed_tensor_readout; Theta/I/kappa baseline; DeltaC omitted in current implementation scope; readout/shunt/cooling "
        "rejected) + accumulate_vdsg_velocity_modifier (continuous zero contribution at tpf_vdsg_coupling == 0)");

  c.tpf_vdsg_coupling = 1e-12;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "tpf_dynamics_mode=direct_tpf; Theta/I/kappa; DeltaC omitted; Xi-directed readout; vdsg_coupling=1.000000e-12");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf metrics; Theta/I/kappa; DeltaC omitted; vdsg_coupling=1.000000e-12");
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "TPFCorePackage::compute_direct_tpf_accelerations (principal-part implementation: field_evaluation -> Theta3D -> principal_Cij -> "
        "Xi_directed_tensor_readout; Theta/I/kappa baseline; DeltaC omitted in current implementation scope; readout/shunt/cooling "
        "rejected) + accumulate_vdsg_velocity_modifier (optional additive VDSG extension)");
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
