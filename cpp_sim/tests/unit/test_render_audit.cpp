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
        "TPF_direct_tpf_tensor_principal_part_DeltaC_omitted_"
        "VDSG_off_provisional_readout_off_shunt_off_cooling_off");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf_metrics_tensor_principal_part_DeltaC_omitted_"
        "VDSG_off_provisional_readout_off_shunt_off_cooling_off");
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "TPFCorePackage::compute_direct_tpf_accelerations (field_evaluation -> Theta3D -> principal_Cij -> "
        "tensor_projection; DeltaC omitted; no weak_field_correspondence alpha helper; readout/shunt/cooling "
        "rejected) + accumulate_vdsg_velocity_modifier (continuous zero contribution at tpf_vdsg_coupling == 0)");

  c.tpf_vdsg_coupling = 1e-12;
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "TPF_direct_tpf_tensor_principal_part_DeltaC_omitted_"
        "VDSG_on_provisional_readout_off_shunt_off_cooling_off");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf_metrics_tensor_principal_part_DeltaC_omitted_"
        "VDSG_on_provisional_readout_off_shunt_off_cooling_off");
  CHECK(galaxy::compute_acceleration_code_path(c) ==
        "TPFCorePackage::compute_direct_tpf_accelerations (field_evaluation -> Theta3D -> principal_Cij -> "
        "tensor_projection; DeltaC omitted; no weak_field_correspondence alpha helper; readout/shunt/cooling "
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
