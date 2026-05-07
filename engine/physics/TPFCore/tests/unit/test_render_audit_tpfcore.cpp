#include "config.hpp"
#include "doctest.h"
#include "render_audit.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

TEST_CASE("render audit labels tpf_xi_theta_v1 as active route") {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";

  const std::string dyn = galaxy::compute_active_dynamics_branch(c);
  const std::string met = galaxy::compute_active_metrics_branch(c);
  const std::string acc = galaxy::compute_acceleration_code_path(c);

  CHECK(dyn.find("tpf_xi_theta_v1") != std::string::npos);
  CHECK(met.find("tpf_xi_theta_v1") != std::string::npos);
  CHECK(acc.find("Xi_total") != std::string::npos);
}

TEST_CASE("render audit marks non-v1 modes as unsupported on this branch") {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";

  CHECK(galaxy::compute_active_dynamics_branch(c).find("unsupported_non_v1") != std::string::npos);
  CHECK(galaxy::compute_active_metrics_branch(c).find("unsupported_non_v1") != std::string::npos);
  CHECK(galaxy::compute_acceleration_code_path(c).find("unsupported_non_v1") != std::string::npos);
}


namespace {
std::string slurp_manifest(const std::string& path) {
  std::ifstream in(path.c_str());
  std::ostringstream buf;
  buf << in.rdbuf();
  return buf.str();
}
}  // namespace

TEST_CASE("render manifest v1 law labels do not fall back to legacy_readout_provisional") {
  char dir_template[] = "/tmp/tpfcore_render_manifest_v1_XXXXXX";
  char* out_dir_c = mkdtemp(dir_template);
  REQUIRE(out_dir_c != nullptr);
  const std::string out_dir(out_dir_c);

  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_kernel_coupling = 0.0;

  galaxy::write_render_manifest(out_dir, c, 1, 1, 1, nullptr);
  const std::string txt = slurp_manifest(out_dir + "/render_manifest.txt");
  const std::string json = slurp_manifest(out_dir + "/render_manifest.json");

  CHECK(txt.find("tpf_core_law_mode	tpf_xi_theta_v1") != std::string::npos);
  CHECK(txt.find("legacy_readout_provisional") == std::string::npos);
  CHECK(json.find("\"tpf_core_law_mode\": \"tpf_xi_theta_v1\"") != std::string::npos);
  CHECK(json.find("legacy_readout_provisional") == std::string::npos);

  std::remove((out_dir + "/render_manifest.txt").c_str());
  std::remove((out_dir + "/render_manifest.json").c_str());
  rmdir(out_dir.c_str());
}
