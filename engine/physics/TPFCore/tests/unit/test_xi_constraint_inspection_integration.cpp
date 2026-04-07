#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "types.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::string slurp(const std::string& path) {
  std::ifstream in(path);
  if (!in) return "";
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

bool file_exists(const std::string& path) {
  std::ifstream in(path);
  return static_cast<bool>(in);
}

galaxy::State one_body_state() {
  galaxy::State s;
  s.resize(1);
  s.x[0] = 10.0;
  s.y[0] = 0.0;
  s.vx[0] = 0.0;
  s.vy[0] = 0.0;
  s.mass[0] = 2.0;
  return s;
}

}  // namespace

TEST_CASE("single-source inspect unchanged when xi exterior inspect flag is off") {
  char dir_template[] = "/tmp/tpf_single_source_inspect_off_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.bh_mass = 5.0;
  c.softening = 0.2;
  c.tpfcore_source_softening = 0.2;
  c.tpfcore_probe_radius_min = 1.0;
  c.tpfcore_probe_radius_max = 2.0;
  c.tpfcore_probe_samples = 5;
  c.tpfcore_dump_theta_profile = true;
  c.tpfcore_dump_invariant_profile = true;
  c.tpf_xi_constraint_exterior_inspect = false;

  galaxy::TPFCorePackage pkg;
  pkg.run_single_source_inspect(c, out_dir);

  CHECK(file_exists(std::string(out_dir) + "/theta_profile.csv"));
  CHECK(file_exists(std::string(out_dir) + "/invariant_profile.csv"));
  CHECK(file_exists(std::string(out_dir) + "/field_summary.txt"));
  CHECK_FALSE(file_exists(std::string(out_dir) + "/xi_constraint_exterior_grid.csv"));
  CHECK_FALSE(file_exists(std::string(out_dir) + "/xi_constraint_exterior_summary.txt"));
  CHECK_FALSE(file_exists(std::string(out_dir) + "/xi_constraint_exterior_radial_profile.csv"));
}

TEST_CASE("single-source inspect xi exterior outputs are produced and report free/pinned/masked stats") {
  char dir_template[] = "/tmp/tpf_single_source_inspect_on_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.bh_mass = 5.0;
  c.softening = 0.2;
  c.tpfcore_source_softening = 0.2;
  c.tpfcore_probe_radius_min = 1.0;
  c.tpfcore_probe_radius_max = 2.0;
  c.tpfcore_probe_samples = 5;
  c.tpf_xi_constraint_exterior_inspect = true;
  c.tpf_xi_constraint_grid_n = 33;
  c.tpf_xi_constraint_half_extent = 3.0;
  c.tpf_xi_constraint_inner_radius = 0.5;
  c.tpf_xi_constraint_max_iters = 20;
  c.tpf_xi_constraint_tol = 1e-6;

  galaxy::TPFCorePackage pkg;
  pkg.run_single_source_inspect(c, out_dir);

  const std::string grid_path = std::string(out_dir) + "/xi_constraint_exterior_grid.csv";
  const std::string summary_path = std::string(out_dir) + "/xi_constraint_exterior_summary.txt";
  const std::string radial_path = std::string(out_dir) + "/xi_constraint_exterior_radial_profile.csv";

  CHECK(file_exists(grid_path));
  CHECK(file_exists(summary_path));
  CHECK(file_exists(radial_path));

  const std::string grid_text = slurp(grid_path);
  CHECK(grid_text.find("is_masked") != std::string::npos);
  CHECK(grid_text.find("is_boundary_pinned") != std::string::npos);
  CHECK(grid_text.find("is_free_cell") != std::string::npos);

  const std::string summary_text = slurp(summary_path);
  CHECK(summary_text.find("number of masked cells") != std::string::npos);
  CHECK(summary_text.find("number of pinned cells") != std::string::npos);
  CHECK(summary_text.find("number of free cells") != std::string::npos);
  CHECK(summary_text.find("initial max residual norm over free cells") != std::string::npos);
  CHECK(summary_text.find("final max residual norm over free cells") != std::string::npos);
  CHECK(summary_text.find("ansatz used as boundary data and initial guess") != std::string::npos);
  CHECK(summary_text.find("not full Eq. (10) validation") != std::string::npos);
  CHECK(summary_text.find("no DeltaC included") != std::string::npos);
  CHECK(summary_text.find("single-source only") != std::string::npos);
}

TEST_CASE("xi exterior inspection flag does not change dynamics routing or accelerations") {
  galaxy::Config base;
  base.tpfcore_enable_provisional_readout = true;
  base.tpfcore_readout_mode = "tensor_radial_projection";
  base.tpf_dynamics_mode = "legacy_readout";
  base.tpf_vdsg_coupling = 0.0;

  galaxy::Config with_flag = base;
  with_flag.tpf_xi_constraint_exterior_inspect = true;
  with_flag.tpf_xi_constraint_grid_n = 33;

  galaxy::TPFCorePackage p0;
  galaxy::TPFCorePackage p1;
  p0.init_from_config(base);
  p1.init_from_config(with_flag);

  galaxy::State s = one_body_state();
  std::vector<double> ax0, ay0, ax1, ay1;
  p0.compute_accelerations(s, 100.0, 1.0, false, ax0, ay0);
  p1.compute_accelerations(s, 100.0, 1.0, false, ax1, ay1);

  REQUIRE(ax0.size() == 1);
  REQUIRE(ax1.size() == 1);
  CHECK(ax1[0] == doctest::Approx(ax0[0]).epsilon(1e-12));
  CHECK(ay1[0] == doctest::Approx(ay0[0]).epsilon(1e-12));
}
