#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "types.hpp"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <limits>
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

double parse_summary_scalar(const std::string& summary_text, const std::string& key_prefix) {
  std::istringstream in(summary_text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.find(key_prefix) == 0) {
      const std::size_t p = line.find(':');
      if (p == std::string::npos) return std::numeric_limits<double>::quiet_NaN();
      return std::atof(line.substr(p + 1).c_str());
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  for (std::size_t i = 0; i < line.size(); ++i) {
    if (line[i] == ',') {
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(line[i]);
    }
  }
  out.push_back(cur);
  return out;
}

double parse_csv_double_or_nan(const std::string& s) {
  if (s.empty() || s == "nan" || s == "NaN" || s == "NAN") {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::atof(s.c_str());
}

int find_csv_col(const std::vector<std::string>& header, const std::string& key) {
  for (std::size_t i = 0; i < header.size(); ++i) {
    if (header[i] == key) return static_cast<int>(i);
  }
  return -1;
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

TEST_CASE("source-field benchmark bonded_pair supports explicit unequal masses with fallback") {
  char dir_template[] = "/tmp/tpf_source_field_benchmark_pair_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "bonded_pair";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_source_benchmark_mass1 = 7.5e11;
  c.tpf_source_benchmark_mass2 = 2.5e11;
  c.tpf_source_benchmark_separation = 8.0;
  c.tpf_source_benchmark_orientation_deg = 90.0;
  c.tpf_source_probe_grid_half_extent = 20.0;
  c.tpf_source_probe_grid_n = 9;
  c.softening = 0.1;

  galaxy::TPFCorePackage pkg;
  pkg.run_source_field_benchmark(c, out_dir);

  std::string csv = slurp(std::string(out_dir) + "/tpf_source_field_probe_grid.csv");
  CHECK(csv.find("source_mass1,source_mass2") != std::string::npos);
  CHECK(csv.find("bonded_pair_centered_origin_explicit_mass1_mass2") != std::string::npos);
  CHECK(csv.find(",7.500000e+11,2.500000e+11,") != std::string::npos);

  c.tpf_source_benchmark_mass1 = 0.0;
  c.tpf_source_benchmark_mass2 = 0.0;
  pkg.run_source_field_benchmark(c, out_dir);
  csv = slurp(std::string(out_dir) + "/tpf_source_field_probe_grid.csv");
  CHECK(csv.find("bonded_pair_centered_origin_equal_mass") != std::string::npos);
  CHECK(csv.find(",5.000000e+11,5.000000e+11,") != std::string::npos);
}

TEST_CASE("source-field benchmark writes residual diagnostic columns and summary file") {
  char dir_template[] = "/tmp/tpf_source_field_benchmark_residual_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_source_probe_grid_half_extent = 20.0;
  c.tpf_source_probe_grid_n = 9;
  c.softening = 0.1;
  c.tpf_source_residual_exclusion_radius = 0.25;

  galaxy::TPFCorePackage pkg;
  pkg.run_source_field_benchmark(c, out_dir);

  const std::string csv_path = std::string(out_dir) + "/tpf_source_field_probe_grid.csv";
  const std::string summary_path = std::string(out_dir) + "/tpf_source_field_residual_summary.txt";

  CHECK(file_exists(csv_path));
  CHECK(file_exists(summary_path));

  const std::string csv = slurp(csv_path);
  CHECK(csv.find("residual_x") != std::string::npos);
  CHECK(csv.find("residual_y") != std::string::npos);
  CHECK(csv.find("residual_norm") != std::string::npos);
  CHECK(csv.find("residual_norm_over_theta_norm") != std::string::npos);
  CHECK(csv.find("excluded_boundary") != std::string::npos);
  CHECK(csv.find("excluded_near_source") != std::string::npos);
  CHECK(csv.find(",1,0\n") != std::string::npos);
  CHECK(csv.find(",0,1\n") != std::string::npos);

  const std::string summary = slurp(summary_path);
  CHECK(summary.find("total grid cells") != std::string::npos);
  CHECK(summary.find("interior free cells used") != std::string::npos);
  CHECK(summary.find("excluded boundary count") != std::string::npos);
  CHECK(summary.find("excluded near-source count") != std::string::npos);
}

TEST_CASE("source-field benchmark near-source exclusion radius flags cells for bonded pair") {
  char dir_template[] = "/tmp/tpf_source_field_benchmark_residual_pair_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "bonded_pair";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_source_benchmark_separation = 8.0;
  c.tpf_source_benchmark_orientation_deg = 90.0;
  c.tpf_source_probe_grid_half_extent = 20.0;
  c.tpf_source_probe_grid_n = 9;
  c.softening = 0.1;
  c.tpf_source_residual_exclusion_radius = 6.0;

  galaxy::TPFCorePackage pkg;
  pkg.run_source_field_benchmark(c, out_dir);

  const std::string csv = slurp(std::string(out_dir) + "/tpf_source_field_probe_grid.csv");
  CHECK(csv.find(",0,1\n") != std::string::npos);
}

TEST_CASE("source-field benchmark unknown shape preserves legacy no-output behavior") {
  char dir_template[] = "/tmp/tpf_source_field_benchmark_unknown_shape_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "unknown_shape";
  c.tpf_source_probe_grid_half_extent = 10.0;
  c.tpf_source_probe_grid_n = 9;
  c.softening = 0.1;

  galaxy::TPFCorePackage pkg;
  CHECK_NOTHROW(pkg.run_source_field_benchmark(c, out_dir));
  CHECK_FALSE(file_exists(std::string(out_dir) + "/tpf_source_field_probe_grid.csv"));
}

TEST_CASE("4d static residual benchmark writes summary and slice artifacts") {
  char dir_template[] = "/tmp/tpf_4d_static_residual_artifacts_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_4d_residual_grid_n = 7;
  c.tpf_4d_residual_grid_half_extent = 3.0;
  c.tpf_4d_residual_source_exclusion_radius = 0.5;
  c.tpf_4d_residual_field_softening = 0.1;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_static_residual_benchmark(c, out_dir);

  const std::string summary_path = std::string(out_dir) + "/tpf_4d_static_residual_summary.txt";
  const std::string csv_path = std::string(out_dir) + "/tpf_4d_static_residual_slice.csv";
  const std::string csv_xy_path = std::string(out_dir) + "/tpf_4d_static_residual_slice_xy.csv";
  const std::string csv_xz_path = std::string(out_dir) + "/tpf_4d_static_residual_slice_xz.csv";
  const std::string csv_yz_path = std::string(out_dir) + "/tpf_4d_static_residual_slice_yz.csv";
  const std::string sources_csv_path = std::string(out_dir) + "/tpf_4d_static_residual_sources.csv";
  const std::string bins_nearest_path = std::string(out_dir) + "/tpf_4d_static_residual_bins_nearest_source.csv";
  const std::string bins_origin_path = std::string(out_dir) + "/tpf_4d_static_residual_bins_origin.csv";
  CHECK(file_exists(summary_path));
  CHECK(file_exists(csv_path));
  CHECK(file_exists(csv_xy_path));
  CHECK(file_exists(csv_xz_path));
  CHECK(file_exists(csv_yz_path));
  CHECK(file_exists(sources_csv_path));
  CHECK(file_exists(bins_nearest_path));
  CHECK(file_exists(bins_origin_path));

  const std::string summary = slurp(summary_path);
  CHECK(summary.find("grid_n:") != std::string::npos);
  CHECK(summary.find("free cells used:") != std::string::npos);
  CHECK(summary.find("mean normalized residual:") != std::string::npos);
  CHECK(summary.find("residual bin count:") != std::string::npos);
  CHECK(summary.find("nearest-source residual bin radius max used:") != std::string::npos);
  CHECK(summary.find("origin residual bin radius max used:") != std::string::npos);
  CHECK(summary.find("residual bins nearest-source csv: tpf_4d_static_residual_bins_nearest_source.csv") !=
        std::string::npos);
  CHECK(summary.find("residual bins origin csv: tpf_4d_static_residual_bins_origin.csv") != std::string::npos);

  const std::string expected_header =
      "source_shape,x,y,z,residual_t,residual_x,residual_y,residual_z,residual_spatial_norm,"
      "residual_4_norm_like,theta_spatial_frobenius_norm,normalized_residual,is_boundary,is_near_source,used_in_summary,"
      "xi_t,xi_x,xi_y,xi_z,xi_spatial_norm,theta_trace_4d,invariant_I_4d";
  const std::string csv = slurp(csv_path);
  const std::string csv_xy = slurp(csv_xy_path);
  const std::string csv_xz = slurp(csv_xz_path);
  const std::string csv_yz = slurp(csv_yz_path);
  CHECK(csv.find(expected_header) != std::string::npos);
  CHECK(csv_xy.find(expected_header) != std::string::npos);
  CHECK(csv_xz.find(expected_header) != std::string::npos);
  CHECK(csv_yz.find(expected_header) != std::string::npos);

  const std::string sources_csv = slurp(sources_csv_path);
  CHECK(sources_csv.find("source_index,mass,x,y,z,source_config_id,source_shape") != std::string::npos);

  const std::string bins_header =
      "bin_index,r_min,r_max,r_mid,cell_count_total,cell_count_used,cell_count_boundary,cell_count_near_source,"
      "mean_normalized_residual,median_normalized_residual,p90_normalized_residual,p95_normalized_residual,"
      "p99_normalized_residual,max_normalized_residual,mean_residual_spatial_norm,median_residual_spatial_norm,"
      "p95_residual_spatial_norm,max_residual_spatial_norm,mean_theta_spatial_frobenius_norm,"
      "median_theta_spatial_frobenius_norm";
  CHECK(slurp(bins_nearest_path).find(bins_header) != std::string::npos);
  CHECK(slurp(bins_origin_path).find(bins_header) != std::string::npos);
}

TEST_CASE("4d static residual benchmark plot script parses under py_compile") {
  const int ret = std::system("python3 -m py_compile ../plot_tpf_4d_static_residual.py");
  CHECK(ret == 0);
}

TEST_CASE("4d static residual benchmark monopole smoke has finite summary and free cells") {
  char dir_template[] = "/tmp/tpf_4d_static_residual_monopole_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 2.0e12;
  c.tpf_4d_residual_grid_n = 9;
  c.tpf_4d_residual_grid_half_extent = 4.0;
  c.tpf_4d_residual_source_exclusion_radius = 0.4;
  c.tpf_4d_residual_field_softening = 0.15;
  c.tpf_4d_residual_bin_count = 16;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_static_residual_benchmark(c, out_dir);

  const std::string summary = slurp(std::string(out_dir) + "/tpf_4d_static_residual_summary.txt");
  const double free_cells = parse_summary_scalar(summary, "free cells used");
  const double mean_norm = parse_summary_scalar(summary, "mean residual spatial norm");
  const double max_norm = parse_summary_scalar(summary, "max residual spatial norm");
  CHECK(free_cells > 0.0);
  CHECK(std::isfinite(mean_norm));
  CHECK(std::isfinite(max_norm));

  const std::string bins_text = slurp(std::string(out_dir) + "/tpf_4d_static_residual_bins_nearest_source.csv");
  std::istringstream bins_in(bins_text);
  std::string line;
  REQUIRE(std::getline(bins_in, line));
  int rows = 0;
  int used_positive_bins = 0;
  while (std::getline(bins_in, line)) {
    if (line.empty()) continue;
    ++rows;
    const std::vector<std::string> cols = split_csv_line(line);
    REQUIRE(cols.size() == 20);
    const int cell_count_used = std::atoi(cols[5].c_str());
    if (cell_count_used > 0) {
      ++used_positive_bins;
      const double median_n = parse_csv_double_or_nan(cols[9]);
      const double p95_n = parse_csv_double_or_nan(cols[11]);
      const double max_n = parse_csv_double_or_nan(cols[13]);
      if (std::isfinite(median_n) && std::isfinite(p95_n)) {
        CHECK(p95_n >= median_n);
      }
      if (std::isfinite(p95_n) && std::isfinite(max_n)) {
        CHECK(max_n >= p95_n);
      }
    } else {
      CHECK(cols[8] == "nan");
      CHECK(cols[9] == "nan");
      CHECK(cols[13] == "nan");
    }
  }
  CHECK(rows == c.tpf_4d_residual_bin_count);
  CHECK(used_positive_bins > 0);
}

TEST_CASE("4d static residual benchmark bonded pair equal-mass fallback and unequal masses") {
  char dir_template[] = "/tmp/tpf_4d_static_residual_pair_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "bonded_pair";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_source_benchmark_mass1 = 0.0;
  c.tpf_source_benchmark_mass2 = 0.0;
  c.tpf_source_benchmark_separation = 8.0;
  c.tpf_source_benchmark_orientation_deg = 30.0;
  c.tpf_4d_residual_grid_n = 7;
  c.tpf_4d_residual_grid_half_extent = 4.0;
  c.tpf_4d_residual_source_exclusion_radius = 0.2;
  c.tpf_4d_residual_field_softening = 0.1;
  c.tpf_4d_residual_bin_count = 12;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_static_residual_benchmark(c, out_dir);
  std::string summary = slurp(std::string(out_dir) + "/tpf_4d_static_residual_summary.txt");
  CHECK(summary.find("source config id: bonded_pair_centered_origin_equal_mass") != std::string::npos);
  CHECK(summary.find("source masses: 5.000000e+11, 5.000000e+11") != std::string::npos);
  CHECK(std::isfinite(parse_summary_scalar(summary, "mean normalized residual")));

  c.tpf_source_benchmark_mass1 = 8.0e11;
  c.tpf_source_benchmark_mass2 = 2.0e11;
  pkg.run_4d_static_residual_benchmark(c, out_dir);
  summary = slurp(std::string(out_dir) + "/tpf_4d_static_residual_summary.txt");
  CHECK(summary.find("source config id: bonded_pair_centered_origin_explicit_mass1_mass2") != std::string::npos);
  CHECK(summary.find("source masses: 8.000000e+11, 2.000000e+11") != std::string::npos);
  CHECK(std::isfinite(parse_summary_scalar(summary, "max normalized residual")));

  const std::string bins_origin_text = slurp(std::string(out_dir) + "/tpf_4d_static_residual_bins_origin.csv");
  std::istringstream bins_origin_in(bins_origin_text);
  std::string line;
  REQUIRE(std::getline(bins_origin_in, line));
  int rows = 0;
  bool saw_excluded = false;
  while (std::getline(bins_origin_in, line)) {
    if (line.empty()) continue;
    ++rows;
    const std::vector<std::string> cols = split_csv_line(line);
    REQUIRE(cols.size() == 20);
    const int total = std::atoi(cols[4].c_str());
    const int used = std::atoi(cols[5].c_str());
    const int boundary = std::atoi(cols[6].c_str());
    const int near_source = std::atoi(cols[7].c_str());
    CHECK(total >= used);
    CHECK(total >= boundary);
    CHECK(total >= near_source);
    if ((boundary > 0 || near_source > 0) && used < total) {
      saw_excluded = true;
    }
  }
  CHECK(rows == c.tpf_4d_residual_bin_count);
  CHECK(saw_excluded);
}

TEST_CASE("4d static residual benchmark rejects invalid config values loudly") {
  char dir_template[] = "/tmp/tpf_4d_static_residual_invalid_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0;
  c.tpf_4d_residual_grid_n = 2;
  galaxy::TPFCorePackage pkg;
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_grid_n = 5;
  c.tpf_4d_residual_grid_half_extent = 0.0;
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_grid_half_extent = 1.0;
  c.tpf_4d_residual_field_softening = -0.1;
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_field_softening = 0.1;
  c.tpf_4d_residual_source_exclusion_radius = -0.1;
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_source_exclusion_radius = 0.0;
  c.tpf_source_benchmark_shape = "unknown_shape";
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_4d_residual_bin_count = 0;
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_bin_count = 4097;
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_bin_count = 16;
  c.tpf_4d_residual_bin_radius_max = std::numeric_limits<double>::infinity();
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_bin_radius_max = 0.0;
  c.tpf_4d_residual_field_softening = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_field_softening = 0.1;
  c.tpf_4d_residual_grid_half_extent = std::numeric_limits<double>::infinity();
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_grid_half_extent = 1.0;
  c.tpf_4d_residual_source_exclusion_radius = std::numeric_limits<double>::infinity();
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_4d_residual_source_exclusion_radius = 0.0;
  c.tpf_source_benchmark_total_mass = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_source_benchmark_total_mass = 1.0;
  c.tpf_source_benchmark_shape = "bonded_pair";
  c.tpf_source_benchmark_mass1 = std::numeric_limits<double>::quiet_NaN();
  c.tpf_source_benchmark_mass2 = 1.0;
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_source_benchmark_mass1 = 1.0;
  c.tpf_source_benchmark_mass2 = std::numeric_limits<double>::infinity();
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_source_benchmark_mass2 = 1.0;
  c.tpf_source_benchmark_separation = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));

  c.tpf_source_benchmark_separation = 1.0;
  c.tpf_source_benchmark_orientation_deg = std::numeric_limits<double>::infinity();
  CHECK_THROWS(pkg.run_4d_static_residual_benchmark(c, out_dir));
}

TEST_CASE("4d static motion readout benchmark writes required artifacts and finite summary stats") {
  char dir_template[] = "/tmp/tpf_4d_static_motion_readout_artifacts_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_4d_motion_probe_grid_n = 9;
  c.tpf_4d_motion_probe_grid_half_extent = 4.0;
  c.tpf_4d_motion_source_exclusion_radius = 0.5;
  c.tpf_4d_motion_field_softening = 0.1;
  c.tpf_4d_motion_kappa = 1.0;
  c.tpf_4d_motion_readout_scale = 1.0;
  c.tpf_4d_motion_bin_count = 12;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_static_motion_readout_benchmark(c, out_dir);

  const std::string summary_path = std::string(out_dir) + "/tpf_4d_static_motion_readout_summary.txt";
  const std::string probe_path = std::string(out_dir) + "/tpf_4d_static_motion_readout_probe_grid.csv";
  const std::string bins_path = std::string(out_dir) + "/tpf_4d_static_motion_readout_bins_origin.csv";
  CHECK(file_exists(summary_path));
  CHECK(file_exists(probe_path));
  CHECK(file_exists(bins_path));

  const std::string summary = slurp(summary_path);
  CHECK(summary.find("kappa_motion:") != std::string::npos);
  CHECK(summary.find("motion_readout_scale:") != std::string::npos);
  CHECK(summary.find("interior cells:") != std::string::npos);
  CHECK(summary.find("used/free probe cells:") != std::string::npos);
  CHECK(summary.find("interior/free probe cells:") == std::string::npos);
  CHECK(summary.find("measured log-log falloff slope for monopole if available:") != std::string::npos);
  CHECK(std::isfinite(parse_summary_scalar(summary, "mean acceleration norm")));
  CHECK(std::isfinite(parse_summary_scalar(summary, "max acceleration norm")));
}

TEST_CASE("4d static motion readout monopole axis direction and symmetry") {
  char dir_template[] = "/tmp/tpf_4d_static_motion_readout_axis_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 2.0e12;
  c.tpf_4d_motion_probe_grid_n = 11;
  c.tpf_4d_motion_probe_grid_half_extent = 5.0;
  c.tpf_4d_motion_source_exclusion_radius = 0.5;
  c.tpf_4d_motion_field_softening = 0.1;
  c.tpf_4d_motion_kappa = 1.0;
  c.tpf_4d_motion_readout_scale = 1.0;
  c.tpf_4d_motion_bin_count = 16;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_static_motion_readout_benchmark(c, out_dir);

  std::ifstream in(std::string(out_dir) + "/tpf_4d_static_motion_readout_probe_grid.csv");
  REQUIRE(static_cast<bool>(in));
  std::string line;
  REQUIRE(std::getline(in, line));
  const std::vector<std::string> header = split_csv_line(line);
  const int x_col = find_csv_col(header, "x");
  const int y_col = find_csv_col(header, "y");
  const int z_col = find_csv_col(header, "z");
  const int ax_col = find_csv_col(header, "a_x");
  const int ay_col = find_csv_col(header, "a_y");
  const int az_col = find_csv_col(header, "a_z");
  const int an_col = find_csv_col(header, "a_norm");
  const int used_col = find_csv_col(header, "used_in_summary");
  REQUIRE(x_col >= 0);
  REQUIRE(used_col >= 0);

  double ax_pos = 0.0, ay_pos = 0.0, az_pos = 0.0;
  double ax_neg = 0.0, ay_neg = 0.0, az_neg = 0.0;
  bool found_pos = false;
  bool found_neg = false;
  double a_xaxis = 0.0, a_yaxis = 0.0, a_zaxis = 0.0;
  bool found_x = false, found_y = false, found_z = false;

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> cols = split_csv_line(line);
    if (cols.size() != header.size()) continue;
    if (std::atoi(cols[used_col].c_str()) == 0) continue;
    const double x = std::atof(cols[x_col].c_str());
    const double y = std::atof(cols[y_col].c_str());
    const double z = std::atof(cols[z_col].c_str());
    const double ax = std::atof(cols[ax_col].c_str());
    const double ay = std::atof(cols[ay_col].c_str());
    const double az = std::atof(cols[az_col].c_str());
    const double an = std::atof(cols[an_col].c_str());
    if (!found_pos && x == 3.0 && y == 0.0 && z == 0.0) {
      found_pos = true;
      ax_pos = ax; ay_pos = ay; az_pos = az;
      a_xaxis = an; found_x = true;
    }
    if (!found_neg && x == -3.0 && y == 0.0 && z == 0.0) {
      found_neg = true;
      ax_neg = ax; ay_neg = ay; az_neg = az;
    }
    if (!found_y && x == 0.0 && y == 3.0 && z == 0.0) {
      a_yaxis = an; found_y = true;
    }
    if (!found_z && x == 0.0 && y == 0.0 && z == 3.0) {
      a_zaxis = an; found_z = true;
    }
  }
  REQUIRE(found_pos);
  REQUIRE(found_neg);
  CHECK(ax_pos < 0.0);
  CHECK(std::fabs(ay_pos) < std::fabs(ax_pos) * 1e-6 + 1e-12);
  CHECK(std::fabs(az_pos) < std::fabs(ax_pos) * 1e-6 + 1e-12);
  CHECK(ax_neg > 0.0);
  CHECK(std::fabs(ay_neg) < std::fabs(ax_neg) * 1e-6 + 1e-12);
  CHECK(std::fabs(az_neg) < std::fabs(ax_neg) * 1e-6 + 1e-12);
  REQUIRE(found_x);
  REQUIRE(found_y);
  REQUIRE(found_z);
  CHECK(a_yaxis == doctest::Approx(a_xaxis).epsilon(0.25));
  CHECK(a_zaxis == doctest::Approx(a_xaxis).epsilon(0.25));
}

TEST_CASE("4d static motion readout monopole falloff slope is finite when available") {
  char dir_template[] = "/tmp/tpf_4d_static_motion_readout_falloff_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_4d_motion_probe_grid_n = 11;
  c.tpf_4d_motion_probe_grid_half_extent = 6.0;
  c.tpf_4d_motion_source_exclusion_radius = 0.75;
  c.tpf_4d_motion_field_softening = 0.1;
  c.tpf_4d_motion_bin_count = 16;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_static_motion_readout_benchmark(c, out_dir);
  const std::string summary = slurp(std::string(out_dir) + "/tpf_4d_static_motion_readout_summary.txt");
  const double slope = parse_summary_scalar(summary, "measured log-log falloff slope for monopole if available");
  CHECK(std::isfinite(slope));
}

TEST_CASE("4d static motion readout bonded pair smoke and source mass recording") {
  char dir_template[] = "/tmp/tpf_4d_static_motion_readout_pair_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "bonded_pair";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_source_benchmark_mass1 = 8.0e11;
  c.tpf_source_benchmark_mass2 = 2.0e11;
  c.tpf_source_benchmark_separation = 8.0;
  c.tpf_source_benchmark_orientation_deg = 45.0;
  c.tpf_4d_motion_probe_grid_n = 9;
  c.tpf_4d_motion_probe_grid_half_extent = 4.0;
  c.tpf_4d_motion_source_exclusion_radius = 0.5;
  c.tpf_4d_motion_bin_count = 12;

  galaxy::TPFCorePackage pkg;
  CHECK_NOTHROW(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
  const std::string summary = slurp(std::string(out_dir) + "/tpf_4d_static_motion_readout_summary.txt");
  CHECK(summary.find("source masses: 8.000000e+11, 2.000000e+11") != std::string::npos);
  CHECK(std::isfinite(parse_summary_scalar(summary, "mean acceleration norm")));
}

TEST_CASE("4d static motion readout benchmark rejects invalid config values loudly") {
  char dir_template[] = "/tmp/tpf_4d_static_motion_readout_invalid_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0;
  galaxy::TPFCorePackage pkg;

  c.tpf_4d_motion_probe_grid_n = 2;
  CHECK_THROWS(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
  c.tpf_4d_motion_probe_grid_n = 5;
  c.tpf_4d_motion_probe_grid_half_extent = 0.0;
  CHECK_THROWS(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
  c.tpf_4d_motion_probe_grid_half_extent = 1.0;
  c.tpf_4d_motion_source_exclusion_radius = -0.1;
  CHECK_THROWS(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
  c.tpf_4d_motion_source_exclusion_radius = 0.1;
  c.tpf_4d_motion_field_softening = -0.1;
  CHECK_THROWS(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
  c.tpf_4d_motion_field_softening = 0.1;
  c.tpf_4d_motion_kappa = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
  c.tpf_4d_motion_kappa = 1.0;
  c.tpf_4d_motion_readout_scale = std::numeric_limits<double>::infinity();
  CHECK_THROWS(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
  c.tpf_4d_motion_readout_scale = 1.0;
  c.tpf_4d_motion_bin_count = 0;
  CHECK_THROWS(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
  c.tpf_4d_motion_bin_count = 4097;
  CHECK_THROWS(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
  c.tpf_4d_motion_bin_count = 16;
  c.tpf_source_benchmark_shape = "unknown";
  CHECK_THROWS(pkg.run_4d_static_motion_readout_benchmark(c, out_dir));
}

TEST_CASE("4d xi motion probe benchmark writes required artifacts and finite valid rows") {
  char dir_template[] = "/tmp/tpf_4d_xi_motion_probe_artifacts_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_4d_xi_motion_dt = 1.0e-3;
  c.tpf_4d_xi_motion_steps = 20;
  c.tpf_4d_xi_motion_probe_layout = "ring";
  c.tpf_4d_xi_motion_probe_count = 12;
  c.tpf_4d_xi_motion_probe_radius = 10.0;
  c.tpf_4d_xi_motion_probe_speed = 1.0;
  c.tpf_4d_xi_motion_integrator = "velocity_verlet";
  c.tpf_4d_xi_motion_dump_every = 1;

  galaxy::TPFCorePackage pkg;
  CHECK_NOTHROW(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));

  const std::string summary_path = std::string(out_dir) + "/tpf_4d_xi_motion_probe_summary.txt";
  const std::string traj_path = std::string(out_dir) + "/tpf_4d_xi_motion_probe_trajectories.csv";
  const std::string init_path = std::string(out_dir) + "/tpf_4d_xi_motion_probe_initial_readout.csv";
  CHECK(file_exists(summary_path));
  CHECK(file_exists(traj_path));
  CHECK(file_exists(init_path));

  const std::string summary = slurp(summary_path);
  CHECK(summary.find("readout name: GravityXiMotionReadout_v1") != std::string::npos);
  CHECK(summary.find("xi kernel readout name: GravityXiKernelDeformation_v1") != std::string::npos);
  CHECK(summary.find("acceleration formula: a=-K_xi*Xi_eff_spatial") != std::string::npos);
  CHECK(summary.find("no compute_accelerations(...) calls: true") != std::string::npos);
  CHECK(summary.find("no compute_direct_tpf_accelerations(...) calls: true") != std::string::npos);
  CHECK(summary.find("note: old additive acceleration VDSG path is not used") != std::string::npos);
  CHECK(summary.find("escaped/invalid probe count:") != std::string::npos);
  CHECK(summary.find("escaped/invalid trajectory row count:") != std::string::npos);

  std::ifstream in(traj_path);
  REQUIRE(static_cast<bool>(in));
  std::string line;
  REQUIRE(std::getline(in, line));
  const std::vector<std::string> header = split_csv_line(line);
  const int step_col = find_csv_col(header, "step");
  const int x_col = find_csv_col(header, "x");
  const int valid_col = find_csv_col(header, "valid");
  CHECK(find_csv_col(header, "xi_t") >= 0);
  CHECK(find_csv_col(header, "xi_x_base") >= 0);
  CHECK(find_csv_col(header, "xi_x_eff") >= 0);
  CHECK(find_csv_col(header, "beta_rel") >= 0);
  CHECK(find_csv_col(header, "gamma_rel") >= 0);
  CHECK(find_csv_col(header, "xi_kernel_factor_raw") >= 0);
  CHECK(find_csv_col(header, "xi_kernel_metric_scale") >= 0);
  CHECK(find_csv_col(header, "xi_kernel_mode") >= 0);
  CHECK(find_csv_col(header, "xi_temporal_mode") >= 0);
  REQUIRE(step_col >= 0);
  REQUIRE(x_col >= 0);
  REQUIRE(valid_col >= 0);

  bool saw_step0 = false;
  bool saw_later = false;
  bool finite_valid_rows = true;
  double x0 = std::numeric_limits<double>::quiet_NaN();
  double x_last = std::numeric_limits<double>::quiet_NaN();
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> cols = split_csv_line(line);
    if (cols.size() != header.size()) continue;
    const int step = std::atoi(cols[step_col].c_str());
    const int valid = std::atoi(cols[valid_col].c_str());
    const double x = std::atof(cols[x_col].c_str());
    if (step == 0) {
      saw_step0 = true;
      if (!std::isfinite(x0)) x0 = x;
    } else {
      saw_later = true;
      x_last = x;
    }
    if (valid == 1 && !std::isfinite(x)) finite_valid_rows = false;
  }
  CHECK(saw_step0);
  CHECK(saw_later);
  CHECK(finite_valid_rows);
  CHECK(std::isfinite(x0));
  CHECK(std::isfinite(x_last));
  CHECK(x_last != doctest::Approx(x0));
}

TEST_CASE("4d xi motion probe benchmark off mode preserves Xi base and acceleration readout") {
  char dir_template[] = "/tmp/tpf_4d_xi_motion_probe_off_regression_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_4d_xi_motion_steps = 1;
  c.tpf_4d_xi_motion_probe_layout = "axis";
  c.tpf_4d_xi_motion_probe_count = 6;
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_motion_readout_scale = 2.0e-12;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_xi_motion_probe_benchmark(c, out_dir);

  std::ifstream in(std::string(out_dir) + "/tpf_4d_xi_motion_probe_initial_readout.csv");
  REQUIRE(static_cast<bool>(in));
  std::string line;
  REQUIRE(std::getline(in, line));
  const std::vector<std::string> header = split_csv_line(line);
  const int xi_x_col = find_csv_col(header, "xi_x");
  const int xi_y_col = find_csv_col(header, "xi_y");
  const int xi_z_col = find_csv_col(header, "xi_z");
  const int xi_x_base_col = find_csv_col(header, "xi_x_base");
  const int xi_y_base_col = find_csv_col(header, "xi_y_base");
  const int xi_z_base_col = find_csv_col(header, "xi_z_base");
  const int xi_x_eff_col = find_csv_col(header, "xi_x_eff");
  const int xi_y_eff_col = find_csv_col(header, "xi_y_eff");
  const int xi_z_eff_col = find_csv_col(header, "xi_z_eff");
  const int ax_col = find_csv_col(header, "ax");
  const int ay_col = find_csv_col(header, "ay");
  const int az_col = find_csv_col(header, "az");
  const int xi_t_col = find_csv_col(header, "xi_t");
  REQUIRE(xi_x_col >= 0);
  REQUIRE(ax_col >= 0);
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> cols = split_csv_line(line);
    if (cols.size() != header.size()) continue;
    const double xi_x = std::atof(cols[xi_x_col].c_str());
    const double xi_y = std::atof(cols[xi_y_col].c_str());
    const double xi_z = std::atof(cols[xi_z_col].c_str());
    CHECK(xi_x == doctest::Approx(std::atof(cols[xi_x_base_col].c_str())));
    CHECK(xi_y == doctest::Approx(std::atof(cols[xi_y_base_col].c_str())));
    CHECK(xi_z == doctest::Approx(std::atof(cols[xi_z_base_col].c_str())));
    CHECK(xi_x == doctest::Approx(std::atof(cols[xi_x_eff_col].c_str())));
    CHECK(xi_y == doctest::Approx(std::atof(cols[xi_y_eff_col].c_str())));
    CHECK(xi_z == doctest::Approx(std::atof(cols[xi_z_eff_col].c_str())));
    CHECK(std::atof(cols[ax_col].c_str()) == doctest::Approx(-c.tpf_4d_xi_motion_readout_scale * xi_x));
    CHECK(std::atof(cols[ay_col].c_str()) == doctest::Approx(-c.tpf_4d_xi_motion_readout_scale * xi_y));
    CHECK(std::atof(cols[az_col].c_str()) == doctest::Approx(-c.tpf_4d_xi_motion_readout_scale * xi_z));
    CHECK(std::atof(cols[xi_t_col].c_str()) == doctest::Approx(0.0));
  }
}

TEST_CASE("4d xi motion probe benchmark scalar_beta scales Xi_eff by 1+factor_raw") {
  char dir_template[] = "/tmp/tpf_4d_xi_motion_probe_scalar_beta_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_4d_xi_motion_steps = 1;
  c.tpf_4d_xi_motion_probe_layout = "axis";
  c.tpf_4d_xi_motion_probe_count = 6;
  c.tpf_4d_xi_motion_probe_speed = 0.1 * 299792458.0;
  c.tpf_4d_xi_kernel_mode = "scalar_beta";
  c.tpf_4d_xi_kernel_coupling = 0.5;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_xi_motion_probe_benchmark(c, out_dir);
  std::ifstream in(std::string(out_dir) + "/tpf_4d_xi_motion_probe_initial_readout.csv");
  REQUIRE(static_cast<bool>(in));
  std::string line;
  REQUIRE(std::getline(in, line));
  const std::vector<std::string> header = split_csv_line(line);
  const int base_col = find_csv_col(header, "xi_x_base");
  const int eff_col = find_csv_col(header, "xi_x_eff");
  const int factor_col = find_csv_col(header, "xi_kernel_factor_raw");
  const int ax_col = find_csv_col(header, "ax");
  const int xicol = find_csv_col(header, "xi_x");
  REQUIRE(base_col >= 0);
  REQUIRE(eff_col >= 0);
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> cols = split_csv_line(line);
    if (cols.size() != header.size()) continue;
    const double base = std::atof(cols[base_col].c_str());
    const double factor = std::atof(cols[factor_col].c_str());
    const double eff = std::atof(cols[eff_col].c_str());
    CHECK(eff == doctest::Approx(base * (1.0 + factor)));
    CHECK(std::atof(cols[xicol].c_str()) == doctest::Approx(eff));
    CHECK(std::atof(cols[ax_col].c_str()) == doctest::Approx(-c.tpf_4d_xi_motion_readout_scale * eff));
  }
}

TEST_CASE("4d xi motion probe benchmark spacetime_metric emits xi_t diagnostics only") {
  char dir_template[] = "/tmp/tpf_4d_xi_motion_probe_spacetime_metric_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_4d_xi_motion_steps = 1;
  c.tpf_4d_xi_motion_probe_layout = "axis";
  c.tpf_4d_xi_motion_probe_count = 6;
  c.tpf_4d_xi_motion_probe_speed = 0.05 * 299792458.0;
  c.tpf_4d_xi_kernel_mode = "spacetime_metric";
  c.tpf_4d_xi_kernel_coupling = 0.75;
  c.tpf_4d_xi_temporal_mode = "norm_scaled";
  c.tpf_4d_xi_temporal_coupling = 1.0;
  c.tpf_4d_xi_source_speed_x = 1.0e5;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_xi_motion_probe_benchmark(c, out_dir);
  std::ifstream in(std::string(out_dir) + "/tpf_4d_xi_motion_probe_initial_readout.csv");
  REQUIRE(static_cast<bool>(in));
  std::string line;
  REQUIRE(std::getline(in, line));
  const std::vector<std::string> header = split_csv_line(line);
  const int xi_t_col = find_csv_col(header, "xi_t");
  const int ax_col = find_csv_col(header, "ax");
  const int xi_x_eff_col = find_csv_col(header, "xi_x_eff");
  REQUIRE(xi_t_col >= 0);

  bool saw_nonzero_xi_t = false;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> cols = split_csv_line(line);
    if (cols.size() != header.size()) continue;
    const double xi_t = std::atof(cols[xi_t_col].c_str());
    const double xi_x_eff = std::atof(cols[xi_x_eff_col].c_str());
    const double ax = std::atof(cols[ax_col].c_str());
    if (std::fabs(xi_t) > 0.0) saw_nonzero_xi_t = true;
    CHECK(ax == doctest::Approx(-c.tpf_4d_xi_motion_readout_scale * xi_x_eff));
  }
  CHECK(saw_nonzero_xi_t);
}

TEST_CASE("4d xi motion probe benchmark monopole axis acceleration points inward with low transverse fraction") {
  char dir_template[] = "/tmp/tpf_4d_xi_motion_probe_direction_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_4d_xi_motion_steps = 1;
  c.tpf_4d_xi_motion_probe_layout = "axis";
  c.tpf_4d_xi_motion_probe_count = 6;
  c.tpf_4d_xi_motion_probe_radius = 10.0;
  c.tpf_4d_xi_motion_probe_speed = 0.0;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_xi_motion_probe_benchmark(c, out_dir);
  std::ifstream in(std::string(out_dir) + "/tpf_4d_xi_motion_probe_initial_readout.csv");
  REQUIRE(static_cast<bool>(in));
  std::string line;
  REQUIRE(std::getline(in, line));
  const std::vector<std::string> header = split_csv_line(line);
  const int x_col = find_csv_col(header, "x");
  const int y_col = find_csv_col(header, "y");
  const int z_col = find_csv_col(header, "z");
  const int ax_col = find_csv_col(header, "ax");
  const int align_col = find_csv_col(header, "radial_alignment_to_origin_inward");
  const int tr_col = find_csv_col(header, "transverse_fraction_origin");
  REQUIRE(ax_col >= 0);

  bool found_xpos = false;
  double ax_xpos = 0.0;
  double mean_align = 0.0;
  double mean_transverse = 0.0;
  int count = 0;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> cols = split_csv_line(line);
    if (cols.size() != header.size()) continue;
    const double x = std::atof(cols[x_col].c_str());
    const double y = std::atof(cols[y_col].c_str());
    const double z = std::atof(cols[z_col].c_str());
    const double ax = std::atof(cols[ax_col].c_str());
    if (!found_xpos && x > 0.0 && std::fabs(y) < 1e-12 && std::fabs(z) < 1e-12) {
      found_xpos = true;
      ax_xpos = ax;
    }
    const double align = parse_csv_double_or_nan(cols[align_col]);
    const double tr = parse_csv_double_or_nan(cols[tr_col]);
    if (std::isfinite(align) && std::isfinite(tr)) {
      mean_align += align;
      mean_transverse += tr;
      ++count;
    }
  }
  REQUIRE(found_xpos);
  CHECK(ax_xpos < 0.0);
  REQUIRE(count > 0);
  mean_align /= static_cast<double>(count);
  mean_transverse /= static_cast<double>(count);
  CHECK(mean_align > 0.95);
  CHECK(mean_transverse < 0.1);
}

TEST_CASE("4d xi motion probe benchmark monopole falloff slope near -2 and near-source guard marks invalid") {
  char dir_template[] = "/tmp/tpf_4d_xi_motion_probe_falloff_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0e12;
  c.tpf_4d_xi_motion_steps = 2;
  c.tpf_4d_xi_motion_probe_layout = "axis";
  c.tpf_4d_xi_motion_probe_count = 6;
  c.tpf_4d_xi_motion_probe_radius = 0.1;
  c.tpf_4d_xi_motion_source_exclusion_radius = 0.5;

  galaxy::TPFCorePackage pkg;
  pkg.run_4d_xi_motion_probe_benchmark(c, out_dir);
  const std::string summary = slurp(std::string(out_dir) + "/tpf_4d_xi_motion_probe_summary.txt");
  const double slope = parse_summary_scalar(summary, "measured Xi acceleration falloff slope from initial probe samples if available");
  if (std::isfinite(slope)) {
    CHECK(slope > -2.3);
    CHECK(slope < -1.7);
  }
  const double invalid_probe_count = parse_summary_scalar(summary, "escaped/invalid probe count");
  const double invalid_row_count = parse_summary_scalar(summary, "escaped/invalid trajectory row count");
  CHECK(invalid_probe_count > 0.0);
  CHECK(invalid_probe_count <= static_cast<double>(c.tpf_4d_xi_motion_probe_count));
  CHECK(invalid_row_count >= invalid_probe_count);
  CHECK(invalid_row_count > invalid_probe_count);
  CHECK(std::isnan(parse_summary_scalar(summary, "minimum radius reached")));
  CHECK(std::isnan(parse_summary_scalar(summary, "maximum radius reached")));
}

TEST_CASE("4d xi motion probe benchmark rejects invalid config values loudly") {
  char dir_template[] = "/tmp/tpf_4d_xi_motion_probe_invalid_XXXXXX";
  char* out_dir = mkdtemp(dir_template);
  REQUIRE(out_dir != nullptr);

  galaxy::Config c;
  c.tpf_source_benchmark_shape = "monopole";
  c.tpf_source_benchmark_total_mass = 1.0;
  galaxy::TPFCorePackage pkg;

  c.tpf_4d_xi_motion_dt = 0.0;
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_motion_dt = 0.001;
  c.tpf_4d_xi_motion_steps = 0;
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_motion_steps = 1;
  c.tpf_4d_xi_motion_probe_count = 0;
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_motion_probe_count = 8;
  c.tpf_4d_xi_motion_probe_layout = "grid";
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_motion_probe_layout = "ring";
  c.tpf_4d_xi_motion_integrator = "rk4";
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_motion_integrator = "velocity_verlet";
  c.tpf_4d_xi_kernel_mode = "bad";
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_kernel_factor_mode = "bad";
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_kernel_factor_mode = "beta_power";
  c.tpf_4d_xi_temporal_mode = "bad";
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_temporal_mode = "off";
  c.tpf_4d_xi_kernel_beta_power = -1.0;
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_kernel_beta_power = 1.0;
  c.tpf_4d_xi_kernel_metric_min = 0.0;
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_kernel_metric_min = 0.1;
  c.tpf_4d_xi_kernel_metric_max = 0.05;
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_kernel_metric_max = 10.0;
  c.tpf_4d_xi_source_speed_x = std::numeric_limits<double>::infinity();
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
  c.tpf_4d_xi_source_speed_x = 0.0;
  c.tpf_source_benchmark_shape = "unknown";
  CHECK_THROWS(pkg.run_4d_xi_motion_probe_benchmark(c, out_dir));
}
