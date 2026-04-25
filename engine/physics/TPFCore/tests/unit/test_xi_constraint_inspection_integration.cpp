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
