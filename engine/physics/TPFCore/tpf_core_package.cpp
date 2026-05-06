/**
 * TPFCore package implementation.
 * Canonical branch route: tpf_xi_theta_v1.
 * - Sources contribute Xi_a.
 * - Xi_total = sum_a Xi_a.
 * - Motion update uses Xi_total only.
 * - Theta is the diagnostic unsymmetrized spatial Jacobian: Theta_ij = d_j Xi_i,total.
 * - Theta must not alter motion in v1.
 */

#include "tpf_core_package.hpp"
#include "../../accel_pipeline_stats.hpp"
#include "../../config.hpp"
#include "../../softening_policy.hpp"
#include "../physics_package.hpp"  /* get_physics_package for Newtonian benchmark and live audits */
#include "derived_tpf_radial.hpp"
#include "field_evaluation.hpp"
#include "direct_tpf_baseline.hpp"
#include "extensions_vdsg.hpp"
#include "provisional_readout.hpp"
#include "regime_diagnostics.hpp"
#include "runtime_package_helpers.hpp"
#include "source_ansatz.hpp"
#include "tpf_core_params.hpp"
#include "tpf_4d_static_residual.hpp"
#include "v11_weak_field_correspondence.hpp"
#include "xi_constraint_exterior_solver.hpp"
#include "source_iteration.hpp"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace galaxy {

/** Single place that maps simulator Config -> TPFCore params. Keeps TPFCore decoupled from Config. */
static tpfcore::TPFCoreParams build_params(const Config& config, const std::string& output_dir) {
  tpfcore::TPFCoreParams p;
  p.output_dir = output_dir;
  p.softening = config.softening;
  p.bh_mass = config.bh_mass;
  p.star_mass = config.star_mass;
  p.enable_star_star_gravity = config.enable_star_star_gravity;
  p.tpfcore_source_softening = config.tpfcore_source_softening;
  p.effective_source_softening = (config.tpfcore_source_softening > 0.0) ? config.tpfcore_source_softening : config.softening;
  p.tpfcore_probe_radius_min = config.tpfcore_probe_radius_min;
  p.tpfcore_probe_radius_max = config.tpfcore_probe_radius_max;
  p.tpfcore_probe_samples = config.tpfcore_probe_samples;
  p.tpfcore_dump_theta_profile = config.tpfcore_dump_theta_profile;
  p.tpfcore_dump_invariant_profile = config.tpfcore_dump_invariant_profile;
  p.tpf_xi_constraint_exterior_inspect = config.tpf_xi_constraint_exterior_inspect;
  p.tpf_xi_constraint_grid_n = config.tpf_xi_constraint_grid_n;
  p.tpf_xi_constraint_half_extent = config.tpf_xi_constraint_half_extent;
  p.tpf_xi_constraint_inner_radius = config.tpf_xi_constraint_inner_radius;
  p.tpf_xi_constraint_max_iters = config.tpf_xi_constraint_max_iters;
  p.tpf_xi_constraint_tol = config.tpf_xi_constraint_tol;
  p.tpfcore_dump_readout_debug = config.tpfcore_dump_readout_debug;
  p.validation_symmetric_separation = config.validation_symmetric_separation;
  return p;
}

TPFCorePackage::TPFCorePackage()
    : tpf_dynamics_mode_("tpf_xi_theta_v1"),
      provisional_readout_(false),
      readout_mode_("tensor_radial_projection"),
      readout_scale_(1.0),
      theta_tt_scale_(1.0),
      theta_tr_scale_(1.0),
      source_softening_(0.0),
      kappa_(1.0e32),
      weak_field_correspondence_alpha_si_(-tpfcore::TPF_G_SI),
      vdsg_coupling_(1.0e-20),
      vdsg_mass_baseline_resolved_kg_(0.0),
      vdsg_mode_("legacy_speed"),
      vdsg_mass_gate_m0_kg_(1.98847e30),
      vdsg_mass_gate_alpha_(1.0),
      vdsg_x_clamp_(0.25),
      vdsg_weak_field_gate_enable_(true),
      vdsg_weak_field_a0_(1.0e-10),
      vdsg_weak_field_power_(1.0),
      vdsg_bounded_amplitude_(0.25),
      simulation_dt_(0.01),
      cooling_fraction_(0.2),
      shunt_enable_(false),
      shunt_fraction_(0.001),
      pipeline_diagnostics_csv_(true),
      xi_motion_readout_scale_(1.0e-12),
      xi_kernel_mode_("off"),
      xi_kernel_coupling_(0.0),
      xi_kernel_beta_power_(1.0),
      xi_kernel_factor_mode_("beta_power"),
      xi_kernel_metric_min_(0.1),
      xi_kernel_metric_max_(10.0),
      xi_temporal_mode_("off"),
      xi_temporal_coupling_(0.0),
      xi_source_speed_x_(0.0),
      xi_source_speed_y_(0.0),
      xi_source_speed_z_(0.0) {}

void TPFCorePackage::init_from_config(const Config& config) {
  tpf_dynamics_mode_ = config.tpf_dynamics_mode;
  if (tpf_dynamics_mode_ != "tpf_xi_theta_v1") {
    throw std::runtime_error("TPFCore on this branch supports only tpf_dynamics_mode=tpf_xi_theta_v1.");
  }
  provisional_readout_ = config.tpfcore_enable_provisional_readout;
  readout_mode_ = config.tpfcore_readout_mode;
  readout_scale_ = config.tpfcore_readout_scale;
  theta_tt_scale_ = config.tpfcore_theta_tt_scale;
  theta_tr_scale_ = config.tpfcore_theta_tr_scale;
  source_softening_ = config.tpfcore_source_softening;  /* 0 => use global softening at runtime */
  weak_field_correspondence_alpha_si_ = config.tpf_weak_field_correspondence_alpha_si;
  kappa_ = config.tpf_kappa;                // external flat key tpf_kappa -> internal direct_tpf paper coupling
  derived_poisson_cfg_.kappa = config.tpf_kappa;  // same incoming key also feeds derived-radial closure ledger
  derived_poisson_cfg_.bins = config.tpf_poisson_bins;
  derived_poisson_cfg_.max_radius = config.tpf_poisson_max_radius;
  derived_poisson_cfg_.galaxy_radius = config.galaxy_radius;
  vdsg_coupling_ = config.tpf_vdsg_coupling;
  vdsg_mass_baseline_resolved_kg_ =
      (config.tpf_vdsg_mass_baseline_kg > 0.0) ? config.tpf_vdsg_mass_baseline_kg : config.star_mass;
  vdsg_mode_ = config.tpf_vdsg_mode;
  if (!tpfcore::is_valid_vdsg_mode(vdsg_mode_)) {
    throw std::runtime_error("invalid tpf_vdsg_mode: " + vdsg_mode_ +
                             "; expected legacy_speed, radial_doppler_rational, radial_doppler_exp, or radial_doppler_bounded");
  }
  vdsg_mass_gate_m0_kg_ = config.tpf_vdsg_mass_gate_m0_kg;
  vdsg_mass_gate_alpha_ = config.tpf_vdsg_mass_gate_alpha;
  vdsg_x_clamp_ = config.tpf_vdsg_x_clamp;
  vdsg_weak_field_gate_enable_ = config.tpf_vdsg_weak_field_gate_enable;
  vdsg_weak_field_a0_ = config.tpf_vdsg_weak_field_a0;
  vdsg_weak_field_power_ = config.tpf_vdsg_weak_field_power;
  vdsg_bounded_amplitude_ = config.tpf_vdsg_bounded_amplitude;
  simulation_dt_ = config.dt;
  cooling_fraction_ = config.tpf_cooling_fraction;
  shunt_enable_ = config.tpf_global_accel_shunt_enable;
  shunt_fraction_ = (config.tpf_global_accel_shunt_fraction > 0.0 && std::isfinite(config.tpf_global_accel_shunt_fraction))
                        ? config.tpf_global_accel_shunt_fraction
                        : 0.001;
  pipeline_diagnostics_csv_ = config.tpf_accel_pipeline_diagnostics_csv;
  xi_motion_readout_scale_ = config.tpf_4d_xi_motion_readout_scale;
  xi_kernel_mode_ = config.tpf_4d_xi_kernel_mode;
  xi_kernel_coupling_ = config.tpf_4d_xi_kernel_coupling;
  xi_kernel_beta_power_ = config.tpf_4d_xi_kernel_beta_power;
  xi_kernel_factor_mode_ = config.tpf_4d_xi_kernel_factor_mode;
  xi_kernel_metric_min_ = config.tpf_4d_xi_kernel_metric_min;
  xi_kernel_metric_max_ = config.tpf_4d_xi_kernel_metric_max;
  xi_temporal_mode_ = config.tpf_4d_xi_temporal_mode;
  xi_temporal_coupling_ = config.tpf_4d_xi_temporal_coupling;
  xi_source_speed_x_ = config.tpf_4d_xi_source_speed_x;
  xi_source_speed_y_ = config.tpf_4d_xi_source_speed_y;
  xi_source_speed_z_ = config.tpf_4d_xi_source_speed_z;
}

namespace {

const double C_SI_LIGHT = 299792458.0;

struct XiKernelRuntimeSource {
  double x = 0.0, y = 0.0, z = 0.0;
  double vx = 0.0, vy = 0.0, vz = 0.0;
  double mass = 0.0;
};

struct BenchmarkSourceSpec2D {
  double x;
  double y;
  double m;
};

struct BenchmarkSourceSpec3D {
  double x;
  double y;
  double z;
  double m;
};

struct ResidualBinAccumulator {
  std::size_t cell_count_total = 0;
  std::size_t cell_count_used = 0;
  std::size_t cell_count_boundary = 0;
  std::size_t cell_count_near_source = 0;
  std::vector<double> normalized_residual;
  std::vector<double> residual_spatial_norm;
  std::vector<double> theta_spatial_frobenius_norm;
};

double nearest_source_distance(const std::vector<BenchmarkSourceSpec3D>& source_specs, double x, double y, double z) {
  if (source_specs.empty()) {
    return 0.0;
  }
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    const double dx = x - source_specs[i].x;
    const double dy = y - source_specs[i].y;
    const double dz = z - source_specs[i].z;
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (r < best) best = r;
  }
  return best;
}

double percentile_sorted_linear(const std::vector<double>& sorted_values, double p) {
  if (sorted_values.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (p <= 0.0) return sorted_values.front();
  if (p >= 1.0) return sorted_values.back();
  const double rank = p * static_cast<double>(sorted_values.size() - 1);
  const std::size_t lo = static_cast<std::size_t>(std::floor(rank));
  const std::size_t hi = static_cast<std::size_t>(std::ceil(rank));
  if (lo == hi) return sorted_values[lo];
  const double w = rank - static_cast<double>(lo);
  return (1.0 - w) * sorted_values[lo] + w * sorted_values[hi];
}

double mean_or_nan(const std::vector<double>& values) {
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  double sum = 0.0;
  for (std::size_t i = 0; i < values.size(); ++i) sum += values[i];
  return sum / static_cast<double>(values.size());
}

double max_or_nan(const std::vector<double>& values) {
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  double out = values[0];
  for (std::size_t i = 1; i < values.size(); ++i) out = std::max(out, values[i]);
  return out;
}

std::size_t compute_bin_index(double r, double radius_max, int bin_count) {
  if (bin_count <= 1 || !(radius_max > 0.0) || !std::isfinite(radius_max)) {
    return 0;
  }
  if (r <= 0.0) return 0;
  if (r >= radius_max) return static_cast<std::size_t>(bin_count - 1);
  const double unit = r / radius_max;
  const std::size_t idx = static_cast<std::size_t>(std::floor(unit * static_cast<double>(bin_count)));
  const std::size_t max_idx = static_cast<std::size_t>(bin_count - 1);
  return std::min(idx, max_idx);
}

double write_residual_bins_csv(const std::string& csv_path,
                               const tpfcore::StaticResidualGridResult& result,
                               const std::vector<BenchmarkSourceSpec3D>& source_specs,
                               int bin_count,
                               double configured_radius_max,
                               bool bin_by_nearest_source) {
  std::vector<double> radius_per_point(result.points.size(), 0.0);
  double observed_radius_max = 0.0;
  for (std::size_t i = 0; i < result.points.size(); ++i) {
    const tpfcore::StaticResidualAtPoint& p = result.points[i];
    const double r = bin_by_nearest_source
                         ? nearest_source_distance(source_specs, p.x, p.y, p.z)
                         : std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    radius_per_point[i] = r;
    observed_radius_max = std::max(observed_radius_max, r);
  }

  double radius_max_used = configured_radius_max;
  if (!(radius_max_used > 0.0)) {
    radius_max_used = observed_radius_max;
  }
  if (!(radius_max_used > 0.0) || !std::isfinite(radius_max_used)) {
    radius_max_used = 1.0;
  }

  std::vector<ResidualBinAccumulator> bins(static_cast<std::size_t>(bin_count));
  for (std::size_t i = 0; i < result.points.size(); ++i) {
    const tpfcore::StaticResidualAtPoint& p = result.points[i];
    const std::size_t bi = compute_bin_index(radius_per_point[i], radius_max_used, bin_count);
    ResidualBinAccumulator& b = bins[bi];
    ++b.cell_count_total;
    if (p.is_boundary) ++b.cell_count_boundary;
    if (p.is_near_source) ++b.cell_count_near_source;
    if (p.used_in_summary) {
      ++b.cell_count_used;
      b.normalized_residual.push_back(p.normalized_residual);
      b.residual_spatial_norm.push_back(p.residual_spatial_norm);
      b.theta_spatial_frobenius_norm.push_back(p.theta_spatial_frobenius_norm);
    }
  }

  std::ofstream csv(csv_path);
  if (!csv) return radius_max_used;
  csv << std::scientific;
  csv << "bin_index,r_min,r_max,r_mid,cell_count_total,cell_count_used,cell_count_boundary,cell_count_near_source,"
         "mean_normalized_residual,median_normalized_residual,p90_normalized_residual,p95_normalized_residual,"
         "p99_normalized_residual,max_normalized_residual,mean_residual_spatial_norm,median_residual_spatial_norm,"
         "p95_residual_spatial_norm,max_residual_spatial_norm,mean_theta_spatial_frobenius_norm,"
         "median_theta_spatial_frobenius_norm\n";

  const double bin_width = radius_max_used / static_cast<double>(bin_count);
  for (int i = 0; i < bin_count; ++i) {
    const ResidualBinAccumulator& b = bins[static_cast<std::size_t>(i)];
    const double r_min = bin_width * static_cast<double>(i);
    const double r_max = (i + 1 == bin_count) ? radius_max_used : bin_width * static_cast<double>(i + 1);
    const double r_mid = 0.5 * (r_min + r_max);
    std::vector<double> n = b.normalized_residual;
    std::vector<double> s = b.residual_spatial_norm;
    std::vector<double> t = b.theta_spatial_frobenius_norm;
    std::sort(n.begin(), n.end());
    std::sort(s.begin(), s.end());
    std::sort(t.begin(), t.end());
    csv << i << "," << r_min << "," << r_max << "," << r_mid << ","
        << b.cell_count_total << "," << b.cell_count_used << "," << b.cell_count_boundary << ","
        << b.cell_count_near_source << ","
        << mean_or_nan(n) << "," << percentile_sorted_linear(n, 0.5) << "," << percentile_sorted_linear(n, 0.9)
        << "," << percentile_sorted_linear(n, 0.95) << "," << percentile_sorted_linear(n, 0.99) << ","
        << max_or_nan(n) << "," << mean_or_nan(s) << "," << percentile_sorted_linear(s, 0.5) << ","
        << percentile_sorted_linear(s, 0.95) << "," << max_or_nan(s) << "," << mean_or_nan(t) << ","
        << percentile_sorted_linear(t, 0.5) << "\n";
  }
  return radius_max_used;
}

std::vector<BenchmarkSourceSpec3D> build_tpf_benchmark_sources_3d(const Config& config, std::string* shape_id_out) {
  const std::string shape = config.tpf_source_benchmark_shape;
  const double total_mass = config.tpf_source_benchmark_total_mass;
  const double configured_mass1 = config.tpf_source_benchmark_mass1;
  const double configured_mass2 = config.tpf_source_benchmark_mass2;
  const double separation = config.tpf_source_benchmark_separation;
  const double orientation_deg = config.tpf_source_benchmark_orientation_deg;

  if (shape == "monopole") {
    if (shape_id_out) *shape_id_out = "monopole_centered_origin";
    return std::vector<BenchmarkSourceSpec3D>{BenchmarkSourceSpec3D{0.0, 0.0, 0.0, total_mass}};
  }
  if (shape == "bonded_pair") {
    const double angle_rad = orientation_deg * (3.14159265358979323846 / 180.0);
    const double dx = 0.5 * separation * std::cos(angle_rad);
    const double dy = 0.5 * separation * std::sin(angle_rad);
    const bool use_unequal_masses = (configured_mass1 > 0.0 && configured_mass2 > 0.0);
    const double source_mass1 = use_unequal_masses ? configured_mass1 : (0.5 * total_mass);
    const double source_mass2 = use_unequal_masses ? configured_mass2 : (0.5 * total_mass);
    if (shape_id_out) {
      *shape_id_out = use_unequal_masses
                          ? "bonded_pair_centered_origin_explicit_mass1_mass2"
                          : "bonded_pair_centered_origin_equal_mass";
    }
    return std::vector<BenchmarkSourceSpec3D>{
        BenchmarkSourceSpec3D{dx, dy, 0.0, source_mass1},
        BenchmarkSourceSpec3D{-dx, -dy, 0.0, source_mass2}};
  }
  throw std::runtime_error("unknown tpf_source_benchmark_shape for benchmark source setup: " + shape);
}

bool try_build_tpf_benchmark_sources_3d(const Config& config,
                                        std::vector<BenchmarkSourceSpec3D>* sources_out,
                                        std::string* shape_id_out) {
  const std::string shape = config.tpf_source_benchmark_shape;
  if (shape != "monopole" && shape != "bonded_pair") {
    return false;
  }
  *sources_out = build_tpf_benchmark_sources_3d(config, shape_id_out);
  return true;
}

void validate_stage4_residual_source_inputs(const Config& config) {
  if (!std::isfinite(config.tpf_source_benchmark_total_mass)) {
    throw std::runtime_error("tpf_source_benchmark_total_mass must be finite");
  }
  if (!std::isfinite(config.tpf_source_benchmark_mass1)) {
    throw std::runtime_error("tpf_source_benchmark_mass1 must be finite");
  }
  if (!std::isfinite(config.tpf_source_benchmark_mass2)) {
    throw std::runtime_error("tpf_source_benchmark_mass2 must be finite");
  }
  if (!std::isfinite(config.tpf_source_benchmark_separation)) {
    throw std::runtime_error("tpf_source_benchmark_separation must be finite");
  }
  if (!std::isfinite(config.tpf_source_benchmark_orientation_deg)) {
    throw std::runtime_error("tpf_source_benchmark_orientation_deg must be finite");
  }
}

void validate_stage4_residual_source_specs(const std::vector<BenchmarkSourceSpec3D>& source_specs) {
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    if (!std::isfinite(source_specs[i].x) ||
        !std::isfinite(source_specs[i].y) ||
        !std::isfinite(source_specs[i].z) ||
        !std::isfinite(source_specs[i].m)) {
      throw std::runtime_error("resolved Stage 4 source positions and masses must be finite");
    }
  }
}

struct MotionReadoutPoint {
  double x = 0.0, y = 0.0, z = 0.0;
  double r_origin = 0.0;
  double r_nearest_source = 0.0;
  double xi_t = 0.0, xi_x = 0.0, xi_y = 0.0, xi_z = 0.0, xi_spatial_norm = 0.0;
  double theta_trace_4d = 0.0;
  double invariant_I_4d = 0.0;
  double c_xx = 0.0, c_xy = 0.0, c_xz = 0.0, c_yx = 0.0, c_yy = 0.0, c_yz = 0.0, c_zx = 0.0, c_zy = 0.0, c_zz = 0.0;
  double a_x = 0.0, a_y = 0.0, a_z = 0.0, a_norm = 0.0;
  double radial_alignment_to_origin_inward = std::numeric_limits<double>::quiet_NaN();
  double transverse_fraction_origin = std::numeric_limits<double>::quiet_NaN();
  bool is_boundary = false;
  bool is_near_source = false;
  bool xi_degenerate = false;
  bool used_in_summary = false;
};

struct XiConstraintExteriorStats {
  size_t n_masked = 0;
  size_t n_pinned = 0;
  size_t n_free = 0;
  double initial_max_residual_free = 0.0;
  double final_max_residual_free = 0.0;
  double initial_max_residual_nonmasked = 0.0;
  double final_max_residual_nonmasked = 0.0;
  double max_delta_xi_free = 0.0;
  double mean_delta_xi_free = 0.0;
};

XiConstraintExteriorStats compute_xi_constraint_exterior_stats(
    const tpfcore::XiConstraintExteriorParams& solve_params,
    const tpfcore::XiConstraintExteriorSolveResult& solve_result) {
  XiConstraintExteriorStats stats;
  const tpfcore::PlanarXiGrid& g = solve_result.grid;
  const std::vector<tpfcore::PlanarConstraintResidual2D> initial_residuals =
      tpfcore::compute_planar_constraint_residual_field(tpfcore::initialize_planar_xi_grid_from_ansatz(solve_params));
  double sum_delta_free = 0.0;

  for (int j = 0; j < g.ny; ++j) {
    for (int i = 0; i < g.nx; ++i) {
      const int idx = g.index(i, j);
      if (!g.is_exterior[idx]) {
        ++stats.n_masked;
        continue;
      }

      if (g.is_pinned[idx]) {
        ++stats.n_pinned;
      } else {
        ++stats.n_free;
      }

      const double initial_norm = initial_residuals[idx].norm;
      const double final_norm = solve_result.residuals[idx].norm;
      stats.initial_max_residual_nonmasked = std::max(stats.initial_max_residual_nonmasked, initial_norm);
      stats.final_max_residual_nonmasked = std::max(stats.final_max_residual_nonmasked, final_norm);

      if (!g.is_pinned[idx]) {
        stats.initial_max_residual_free = std::max(stats.initial_max_residual_free, initial_norm);
        stats.final_max_residual_free = std::max(stats.final_max_residual_free, final_norm);
        const tpfcore::Xi2D xi_ansatz =
            tpfcore::provisional_point_source_field(0.0, 0.0, solve_params.source_mass, g.x_at(i), g.y_at(j), solve_params.softening).xi;
        const double dmag = std::hypot(g.xi_x[idx] - xi_ansatz.x, g.xi_y[idx] - xi_ansatz.y);
        stats.max_delta_xi_free = std::max(stats.max_delta_xi_free, dmag);
        sum_delta_free += dmag;
      }
    }
  }

  if (stats.n_free > 0) {
    stats.mean_delta_xi_free = sum_delta_free / static_cast<double>(stats.n_free);
  }
  return stats;
}

struct TPFCoreRegistrar {
  TPFCoreRegistrar() {
    register_physics_package_factory(
        "TPFCore",
        []() -> std::unique_ptr<PhysicsPackage> {
          return std::unique_ptr<PhysicsPackage>(new TPFCorePackage());
        });
  }
};

TPFCoreRegistrar s_tpfcore_registrar;

}  // namespace

void TPFCorePackage::eval_accel_pipeline(const State& state,
                                         double bh_mass,
                                         double softening,
                                         bool star_star,
                                         std::vector<double>& ax,
                                         std::vector<double>& ay,
                                         AccelPipelineStats* stats_out) const {
  const int n = state.n();
  ax.assign(n, 0.0);
  ay.assign(n, 0.0);

  const double eps = (source_softening_ > 0.0) ? source_softening_ : softening;
  if (tpfcore::is_derived_tpf_radial_readout_mode(readout_mode_)) {
    tpfcore::TpfRadialGravityProfile profile =
        tpfcore::build_tpf_gravity_profile(state, bh_mass, derived_poisson_cfg_, eps);
    for (int i = 0; i < n; ++i) {
      tpfcore::compute_provisional_readout_acceleration(
          state, i, bh_mass, star_star, softening, source_softening_,
          readout_mode_, readout_scale_,
          theta_tt_scale_, theta_tr_scale_, ax[i], ay[i],
          &derived_poisson_cfg_, &profile);
    }
  } else {
    for (int i = 0; i < n; ++i) {
      tpfcore::compute_provisional_readout_acceleration(
          state, i, bh_mass, star_star, softening, source_softening_,
          readout_mode_, readout_scale_,
          theta_tt_scale_, theta_tr_scale_, ax[i], ay[i],
          nullptr, nullptr);
    }
  }

  double sum_b = 0.0;
  for (int i = 0; i < n; ++i) {
    sum_b += std::hypot(ax[i], ay[i]);
  }
  const double mean_b = (n > 0) ? (sum_b / static_cast<double>(n)) : 0.0;

  std::vector<double> dax, day;
  tpfcore::VdsgDiagnosticsSummary vdsg_diag;
  tpfcore::accumulate_vdsg_velocity_modifier(state, bh_mass, softening, star_star, vdsg_coupling_,
                                             vdsg_mass_baseline_resolved_kg_, vdsg_mode_, vdsg_mass_gate_m0_kg_,
                                             vdsg_mass_gate_alpha_, vdsg_x_clamp_, vdsg_weak_field_gate_enable_,
                                             vdsg_weak_field_a0_, vdsg_weak_field_power_, vdsg_bounded_amplitude_, dax, day,
                                             &vdsg_diag);
  double sum_v = 0.0;
  for (int i = 0; i < n; ++i) {
    sum_v += std::hypot(dax[i], day[i]);
  }
  const double mean_v = (n > 0) ? (sum_v / static_cast<double>(n)) : 0.0;

  for (int i = 0; i < n; ++i) {
    ax[i] += dax[i];
    ay[i] += day[i];
  }

  const unsigned shunt_n =
      tpfcore::apply_global_accel_magnitude_shunt(state, simulation_dt_, shunt_enable_, shunt_fraction_, ax, ay);

  if (stats_out) {
    stats_out->valid = true;
    stats_out->mean_baseline_mag = mean_b;
    stats_out->mean_vdsg_mag = mean_v;
    stats_out->vdsg_pairs_evaluated = vdsg_diag.pairs_evaluated;
    stats_out->vdsg_min_beta_rad = vdsg_diag.min_beta_rad;
    stats_out->vdsg_max_beta_rad = vdsg_diag.max_beta_rad;
    stats_out->vdsg_mean_abs_beta_rad = (vdsg_diag.pairs_evaluated > 0) ? (vdsg_diag.sum_abs_beta_rad / vdsg_diag.pairs_evaluated) : 0.0;
    stats_out->vdsg_over_baseline_ratio = (mean_b > 1e-300) ? (mean_v / mean_b) : 0.0;
    stats_out->shunt_events_last_step = shunt_n;
    stats_out->frac_capped_last_step =
        (n > 0) ? (static_cast<double>(shunt_n) / static_cast<double>(n)) : 0.0;
    stats_out->shunt_enabled = shunt_enable_;
    stats_out->shunt_fraction = shunt_fraction_;
  }
}

void TPFCorePackage::compute_direct_tpf_accelerations(const State& state,
                                                      double bh_mass,
                                                      double softening,
                                                      bool star_star,
                                                      std::vector<double>& ax,
                                                      std::vector<double>& ay) const {
  const int n = state.n();
  ax.assign(n, 0.0);
  ay.assign(n, 0.0);
  const double eps = (source_softening_ > 0.0) ? source_softening_ : softening;
  const double kappa = kappa_;
  const double lambda = tpfcore::LAMBDA_4D;

  for (int i = 0; i < n; ++i) {
    // Explicit upstream boundary: Eq. (10) direct_tpf principal-part baseline is applied
    // over the current provisional field ansatz from evaluate_provisional_field_multi_source.
    const tpfcore::FieldAtPoint field =
        tpfcore::evaluate_provisional_field_multi_source(state, i, bh_mass, star_star, eps);
    const tpfcore::DirectTpfBaselineAccelerationResult baseline =
        tpfcore::compute_direct_tpf_baseline_acceleration(field, kappa, lambda);
    ax[i] = baseline.ax;
    ay[i] = baseline.ay;
  }
}

void TPFCorePackage::compute_v11_weak_field_truncation_accelerations(const State& state,
                                                                     double bh_mass,
                                                                     double softening,
                                                                     bool star_star,
                                                                     std::vector<double>& ax,
                                                                     std::vector<double>& ay) const {
  /**
   * Canonical weak-field/geodesic correspondence chain (Paper A draft appendix + TPF_FOUNDATIONS):
   *   matter -> Xi -> rho_Xi -> Poisson -> Phi -> a = -grad(Phi).
   *
   * Runtime on this route uses the analytic closed-form point-source shortcut of that chain.
   * It does NOT numerically compute rho_Xi, solve Poisson on a grid, or differentiate Phi on a grid.
   *
   * Singular correspondence limit for a point source:
   *   rho_Xi = M delta^3(x).
   *
   * Softened implementation uses:
   *   rho_eps(r) = 3 M eps^2 / [4 pi (r^2 + eps^2)^(5/2)],
   * which integrates to M over R^3.
   *
   * Theta_ij Theta^ij, I, and Q_wf are nonlinear ledger/diagnostic quantities.
   * They are NOT the leading density source used by geodesic_correspondence acceleration.
   */
  compute_v11_weak_field_correspondence_accelerations(state, bh_mass, softening, star_star,
                                                      weak_field_correspondence_alpha_si_, ax, ay);
}

void TPFCorePackage::apply_vdsg_additive_extension(const State& state,
                                                   double bh_mass,
                                                   double softening,
                                                   bool star_star,
                                                   std::vector<double>& ax,
                                                   std::vector<double>& ay) const {
  std::vector<double> dax, day;
  tpfcore::VdsgDiagnosticsSummary vdsg_diag;
  tpfcore::accumulate_vdsg_velocity_modifier(state, bh_mass, softening, star_star, vdsg_coupling_,
                                             vdsg_mass_baseline_resolved_kg_, vdsg_mode_, vdsg_mass_gate_m0_kg_,
                                             vdsg_mass_gate_alpha_, vdsg_x_clamp_, vdsg_weak_field_gate_enable_,
                                             vdsg_weak_field_a0_, vdsg_weak_field_power_, vdsg_bounded_amplitude_, dax, day,
                                             &vdsg_diag);
  for (int i = 0; i < state.n(); ++i) {
    ax[i] += dax[i];
    ay[i] += day[i];
  }
}

bool TPFCorePackage::xi_kernel_deformation_active() const {
  return (xi_kernel_mode_ != "off" && xi_kernel_coupling_ != 0.0);
}

void TPFCorePackage::validate_xi_kernel_runtime_config() const {
  if (!std::isfinite(xi_motion_readout_scale_)) throw std::runtime_error("tpf_4d_xi_motion_readout_scale must be finite");
  if (xi_kernel_mode_ != "off" && xi_kernel_mode_ != "scalar_beta" && xi_kernel_mode_ != "metric_radial" &&
      xi_kernel_mode_ != "metric_velocity" && xi_kernel_mode_ != "metric_transverse_wake" &&
      xi_kernel_mode_ != "metric_transverse_continuous" && xi_kernel_mode_ != "spacetime_metric") {
    throw std::runtime_error(
        "tpf_4d_xi_kernel_mode must be one of: off, scalar_beta, metric_radial, metric_velocity, metric_transverse_wake, metric_transverse_continuous, spacetime_metric");
  }
  if (xi_kernel_factor_mode_ != "beta_power" && xi_kernel_factor_mode_ != "gamma_minus_one") {
    throw std::runtime_error("tpf_4d_xi_kernel_factor_mode must be beta_power or gamma_minus_one");
  }
  if (xi_temporal_mode_ != "off" && xi_temporal_mode_ != "norm_scaled") {
    throw std::runtime_error("tpf_4d_xi_temporal_mode must be off or norm_scaled");
  }
  if (!std::isfinite(xi_kernel_coupling_)) throw std::runtime_error("tpf_4d_xi_kernel_coupling must be finite");
  if (!std::isfinite(xi_kernel_beta_power_) || xi_kernel_beta_power_ < 0.0) {
    throw std::runtime_error("tpf_4d_xi_kernel_beta_power must be finite and >= 0");
  }
  if (!std::isfinite(xi_kernel_metric_min_) || xi_kernel_metric_min_ <= 0.0) {
    throw std::runtime_error("tpf_4d_xi_kernel_metric_min must be finite and > 0");
  }
  if (!std::isfinite(xi_kernel_metric_max_) || xi_kernel_metric_max_ < xi_kernel_metric_min_) {
    throw std::runtime_error("tpf_4d_xi_kernel_metric_max must be finite and >= tpf_4d_xi_kernel_metric_min");
  }
  if (!std::isfinite(xi_temporal_coupling_)) throw std::runtime_error("tpf_4d_xi_temporal_coupling must be finite");
  if (!std::isfinite(xi_source_speed_x_) || !std::isfinite(xi_source_speed_y_) || !std::isfinite(xi_source_speed_z_)) {
    throw std::runtime_error("tpf_4d_xi_source_speed_x/y/z must all be finite");
  }
}

void TPFCorePackage::compute_xi_kernel_deformed_accelerations(const State& state,
                                                              double bh_mass,
                                                              double softening,
                                                              bool star_star,
                                                              std::vector<double>& ax,
                                                              std::vector<double>& ay) const {
  validate_xi_kernel_runtime_config();
  const int n = state.n();
  ax.assign(static_cast<std::size_t>(n), 0.0);
  ay.assign(static_cast<std::size_t>(n), 0.0);
  const double eps = (source_softening_ > 0.0) ? source_softening_ : softening;
  const std::uint64_t bh_pairs = (bh_mass > 0.0) ? static_cast<std::uint64_t>(std::max(0, n)) : 0;
  const std::uint64_t star_pairs = star_star ? static_cast<std::uint64_t>(std::max(0, n)) * static_cast<std::uint64_t>(std::max(0, n - 1)) : 0;
  const std::uint64_t total_pairs = bh_pairs + star_pairs;
  xi_runtime_counters_.xi_last_call_pair_evaluations = total_pairs;
  xi_runtime_counters_.xi_total_pair_evaluations += total_pairs;

  std::vector<XiKernelRuntimeSource> star_sources;
  if (star_star) {
    star_sources.reserve(static_cast<std::size_t>(std::max(0, n)));
    for (int j = 0; j < n; ++j) {
      XiKernelRuntimeSource s;
      s.x = state.x[j];
      s.y = state.y[j];
      s.vx = state.vx[j];
      s.vy = state.vy[j];
      s.mass = state.mass[j];
      star_sources.push_back(s);
    }
  }

  auto accumulate_source = [&](int i, const XiKernelRuntimeSource& src, double& ax_i, double& ay_i) {
    const double dx = state.x[i] - src.x;
    const double dy = state.y[i] - src.y;
    const double dz = 0.0;
    const double vx_rel = state.vx[i] - src.vx;
    const double vy_rel = state.vy[i] - src.vy;
    const double vz_rel = 0.0 - src.vz;
    const double v_rel_norm = std::sqrt(vx_rel * vx_rel + vy_rel * vy_rel + vz_rel * vz_rel);
    const double beta_rel = v_rel_norm / C_SI_LIGHT;
    if (!std::isfinite(beta_rel)) throw std::runtime_error("non-finite beta_rel in Xi kernel deformation");
    if (beta_rel >= 1.0) throw std::runtime_error("beta_rel must be < 1.0 in Xi kernel deformation");
    const tpfcore::XiWakeKinematics wake =
        tpfcore::compute_xi_wake_kinematics(dx, dy, dz, vx_rel, vy_rel, vz_rel, C_SI_LIGHT,
                                            xi_kernel_mode_ == "metric_transverse_wake");
    const double beta_effective = wake.beta_effective;
    const double beta_for_gamma = std::min(beta_effective, 1.0 - 1.0e-12);
    const double gamma_rel = 1.0 / std::sqrt(1.0 - beta_for_gamma * beta_for_gamma);
    double factor_raw = 0.0;
    if (xi_kernel_factor_mode_ == "beta_power") factor_raw = xi_kernel_coupling_ * std::pow(beta_effective, xi_kernel_beta_power_);
    else factor_raw = xi_kernel_coupling_ * (gamma_rel - 1.0);
    double metric_scale = 1.0;
    if (xi_kernel_mode_ == "metric_transverse_wake" || xi_kernel_mode_ == "metric_transverse_continuous") {
      constexpr double k_factor_floor_eps = 1.0e-12;
      const double denom_floor = -1.0 + k_factor_floor_eps;
      const double guarded_factor = std::max(factor_raw, denom_floor);
      metric_scale = std::max(xi_kernel_metric_min_,
                              std::min(xi_kernel_metric_max_, 1.0 / (1.0 + guarded_factor)));
    } else {
      metric_scale = std::max(xi_kernel_metric_min_, std::min(xi_kernel_metric_max_, 1.0 + factor_raw));
    }
    constexpr double k_runtime_pair_tiny = 1.0e-30;
    const double r2 = dx * dx + dy * dy + dz * dz;
    if (!std::isfinite(r2) || r2 <= k_runtime_pair_tiny) return;
    const double r = std::sqrt(r2);
    const double inv_r3 = 1.0 / (r2 * r);
    double r_sq_basis = r2;
    double xi_sx = src.mass * dx * inv_r3;
    double xi_sy = src.mass * dy * inv_r3;
    if (xi_kernel_mode_ == "off") {
    } else if (xi_kernel_mode_ == "scalar_beta") {
      const double scalar = 1.0 + factor_raw;
      xi_sx *= scalar;
      xi_sy *= scalar;
    } else {
      double nx = 0.0, ny = 0.0, nz = 0.0;
      double n_norm = 0.0;
      if (xi_kernel_mode_ == "metric_radial") {
        n_norm = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (n_norm <= 1.0e-30) throw std::runtime_error("near-source invalid radial direction in metric_radial Xi kernel mode");
        nx = dx / n_norm; ny = dy / n_norm; nz = dz / n_norm;
      } else if (xi_kernel_mode_ == "metric_transverse_wake" || xi_kernel_mode_ == "metric_transverse_continuous") {
        n_norm = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (n_norm <= k_runtime_pair_tiny) return;
        nx = dx / n_norm;
        ny = dy / n_norm;
        nz = dz / n_norm;
      } else {
        n_norm = v_rel_norm;
        if (n_norm > 1.0e-30) {
          nx = vx_rel / n_norm; ny = vy_rel / n_norm; nz = vz_rel / n_norm;
        }
      }
      const double alpha = (n_norm > 1.0e-30) ? (metric_scale - 1.0) : 0.0;
      const double nd = nx * dx + ny * dy + nz * dz;
      const double gx = dx + alpha * nx * nd;
      const double gy = dy + alpha * ny * nd;
      const double gz = dz + alpha * nz * nd;
      const double r_eff2 = dx * gx + dy * gy + dz * gz;
      if (std::isfinite(r_eff2) && r_eff2 > k_runtime_pair_tiny) {
        const double r_eff = std::sqrt(r_eff2);
        const double inv_r_eff3 = 1.0 / (r_eff2 * r_eff);
        r_sq_basis = r_eff2;
        xi_sx = src.mass * gx * inv_r_eff3;
        xi_sy = src.mass * gy * inv_r_eff3;
      }
    }
    double dax = -xi_motion_readout_scale_ * xi_sx;
    double day = -xi_motion_readout_scale_ * xi_sy;
    const double softening_scale = plummer_softening_scale(r_sq_basis, eps);
    dax *= softening_scale;
    day *= softening_scale;
    ax_i += dax;
    ay_i += day;
  };

  XiKernelRuntimeSource bh_src;
  bh_src.mass = bh_mass;
  bh_src.vx = xi_source_speed_x_;
  bh_src.vy = xi_source_speed_y_;
  bh_src.vz = xi_source_speed_z_;

  for (int i = 0; i < n; ++i) {
    double ax_i = 0.0;
    double ay_i = 0.0;
    if (bh_mass > 0.0) {
      accumulate_source(i, bh_src, ax_i, ay_i);
    }
    if (star_star) {
      for (int j = 0; j < n; ++j) {
        if (j == i) continue;
        accumulate_source(i, star_sources[static_cast<std::size_t>(j)], ax_i, ay_i);
      }
    }
    ax[i] = ax_i;
    ay[i] = ay_i;

    if (bh_mass > 0.0 && !star_star && xi_motion_readout_scale_ > 0.0) {
      const double r = std::sqrt(state.x[i] * state.x[i] + state.y[i] * state.y[i]);
      const double a_norm = std::sqrt(ax[i] * ax[i] + ay[i] * ay[i]);
      if (r > 1.0e-30 && a_norm > 1.0e-30) {
        const double inward_x = -state.x[i] / r;
        const double inward_y = -state.y[i] / r;
        const double radial_alignment_to_bh_inward = (ax[i] * inward_x + ay[i] * inward_y) / a_norm;
        if (!std::isfinite(radial_alignment_to_bh_inward) || radial_alignment_to_bh_inward < -1.0e-12) {
          throw std::runtime_error("xi_kernel_deformed BH-only sign audit failed: acceleration points outward");
        }
      }
    }
  }
}

void TPFCorePackage::compute_accelerations(const State& state,
                                            double bh_mass,
                                            double softening,
                                            bool star_star,
                                            std::vector<double>& ax,
                                            std::vector<double>& ay) const {
  if (tpf_dynamics_mode_ == "tpf_xi_theta_v1") {
    xi_runtime_counters_.xi_last_call_pair_evaluations = 0;
    compute_xi_kernel_deformed_accelerations(state, bh_mass, softening, star_star, ax, ay);
    last_pipeline_ = AccelPipelineStats{};
    return;
  }
  throw std::runtime_error("Unsupported tpf_dynamics_mode for this branch; expected tpf_xi_theta_v1.");

  const int n = state.n();

  static bool has_branch_audited = false;
  if (!has_branch_audited && n > 0) {
    has_branch_audited = true;
    const bool vdsg_active = (vdsg_coupling_ != 0.0);
    const bool readout_enabled = provisional_readout_;

    std::cerr << "\n========== ACTIVE PHYSICS LEDGER (one-time branch audit) ==========\n";
    std::cerr << "TPF readout baseline: provisional readout closures per tpfcore_readout_mode "
                 "(tensor / derived radial / experimental)\n";
    std::cerr << "VDSG velocity modifier: " << (vdsg_active ? "ACTIVE" : "INACTIVE")
              << "  (tpf_vdsg_coupling " << (vdsg_active ? "!=" : "==") << " 0; additive excess "
                 "a_N*(doppler_scale-1), doppler_scale=1+lambda_eff*|v_rel|/c)\n";
    std::cerr << "tpf_dynamics_mode: legacy_readout (provisional readout + optional VDSG)\n";
    std::cerr << "provisional_readout GATE: " << (readout_enabled ? "ENABLED" : "DISABLED")
              << "  (tpfcore_enable_provisional_readout; required for legacy_readout accelerations)\n";
    std::cerr << "Global |a| shunt (velocity cap): "
              << (shunt_enable_ ? "ENABLED" : "OFF")
              << "  (tpf_global_accel_shunt_enable; independent of tpf_vdsg_coupling; clean λ=0 baseline when off)\n";
    if (shunt_enable_) {
      std::cerr << "  shunt cap = tpf_global_accel_shunt_fraction * |v|/dt  (fraction=" << shunt_fraction_
                << ")\n";
    }
    std::cerr << std::scientific << std::setprecision(16);
    std::cerr << "tpf_vdsg_coupling = " << vdsg_coupling_ << "\n";
    std::cerr << "tpf_kappa (external config key; mapped to internal direct_tpf kappa) = " << kappa_ << "\n";
    std::cerr << "tpfcore_readout_mode = " << readout_mode_ << "\n";
    std::cerr << "VDSG mass normalization (PROVISIONAL / exploratory heuristic): "
                 "lambda_eff = lambda * log10(max(M_ref,eps_kg)) / log10(max(M_src,eps_kg)); "
                 "M_ref = tpf_vdsg_mass_baseline_kg if >0 else star_mass; "
                 "if M_ref<=0 => lambda_eff=lambda (identity). Single tpf_vdsg_coupling only.\n";
    std::cerr << "tpf_vdsg_mass_baseline_resolved_kg (runtime M_ref) = " << vdsg_mass_baseline_resolved_kg_
              << "\n";
    if (vdsg_active && vdsg_mass_baseline_resolved_kg_ > 0.0) {
      const double le_bh =
          tpfcore::vdsg_effective_coupling(vdsg_coupling_, bh_mass, vdsg_mass_baseline_resolved_kg_);
      const int idx_star = (n > 1) ? 1 : 0;
      const double m_st = state.mass[idx_star];
      const double le_st =
          tpfcore::vdsg_effective_coupling(vdsg_coupling_, m_st, vdsg_mass_baseline_resolved_kg_);
      std::cerr << "lambda_eff sample (gravitational source = BH, M_src = M_BH) = " << le_bh << "\n";
      std::cerr << "lambda_eff sample (gravitational source = star, M_src = mass[" << idx_star << "]) = " << le_st
                << "\n";
    } else if (vdsg_active) {
      std::cerr << "lambda_eff mass scaling: OFF (M_ref <= 0; set tpf_vdsg_mass_baseline_kg or star_mass > 0).\n";
    }

    const double eps_soft = (source_softening_ > 0.0) ? source_softening_ : softening;
    const double G_si = tpfcore::TPF_G_SI;

    if (bh_mass > 0.0) {
      const double xi = state.x[0];
      const double yi = state.y[0];
      const double r_sq = xi * xi + yi * yi;
      const double eps_sq = softening * softening;
      const double denom = r_sq + eps_sq;
      const double r_mag = std::sqrt(denom);
      if (r_mag > 1e-300) {
        const double v_mag =
            std::sqrt(state.vx[0] * state.vx[0] + state.vy[0] * state.vy[0] + 1e-300);
        const double beta_sample = v_mag / C_SI_LIGHT;
        const double lambda_eff_bh =
            tpfcore::vdsg_effective_coupling(vdsg_coupling_, bh_mass, vdsg_mass_baseline_resolved_kg_);
        const double doppler_scale = 1.0 + lambda_eff_bh * beta_sample;
        const double a_newtonian = G_si * bh_mass / denom;
        const double a_VDSG_modifier = a_newtonian * (doppler_scale - 1.0);

        std::cerr << "--- Sample: Star 0 vs BH (SI magnitudes) ---\n";
        std::cerr << "a_newtonian = " << a_newtonian << "  (G*M_BH/(r^2+eps^2))\n";
        std::cerr << "a_VDSG_modifier = " << a_VDSG_modifier
                  << "  (= a_newtonian * (doppler_scale - 1); doppler_scale = 1 + lambda_eff*|v|/c = "
                  << doppler_scale;
        if (vdsg_active)
          std::cerr << "; lambda_eff(BH source) = " << lambda_eff_bh;
        std::cerr << ")\n";

        double ax_r = 0.0;
        double ay_r = 0.0;
        if (tpfcore::is_derived_tpf_radial_readout_mode(readout_mode_)) {
          tpfcore::TpfRadialGravityProfile profile =
              tpfcore::build_tpf_gravity_profile(state, bh_mass, derived_poisson_cfg_, eps_soft);
          tpfcore::compute_provisional_readout_acceleration(
              state, 0, bh_mass, star_star, softening, source_softening_, readout_mode_, readout_scale_,
              theta_tt_scale_, theta_tr_scale_, ax_r, ay_r, &derived_poisson_cfg_, &profile);
        } else {
          tpfcore::compute_provisional_readout_acceleration(
              state, 0, bh_mass, star_star, softening, source_softening_, readout_mode_, readout_scale_,
              theta_tt_scale_, theta_tr_scale_, ax_r, ay_r, nullptr, nullptr);
        }
        const double a_tpf_readout = std::hypot(ax_r, ay_r);
        std::cerr << "a_readout_baseline (sample |a|) = " << a_tpf_readout << "\n";
      } else {
        std::cerr << "--- Sample: Star 0 vs BH skipped (r_mag too small) ---\n";
      }
    } else {
      std::cerr << "--- Sample: Star 0 vs BH skipped (bh_mass <= 0) ---\n";
    }
    std::cerr << "===================================================================\n\n" << std::flush;
  }

  eval_accel_pipeline(state, bh_mass, softening, star_star, ax, ay, &last_pipeline_);
}

void TPFCorePackage::write_readout_debug(const std::vector<Snapshot>& snapshots,
                                         const Config& config,
                                         const std::string& output_dir) const {
  tpfcore::TPFCoreParams params = build_params(config, output_dir);
  if (config.tpf_dynamics_mode == "xi_kernel_deformed" && !config.tpf_xi_kernel_dump_field_diagnostics) return;
  if (!provisional_readout_ || !params.tpfcore_dump_readout_debug || snapshots.empty())
    return;
  tpfcore::write_readout_debug_csv(snapshots, params.output_dir,
                                   params.softening, params.bh_mass, params.enable_star_star_gravity,
                                   source_softening_,
                                   readout_mode_, readout_scale_, theta_tt_scale_, theta_tr_scale_,
                                   derived_poisson_cfg_);
}

namespace {

/* Same radial extraction as force_compare: |ax*rx + ay*ry| with r_eff = sqrt(r^2+eps^2), (rx,ry) = (x,y)/r_eff. */
double newton_radial_magnitude_on_axis(double r, double ax, double ay, double eps_sq) {
  (void)ay;
  double r_eff_sq = r * r + eps_sq;
  double r_eff = std::sqrt(r_eff_sq);
  if (r_eff < 1e-30) return 0.0;
  double rx = r / r_eff;
  return std::abs(ax * rx);
}

}  // namespace

void TPFCorePackage::run_weak_field_calibration(const Config& config, const std::string& output_dir) {
  if (!provisional_readout_) return;

  PhysicsPackage* newton = get_physics_package("Newtonian");
  if (!newton) return;

  tpfcore::TPFCoreParams params = build_params(config, output_dir);
  double r_min = params.tpfcore_probe_radius_min;
  double r_max = params.tpfcore_probe_radius_max;
  int n_samples = params.tpfcore_probe_samples;
  double M = params.bh_mass;
  double softening = params.softening;
  const double eps_sq = softening * softening;
  if (r_min >= r_max || n_samples < 2) return;

  /* Closure candidates to compare (always run both and write comparison). */
  static const char* const CALIBRATION_MODES[] = {"tr_coherence_readout", "experimental_radial_r_scaling"};
  const int num_modes = static_cast<int>(sizeof(CALIBRATION_MODES) / sizeof(CALIBRATION_MODES[0]));

  struct ModeResult {
    std::string mode;
    std::vector<double> r_vals, a_tpf_vals, ratio_vals;
    double K_eff;
    double ratio_min, ratio_max, ratio_spread;
    bool one_constant_sufficient;
  };
  std::vector<ModeResult> results;
  results.reserve(num_modes);

  State state;
  state.resize(1);
  state.mass[0] = params.star_mass;

  /* Precompute Newtonian benchmark at each radius (same for all modes; from simulator path). */
  std::vector<double> r_vals_common, a_newton_vals;
  r_vals_common.reserve(n_samples);
  a_newton_vals.reserve(n_samples);
  for (int i = 0; i < n_samples; ++i) {
    double frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
    double r = r_min + frac * (r_max - r_min);
    state.x[0] = r;
    state.y[0] = 0.0;
    state.vx[0] = state.vy[0] = 0.0;
    std::vector<double> ax_n, ay_n;
    newton->compute_accelerations(state, M, softening, false, ax_n, ay_n);
    double a_newton = newton_radial_magnitude_on_axis(r, ax_n[0], ay_n[0], eps_sq);
    r_vals_common.push_back(r);
    a_newton_vals.push_back(a_newton);
  }

  std::string current_mode = readout_mode_;
  double current_scale = readout_scale_;

  for (int m = 0; m < num_modes; ++m) {
    readout_mode_ = CALIBRATION_MODES[m];
    readout_scale_ = 1.0;
    ModeResult res;
    res.mode = readout_mode_;
    res.r_vals = r_vals_common;
    res.a_tpf_vals.resize(n_samples);
    res.ratio_vals.resize(n_samples);

    for (int i = 0; i < n_samples; ++i) {
      double r = res.r_vals[i];
      state.x[0] = r;
      state.y[0] = 0.0;
      state.vx[0] = state.vy[0] = 0.0;
      std::vector<double> ax, ay;
      compute_accelerations(state, M, softening, false, ax, ay);
      double a_tpf = std::abs(ax[0]);
      double a_newton = a_newton_vals[i];
      res.a_tpf_vals[i] = a_tpf;
      res.ratio_vals[i] = (a_newton > 1e-300) ? (a_tpf / a_newton) : 0.0;
    }

    double sum_prod = 0.0, sum_tpf_sq = 0.0;
    for (int k = 0; k < n_samples; ++k) {
      sum_prod += res.a_tpf_vals[k] * a_newton_vals[k];
      sum_tpf_sq += res.a_tpf_vals[k] * res.a_tpf_vals[k];
    }
    res.K_eff = (sum_tpf_sq > 1e-300) ? (sum_prod / sum_tpf_sq) : 0.0;

    res.ratio_min = res.ratio_vals[0];
    res.ratio_max = res.ratio_vals[0];
    for (double q : res.ratio_vals) {
      if (q < res.ratio_min) res.ratio_min = q;
      if (q > res.ratio_max) res.ratio_max = q;
    }
    res.ratio_spread = res.ratio_max - res.ratio_min;
    res.one_constant_sufficient = (res.ratio_spread <= 0.3);
    results.push_back(res);
  }

  readout_mode_ = current_mode;
  readout_scale_ = current_scale;

  /* Per-mode CSV/txt: use config mode if it was one of the two, else first mode. */
  int current_idx = 0;
  for (int m = 0; m < num_modes; ++m) {
    if (results[m].mode == current_mode) { current_idx = m; break; }
  }

  std::string csv_path = params.output_dir + "/tpf_weak_field_calibration.csv";
  std::ofstream csv(csv_path);
  if (csv) {
    csv << "radius,a_tpf,a_newton,ratio,mode\n";
    const ModeResult& res = results[current_idx];
    for (size_t k = 0; k < res.r_vals.size(); ++k)
      csv << std::scientific << res.r_vals[k] << "," << res.a_tpf_vals[k] << "," << a_newton_vals[k] << "," << res.ratio_vals[k] << "," << res.mode << "\n";
  }

  std::string txt_path = params.output_dir + "/tpf_weak_field_calibration.txt";
  std::ofstream txt(txt_path);
  if (txt) {
    const ModeResult& cur = results[current_idx];
    txt << "TPFCore correspondence calibration (diagnostic / provisional readout sector)\n";
    txt << "Newtonian benchmark: from Newtonian package (same path as simulator and force_compare).\n";
    txt << "Radial extraction: |a_rad| = |ax*rx + ay*ry| with r_eff = sqrt(r^2+eps^2), (rx,ry) = (r,0)/r_eff on x-axis.\n";
    txt << "Previous scalar formula M/(r^2+eps^2)^(3/2) was wrong by factor r; value 0.2046442 is INVALIDATED.\n\n";
    txt << "Source mass M = " << std::scientific << M << ", softening = " << softening << "\n";
    txt << "Probe: +x axis, r in [" << r_min << ", " << r_max << "], n = " << n_samples << "\n\n";
    txt << "Current mode: " << cur.mode << "\n";
    txt << "  K_eff = " << std::scientific << cur.K_eff << "\n";
    txt << "  Ratio a_tpf/a_newton: min = " << cur.ratio_min << ", max = " << cur.ratio_max << ", spread = " << cur.ratio_spread << "\n";
    txt << "  One constant sufficient: " << (cur.one_constant_sufficient ? "yes" : "no") << "\n\n";
    const char* interpretation = cur.ratio_spread > 0.3 ? "radius-inconsistent (one constant scale factor not sufficient)"
      : (cur.ratio_max < 0.9 ? "underpowered relative to Newtonian"
         : (cur.ratio_min > 1.1 ? "overpowered relative to Newtonian"
            : (cur.ratio_spread < 0.2 ? "roughly matched by one constant scale factor (low spread)"
               : "roughly matched by one constant scale factor (moderate spread)")));
    txt << "Interpretation: " << interpretation << "\n\n";
    txt << "--- Table ---\nradius\ta_tpf\ta_newton\tratio\n";
    for (size_t k = 0; k < cur.r_vals.size(); ++k)
      txt << std::scientific << cur.r_vals[k] << "\t" << cur.a_tpf_vals[k] << "\t" << a_newton_vals[k] << "\t" << cur.ratio_vals[k] << "\n";
    txt << "\nSee tpf_weak_field_calibration_comparison.txt for closure/correspondence comparison details.\n";
  }

  std::string comp_path = params.output_dir + "/tpf_weak_field_calibration_comparison.txt";
  std::ofstream comp(comp_path);
  if (comp) {
    comp << "TPFCore correspondence calibration: comparison of closure candidates\n";
    comp << "Benchmark: Newtonian package, same radial extraction as simulator/force_compare.\n";
    comp << "Previous 0.2046442 was invalidated by benchmark formula mismatch (factor r).\n\n";
    comp << "--- Per-mode results ---\n\n";
    for (size_t m = 0; m < results.size(); ++m) {
      const ModeResult& res = results[m];
      comp << "Mode: " << res.mode << "\n";
      comp << "  K_eff = " << std::scientific << res.K_eff << "\n";
      comp << "  ratio min = " << res.ratio_min << ", max = " << res.ratio_max << ", spread = " << res.ratio_spread << "\n";
      comp << "  One constant sufficient: " << (res.one_constant_sufficient ? "yes" : "no") << "\n\n";
    }
    comp << "--- Correspondence-shape verdict ---\n";
    bool tr_ok = results[0].one_constant_sufficient;
    bool exp_ok = results[1].one_constant_sufficient;
    double spread_tr = results[0].ratio_spread;
    double spread_exp = results[1].ratio_spread;
    if (tr_ok && !exp_ok)
      comp << "  tr_coherence_readout has better correspondence/Newtonian shape match (one constant sufficient; experimental_radial_r_scaling ratio varies more with r).\n";
    else if (!tr_ok && exp_ok)
      comp << "  experimental_radial_r_scaling has better correspondence/Newtonian shape match (one constant sufficient; tr_coherence_readout ratio varies more with r).\n";
    else if (tr_ok && exp_ok)
      comp << "  Both closures have acceptable shape (one constant sufficient). Lower spread wins; spread tr=" << std::scientific << spread_tr << " exp=" << spread_exp << ".\n";
    else {
      comp << "  Neither closure has one constant sufficient over the probe range; lower ratio spread is better correspondence shape.\n";
      if (spread_tr < spread_exp)
        comp << "  tr_coherence_readout is the better correspondence/Newtonian shape match (spread " << std::scientific << spread_tr << " < " << spread_exp << ").\n";
      else if (spread_exp < spread_tr)
        comp << "  experimental_radial_r_scaling is the better correspondence/Newtonian shape match (spread " << std::scientific << spread_exp << " < " << spread_tr << ").\n";
      else
        comp << "  Tied on spread.\n";
    }
  }

  std::string comp_csv_path = params.output_dir + "/tpf_weak_field_calibration_comparison.csv";
  std::ofstream comp_csv(comp_csv_path);
  if (comp_csv) {
    comp_csv << "mode,K_eff,ratio_min,ratio_max,ratio_spread,one_constant_sufficient\n";
    for (size_t m = 0; m < results.size(); ++m) {
      const ModeResult& res = results[m];
      comp_csv << res.mode << "," << std::scientific << res.K_eff << "," << res.ratio_min << "," << res.ratio_max << "," << res.ratio_spread << "," << (res.one_constant_sufficient ? "1" : "0") << "\n";
    }
  }
}

void TPFCorePackage::run_single_source_inspect(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  TPFCoreParams params = build_params(config, output_dir);
  double r_min = params.tpfcore_probe_radius_min;
  double r_max = params.tpfcore_probe_radius_max;
  int n_samples = params.tpfcore_probe_samples;
  double eps = params.effective_source_softening;
  double m = params.bh_mass;

  if (r_min >= r_max || n_samples < 2) return;

  std::vector<double> radii, x_v, y_v, xi_x_v, xi_y_v;
  std::vector<double> theta_xx_v, theta_xy_v, theta_yy_v, theta_trace_v, invariant_I_v;
  std::vector<double> residual_x_v, residual_y_v, residual_norm_v;

  double max_theta_xy_abs = 0.0;
  bool theta_xx_vs_yy_differ = false;
  double max_residual_x_abs = 0.0, max_residual_y_abs = 0.0, max_residual_norm = 0.0;
  double max_theta_norm = 0.0;
  double max_invariant_I = 0.0;

  for (int i = 0; i < n_samples; ++i) {
    double frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
    double r = r_min + frac * (r_max - r_min);
    double x = r, y = 0.0;

    FieldAtPoint field = evaluate_provisional_field_single_source(0, 0, m, x, y, eps);

    radii.push_back(r);
    x_v.push_back(x);
    y_v.push_back(y);
    xi_x_v.push_back(field.xi.x);
    xi_y_v.push_back(field.xi.y);
    theta_xx_v.push_back(field.theta.xx);
    theta_xy_v.push_back(field.theta.xy);
    theta_yy_v.push_back(field.theta.yy);
    theta_trace_v.push_back(field.theta_trace);
    invariant_I_v.push_back(field.invariant_I);
    residual_x_v.push_back(field.residual.x);
    residual_y_v.push_back(field.residual.y);
    residual_norm_v.push_back(field.residual.norm());

    double txy = std::abs(field.theta.xy);
    if (txy > max_theta_xy_abs) max_theta_xy_abs = txy;
    if (std::abs(field.theta.xx - field.theta.yy) > 1e-14) theta_xx_vs_yy_differ = true;
    double rax = std::abs(field.residual.x), ray = std::abs(field.residual.y), rn = field.residual.norm();
    if (rax > max_residual_x_abs) max_residual_x_abs = rax;
    if (ray > max_residual_y_abs) max_residual_y_abs = ray;
    if (rn > max_residual_norm) max_residual_norm = rn;
    double tn = tpfcore::theta_frobenius_norm(field.theta);
    if (tn > max_theta_norm) max_theta_norm = tn;
    if (std::abs(field.invariant_I) > max_invariant_I) max_invariant_I = std::abs(field.invariant_I);
  }

  bool residual_y_near_zero = (max_residual_y_abs < 1e-10);

  if (params.tpfcore_dump_theta_profile) {
    std::ofstream f(params.output_dir + "/theta_profile.csv");
    if (f) {
      f << "radius,x,y,xi_x,xi_y,theta_xx,theta_xy,theta_yy,theta_trace,invariant_I,residual_x,residual_y,residual_norm\n";
      for (size_t i = 0; i < radii.size(); ++i) {
        f << std::scientific << radii[i] << "," << x_v[i] << "," << y_v[i] << ","
          << xi_x_v[i] << "," << xi_y_v[i] << ","
          << theta_xx_v[i] << "," << theta_xy_v[i] << "," << theta_yy_v[i] << ","
          << theta_trace_v[i] << "," << invariant_I_v[i] << ","
          << residual_x_v[i] << "," << residual_y_v[i] << "," << residual_norm_v[i] << "\n";
      }
    }
  }

  if (params.tpfcore_dump_invariant_profile) {
    std::ofstream f(params.output_dir + "/invariant_profile.csv");
    if (f) {
      f << "radius,invariant_I,theta_trace,residual_norm\n";
      for (size_t i = 0; i < radii.size(); ++i) {
        f << std::scientific << radii[i] << "," << invariant_I_v[i] << "," << theta_trace_v[i]
          << "," << residual_norm_v[i] << "\n";
      }
    }
  }

  {
    std::ofstream f(params.output_dir + "/field_summary.txt");
    if (f) {
      f << "TPFCore single-source inspection\n";
      f << "--- Parameter classification for this run ---\n";
      f << "  fixed_theory:        lambda=" << LAMBDA_4D << " (4D; not tunable)\n";
      f << "  numerical_reg:       source softening eps=" << std::scientific << eps << "\n";
      f << "  inspection:          probe r in [" << r_min << ", " << r_max << "], n=" << n_samples << "\n";
      f << "---\n";
      f << "Ansatz: 3D Phi=-M/R with R^2=dx^2+dy^2+eps^2 on z=0; Xi=(partial_x Phi, partial_y Phi); Theta=Hess_3D(Phi) (xz=yz=0 on plane)\n";
      f << "Source: (0,0), mass=" << m << "\n";
      f << "Probe: +x axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << "\n";
      f << "Source softening eps=" << eps << " (numerical regularization)\n";
      f << "Provisional point-source ansatz: yes\n";
      f << "Lambda: " << LAMBDA_4D << " (fixed, 4D)\n";
      f << "\n--- Regime diagnostics (reporting only; no equation change) ---\n";
      f << "  Field intensity (Theta Frobenius): max = " << std::scientific << max_theta_norm << "\n";
      f << "  Invariant I: max|I| = " << max_invariant_I << "\n";
      f << "  Residual: available (analytic, single-source)\n";
      f << "  Max residual norm = " << max_residual_norm << "\n";
      f << "  Regime (at max-intensity probe): " << tpfcore::regime_label_from_theta_norm(max_theta_norm) << "\n";
      f << "  (Thresholds: low-intensity < " << tpfcore::THETA_NORM_LOW_MAX
        << ", transitional < " << tpfcore::THETA_NORM_TRANSITIONAL_MAX << ", else high-intensity; heuristic.)\n";
      f << "---\n";
      f << "\nXi/Theta configuration-equation residual (full 3D spatial divergence; i=x,y,z; nu=x,y):\n";
      f << "  max|residual_x|=" << std::scientific << max_residual_x_abs << "\n";
      f << "  max|residual_y|=" << max_residual_y_abs << "\n";
      f << "  max residual_norm=" << max_residual_norm << "\n";
      f << "  (Vanishes as eps->0 away from source; O(eps^2) with softening.)\n";
      f << "  residual_y near zero on +x axis (y=0 symmetry): " << (residual_y_near_zero ? "yes [OK]\n" : "no [check]\n");
      f << "\nSymmetry expectations on +x axis:\n";
      f << "  theta_xy should be zero (y=0 symmetry): max|theta_xy|=" << max_theta_xy_abs;
      f << (max_theta_xy_abs < 1e-10 ? " [OK]\n" : " [check]\n");
      f << "  theta_xx and theta_yy should differ (radial vs transverse): " << (theta_xx_vs_yy_differ ? "yes [OK]\n" : "no [check]\n");
      f << "  invariant_I should decay smoothly with radius\n";
    }
  }

  if (params.tpf_xi_constraint_exterior_inspect) {
    XiConstraintExteriorParams solve_params;
    solve_params.n = params.tpf_xi_constraint_grid_n;
    solve_params.L = params.tpf_xi_constraint_half_extent;
    solve_params.r_inner = params.tpf_xi_constraint_inner_radius;
    solve_params.source_mass = m;
    solve_params.softening = eps;
    solve_params.max_iterations = params.tpf_xi_constraint_max_iters;
    solve_params.tolerance = params.tpf_xi_constraint_tol;

    XiConstraintExteriorSolveResult solve_result = solve_xi_constraint_exterior(solve_params);
    XiConstraintExteriorStats stats = compute_xi_constraint_exterior_stats(solve_params, solve_result);

    std::ofstream grid_csv(params.output_dir + "/xi_constraint_exterior_grid.csv");
    if (grid_csv) {
      grid_csv << "x,y,is_masked,is_boundary_pinned,is_free_cell,xi_ansatz_x,xi_ansatz_y,xi_solved_x,xi_solved_y,delta_xi_x,delta_xi_y,residual_x,residual_y,residual_norm\n";
      for (int j = 0; j < solve_result.grid.ny; ++j) {
        for (int i = 0; i < solve_result.grid.nx; ++i) {
          const int idx = solve_result.grid.index(i, j);
          const double x = solve_result.grid.x_at(i);
          const double y = solve_result.grid.y_at(j);
          const bool is_masked = (solve_result.grid.is_exterior[idx] == 0);
          const bool is_pinned = (solve_result.grid.is_exterior[idx] != 0) && (solve_result.grid.is_pinned[idx] != 0);
          const bool is_free = (solve_result.grid.is_exterior[idx] != 0) && (solve_result.grid.is_pinned[idx] == 0);
          const Xi2D xi_ansatz = provisional_point_source_field(0.0, 0.0, m, x, y, eps).xi;
          const double xi_solved_x = solve_result.grid.xi_x[idx];
          const double xi_solved_y = solve_result.grid.xi_y[idx];
          const double delta_x = xi_solved_x - xi_ansatz.x;
          const double delta_y = xi_solved_y - xi_ansatz.y;
          const PlanarConstraintResidual2D res = solve_result.residuals[idx];
          grid_csv << std::scientific << x << "," << y << ","
                   << (is_masked ? 1 : 0) << "," << (is_pinned ? 1 : 0) << "," << (is_free ? 1 : 0) << ","
                   << xi_ansatz.x << "," << xi_ansatz.y << "," << xi_solved_x << "," << xi_solved_y << ","
                   << delta_x << "," << delta_y << "," << res.Rx << "," << res.Ry << "," << res.norm << "\n";
        }
      }
    }

    std::ofstream summary(params.output_dir + "/xi_constraint_exterior_summary.txt");
    if (summary) {
      summary << "TPFCore Xi constraint exterior inspection summary\n";
      summary << "solver type: experimental planar Xi solve\n";
      summary << "grid size: " << solve_params.n << " x " << solve_params.n << "\n";
      summary << "half extent: " << std::scientific << solve_params.L << "\n";
      summary << "inner radius: " << solve_params.r_inner << "\n";
      summary << "iteration count: " << solve_result.iterations << "\n";
      summary << "convergence reached: " << (solve_result.converged ? "yes" : "no") << "\n";
      summary << "number of masked cells: " << stats.n_masked << "\n";
      summary << "number of pinned cells: " << stats.n_pinned << "\n";
      summary << "number of free cells: " << stats.n_free << "\n";
      summary << "initial max residual norm over free cells: " << stats.initial_max_residual_free << "\n";
      summary << "final max residual norm over free cells: " << stats.final_max_residual_free << "\n";
      summary << "initial max residual norm over all non-masked cells: " << stats.initial_max_residual_nonmasked << "\n";
      summary << "final max residual norm over all non-masked cells: " << stats.final_max_residual_nonmasked << "\n";
      summary << "max |delta Xi| over free cells: " << stats.max_delta_xi_free << "\n";
      summary << "mean |delta Xi| over free cells: " << stats.mean_delta_xi_free << "\n";
      summary << "experimental planar configuration-equation exterior solve\n";
      summary << "ansatz used as boundary data and initial guess\n";
      summary << "ansatz boundary comparison\n";
      summary << "free-cell residual is the headline metric in this report\n";
      summary << "not full Eq. (10) validation\n";
      summary << "no DeltaC included\n";
      summary << "single-source only\n";
    }

    std::ofstream radial(params.output_dir + "/xi_constraint_exterior_radial_profile.csv");
    if (radial) {
      radial << "radius,xi_ansatz_mag,xi_solved_mag,delta_xi_mag\n";
      for (int j = 0; j < solve_result.grid.ny; ++j) {
        for (int i = 0; i < solve_result.grid.nx; ++i) {
          const int idx = solve_result.grid.index(i, j);
          if (!solve_result.grid.is_exterior[idx] || solve_result.grid.is_pinned[idx]) continue;
          const double x = solve_result.grid.x_at(i);
          const double y = solve_result.grid.y_at(j);
          const Xi2D xi_ansatz = provisional_point_source_field(0.0, 0.0, m, x, y, eps).xi;
          const double ansatz_mag = std::hypot(xi_ansatz.x, xi_ansatz.y);
          const double solved_mag = std::hypot(solve_result.grid.xi_x[idx], solve_result.grid.xi_y[idx]);
          radial << std::scientific << std::hypot(x, y) << "," << ansatz_mag << "," << solved_mag
                 << "," << std::fabs(solved_mag - ansatz_mag) << "\n";
        }
      }
    }
  }
}

void TPFCorePackage::run_symmetric_pair_inspect(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  TPFCoreParams params = build_params(config, output_dir);
  double d = params.validation_symmetric_separation;
  double m = params.star_mass;
  double r_min = params.tpfcore_probe_radius_min;
  double r_max = params.tpfcore_probe_radius_max;
  int n_samples = params.tpfcore_probe_samples;
  double eps = params.effective_source_softening;

  if (r_min >= r_max || n_samples < 2) return;

  std::vector<double> x_v, y_v, xi_x_v, xi_y_v;
  std::vector<double> theta_xx_v, theta_xy_v, theta_yy_v, theta_trace_v, invariant_I_v;
  std::vector<double> residual_x_v, residual_y_v, residual_norm_v;
  std::vector<std::string> axis_v;

  double max_residual_x_abs = 0.0, max_residual_y_abs = 0.0, max_residual_norm = 0.0;
  double max_residual_x_abs_x_axis = 0.0, max_residual_y_abs_x_axis = 0.0;
  double max_residual_x_abs_y_axis = 0.0, max_residual_y_abs_y_axis = 0.0;
  double max_theta_norm = 0.0;
  double max_invariant_I = 0.0;

  auto append = [&](const char* ax, double (*px)(double), double (*py)(double)) {
    for (int i = 0; i < n_samples; ++i) {
      double frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
      double r = r_min + frac * (r_max - r_min);
      double x = px(r), y = py(r);

      FieldAtPoint f1 = evaluate_provisional_field_single_source(d, 0, m, x, y, eps);
      FieldAtPoint f2 = evaluate_provisional_field_single_source(-d, 0, m, x, y, eps);
      FieldAtPoint field = add_provisional_fields(f1, f2);

      axis_v.push_back(ax);
      x_v.push_back(x);
      y_v.push_back(y);
      xi_x_v.push_back(field.xi.x);
      xi_y_v.push_back(field.xi.y);
      theta_xx_v.push_back(field.theta.xx);
      theta_xy_v.push_back(field.theta.xy);
      theta_yy_v.push_back(field.theta.yy);
      theta_trace_v.push_back(field.theta_trace);
      invariant_I_v.push_back(field.invariant_I);
      residual_x_v.push_back(field.residual.x);
      residual_y_v.push_back(field.residual.y);
      residual_norm_v.push_back(field.residual.norm());

      double rax = std::abs(field.residual.x), ray = std::abs(field.residual.y), rn = field.residual.norm();
      if (rax > max_residual_x_abs) max_residual_x_abs = rax;
      if (ray > max_residual_y_abs) max_residual_y_abs = ray;
      if (rn > max_residual_norm) max_residual_norm = rn;
      if (ax[0] == 'x') {
        if (rax > max_residual_x_abs_x_axis) max_residual_x_abs_x_axis = rax;
        if (ray > max_residual_y_abs_x_axis) max_residual_y_abs_x_axis = ray;
      } else {
        if (rax > max_residual_x_abs_y_axis) max_residual_x_abs_y_axis = rax;
        if (ray > max_residual_y_abs_y_axis) max_residual_y_abs_y_axis = ray;
      }
      double tn = tpfcore::theta_frobenius_norm(field.theta);
      if (tn > max_theta_norm) max_theta_norm = tn;
      if (std::abs(field.invariant_I) > max_invariant_I) max_invariant_I = std::abs(field.invariant_I);
    }
  };

  append("x", [](double r) { return r; }, [](double) { return 0.0; });
  append("y", [](double) { return 0.0; }, [](double r) { return r; });

  bool residual_y_near_zero_x_axis = (max_residual_y_abs_x_axis < 1e-10);
  bool residual_x_near_zero_y_axis = (max_residual_x_abs_y_axis < 1e-10);

  if (params.tpfcore_dump_theta_profile) {
    std::ofstream f(params.output_dir + "/theta_profile.csv");
    if (f) {
      f << "axis,radius,x,y,xi_x,xi_y,theta_xx,theta_xy,theta_yy,theta_trace,invariant_I,residual_x,residual_y,residual_norm\n";
      for (size_t i = 0; i < axis_v.size(); ++i) {
        double r = (axis_v[i] == "x") ? x_v[i] : y_v[i];
        f << axis_v[i] << "," << std::scientific << r << "," << x_v[i] << "," << y_v[i] << ","
          << xi_x_v[i] << "," << xi_y_v[i] << ","
          << theta_xx_v[i] << "," << theta_xy_v[i] << "," << theta_yy_v[i] << ","
          << theta_trace_v[i] << "," << invariant_I_v[i] << ","
          << residual_x_v[i] << "," << residual_y_v[i] << "," << residual_norm_v[i] << "\n";
      }
    }
  }

  if (params.tpfcore_dump_invariant_profile) {
    std::ofstream f(params.output_dir + "/invariant_profile.csv");
    if (f) {
      f << "axis,radius,invariant_I,theta_trace,residual_norm\n";
      for (size_t i = 0; i < axis_v.size(); ++i) {
        double r = (axis_v[i] == "x") ? x_v[i] : y_v[i];
        f << axis_v[i] << "," << std::scientific << r << "," << invariant_I_v[i] << "," << theta_trace_v[i]
          << "," << residual_norm_v[i] << "\n";
      }
    }
  }

  {
    std::ofstream f(params.output_dir + "/field_summary.txt");
    if (f) {
      f << "TPFCore symmetric-pair inspection\n";
      f << "--- Parameter classification for this run ---\n";
      f << "  fixed_theory:        lambda=" << LAMBDA_4D << " (4D; not tunable)\n";
      f << "  numerical_reg:       source softening eps=" << std::scientific << eps << "\n";
      f << "  inspection:          probe r in [" << r_min << ", " << r_max << "], n=" << n_samples << " per axis\n";
      f << "---\n";
      f << "Ansatz: 3D Phi=-M/R on z=0 (R^2=dx^2+dy^2+eps^2); Theta=Hess_3D(Phi)\n";
      f << "Source positions: (" << d << ",0) and (-" << d << ",0)\n";
      f << "Source masses: " << m << " each\n";
      f << "Probe geometry: +x axis and +y axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << " per axis\n";
      f << "Source softening eps=" << eps << " (numerical regularization)\n";
      f << "Provisional point-source ansatz: yes\n";
      f << "Lambda: " << LAMBDA_4D << " (fixed, 4D)\n";
      f << "\n--- Regime diagnostics (reporting only; no equation change) ---\n";
      f << "  Field intensity (Theta Frobenius): max = " << std::scientific << max_theta_norm << "\n";
      f << "  Invariant I: max|I| = " << max_invariant_I << "\n";
      f << "  Residual: available (analytic, symmetric pair)\n";
      f << "  Max residual norm = " << max_residual_norm << "\n";
      f << "  Regime (at max-intensity probe): " << tpfcore::regime_label_from_theta_norm(max_theta_norm) << "\n";
      f << "  (Thresholds: low-intensity < " << tpfcore::THETA_NORM_LOW_MAX
        << ", transitional < " << tpfcore::THETA_NORM_TRANSITIONAL_MAX << ", else high-intensity; heuristic.)\n";
      f << "---\n";
      f << "\nXi/Theta configuration-equation residual (R_nu = partial_i(Theta_i_nu - lambda*delta_i_nu*Theta)):\n";
      f << "  max|residual_x|=" << std::scientific << max_residual_x_abs << "\n";
      f << "  max|residual_y|=" << max_residual_y_abs << "\n";
      f << "  max residual_norm=" << max_residual_norm << "\n";
      f << "  residual_y near zero on +x axis (y=0 symmetry): " << (residual_y_near_zero_x_axis ? "yes [OK]\n" : "no [check]\n");
      f << "  residual_x near zero on +y axis (x=0 symmetry): " << (residual_x_near_zero_y_axis ? "yes [OK]\n" : "no [check]\n");
      f << "Output files: theta_profile.csv, invariant_profile.csv (axis column: x or y)\n";
    }
  }
}

void TPFCorePackage::run_source_field_benchmark(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  const std::string shape = config.tpf_source_benchmark_shape;
  const double orientation_deg = config.tpf_source_benchmark_orientation_deg;
  const double half_extent = config.tpf_source_probe_grid_half_extent;
  const int n = config.tpf_source_probe_grid_n;
  const double eps = (config.tpfcore_source_softening > 0.0) ? config.tpfcore_source_softening : config.softening;
  const double exclusion_radius = (config.tpf_source_residual_exclusion_radius > 0.0)
      ? config.tpf_source_residual_exclusion_radius
      : (2.0 * eps);
  const double lambda = LAMBDA_4D;

  if (n < 2 || half_extent <= 0.0) return;

  std::string source_config_id;
  std::vector<BenchmarkSourceSpec3D> sources3d;
  if (!try_build_tpf_benchmark_sources_3d(config, &sources3d, &source_config_id)) {
    return;
  }
  std::vector<BenchmarkSourceSpec2D> sources;
  sources.reserve(sources3d.size());
  for (std::size_t i = 0; i < sources3d.size(); ++i) {
    sources.push_back(BenchmarkSourceSpec2D{sources3d[i].x, sources3d[i].y, sources3d[i].m});
  }

  std::ofstream csv(output_dir + "/tpf_source_field_probe_grid.csv");
  if (!csv) return;

  const double source_mass1 = sources.empty() ? 0.0 : sources[0].m;
  const double source_mass2 = (sources.size() >= 2) ? sources[1].m : 0.0;

  const size_t n_cells = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> x_v(n_cells), y_v(n_cells);
  std::vector<double> xi_x_v(n_cells), xi_y_v(n_cells), xi_norm_v(n_cells);
  std::vector<double> theta_xx_v(n_cells), theta_xy_v(n_cells), theta_yy_v(n_cells), theta_zz_v(n_cells);
  std::vector<double> theta_trace_v(n_cells), invariant_I_v(n_cells), theta_norm_v(n_cells);
  std::vector<double> residual_x_v(n_cells, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> residual_y_v(n_cells, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> residual_norm_v(n_cells, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> residual_norm_over_theta_norm_v(n_cells, std::numeric_limits<double>::quiet_NaN());
  std::vector<int> excluded_boundary_v(n_cells, 0);
  std::vector<int> excluded_near_source_v(n_cells, 0);

  const auto idx_2d = [n](int ix, int iy) -> size_t {
    return static_cast<size_t>(iy) * static_cast<size_t>(n) + static_cast<size_t>(ix);
  };

  const double dx = (2.0 * half_extent) / static_cast<double>(n - 1);

  for (int iy = 0; iy < n; ++iy) {
    const double fy = static_cast<double>(iy) / static_cast<double>(n - 1);
    const double y = -half_extent + 2.0 * half_extent * fy;
    for (int ix = 0; ix < n; ++ix) {
      const double fx = static_cast<double>(ix) / static_cast<double>(n - 1);
      const double x = -half_extent + 2.0 * half_extent * fx;
      const size_t idx = idx_2d(ix, iy);

      FieldAtPoint combined = evaluate_provisional_field_single_source(
          sources[0].x, sources[0].y, sources[0].m, x, y, eps);
      for (size_t k = 1; k < sources.size(); ++k) {
        FieldAtPoint fk = evaluate_provisional_field_single_source(
            sources[k].x, sources[k].y, sources[k].m, x, y, eps);
        combined = add_provisional_fields(combined, fk);
      }

      x_v[idx] = x;
      y_v[idx] = y;
      xi_x_v[idx] = combined.xi.x;
      xi_y_v[idx] = combined.xi.y;
      xi_norm_v[idx] = std::hypot(combined.xi.x, combined.xi.y);
      theta_xx_v[idx] = combined.theta.xx;
      theta_xy_v[idx] = combined.theta.xy;
      theta_yy_v[idx] = combined.theta.yy;
      theta_zz_v[idx] = combined.theta.zz;
      theta_trace_v[idx] = combined.theta_trace;
      invariant_I_v[idx] = combined.invariant_I;
      theta_norm_v[idx] = theta_frobenius_norm(combined.theta);
      excluded_boundary_v[idx] = (ix == 0 || iy == 0 || ix == (n - 1) || iy == (n - 1)) ? 1 : 0;

      bool near_source = false;
      for (const BenchmarkSourceSpec2D& source : sources) {
        if (std::hypot(x - source.x, y - source.y) <= exclusion_radius) {
          near_source = true;
          break;
        }
      }
      excluded_near_source_v[idx] = near_source ? 1 : 0;
    }
  }

  const double inv_2dx = 1.0 / (2.0 * dx);
  for (int iy = 1; iy < n - 1; ++iy) {
    for (int ix = 1; ix < n - 1; ++ix) {
      const size_t idx = idx_2d(ix, iy);
      const size_t idx_xp = idx_2d(ix + 1, iy);
      const size_t idx_xm = idx_2d(ix - 1, iy);
      const size_t idx_yp = idx_2d(ix, iy + 1);
      const size_t idx_ym = idx_2d(ix, iy - 1);

      const double a_xx_p = theta_xx_v[idx_xp] - lambda * theta_trace_v[idx_xp];
      const double a_xx_m = theta_xx_v[idx_xm] - lambda * theta_trace_v[idx_xm];
      const double a_xy_p = theta_xy_v[idx_xp];
      const double a_xy_m = theta_xy_v[idx_xm];
      const double a_yx_p = theta_xy_v[idx_yp];
      const double a_yx_m = theta_xy_v[idx_ym];
      const double a_yy_p = theta_yy_v[idx_yp] - lambda * theta_trace_v[idx_yp];
      const double a_yy_m = theta_yy_v[idx_ym] - lambda * theta_trace_v[idx_ym];

      const double residual_x = (a_xx_p - a_xx_m) * inv_2dx + (a_yx_p - a_yx_m) * inv_2dx;
      const double residual_y = (a_xy_p - a_xy_m) * inv_2dx + (a_yy_p - a_yy_m) * inv_2dx;
      const double residual_norm = std::hypot(residual_x, residual_y);
      const double theta_norm = theta_norm_v[idx];
      const double residual_norm_over_theta_norm = (theta_norm > 0.0)
          ? (residual_norm / theta_norm)
          : std::numeric_limits<double>::quiet_NaN();

      residual_x_v[idx] = residual_x;
      residual_y_v[idx] = residual_y;
      residual_norm_v[idx] = residual_norm;
      residual_norm_over_theta_norm_v[idx] = residual_norm_over_theta_norm;
    }
  }

  size_t excluded_boundary_count = 0;
  size_t excluded_near_source_count = 0;
  size_t free_cell_count = 0;
  double max_residual_norm = 0.0;
  double sum_residual_norm = 0.0;
  double max_residual_ratio = 0.0;
  double sum_residual_ratio = 0.0;
  size_t residual_ratio_count = 0;
  std::vector<double> residual_norm_free;
  residual_norm_free.reserve(n_cells);

  for (size_t idx = 0; idx < n_cells; ++idx) {
    const bool excluded_boundary = excluded_boundary_v[idx] != 0;
    const bool excluded_near_source = excluded_near_source_v[idx] != 0;
    if (excluded_boundary) ++excluded_boundary_count;
    if (excluded_near_source) ++excluded_near_source_count;
    if (excluded_boundary || excluded_near_source) continue;
    if (!std::isfinite(residual_norm_v[idx])) continue;
    ++free_cell_count;
    const double rn = residual_norm_v[idx];
    const double rr = residual_norm_over_theta_norm_v[idx];
    max_residual_norm = std::max(max_residual_norm, rn);
    sum_residual_norm += rn;
    if (std::isfinite(rr)) {
      max_residual_ratio = std::max(max_residual_ratio, rr);
      sum_residual_ratio += rr;
      ++residual_ratio_count;
    }
    residual_norm_free.push_back(rn);
  }

  double median_residual_norm = std::numeric_limits<double>::quiet_NaN();
  if (!residual_norm_free.empty()) {
    std::sort(residual_norm_free.begin(), residual_norm_free.end());
    const size_t mid = residual_norm_free.size() / 2;
    if ((residual_norm_free.size() % 2) == 0) {
      median_residual_norm = 0.5 * (residual_norm_free[mid - 1] + residual_norm_free[mid]);
    } else {
      median_residual_norm = residual_norm_free[mid];
    }
  }

  const double mean_residual_norm = (free_cell_count > 0)
      ? (sum_residual_norm / static_cast<double>(free_cell_count))
      : std::numeric_limits<double>::quiet_NaN();
  const double mean_residual_ratio = (residual_ratio_count > 0)
      ? (sum_residual_ratio / static_cast<double>(residual_ratio_count))
      : std::numeric_limits<double>::quiet_NaN();

  std::ofstream summary(output_dir + "/tpf_source_field_residual_summary.txt");
  if (summary) {
    summary << std::scientific;
    summary << "total grid cells: " << n_cells << "\n";
    summary << "interior free cells used: " << free_cell_count << "\n";
    summary << "excluded boundary count: " << excluded_boundary_count << "\n";
    summary << "excluded near-source count: " << excluded_near_source_count << "\n";
    summary << "max residual norm over free cells: " << max_residual_norm << "\n";
    summary << "mean residual norm over free cells: " << mean_residual_norm << "\n";
    summary << "median residual norm over free cells: " << median_residual_norm << "\n";
    summary << "max residual_norm_over_theta_norm over free cells: " << max_residual_ratio << "\n";
    summary << "mean residual_norm_over_theta_norm over free cells: " << mean_residual_ratio << "\n";
  }

  csv << "source_shape,source_config_id,source_orientation_deg,source_mass1,source_mass2,x,y,r,xi_x,xi_y,xi_norm,"
         "theta_xx,theta_xy,theta_yy,theta_zz,theta_trace,invariant_I,residual_x,residual_y,residual_norm,"
         "residual_norm_over_theta_norm,excluded_boundary,excluded_near_source\n";
  for (size_t idx = 0; idx < n_cells; ++idx) {
    const double r = std::hypot(x_v[idx], y_v[idx]);
    csv << shape << "," << source_config_id << "," << std::scientific << orientation_deg << ","
        << source_mass1 << "," << source_mass2 << ","
        << x_v[idx] << "," << y_v[idx] << "," << r << ","
        << xi_x_v[idx] << "," << xi_y_v[idx] << "," << xi_norm_v[idx] << ","
        << theta_xx_v[idx] << "," << theta_xy_v[idx] << "," << theta_yy_v[idx] << ","
        << theta_zz_v[idx] << "," << theta_trace_v[idx] << "," << invariant_I_v[idx] << ","
        << residual_x_v[idx] << "," << residual_y_v[idx] << "," << residual_norm_v[idx] << ","
        << residual_norm_over_theta_norm_v[idx] << ","
        << excluded_boundary_v[idx] << "," << excluded_near_source_v[idx] << "\n";
  }
}

void TPFCorePackage::run_4d_static_residual_benchmark(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  const int grid_n = config.tpf_4d_residual_grid_n;
  const int bin_count = config.tpf_4d_residual_bin_count;
  const double half_extent = config.tpf_4d_residual_grid_half_extent;
  const double source_exclusion_radius = config.tpf_4d_residual_source_exclusion_radius;
  const double field_softening = config.tpf_4d_residual_field_softening;
  const double configured_bin_radius_max = config.tpf_4d_residual_bin_radius_max;

  if (grid_n < 3) {
    throw std::runtime_error("tpf_4d_residual_grid_n must be >= 3");
  }
  if (!std::isfinite(half_extent) || !(half_extent > 0.0)) {
    throw std::runtime_error("tpf_4d_residual_grid_half_extent must be finite and > 0");
  }
  if (!std::isfinite(source_exclusion_radius) || source_exclusion_radius < 0.0) {
    throw std::runtime_error("tpf_4d_residual_source_exclusion_radius must be finite and >= 0");
  }
  if (!std::isfinite(field_softening) || field_softening < 0.0) {
    throw std::runtime_error("tpf_4d_residual_field_softening must be finite and >= 0");
  }
  if (bin_count <= 0 || bin_count > 4096) {
    throw std::runtime_error("tpf_4d_residual_bin_count must be in [1, 4096]");
  }
  if (!std::isfinite(configured_bin_radius_max)) {
    throw std::runtime_error("tpf_4d_residual_bin_radius_max must be finite");
  }

  const double spacing = (2.0 * half_extent) / static_cast<double>(grid_n - 1);
  if (!std::isfinite(spacing) || !(spacing > 0.0)) {
    throw std::runtime_error("tpf_4d_residual spacing must be finite and > 0");
  }

  validate_stage4_residual_source_inputs(config);
  std::string source_config_id;
  const std::vector<BenchmarkSourceSpec3D> source_specs = build_tpf_benchmark_sources_3d(config, &source_config_id);
  validate_stage4_residual_source_specs(source_specs);
  std::vector<StaticSourcePoint4D> sources;
  sources.reserve(source_specs.size());
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    sources.push_back(StaticSourcePoint4D{source_specs[i].m, source_specs[i].x, source_specs[i].y, source_specs[i].z});
  }

  StaticResidualGridConfig grid_cfg{};
  grid_cfg.nx = static_cast<std::size_t>(grid_n);
  grid_cfg.ny = static_cast<std::size_t>(grid_n);
  grid_cfg.nz = static_cast<std::size_t>(grid_n);
  grid_cfg.spacing = spacing;
  grid_cfg.origin_x = -half_extent;
  grid_cfg.origin_y = -half_extent;
  grid_cfg.origin_z = -half_extent;
  grid_cfg.field_softening_eps = field_softening;
  grid_cfg.source_exclusion_radius = source_exclusion_radius;

  const StaticResidualGridResult result = evaluate_static_configuration_residual_4d(sources, grid_cfg);
  const double nearest_source_bin_radius_max_used = write_residual_bins_csv(
      output_dir + "/tpf_4d_static_residual_bins_nearest_source.csv",
      result,
      source_specs,
      bin_count,
      configured_bin_radius_max,
      true);
  const double origin_bin_radius_max_used = write_residual_bins_csv(
      output_dir + "/tpf_4d_static_residual_bins_origin.csv",
      result,
      source_specs,
      bin_count,
      configured_bin_radius_max,
      false);

  std::ofstream summary(output_dir + "/tpf_4d_static_residual_summary.txt");
  if (summary) {
    summary << std::scientific;
    summary << "source shape: " << config.tpf_source_benchmark_shape << "\n";
    summary << "source config id: " << source_config_id << "\n";
    summary << "source masses: ";
    for (std::size_t i = 0; i < source_specs.size(); ++i) {
      if (i > 0) summary << ", ";
      summary << source_specs[i].m;
    }
    summary << "\n";
    summary << "source positions: ";
    for (std::size_t i = 0; i < source_specs.size(); ++i) {
      if (i > 0) summary << "; ";
      summary << "(" << source_specs[i].x << "," << source_specs[i].y << "," << source_specs[i].z << ")";
    }
    summary << "\n";
    summary << "grid_n: " << grid_n << "\n";
    summary << "half_extent: " << half_extent << "\n";
    summary << "spacing: " << spacing << "\n";
    summary << "softening: " << field_softening << "\n";
    summary << "exclusion radius: " << source_exclusion_radius << "\n";
    summary << "total grid cells: " << result.summary.total_grid_cells << "\n";
    summary << "interior cells: " << result.summary.interior_cells << "\n";
    summary << "excluded boundary count: " << result.summary.excluded_boundary_count << "\n";
    summary << "excluded near-source count: " << result.summary.excluded_near_source_count << "\n";
    summary << "free cells used: " << result.summary.free_cells_used << "\n";
    summary << "z derivative samples: " << result.summary.z_derivative_samples << "\n";
    summary << "max residual spatial norm: " << result.summary.max_residual_spatial_norm << "\n";
    summary << "mean residual spatial norm: " << result.summary.mean_residual_spatial_norm << "\n";
    summary << "median residual spatial norm: " << result.summary.median_residual_spatial_norm << "\n";
    summary << "max normalized residual: " << result.summary.max_normalized_residual << "\n";
    summary << "mean normalized residual: " << result.summary.mean_normalized_residual << "\n";
    summary << "residual bin count: " << bin_count << "\n";
    summary << "nearest-source residual bin radius max used: " << nearest_source_bin_radius_max_used << "\n";
    summary << "origin residual bin radius max used: " << origin_bin_radius_max_used << "\n";
    summary << "residual bins nearest-source csv: tpf_4d_static_residual_bins_nearest_source.csv\n";
    summary << "residual bins origin csv: tpf_4d_static_residual_bins_origin.csv\n";
  }

  const std::size_t slice_i = static_cast<std::size_t>(
      std::round((0.0 - result.config.origin_x) / result.config.spacing));
  const std::size_t slice_j = static_cast<std::size_t>(
      std::round((0.0 - result.config.origin_y) / result.config.spacing));
  const std::size_t slice_k = static_cast<std::size_t>(
      std::round((0.0 - result.config.origin_z) / result.config.spacing));
  const std::size_t clamped_slice_i = std::min(result.config.nx - 1, slice_i);
  const std::size_t clamped_slice_j = std::min(result.config.ny - 1, slice_j);
  const std::size_t clamped_slice_k = std::min(result.config.nz - 1, slice_k);

  const char* header =
      "source_shape,x,y,z,residual_t,residual_x,residual_y,residual_z,residual_spatial_norm,"
      "residual_4_norm_like,theta_spatial_frobenius_norm,normalized_residual,is_boundary,is_near_source,used_in_summary,"
      "xi_t,xi_x,xi_y,xi_z,xi_spatial_norm,theta_trace_4d,invariant_I_4d\n";

  auto write_view_plane_csv =
      [&](const std::string& csv_path, const std::size_t fixed_axis_idx, const char fixed_axis) {
        std::ofstream csv(csv_path);
        if (!csv) return;
        csv << std::scientific;
        csv << header;
        for (std::size_t idx = 0; idx < result.points.size(); ++idx) {
          const StaticResidualAtPoint& p = result.points[idx];
          bool in_plane = false;
          if (fixed_axis == 'x') in_plane = (p.i == fixed_axis_idx);
          if (fixed_axis == 'y') in_plane = (p.j == fixed_axis_idx);
          if (fixed_axis == 'z') in_plane = (p.k == fixed_axis_idx);
          if (!in_plane) continue;
          const StaticProbePoint4D probe{p.x, p.y, p.z};
          const Field4DAtPoint field = evaluate_static_sources_field_4d(sources, probe, result.config.field_softening_eps);
          const double xi_spatial_norm = std::sqrt(field.xi.x * field.xi.x + field.xi.y * field.xi.y + field.xi.z * field.xi.z);
          csv << config.tpf_source_benchmark_shape << ","
              << p.x << "," << p.y << "," << p.z << ","
              << p.residual_t << "," << p.residual_x << "," << p.residual_y << "," << p.residual_z << ","
              << p.residual_spatial_norm << "," << p.residual_4_norm_like << ","
              << p.theta_spatial_frobenius_norm << "," << p.normalized_residual << ","
              << (p.is_boundary ? 1 : 0) << "," << (p.is_near_source ? 1 : 0) << "," << (p.used_in_summary ? 1 : 0) << ","
              << field.xi.t << "," << field.xi.x << "," << field.xi.y << "," << field.xi.z << ","
              << xi_spatial_norm << "," << field.theta_trace_4d << "," << field.invariant_I_4d << "\n";
        }
      };

  write_view_plane_csv(output_dir + "/tpf_4d_static_residual_slice_xy.csv", clamped_slice_k, 'z');
  write_view_plane_csv(output_dir + "/tpf_4d_static_residual_slice_xz.csv", clamped_slice_j, 'y');
  write_view_plane_csv(output_dir + "/tpf_4d_static_residual_slice_yz.csv", clamped_slice_i, 'x');
  write_view_plane_csv(output_dir + "/tpf_4d_static_residual_slice.csv", clamped_slice_k, 'z');

  std::ofstream sources_csv(output_dir + "/tpf_4d_static_residual_sources.csv");
  if (!sources_csv) return;
  sources_csv << std::scientific;
  sources_csv << "source_index,mass,x,y,z,source_config_id,source_shape\n";
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    sources_csv << i << "," << source_specs[i].m << "," << source_specs[i].x << "," << source_specs[i].y << ","
                << source_specs[i].z << "," << source_config_id << "," << config.tpf_source_benchmark_shape << "\n";
  }
}

void TPFCorePackage::run_4d_static_motion_readout_benchmark(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  const int grid_n = config.tpf_4d_motion_probe_grid_n;
  const int bin_count = config.tpf_4d_motion_bin_count;
  const double half_extent = config.tpf_4d_motion_probe_grid_half_extent;
  const double source_exclusion_radius = config.tpf_4d_motion_source_exclusion_radius;
  const double field_softening = config.tpf_4d_motion_field_softening;
  const double kappa_motion = config.tpf_4d_motion_kappa;
  const double motion_readout_scale = config.tpf_4d_motion_readout_scale;

  if (grid_n < 3) throw std::runtime_error("tpf_4d_motion_probe_grid_n must be >= 3");
  if (!std::isfinite(half_extent) || !(half_extent > 0.0)) throw std::runtime_error("tpf_4d_motion_probe_grid_half_extent must be finite and > 0");
  if (!std::isfinite(source_exclusion_radius) || source_exclusion_radius < 0.0) throw std::runtime_error("tpf_4d_motion_source_exclusion_radius must be finite and >= 0");
  if (!std::isfinite(field_softening) || field_softening < 0.0) throw std::runtime_error("tpf_4d_motion_field_softening must be finite and >= 0");
  if (!std::isfinite(kappa_motion)) throw std::runtime_error("tpf_4d_motion_kappa must be finite");
  if (!std::isfinite(motion_readout_scale)) throw std::runtime_error("tpf_4d_motion_readout_scale must be finite");
  if (bin_count < 1 || bin_count > 4096) throw std::runtime_error("tpf_4d_motion_bin_count must be in [1, 4096]");

  validate_stage4_residual_source_inputs(config);
  std::string source_config_id;
  const std::vector<BenchmarkSourceSpec3D> source_specs = build_tpf_benchmark_sources_3d(config, &source_config_id);
  validate_stage4_residual_source_specs(source_specs);
  std::vector<StaticSourcePoint4D> sources;
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    sources.push_back(StaticSourcePoint4D{source_specs[i].m, source_specs[i].x, source_specs[i].y, source_specs[i].z});
  }

  const double spacing = (2.0 * half_extent) / static_cast<double>(grid_n - 1);
  if (!std::isfinite(spacing) || !(spacing > 0.0)) {
    throw std::runtime_error("tpf_4d_motion spacing must be finite and > 0");
  }

  const double xi_eps = 1e-30;
  const std::size_t total_cells = static_cast<std::size_t>(grid_n) * static_cast<std::size_t>(grid_n) * static_cast<std::size_t>(grid_n);
  std::vector<MotionReadoutPoint> points;
  points.reserve(total_cells);
  std::size_t boundary_count = 0;
  std::size_t near_source_count = 0;
  std::size_t interior_cells = (grid_n >= 3) ? static_cast<std::size_t>(grid_n - 2) * static_cast<std::size_t>(grid_n - 2) * static_cast<std::size_t>(grid_n - 2) : 0;
  std::size_t degenerate_count = 0;

  for (int k = 0; k < grid_n; ++k) {
    for (int j = 0; j < grid_n; ++j) {
      for (int i = 0; i < grid_n; ++i) {
        MotionReadoutPoint p{};
        p.x = -half_extent + spacing * static_cast<double>(i);
        p.y = -half_extent + spacing * static_cast<double>(j);
        p.z = -half_extent + spacing * static_cast<double>(k);
        p.r_origin = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        p.is_boundary = (i == 0 || j == 0 || k == 0 || i + 1 == grid_n || j + 1 == grid_n || k + 1 == grid_n);
        if (p.is_boundary) ++boundary_count;
        p.r_nearest_source = nearest_source_distance(source_specs, p.x, p.y, p.z);
        p.is_near_source = (source_exclusion_radius > 0.0 && p.r_nearest_source <= source_exclusion_radius);
        if (!p.is_boundary && p.is_near_source) ++near_source_count;

        const StaticProbePoint4D probe{p.x, p.y, p.z};
        const Field4DAtPoint field = evaluate_static_sources_field_4d(sources, probe, field_softening);
        p.xi_t = field.xi.t;
        p.xi_x = field.xi.x;
        p.xi_y = field.xi.y;
        p.xi_z = field.xi.z;
        p.theta_trace_4d = field.theta_trace_4d;
        p.invariant_I_4d = field.invariant_I_4d;
        p.xi_spatial_norm = std::sqrt(p.xi_x * p.xi_x + p.xi_y * p.xi_y + p.xi_z * p.xi_z);

        const double theta_xx = field.theta.xx;
        const double theta_xy = field.theta.xy;
        const double theta_xz = field.theta.xz;
        const double theta_yx = field.theta.yx;
        const double theta_yy = field.theta.yy;
        const double theta_yz = field.theta.yz;
        const double theta_zx = field.theta.zx;
        const double theta_zy = field.theta.zy;
        const double theta_zz = field.theta.zz;
        const double theta_tr = field.theta_trace_4d;
        const double inv_I = field.invariant_I_4d;

        const double t_xx = theta_xx * theta_xx + theta_xy * theta_xy + theta_xz * theta_xz;
        const double t_xy = theta_xx * theta_yx + theta_xy * theta_yy + theta_xz * theta_yz;
        const double t_xz = theta_xx * theta_zx + theta_xy * theta_zy + theta_xz * theta_zz;
        const double t_yx = theta_yx * theta_xx + theta_yy * theta_xy + theta_yz * theta_xz;
        const double t_yy = theta_yx * theta_yx + theta_yy * theta_yy + theta_yz * theta_yz;
        const double t_yz = theta_yx * theta_zx + theta_yy * theta_zy + theta_yz * theta_zz;
        const double t_zx = theta_zx * theta_xx + theta_zy * theta_xy + theta_zz * theta_xz;
        const double t_zy = theta_zx * theta_yx + theta_zy * theta_yy + theta_zz * theta_yz;
        const double t_zz = theta_zx * theta_zx + theta_zy * theta_zy + theta_zz * theta_zz;

        p.c_xx = kappa_motion * (t_xx - LAMBDA_4D * theta_tr * theta_xx - 0.5 * inv_I);
        p.c_xy = kappa_motion * (t_xy - LAMBDA_4D * theta_tr * theta_xy);
        p.c_xz = kappa_motion * (t_xz - LAMBDA_4D * theta_tr * theta_xz);
        p.c_yx = kappa_motion * (t_yx - LAMBDA_4D * theta_tr * theta_yx);
        p.c_yy = kappa_motion * (t_yy - LAMBDA_4D * theta_tr * theta_yy - 0.5 * inv_I);
        p.c_yz = kappa_motion * (t_yz - LAMBDA_4D * theta_tr * theta_yz);
        p.c_zx = kappa_motion * (t_zx - LAMBDA_4D * theta_tr * theta_zx);
        p.c_zy = kappa_motion * (t_zy - LAMBDA_4D * theta_tr * theta_zy);
        p.c_zz = kappa_motion * (t_zz - LAMBDA_4D * theta_tr * theta_zz - 0.5 * inv_I);

        if (p.xi_spatial_norm <= xi_eps) {
          p.xi_degenerate = true;
          ++degenerate_count;
          p.a_x = p.a_y = p.a_z = 0.0;
        } else {
          const double hx = p.xi_x / p.xi_spatial_norm;
          const double hy = p.xi_y / p.xi_spatial_norm;
          const double hz = p.xi_z / p.xi_spatial_norm;
          p.a_x = -motion_readout_scale * (p.c_xx * hx + p.c_xy * hy + p.c_xz * hz);
          p.a_y = -motion_readout_scale * (p.c_yx * hx + p.c_yy * hy + p.c_yz * hz);
          p.a_z = -motion_readout_scale * (p.c_zx * hx + p.c_zy * hy + p.c_zz * hz);
        }
        p.a_norm = std::sqrt(p.a_x * p.a_x + p.a_y * p.a_y + p.a_z * p.a_z);

        if (p.r_origin > 1e-30 && p.a_norm > 1e-30) {
          const double inward_x = -p.x / p.r_origin;
          const double inward_y = -p.y / p.r_origin;
          const double inward_z = -p.z / p.r_origin;
          const double radial = (p.a_x * inward_x + p.a_y * inward_y + p.a_z * inward_z) / p.a_norm;
          p.radial_alignment_to_origin_inward = radial;
          p.transverse_fraction_origin = std::sqrt(std::max(0.0, 1.0 - std::min(1.0, radial * radial)));
        }
        p.used_in_summary = (!p.is_boundary && !p.is_near_source && !p.xi_degenerate);
        points.push_back(p);
      }
    }
  }

  std::vector<double> used_norms, used_align, used_transverse;
  double sum_norm = 0.0, sum_align = 0.0, sum_transverse = 0.0;
  double max_norm = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (!points[i].used_in_summary) continue;
    used_norms.push_back(points[i].a_norm);
    sum_norm += points[i].a_norm;
    max_norm = std::max(max_norm, points[i].a_norm);
    if (std::isfinite(points[i].radial_alignment_to_origin_inward)) {
      used_align.push_back(points[i].radial_alignment_to_origin_inward);
      sum_align += points[i].radial_alignment_to_origin_inward;
    }
    if (std::isfinite(points[i].transverse_fraction_origin)) {
      used_transverse.push_back(points[i].transverse_fraction_origin);
      sum_transverse += points[i].transverse_fraction_origin;
    }
  }
  std::sort(used_norms.begin(), used_norms.end());
  std::sort(used_align.begin(), used_align.end());
  std::sort(used_transverse.begin(), used_transverse.end());
  const double mean_norm = used_norms.empty() ? 0.0 : (sum_norm / static_cast<double>(used_norms.size()));
  const double median_norm = percentile_sorted_linear(used_norms, 0.5);
  const double mean_align = used_align.empty() ? 0.0 : (sum_align / static_cast<double>(used_align.size()));
  const double median_align = percentile_sorted_linear(used_align, 0.5);
  const double mean_transverse = used_transverse.empty() ? 0.0 : (sum_transverse / static_cast<double>(used_transverse.size()));

  std::ofstream probe_csv(output_dir + "/tpf_4d_static_motion_readout_probe_grid.csv");
  probe_csv << std::scientific;
  probe_csv << "source_shape,x,y,z,r_origin,r_nearest_source,xi_t,xi_x,xi_y,xi_z,xi_spatial_norm,theta_trace_4d,invariant_I_4d,"
               "c_xx,c_xy,c_xz,c_yx,c_yy,c_yz,c_zx,c_zy,c_zz,a_x,a_y,a_z,a_norm,radial_alignment_to_origin_inward,"
               "transverse_fraction_origin,is_boundary,is_near_source,xi_degenerate,used_in_summary\n";
  for (std::size_t i = 0; i < points.size(); ++i) {
    const MotionReadoutPoint& p = points[i];
    probe_csv << config.tpf_source_benchmark_shape << "," << p.x << "," << p.y << "," << p.z << "," << p.r_origin << ","
              << p.r_nearest_source << "," << p.xi_t << "," << p.xi_x << "," << p.xi_y << "," << p.xi_z << ","
              << p.xi_spatial_norm << "," << p.theta_trace_4d << "," << p.invariant_I_4d << ","
              << p.c_xx << "," << p.c_xy << "," << p.c_xz << "," << p.c_yx << "," << p.c_yy << "," << p.c_yz << ","
              << p.c_zx << "," << p.c_zy << "," << p.c_zz << "," << p.a_x << "," << p.a_y << "," << p.a_z << "," << p.a_norm
              << "," << p.radial_alignment_to_origin_inward << "," << p.transverse_fraction_origin << ","
              << (p.is_boundary ? 1 : 0) << "," << (p.is_near_source ? 1 : 0) << "," << (p.xi_degenerate ? 1 : 0) << ","
              << (p.used_in_summary ? 1 : 0) << "\n";
  }

  struct BinAccum {
    int total = 0;
    int used = 0;
    std::vector<double> a_norm, align, transverse;
  };
  const double radius_max = std::sqrt(3.0) * half_extent;
  const double bin_width = radius_max / static_cast<double>(bin_count);
  std::vector<BinAccum> bins(static_cast<std::size_t>(bin_count));
  for (std::size_t i = 0; i < points.size(); ++i) {
    int idx = (bin_width > 0.0) ? static_cast<int>(points[i].r_origin / bin_width) : 0;
    if (idx < 0) idx = 0;
    if (idx >= bin_count) idx = bin_count - 1;
    BinAccum& b = bins[static_cast<std::size_t>(idx)];
    ++b.total;
    if (!points[i].used_in_summary) continue;
    ++b.used;
    b.a_norm.push_back(points[i].a_norm);
    if (std::isfinite(points[i].radial_alignment_to_origin_inward)) b.align.push_back(points[i].radial_alignment_to_origin_inward);
    if (std::isfinite(points[i].transverse_fraction_origin)) b.transverse.push_back(points[i].transverse_fraction_origin);
  }

  std::ofstream bins_csv(output_dir + "/tpf_4d_static_motion_readout_bins_origin.csv");
  bins_csv << std::scientific;
  bins_csv << "bin_index,r_min,r_max,r_mid,cell_count_total,cell_count_used,mean_a_norm,median_a_norm,p95_a_norm,max_a_norm,"
              "mean_radial_alignment,median_radial_alignment,mean_transverse_fraction,median_transverse_fraction\n";
  for (int i = 0; i < bin_count; ++i) {
    BinAccum& b = bins[static_cast<std::size_t>(i)];
    std::sort(b.a_norm.begin(), b.a_norm.end());
    std::sort(b.align.begin(), b.align.end());
    std::sort(b.transverse.begin(), b.transverse.end());
    const double r_min = bin_width * static_cast<double>(i);
    const double r_max = (i + 1 == bin_count) ? radius_max : (bin_width * static_cast<double>(i + 1));
    const double r_mid = 0.5 * (r_min + r_max);
    bins_csv << i << "," << r_min << "," << r_max << "," << r_mid << "," << b.total << "," << b.used << ","
             << mean_or_nan(b.a_norm) << "," << percentile_sorted_linear(b.a_norm, 0.5) << ","
             << percentile_sorted_linear(b.a_norm, 0.95) << "," << max_or_nan(b.a_norm) << ","
             << mean_or_nan(b.align) << "," << percentile_sorted_linear(b.align, 0.5) << ","
             << mean_or_nan(b.transverse) << "," << percentile_sorted_linear(b.transverse, 0.5) << "\n";
  }

  double falloff_slope = std::numeric_limits<double>::quiet_NaN();
  if (config.tpf_source_benchmark_shape == "monopole") {
    std::vector<double> xs;
    std::vector<double> ys;
    for (int i = 0; i < bin_count; ++i) {
      const BinAccum& b = bins[static_cast<std::size_t>(i)];
      if (b.used <= 0) continue;
      const double r_mid = (static_cast<double>(i) + 0.5) * bin_width;
      const double mean_a = mean_or_nan(b.a_norm);
      if (r_mid > 0.0 && std::isfinite(mean_a) && mean_a > 0.0) {
        xs.push_back(std::log(r_mid));
        ys.push_back(std::log(mean_a));
      }
    }
    if (xs.size() >= 2) {
      double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
      for (std::size_t i = 0; i < xs.size(); ++i) {
        sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i];
      }
      const double n = static_cast<double>(xs.size());
      const double denom = n * sxx - sx * sx;
      if (std::fabs(denom) > 1e-30) {
        falloff_slope = (n * sxy - sx * sy) / denom;
      }
    }
  }

  std::ofstream summary(output_dir + "/tpf_4d_static_motion_readout_summary.txt");
  summary << std::scientific;
  summary << "source shape: " << config.tpf_source_benchmark_shape << "\n";
  summary << "source config id: " << source_config_id << "\n";
  summary << "source masses: ";
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    if (i > 0) summary << ", ";
    summary << source_specs[i].m;
  }
  summary << "\nsource positions: ";
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    if (i > 0) summary << "; ";
    summary << "(" << source_specs[i].x << "," << source_specs[i].y << "," << source_specs[i].z << ")";
  }
  summary << "\n";
  summary << "grid_n: " << grid_n << "\n";
  summary << "half_extent: " << half_extent << "\n";
  summary << "spacing: " << spacing << "\n";
  summary << "field softening: " << field_softening << "\n";
  summary << "exclusion radius: " << source_exclusion_radius << "\n";
  summary << "kappa_motion: " << kappa_motion << "\n";
  summary << "motion_readout_scale: " << motion_readout_scale << "\n";
  summary << "total grid cells: " << total_cells << "\n";
  summary << "interior cells: " << interior_cells << "\n";
  summary << "used/free probe cells: " << used_norms.size() << "\n";
  summary << "excluded boundary count: " << boundary_count << "\n";
  summary << "excluded near-source count: " << near_source_count << "\n";
  summary << "degenerate Xi count: " << degenerate_count << "\n";
  summary << "mean acceleration norm: " << mean_norm << "\n";
  summary << "median acceleration norm: " << median_norm << "\n";
  summary << "max acceleration norm: " << max_norm << "\n";
  summary << "mean radial alignment: " << mean_align << "\n";
  summary << "median radial alignment: " << median_align << "\n";
  summary << "radial alignment sample count: " << used_align.size() << "\n";
  summary << "mean transverse fraction: " << mean_transverse << "\n";
  summary << "transverse fraction sample count: " << used_transverse.size() << "\n";
  summary << "measured log-log falloff slope for monopole if available: " << falloff_slope << "\n";
  summary << "falloff slope note: slope is measured from this readout and is not forced to Newtonian.\n";
  summary << "bins origin csv: tpf_4d_static_motion_readout_bins_origin.csv\n";
  summary << "probe grid csv: tpf_4d_static_motion_readout_probe_grid.csv\n";
}

void TPFCorePackage::run_4d_xi_motion_probe_benchmark(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  const double dt = config.tpf_4d_xi_motion_dt;
  const int steps = config.tpf_4d_xi_motion_steps;
  const double readout_scale = config.tpf_4d_xi_motion_readout_scale;
  const double field_softening = config.tpf_4d_xi_motion_field_softening;
  const double source_exclusion_radius = config.tpf_4d_xi_motion_source_exclusion_radius;
  const std::string layout = config.tpf_4d_xi_motion_probe_layout;
  const int probe_count_cfg = config.tpf_4d_xi_motion_probe_count;
  const double probe_radius = config.tpf_4d_xi_motion_probe_radius;
  const double probe_speed = config.tpf_4d_xi_motion_probe_speed;
  const std::string integrator = config.tpf_4d_xi_motion_integrator;
  const int dump_every = config.tpf_4d_xi_motion_dump_every;
  const std::string kernel_mode = config.tpf_4d_xi_kernel_mode;
  const double kernel_coupling = config.tpf_4d_xi_kernel_coupling;
  const double kernel_beta_power = config.tpf_4d_xi_kernel_beta_power;
  const std::string kernel_factor_mode = config.tpf_4d_xi_kernel_factor_mode;
  const double kernel_metric_min = config.tpf_4d_xi_kernel_metric_min;
  const double kernel_metric_max = config.tpf_4d_xi_kernel_metric_max;
  const std::string temporal_mode = config.tpf_4d_xi_temporal_mode;
  const double temporal_coupling = config.tpf_4d_xi_temporal_coupling;
  const double source_speed_x = config.tpf_4d_xi_source_speed_x;
  const double source_speed_y = config.tpf_4d_xi_source_speed_y;
  const double source_speed_z = config.tpf_4d_xi_source_speed_z;

  if (!std::isfinite(dt) || dt <= 0.0) throw std::runtime_error("tpf_4d_xi_motion_dt must be finite and > 0");
  if (steps < 1) throw std::runtime_error("tpf_4d_xi_motion_steps must be >= 1");
  if (!std::isfinite(readout_scale)) throw std::runtime_error("tpf_4d_xi_motion_readout_scale must be finite");
  if (!std::isfinite(field_softening) || field_softening < 0.0) throw std::runtime_error("tpf_4d_xi_motion_field_softening must be finite and >= 0");
  if (!std::isfinite(source_exclusion_radius) || source_exclusion_radius < 0.0) throw std::runtime_error("tpf_4d_xi_motion_source_exclusion_radius must be finite and >= 0");
  if (probe_count_cfg < 1 || probe_count_cfg > 4096) throw std::runtime_error("tpf_4d_xi_motion_probe_count must be in [1, 4096]");
  if (!std::isfinite(probe_radius) || probe_radius < 0.0) throw std::runtime_error("tpf_4d_xi_motion_probe_radius must be finite and >= 0");
  if (!std::isfinite(probe_speed)) throw std::runtime_error("tpf_4d_xi_motion_probe_speed must be finite");
  if (layout != "ring" && layout != "axis") throw std::runtime_error("tpf_4d_xi_motion_probe_layout must be ring or axis");
  if (integrator != "velocity_verlet" && integrator != "semi_implicit_euler") throw std::runtime_error("tpf_4d_xi_motion_integrator must be velocity_verlet or semi_implicit_euler");
  if (dump_every < 1) throw std::runtime_error("tpf_4d_xi_motion_dump_every must be >= 1");
  if (kernel_mode != "off" && kernel_mode != "scalar_beta" && kernel_mode != "metric_radial" && kernel_mode != "metric_velocity" &&
      kernel_mode != "metric_transverse_wake" && kernel_mode != "metric_transverse_continuous" && kernel_mode != "spacetime_metric") {
    throw std::runtime_error(
        "tpf_4d_xi_kernel_mode must be one of: off, scalar_beta, metric_radial, metric_velocity, metric_transverse_wake, metric_transverse_continuous, spacetime_metric");
  }
  if (kernel_factor_mode != "beta_power" && kernel_factor_mode != "gamma_minus_one") {
    throw std::runtime_error("tpf_4d_xi_kernel_factor_mode must be beta_power or gamma_minus_one");
  }
  if (temporal_mode != "off" && temporal_mode != "norm_scaled") {
    throw std::runtime_error("tpf_4d_xi_temporal_mode must be off or norm_scaled");
  }
  if (!std::isfinite(kernel_coupling)) throw std::runtime_error("tpf_4d_xi_kernel_coupling must be finite");
  const bool xi_kernel_deformation_active = (kernel_mode != "off" && kernel_coupling != 0.0);
  if (xi_kernel_deformation_active && integrator == "velocity_verlet") {
    throw std::runtime_error("velocity_verlet is disabled for active velocity-dependent Xi kernel modes; use semi_implicit_euler");
  }
  if (!std::isfinite(kernel_beta_power) || kernel_beta_power < 0.0) {
    throw std::runtime_error("tpf_4d_xi_kernel_beta_power must be finite and >= 0");
  }
  if (!std::isfinite(kernel_metric_min) || kernel_metric_min <= 0.0) {
    throw std::runtime_error("tpf_4d_xi_kernel_metric_min must be finite and > 0");
  }
  if (!std::isfinite(kernel_metric_max) || kernel_metric_max < kernel_metric_min) {
    throw std::runtime_error("tpf_4d_xi_kernel_metric_max must be finite and >= tpf_4d_xi_kernel_metric_min");
  }
  if (!std::isfinite(temporal_coupling)) throw std::runtime_error("tpf_4d_xi_temporal_coupling must be finite");
  if (!std::isfinite(source_speed_x) || !std::isfinite(source_speed_y) || !std::isfinite(source_speed_z)) {
    throw std::runtime_error("tpf_4d_xi_source_speed_x/y/z must all be finite");
  }
  validate_stage4_residual_source_inputs(config);

  std::string source_config_id;
  const std::vector<BenchmarkSourceSpec3D> source_specs = build_tpf_benchmark_sources_3d(config, &source_config_id);
  validate_stage4_residual_source_specs(source_specs);
  std::vector<StaticSourcePoint4D> static_sources;
  static_sources.reserve(source_specs.size());
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    static_sources.push_back(StaticSourcePoint4D{source_specs[i].m, source_specs[i].x, source_specs[i].y, source_specs[i].z});
  }

  struct ProbeState {
    int id = 0;
    double x = 0.0, y = 0.0, z = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double ax = 0.0, ay = 0.0, az = 0.0;
    double xi_x = 0.0, xi_y = 0.0, xi_z = 0.0;
    double xi_x_base = 0.0, xi_y_base = 0.0, xi_z_base = 0.0;
    double xi_x_eff = 0.0, xi_y_eff = 0.0, xi_z_eff = 0.0;
    double xi_spatial_norm = 0.0;
    double xi_t = 0.0;
    double beta_rel = 0.0;
    double gamma_rel = 1.0;
    double v_radial = 0.0;
    double v_transverse = 0.0;
    double beta_pass = 0.0;
    double wake_gate = 0.0;
    double beta_effective = 0.0;
    double v_radial_min = 0.0;
    double v_radial_max = 0.0;
    double wake_gate_min = 0.0;
    double wake_gate_max = 0.0;
    double xi_kernel_factor_raw = 0.0;
    double xi_kernel_metric_scale = 1.0;
    double theta_trace_4d = 0.0;
    double invariant_I_4d = 0.0;
    bool is_near_source = false;
    bool escaped = false;
    bool valid = true;
  };

  struct KernelStats {
    double sum_beta_rel = 0.0;
    double max_beta_rel = 0.0;
    double sum_gamma_rel = 0.0;
    double max_gamma_rel = 1.0;
    double sum_factor_raw = 0.0;
    double max_abs_factor_raw = 0.0;
    double sum_metric_scale = 0.0;
    double max_metric_scale = 1.0;
    double sum_xi_t = 0.0;
    double max_abs_xi_t = 0.0;
    std::size_t sample_count = 0;
  } kernel_stats;

  auto eval_probe = [&](ProbeState& p, bool track_stats) {
    p.is_near_source = (source_exclusion_radius > 0.0 &&
                        nearest_source_distance(source_specs, p.x, p.y, p.z) <= source_exclusion_radius);
    p.valid = std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) &&
              std::isfinite(p.vx) && std::isfinite(p.vy) && std::isfinite(p.vz) && !p.is_near_source;
    if (!p.valid) {
      p.escaped = true;
      p.ax = p.ay = p.az = 0.0;
      p.xi_x = p.xi_y = p.xi_z = 0.0;
      p.xi_x_base = p.xi_y_base = p.xi_z_base = 0.0;
      p.xi_x_eff = p.xi_y_eff = p.xi_z_eff = 0.0;
      p.xi_spatial_norm = 0.0;
      p.xi_t = 0.0;
      p.beta_rel = 0.0;
      p.gamma_rel = 1.0;
      p.v_radial = 0.0;
      p.v_transverse = 0.0;
      p.beta_pass = 0.0;
      p.wake_gate = 0.0;
      p.beta_effective = 0.0;
      p.v_radial_min = 0.0;
      p.v_radial_max = 0.0;
      p.wake_gate_min = 0.0;
      p.wake_gate_max = 0.0;
      p.xi_kernel_factor_raw = 0.0;
      p.xi_kernel_metric_scale = 1.0;
      p.theta_trace_4d = std::numeric_limits<double>::quiet_NaN();
      p.invariant_I_4d = std::numeric_limits<double>::quiet_NaN();
      return;
    }

    const double vx_rel = p.vx - source_speed_x;
    const double vy_rel = p.vy - source_speed_y;
    const double vz_rel = p.vz - source_speed_z;
    const double v_rel_norm = std::sqrt(vx_rel * vx_rel + vy_rel * vy_rel + vz_rel * vz_rel);
    const double beta_rel = v_rel_norm / C_SI_LIGHT;
    if (!std::isfinite(beta_rel)) throw std::runtime_error("non-finite beta_rel in Xi kernel deformation");
    if (beta_rel >= 1.0) throw std::runtime_error("beta_rel must be < 1.0 in Xi kernel deformation");
    const double beta_for_gamma = std::min(beta_rel, 1.0 - 1.0e-12);
    const double gamma_rel = 1.0 / std::sqrt(1.0 - beta_for_gamma * beta_for_gamma);
    if (!std::isfinite(gamma_rel)) throw std::runtime_error("non-finite gamma_rel in Xi kernel deformation");

    double xi_base_x = 0.0, xi_base_y = 0.0, xi_base_z = 0.0;
    double xi_eff_x = 0.0, xi_eff_y = 0.0, xi_eff_z = 0.0;
    double wake_weight = 0.0;
    double v_radial_wsum = 0.0, v_transverse_wsum = 0.0, beta_pass_wsum = 0.0, wake_gate_wsum = 0.0, beta_effective_wsum = 0.0;
    double factor_raw_wsum = 0.0, metric_scale_wsum = 0.0;
    double v_radial_min = std::numeric_limits<double>::infinity();
    double v_radial_max = -std::numeric_limits<double>::infinity();
    double wake_gate_min = std::numeric_limits<double>::infinity();
    double wake_gate_max = -std::numeric_limits<double>::infinity();
    std::size_t wake_sample_count = 0;

    for (std::size_t si = 0; si < source_specs.size(); ++si) {
      const double dx = p.x - source_specs[si].x;
      const double dy = p.y - source_specs[si].y;
      const double dz = p.z - source_specs[si].z;
      const tpfcore::XiWakeKinematics wake =
          tpfcore::compute_xi_wake_kinematics(dx, dy, dz, vx_rel, vy_rel, vz_rel, C_SI_LIGHT,
                                              kernel_mode == "metric_transverse_wake");
      const double beta_for_gamma_source = std::min(wake.beta_effective, 1.0 - 1.0e-12);
      const double gamma_source = 1.0 / std::sqrt(1.0 - beta_for_gamma_source * beta_for_gamma_source);
      double factor_raw = 0.0;
      if (kernel_factor_mode == "beta_power") {
        factor_raw = kernel_coupling * std::pow(wake.beta_effective, kernel_beta_power);
      } else {
        factor_raw = kernel_coupling * (gamma_source - 1.0);
      }
      if (!std::isfinite(factor_raw)) throw std::runtime_error("non-finite factor_raw in Xi kernel deformation");
      double metric_scale = 1.0;
      if (kernel_mode == "metric_transverse_wake" || kernel_mode == "metric_transverse_continuous") {
        constexpr double k_factor_floor_eps = 1.0e-12;
        const double denom_floor = -1.0 + k_factor_floor_eps;
        const double guarded_factor = std::max(factor_raw, denom_floor);
        metric_scale = std::max(kernel_metric_min,
                                std::min(kernel_metric_max, 1.0 / (1.0 + guarded_factor)));
      } else {
        metric_scale = std::max(kernel_metric_min, std::min(kernel_metric_max, 1.0 + factor_raw));
      }
      const double w = std::fabs(source_specs[si].m);
      wake_weight += w;
      v_radial_wsum += w * wake.v_radial;
      v_transverse_wsum += w * wake.v_transverse;
      beta_pass_wsum += w * wake.beta_pass;
      wake_gate_wsum += w * wake.wake_gate;
      beta_effective_wsum += w * wake.beta_effective;
      factor_raw_wsum += w * factor_raw;
      metric_scale_wsum += w * metric_scale;
      v_radial_min = std::min(v_radial_min, wake.v_radial);
      v_radial_max = std::max(v_radial_max, wake.v_radial);
      wake_gate_min = std::min(wake_gate_min, wake.wake_gate);
      wake_gate_max = std::max(wake_gate_max, wake.wake_gate);
      ++wake_sample_count;

      constexpr double k_probe_pair_tiny = 1.0e-30;
      const double r2 = dx * dx + dy * dy + dz * dz + field_softening * field_softening;
      if (!std::isfinite(r2) || r2 <= k_probe_pair_tiny) continue;
      const double r = std::sqrt(r2);
      const double inv_r3 = 1.0 / (r2 * r);
      const double xi_sx = source_specs[si].m * dx * inv_r3;
      const double xi_sy = source_specs[si].m * dy * inv_r3;
      const double xi_sz = source_specs[si].m * dz * inv_r3;
      xi_base_x += xi_sx;
      xi_base_y += xi_sy;
      xi_base_z += xi_sz;

      if (kernel_mode == "off") {
        xi_eff_x += xi_sx;
        xi_eff_y += xi_sy;
        xi_eff_z += xi_sz;
      } else if (kernel_mode == "scalar_beta") {
        const double scalar = 1.0 + factor_raw;
        xi_eff_x += xi_sx * scalar;
        xi_eff_y += xi_sy * scalar;
        xi_eff_z += xi_sz * scalar;
      } else {
        double nx = 0.0, ny = 0.0, nz = 0.0;
        double n_norm = 0.0;
        if (kernel_mode == "metric_radial") {
          n_norm = std::sqrt(dx * dx + dy * dy + dz * dz);
          if (n_norm <= 1.0e-30) throw std::runtime_error("near-source invalid radial direction in metric_radial Xi kernel mode");
          nx = dx / n_norm;
          ny = dy / n_norm;
          nz = dz / n_norm;
        } else if (kernel_mode == "metric_transverse_wake" || kernel_mode == "metric_transverse_continuous") {
          n_norm = std::sqrt(dx * dx + dy * dy + dz * dz);
          if (n_norm <= k_probe_pair_tiny) continue;
          nx = dx / n_norm;
          ny = dy / n_norm;
          nz = dz / n_norm;
        } else {
          n_norm = v_rel_norm;
          if (n_norm > 1.0e-30) {
            nx = vx_rel / n_norm;
            ny = vy_rel / n_norm;
            nz = vz_rel / n_norm;
          }
        }

        const double alpha = (n_norm > 1.0e-30) ? (metric_scale - 1.0) : 0.0;
        const double nd = nx * dx + ny * dy + nz * dz;
        const double gx = dx + alpha * nx * nd;
        const double gy = dy + alpha * ny * nd;
        const double gz = dz + alpha * nz * nd;
        const double r_eff2 = dx * gx + dy * gy + dz * gz + field_softening * field_softening;
        if (std::isfinite(r_eff2) && r_eff2 > k_probe_pair_tiny) {
          const double r_eff = std::sqrt(r_eff2);
          const double inv_r_eff3 = 1.0 / (r_eff2 * r_eff);
          xi_eff_x += source_specs[si].m * gx * inv_r_eff3;
          xi_eff_y += source_specs[si].m * gy * inv_r_eff3;
          xi_eff_z += source_specs[si].m * gz * inv_r_eff3;
        } else {
          xi_eff_x += xi_sx;
          xi_eff_y += xi_sy;
          xi_eff_z += xi_sz;
        }
      }
    }

    const StaticProbePoint4D probe_field{p.x, p.y, p.z};
    const Field4DAtPoint field_diag = evaluate_static_sources_field_4d(static_sources, probe_field, field_softening);

    const double xi_eff_norm = std::sqrt(xi_eff_x * xi_eff_x + xi_eff_y * xi_eff_y + xi_eff_z * xi_eff_z);
    const double inv_wake_weight = (wake_weight > 0.0) ? (1.0 / wake_weight) : 0.0;
    const double mean_factor_raw = factor_raw_wsum * inv_wake_weight;
    double xi_t = 0.0;
    if (temporal_mode == "norm_scaled" && kernel_mode == "spacetime_metric") {
      xi_t = temporal_coupling * mean_factor_raw * xi_eff_norm;
    }

    if (!std::isfinite(xi_eff_x) || !std::isfinite(xi_eff_y) || !std::isfinite(xi_eff_z) || !std::isfinite(xi_t)) {
      throw std::runtime_error("non-finite Xi kernel deformation outputs");
    }

    p.xi_x_base = xi_base_x;
    p.xi_y_base = xi_base_y;
    p.xi_z_base = xi_base_z;
    p.xi_x_eff = xi_eff_x;
    p.xi_y_eff = xi_eff_y;
    p.xi_z_eff = xi_eff_z;
    p.xi_x = xi_eff_x;
    p.xi_y = xi_eff_y;
    p.xi_z = xi_eff_z;
    p.xi_spatial_norm = xi_eff_norm;
    p.xi_t = xi_t;
    p.beta_rel = beta_rel;
    p.gamma_rel = gamma_rel;
    p.v_radial = v_radial_wsum * inv_wake_weight;
    p.v_transverse = v_transverse_wsum * inv_wake_weight;
    p.beta_pass = beta_pass_wsum * inv_wake_weight;
    p.wake_gate = wake_gate_wsum * inv_wake_weight;
    p.beta_effective = beta_effective_wsum * inv_wake_weight;
    p.v_radial_min = (wake_sample_count > 0) ? v_radial_min : 0.0;
    p.v_radial_max = (wake_sample_count > 0) ? v_radial_max : 0.0;
    p.wake_gate_min = (wake_sample_count > 0) ? wake_gate_min : 0.0;
    p.wake_gate_max = (wake_sample_count > 0) ? wake_gate_max : 0.0;
    p.xi_kernel_factor_raw = (kernel_mode == "off") ? 0.0 : mean_factor_raw;
    p.xi_kernel_metric_scale = (kernel_mode == "metric_radial" || kernel_mode == "metric_velocity" ||
                                kernel_mode == "metric_transverse_wake" || kernel_mode == "metric_transverse_continuous" || kernel_mode == "spacetime_metric")
                                   ? (metric_scale_wsum * inv_wake_weight)
                                   : 1.0;
    p.theta_trace_4d = field_diag.theta_trace_4d;
    p.invariant_I_4d = field_diag.invariant_I_4d;
    p.ax = -readout_scale * p.xi_x_eff;
    p.ay = -readout_scale * p.xi_y_eff;
    p.az = -readout_scale * p.xi_z_eff;
    p.valid = p.valid && std::isfinite(p.ax) && std::isfinite(p.ay) && std::isfinite(p.az);
    if (!p.valid) p.escaped = true;

    if (track_stats) {
      kernel_stats.sum_beta_rel += p.beta_rel;
      kernel_stats.max_beta_rel = std::max(kernel_stats.max_beta_rel, p.beta_rel);
      kernel_stats.sum_gamma_rel += p.gamma_rel;
      kernel_stats.max_gamma_rel = std::max(kernel_stats.max_gamma_rel, p.gamma_rel);
      kernel_stats.sum_factor_raw += p.xi_kernel_factor_raw;
      kernel_stats.max_abs_factor_raw = std::max(kernel_stats.max_abs_factor_raw, std::fabs(p.xi_kernel_factor_raw));
      kernel_stats.sum_metric_scale += p.xi_kernel_metric_scale;
      kernel_stats.max_metric_scale = std::max(kernel_stats.max_metric_scale, p.xi_kernel_metric_scale);
      kernel_stats.sum_xi_t += p.xi_t;
      kernel_stats.max_abs_xi_t = std::max(kernel_stats.max_abs_xi_t, std::fabs(p.xi_t));
      ++kernel_stats.sample_count;
    }
  };

  std::vector<ProbeState> probes;
  if (layout == "ring") {
    probes.resize(static_cast<std::size_t>(probe_count_cfg));
    for (int i = 0; i < probe_count_cfg; ++i) {
      const double t = (2.0 * 3.14159265358979323846 * static_cast<double>(i)) / static_cast<double>(probe_count_cfg);
      ProbeState p;
      p.id = i;
      p.x = probe_radius * std::cos(t);
      p.y = probe_radius * std::sin(t);
      p.z = 0.0;
      p.vx = -probe_speed * std::sin(t);
      p.vy = probe_speed * std::cos(t);
      p.vz = 0.0;
      probes[static_cast<std::size_t>(i)] = p;
    }
  } else {
    const int n = std::max(6, probe_count_cfg);
    probes.resize(static_cast<std::size_t>(n));
    const double coords[6][3] = {
        {probe_radius, 0.0, 0.0}, {-probe_radius, 0.0, 0.0},
        {0.0, probe_radius, 0.0}, {0.0, -probe_radius, 0.0},
        {0.0, 0.0, probe_radius}, {0.0, 0.0, -probe_radius}};
    for (int i = 0; i < n; ++i) {
      ProbeState p;
      p.id = i;
      p.x = coords[i % 6][0];
      p.y = coords[i % 6][1];
      p.z = coords[i % 6][2];
      probes[static_cast<std::size_t>(i)] = p;
    }
  }

  std::ofstream traj(output_dir + "/tpf_4d_xi_motion_probe_trajectories.csv");
  std::ofstream init_csv(output_dir + "/tpf_4d_xi_motion_probe_initial_readout.csv");
  if (!traj || !init_csv) throw std::runtime_error("failed to open Xi motion benchmark CSV outputs");
  traj << std::scientific;
  init_csv << std::scientific;
  traj << "step,time,probe_id,x,y,z,vx,vy,vz,ax,ay,az,a_norm,xi_x,xi_y,xi_z,xi_spatial_norm,xi_t,xi_x_base,xi_y_base,xi_z_base,"
          "xi_x_eff,xi_y_eff,xi_z_eff,beta_rel,gamma_rel,v_radial,v_transverse,beta_pass,wake_gate,beta_effective,"
          "v_radial_min,v_radial_max,wake_gate_min,wake_gate_max,"
          "xi_kernel_factor_raw,xi_kernel_metric_scale,xi_kernel_mode,xi_temporal_mode,"
          "theta_trace_4d,invariant_I_4d,r_origin,radial_alignment_to_origin_inward,transverse_fraction_origin,is_near_source,escaped,valid\n";
  init_csv << "step,time,probe_id,x,y,z,vx,vy,vz,ax,ay,az,a_norm,xi_x,xi_y,xi_z,xi_spatial_norm,xi_t,xi_x_base,xi_y_base,xi_z_base,"
              "xi_x_eff,xi_y_eff,xi_z_eff,beta_rel,gamma_rel,v_radial,v_transverse,beta_pass,wake_gate,beta_effective,"
              "v_radial_min,v_radial_max,wake_gate_min,wake_gate_max,"
              "xi_kernel_factor_raw,xi_kernel_metric_scale,xi_kernel_mode,xi_temporal_mode,"
              "theta_trace_4d,invariant_I_4d,r_origin,radial_alignment_to_origin_inward,transverse_fraction_origin,is_near_source,escaped,valid\n";

  std::size_t rows_written = 0;
  std::size_t invalid_row_count = 0;
  double min_radius = std::numeric_limits<double>::infinity();
  double max_radius = 0.0;
  double sum_align = 0.0, sum_transverse = 0.0;
  std::size_t align_count = 0, trans_count = 0;

  auto dump_row = [&](int step_idx, double tcur, const ProbeState& p, std::ostream& out) {
    const double a_norm = std::sqrt(p.ax * p.ax + p.ay * p.ay + p.az * p.az);
    const double r_origin = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    double align = std::numeric_limits<double>::quiet_NaN();
    double transverse = std::numeric_limits<double>::quiet_NaN();
    if (r_origin > 1e-30 && a_norm > 1e-30) {
      const double inx = -p.x / r_origin;
      const double iny = -p.y / r_origin;
      const double inz = -p.z / r_origin;
      align = (p.ax * inx + p.ay * iny + p.az * inz) / a_norm;
      const double clipped = std::max(-1.0, std::min(1.0, align));
      transverse = std::sqrt(std::max(0.0, 1.0 - clipped * clipped));
    }
    out << step_idx << "," << tcur << "," << p.id << ","
        << p.x << "," << p.y << "," << p.z << ","
        << p.vx << "," << p.vy << "," << p.vz << ","
        << p.ax << "," << p.ay << "," << p.az << "," << a_norm << ","
        << p.xi_x << "," << p.xi_y << "," << p.xi_z << "," << p.xi_spatial_norm << ","
        << p.xi_t << ","
        << p.xi_x_base << "," << p.xi_y_base << "," << p.xi_z_base << ","
        << p.xi_x_eff << "," << p.xi_y_eff << "," << p.xi_z_eff << ","
        << p.beta_rel << "," << p.gamma_rel << ","
        << p.v_radial << "," << p.v_transverse << "," << p.beta_pass << "," << p.wake_gate << "," << p.beta_effective << ","
        << p.v_radial_min << "," << p.v_radial_max << "," << p.wake_gate_min << "," << p.wake_gate_max << ","
        << p.xi_kernel_factor_raw << "," << p.xi_kernel_metric_scale << ","
        << kernel_mode << "," << temporal_mode << ","
        << p.theta_trace_4d << "," << p.invariant_I_4d << ","
        << r_origin << "," << align << "," << transverse << ","
        << (p.is_near_source ? 1 : 0) << "," << (p.escaped ? 1 : 0) << "," << (p.valid ? 1 : 0) << "\n";
  };

  for (std::size_t i = 0; i < probes.size(); ++i) eval_probe(probes[i], true);
  for (std::size_t i = 0; i < probes.size(); ++i) {
    dump_row(0, 0.0, probes[i], traj);
    dump_row(0, 0.0, probes[i], init_csv);
    ++rows_written;
    const double r_origin = std::sqrt(probes[i].x * probes[i].x + probes[i].y * probes[i].y + probes[i].z * probes[i].z);
    if (probes[i].valid) {
      min_radius = std::min(min_radius, r_origin);
      max_radius = std::max(max_radius, r_origin);
    } else {
      ++invalid_row_count;
    }
  }

  std::vector<double> ax_old(probes.size(), 0.0), ay_old(probes.size(), 0.0), az_old(probes.size(), 0.0);
  std::vector<double> ax_new(probes.size(), 0.0), ay_new(probes.size(), 0.0), az_new(probes.size(), 0.0);
  for (int step = 1; step <= steps; ++step) {
    if (integrator == "velocity_verlet") {
      for (std::size_t i = 0; i < probes.size(); ++i) {
        ProbeState& p = probes[i];
        ax_old[i] = p.ax;
        ay_old[i] = p.ay;
        az_old[i] = p.az;
        if (!p.valid) continue;
        p.x += p.vx * dt + 0.5 * p.ax * dt * dt;
        p.y += p.vy * dt + 0.5 * p.ay * dt * dt;
        p.z += p.vz * dt + 0.5 * p.az * dt * dt;
      }
      for (std::size_t i = 0; i < probes.size(); ++i) {
        ProbeState& p = probes[i];
        eval_probe(p, true);
        ax_new[i] = p.ax;
        ay_new[i] = p.ay;
        az_new[i] = p.az;
      }
      for (std::size_t i = 0; i < probes.size(); ++i) {
        ProbeState& p = probes[i];
        if (!p.valid) continue;
        p.vx += 0.5 * (ax_old[i] + ax_new[i]) * dt;
        p.vy += 0.5 * (ay_old[i] + ay_new[i]) * dt;
        p.vz += 0.5 * (az_old[i] + az_new[i]) * dt;
      }
    } else {
      for (std::size_t i = 0; i < probes.size(); ++i) {
        ProbeState& p = probes[i];
        eval_probe(p, true);
        if (!p.valid) continue;
        p.vx += p.ax * dt;
        p.vy += p.ay * dt;
        p.vz += p.az * dt;
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
      }
      for (std::size_t i = 0; i < probes.size(); ++i) eval_probe(probes[i], true);
    }

    if ((step % dump_every) == 0) {
      const double tcur = dt * static_cast<double>(step);
      for (std::size_t i = 0; i < probes.size(); ++i) {
        const ProbeState& p = probes[i];
        dump_row(step, tcur, p, traj);
        ++rows_written;
        const double a_norm = std::sqrt(p.ax * p.ax + p.ay * p.ay + p.az * p.az);
        const double r_origin = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        if (!p.valid) {
          ++invalid_row_count;
          continue;
        }
        min_radius = std::min(min_radius, r_origin);
        max_radius = std::max(max_radius, r_origin);
        if (r_origin > 1e-30 && a_norm > 1e-30) {
          const double inx = -p.x / r_origin;
          const double iny = -p.y / r_origin;
          const double inz = -p.z / r_origin;
          const double align = (p.ax * inx + p.ay * iny + p.az * inz) / a_norm;
          const double clipped = std::max(-1.0, std::min(1.0, align));
          const double transverse = std::sqrt(std::max(0.0, 1.0 - clipped * clipped));
          sum_align += align;
          sum_transverse += transverse;
          ++align_count;
          ++trans_count;
        }
      }
    }
  }

  std::size_t invalid_probe_count = 0;
  for (std::size_t i = 0; i < probes.size(); ++i) {
    if (probes[i].escaped || !probes[i].valid) ++invalid_probe_count;
  }
  if (!std::isfinite(min_radius) || min_radius == std::numeric_limits<double>::infinity()) {
    min_radius = std::numeric_limits<double>::quiet_NaN();
    max_radius = std::numeric_limits<double>::quiet_NaN();
  }

  double slope = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> slope_log_r;
  std::vector<double> slope_log_a;
  if (config.tpf_source_benchmark_shape == "monopole") {
    const double r_base = std::max(1.0, probe_radius);
    for (int i = 0; i < 12; ++i) {
      const double r = r_base * (0.7 + 0.3 * static_cast<double>(i));
      if (!(r > 0.0)) continue;
      if (source_exclusion_radius > 0.0 && r <= source_exclusion_radius) continue;
      ProbeState p;
      p.x = r;
      p.y = 0.0;
      p.z = 0.0;
      p.vx = p.vy = p.vz = 0.0;
      eval_probe(p, false);
      const double a = std::sqrt(p.ax * p.ax + p.ay * p.ay + p.az * p.az);
      if (std::isfinite(a) && a > 0.0) {
        slope_log_r.push_back(std::log(r));
        slope_log_a.push_back(std::log(a));
      }
    }
  }
  if (slope_log_r.size() >= 2) {
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (std::size_t i = 0; i < slope_log_r.size(); ++i) {
      sx += slope_log_r[i];
      sy += slope_log_a[i];
      sxx += slope_log_r[i] * slope_log_r[i];
      sxy += slope_log_r[i] * slope_log_a[i];
    }
    const double n = static_cast<double>(slope_log_r.size());
    const double d = n * sxx - sx * sx;
    if (std::fabs(d) > 1e-30) slope = (n * sxy - sx * sy) / d;
  }

  std::ofstream summary(output_dir + "/tpf_4d_xi_motion_probe_summary.txt");
  if (!summary) throw std::runtime_error("failed to open Xi motion benchmark summary");
  summary << std::scientific;
  summary << "mode: tpf_4d_xi_motion_probe_benchmark\n";
  summary << "readout name: GravityXiMotionReadout_v1\n";
  summary << "xi kernel readout name: GravityXiKernelDeformation_v1\n";
  summary << "source shape: " << config.tpf_source_benchmark_shape << "\n";
  summary << "source config id: " << source_config_id << "\n";
  summary << "source masses: ";
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    if (i > 0) summary << ", ";
    summary << source_specs[i].m;
  }
  summary << "\nsource positions: ";
  for (std::size_t i = 0; i < source_specs.size(); ++i) {
    if (i > 0) summary << "; ";
    summary << "(" << source_specs[i].x << "," << source_specs[i].y << "," << source_specs[i].z << ")";
  }
  summary << "\n";
  summary << "tpf_4d_xi_kernel_mode: " << kernel_mode << "\n";
  summary << "tpf_4d_xi_kernel_coupling: " << kernel_coupling << "\n";
  summary << "tpf_4d_xi_kernel_beta_power: " << kernel_beta_power << "\n";
  summary << "tpf_4d_xi_kernel_factor_mode: " << kernel_factor_mode << "\n";
  summary << "tpf_4d_xi_kernel_metric_min: " << kernel_metric_min << "\n";
  summary << "tpf_4d_xi_kernel_metric_max: " << kernel_metric_max << "\n";
  summary << "tpf_4d_xi_temporal_mode: " << temporal_mode << "\n";
  summary << "tpf_4d_xi_temporal_coupling: " << temporal_coupling << "\n";
  summary << "configured source velocity vector: (" << source_speed_x << ", " << source_speed_y << ", " << source_speed_z << ")\n";
  summary << "dt: " << dt << "\n";
  summary << "steps: " << steps << "\n";
  summary << "integrator: " << integrator << "\n";
  summary << "readout scale: " << readout_scale << "\n";
  summary << "field softening: " << field_softening << "\n";
  summary << "source exclusion radius: " << source_exclusion_radius << "\n";
  summary << "probe layout: " << layout << "\n";
  summary << "probe count: " << probes.size() << "\n";
  summary << "probe radius: " << probe_radius << "\n";
  summary << "probe speed: " << probe_speed << "\n";
  summary << "total trajectory rows written: " << rows_written << "\n";
  summary << "escaped/invalid probe count: " << invalid_probe_count << "\n";
  summary << "escaped/invalid trajectory row count: " << invalid_row_count << "\n";
  summary << "minimum radius reached: " << min_radius << "\n";
  summary << "maximum radius reached: " << max_radius << "\n";
  summary << "mean radial alignment of acceleration with inward source/origin direction: "
          << (align_count > 0 ? sum_align / static_cast<double>(align_count) : std::numeric_limits<double>::quiet_NaN()) << "\n";
  summary << "mean transverse acceleration fraction: "
          << (trans_count > 0 ? sum_transverse / static_cast<double>(trans_count) : std::numeric_limits<double>::quiet_NaN()) << "\n";
  summary << "measured Xi acceleration falloff slope from initial probe samples if available: " << slope << "\n";
  const double n_kernel = (kernel_stats.sample_count > 0) ? static_cast<double>(kernel_stats.sample_count) : 1.0;
  summary << "mean beta_rel: " << (kernel_stats.sample_count > 0 ? kernel_stats.sum_beta_rel / n_kernel : 0.0) << "\n";
  summary << "max beta_rel: " << kernel_stats.max_beta_rel << "\n";
  summary << "mean gamma_rel: " << (kernel_stats.sample_count > 0 ? kernel_stats.sum_gamma_rel / n_kernel : 1.0) << "\n";
  summary << "max gamma_rel: " << kernel_stats.max_gamma_rel << "\n";
  summary << "mean factor_raw: " << (kernel_stats.sample_count > 0 ? kernel_stats.sum_factor_raw / n_kernel : 0.0) << "\n";
  summary << "max abs factor_raw: " << kernel_stats.max_abs_factor_raw << "\n";
  summary << "mean factor_clamped_or_metric_scale: " << (kernel_stats.sample_count > 0 ? kernel_stats.sum_metric_scale / n_kernel : 1.0) << "\n";
  summary << "max metric_scale: " << kernel_stats.max_metric_scale << "\n";
  summary << "mean xi_t: " << (kernel_stats.sample_count > 0 ? kernel_stats.sum_xi_t / n_kernel : 0.0) << "\n";
  summary << "max abs xi_t: " << kernel_stats.max_abs_xi_t << "\n";
  summary << "acceleration formula: a=-K_xi*Xi_eff_spatial\n";
  summary << "xi kernel evaluator: per-source softened Xi kernel with optional GravityXiKernelDeformation_v1\n";
  summary << "theta diagnostic evaluator: evaluate_static_sources_field_4d(...)\n";
  summary << "source behavior: fixed sources for Stage 7B\n";
  summary << "probe behavior: moving probes\n";
  summary << "no Newtonian acceleration calls: true\n";
  summary << "no compute_accelerations(...) calls: true\n";
  summary << "no compute_direct_tpf_accelerations(...) calls: true\n";
  summary << "no principal C tensor acceleration: true\n";
  summary << "no VDSG/shunt/cooling: true\n";
  summary << "no production dynamics route: true\n";
  summary << "note: Xi-kernel deformation is applied before acceleration readout\n";
  summary << "note: acceleration remains a=-K_xi*Xi_eff_spatial\n";
  summary << "note: old additive acceleration VDSG path is not used\n";
  summary << "description: dynamic probe-motion benchmark using Xi-direct acceleration readout from fixed-source 4D field evaluation.\n";
}

void TPFCorePackage::write_regime_diagnostics(const std::vector<Snapshot>& snapshots,
                                              const Config& config,
                                              const std::string& output_dir) const {
  using namespace tpfcore;
  if (snapshots.empty()) return;
  if (config.tpf_dynamics_mode == "xi_kernel_deformed" && !config.tpf_xi_kernel_dump_field_diagnostics) return;

  TPFCoreParams params = build_params(config, output_dir);
  double eps = params.effective_source_softening;
  double bh_mass = params.bh_mass;
  bool star_star = params.enable_star_star_gravity;

  double sum_theta_norm = 0.0, sum_I = 0.0;
  double min_theta_norm = 1e300, max_theta_norm = -1e300;
  double min_I = 1e300, max_I = -1e300;
  size_t n_samples = 0;
  size_t count_low = 0, count_transitional = 0, count_high = 0;

  for (const auto& snap : snapshots) {
    const State& s = snap.state;
    for (int i = 0; i < s.n(); ++i) {
      FieldAtPoint field = evaluate_provisional_field_multi_source(s, i, bh_mass, star_star, eps);
      ++xi_runtime_counters_.theta_evaluations;
      ++xi_runtime_counters_.invariant_I_evaluations;
      double tn = theta_frobenius_norm(field.theta);
      double I = field.invariant_I;

      sum_theta_norm += tn;
      sum_I += I;
      if (tn < min_theta_norm) min_theta_norm = tn;
      if (tn > max_theta_norm) max_theta_norm = tn;
      if (I < min_I) min_I = I;
      if (I > max_I) max_I = I;
      ++n_samples;

      if (tn < THETA_NORM_LOW_MAX) ++count_low;
      else if (tn < THETA_NORM_TRANSITIONAL_MAX) ++count_transitional;
      else ++count_high;
    }
  }

  if (n_samples == 0) return;

  double mean_theta_norm = sum_theta_norm / n_samples;
  double mean_I = sum_I / n_samples;

  std::ofstream f(params.output_dir + "/tpf_regime_diagnostics.txt");
  if (!f) return;

  f << "TPFCore regime diagnostics (dynamical run)\n";
  f << "Reporting only; no equation change. Same TPF law across regimes.\n";
  f << "Residual: not available (multi-source superposition; no analytic residual in this path).\n";
  f << "\n--- Field intensity (Theta Frobenius) ---\n";
  f << "  min = " << std::scientific << min_theta_norm << ", max = " << max_theta_norm
    << ", mean = " << mean_theta_norm << "\n";
  f << "\n--- Invariant I ---\n";
  f << "  min = " << min_I << ", max = " << max_I << ", mean = " << mean_I << "\n";
  f << "\n--- Regime distribution (heuristic thresholds) ---\n";
  f << std::fixed << std::setprecision(1);
  f << "  low-intensity: " << count_low << " (" << (100.0 * count_low / n_samples) << "%)\n";
  f << "  transitional: " << count_transitional << " (" << (100.0 * count_transitional / n_samples) << "%)\n";
  f << "  high-intensity (provisional ansatz caution): " << count_high << " (" << (100.0 * count_high / n_samples) << "%)\n";
  f << std::scientific;
  f << "  Thresholds: low < " << THETA_NORM_LOW_MAX << ", transitional < " << THETA_NORM_TRANSITIONAL_MAX << "\n";
  f << "\nTotal sample points: " << n_samples << " (particles x snapshots)\n";
}

TPFCorePackage::TrajectorySummary TPFCorePackage::compute_trajectory_summary(const std::vector<Snapshot>& snapshots) const {
  TrajectorySummary out;
  if (snapshots.empty()) return out;
  const int n_part = snapshots[0].state.n();
  if (n_part != 1) return out;

  std::vector<double> radii;
  std::vector<double> angles;
  radii.reserve(snapshots.size());
  angles.reserve(snapshots.size());
  for (const auto& snap : snapshots) {
    const State& s = snap.state;
    double x = s.x[0], y = s.y[0];
    radii.push_back(std::sqrt(x * x + y * y));
    angles.push_back(std::atan2(y, x));
  }

  double r_initial = radii.front();
  double r_final = radii.back();
  double r_min = r_initial, r_max = r_initial;
  for (double r : radii) {
    if (r < r_min) r_min = r;
    if (r > r_max) r_max = r;
  }
  double radial_drift = r_final - r_initial;

  double total_angle_sweep = 0.0;
  double prev = angles[0];
  for (size_t k = 1; k < angles.size(); ++k) {
    double a = angles[k];
    double d = a - prev;
    if (d > 3.141592653589793) d -= 2.0 * 3.141592653589793;
    if (d < -3.141592653589793) d += 2.0 * 3.141592653589793;
    total_angle_sweep += d;
    prev = a;
  }
  double revolutions = total_angle_sweep / (2.0 * 3.141592653589793);

  const double PLUNGE_R_MIN_FRAC = 0.25;
  const double PLUNGE_FINAL_FRAC = 0.5;
  const double ESCAPE_R_MAX_FRAC = 2.5;
  const double ESCAPE_FINAL_FRAC = 1.5;
  const double BOUNDED_R_MIN_FRAC = 0.2;
  const double BOUNDED_R_MAX_FRAC = 3.0;
  const double NEAR_CIRCULAR_BAND = 0.25;

  bool radius_stays_bounded = (r_min >= BOUNDED_R_MIN_FRAC * r_initial && r_max <= BOUNDED_R_MAX_FRAC * r_initial);

  const char* trajectory_class = "strongly drifting / unclear";
  if (r_min < PLUNGE_R_MIN_FRAC * r_initial && r_final < PLUNGE_FINAL_FRAC * r_initial)
    trajectory_class = "plunge";
  else if (r_max > ESCAPE_R_MAX_FRAC * r_initial && r_final > ESCAPE_FINAL_FRAC * r_initial)
    trajectory_class = "escape";
  else if (radius_stays_bounded) {
    if ((r_max - r_min) / r_initial <= NEAR_CIRCULAR_BAND)
      trajectory_class = "near-circular-candidate";
    else
      trajectory_class = "bounded-candidate";
  }

  out.valid = true;
  out.r_initial = r_initial;
  out.r_final = r_final;
  out.r_min = r_min;
  out.r_max = r_max;
  out.radial_drift = radial_drift;
  out.revolutions = revolutions;
  out.trajectory_class = trajectory_class;
  return out;
}

TPFCorePackage::RegimeSummary TPFCorePackage::compute_regime_summary(const std::vector<Snapshot>& snapshots,
                                                                       const Config& config,
                                                                       const std::string& output_dir) const {
  using namespace tpfcore;
  RegimeSummary out;
  if (snapshots.empty()) return out;
  if (config.tpf_dynamics_mode == "xi_kernel_deformed" && !config.tpf_xi_kernel_dump_field_diagnostics) return out;

  TPFCoreParams params = build_params(config, output_dir);
  double eps = params.effective_source_softening;
  double bh_mass = params.bh_mass;
  bool star_star = params.enable_star_star_gravity;

  double sum_theta_norm = 0.0;
  double min_theta_norm = 1e300, max_theta_norm = -1e300;
  size_t n_samples = 0;
  size_t count_low = 0, count_transitional = 0, count_high = 0;

  for (const auto& snap : snapshots) {
    const State& s = snap.state;
    for (int i = 0; i < s.n(); ++i) {
      FieldAtPoint field = evaluate_provisional_field_multi_source(s, i, bh_mass, star_star, eps);
      ++xi_runtime_counters_.theta_evaluations;
      double tn = theta_frobenius_norm(field.theta);
      sum_theta_norm += tn;
      if (tn < min_theta_norm) min_theta_norm = tn;
      if (tn > max_theta_norm) max_theta_norm = tn;
      ++n_samples;
      if (tn < THETA_NORM_LOW_MAX) ++count_low;
      else if (tn < THETA_NORM_TRANSITIONAL_MAX) ++count_transitional;
      else ++count_high;
    }
  }

  if (n_samples == 0) return out;
  out.valid = true;
  out.mean_theta_norm = sum_theta_norm / n_samples;
  out.max_theta_norm = max_theta_norm;
  out.min_theta_norm = min_theta_norm;
  out.n_samples = n_samples;
  out.frac_low = static_cast<double>(count_low) / n_samples;
  out.frac_transitional = static_cast<double>(count_transitional) / n_samples;
  out.frac_high = static_cast<double>(count_high) / n_samples;
  return out;
}

void TPFCorePackage::write_trajectory_diagnostics(const std::vector<Snapshot>& snapshots,
                                                  const Config& config,
                                                  const std::string& output_dir) const {
  if (snapshots.empty()) return;

  tpfcore::TPFCoreParams params = build_params(config, output_dir);
  const int n_part = snapshots[0].state.n();
  std::ofstream f(params.output_dir + "/tpf_trajectory_diagnostics.txt");
  if (!f) return;

  f << "TPFCore trajectory diagnostics (dynamical run)\n";
  f << "Analysis/reporting only; does not alter integrator or motion law.\n";
  f << "Conservative geometric/time-series heuristics; labeled diagnostic/provisional.\n";

  TrajectorySummary sum = compute_trajectory_summary(snapshots);
  if (!sum.valid) {
    f << "\nTrajectory classification is only computed for single-body dynamical runs (e.g. bh_orbit_validation).\n";
    f << "This run has " << n_part << " particles. No trajectory class assigned.\n";
    return;
  }

  const double BOUNDED_R_MIN_FRAC = 0.2, BOUNDED_R_MAX_FRAC = 3.0;
  bool radius_stays_bounded = (sum.r_min >= BOUNDED_R_MIN_FRAC * sum.r_initial && sum.r_max <= BOUNDED_R_MAX_FRAC * sum.r_initial);

  f << "\n--- Tracked body (particle 0) ---\n";
  f << "  initial_radius = " << std::scientific << sum.r_initial << "\n";
  f << "  final_radius = " << sum.r_final << "\n";
  f << "  min_radius_over_run = " << sum.r_min << "\n";
  f << "  max_radius_over_run = " << sum.r_max << "\n";
  f << "  radial_drift (final - initial) = " << sum.radial_drift << "\n";
  f << "  radius_stays_within_bounded_band = " << (radius_stays_bounded ? "yes" : "no");
  f << " (heuristic: r_min >= " << BOUNDED_R_MIN_FRAC << "*r_initial, r_max <= " << BOUNDED_R_MAX_FRAC << "*r_initial)\n";
  f << "  approximate_revolutions = " << std::fixed << std::setprecision(2) << sum.revolutions;
  f << " (from snapshot angle unwrap; approximate)\n";
  f << std::scientific;
  f << "\n--- Trajectory classification (diagnostic / provisional) ---\n";
  f << "  class = " << sum.trajectory_class << "\n";
  f << "  (Heuristics: plunge / escape / bounded-candidate / near-circular-candidate / strongly drifting / unclear.\n";
  f << "   Does not validate the provisional motion law.)\n";
}

void TPFCorePackage::write_closure_diagnostics(const std::vector<Snapshot>& snapshots,
                                               const Config& config,
                                               const std::string& output_dir) const {
  const bool closure_diag_mode =
      (tpfcore::is_derived_tpf_radial_readout_mode(readout_mode_) || readout_mode_ == "experimental_radial_r_scaling");
  if (!provisional_readout_ || !closure_diag_mode || snapshots.empty())
    return;
  const int n_part = snapshots[0].state.n();
  if (n_part != 1) return;

  tpfcore::TPFCoreParams params = build_params(config, output_dir);
  double eps = params.effective_source_softening;
  double softening = params.softening;
  double bh_mass = params.bh_mass;
  bool star_star = params.enable_star_star_gravity;

  std::ofstream csv(params.output_dir + "/tpf_closure_diagnostics.csv");
  if (!csv) return;
  csv << "time,r,theta_rr,theta_tt,theta_tr,radial_closure,tangential_closure,a_inward,v_radial,v_tangential,sign_radial_acc_vs_radial_vel\n";

  size_t n_inward_radial = 0;
  double sum_abs_tangential = 0.0;
  size_t n_outward_drift = 0;
  size_t n_outward_drift_with_inward_pull = 0;
  size_t n_same_sign = 0, n_opposite_sign = 0;
  size_t n_rows = 0;

  for (const auto& snap : snapshots) {
    const State& s = snap.state;
    double t = snap.time;
    double x = s.x[0], y = s.y[0], vx = s.vx[0], vy = s.vy[0];
    double r2 = x * x + y * y + eps * eps;
    double r = std::sqrt(r2);
    if (r < 1e-30) continue;

    double rx = x / r, ry = y / r;
    double ax = 0.0, ay = 0.0;
    tpfcore::ReadoutDiagnostics diag;
    tpfcore::TpfRadialGravityProfile prof;
    const tpfcore::TpfRadialGravityProfile* prof_ptr = nullptr;
    if (tpfcore::is_derived_tpf_radial_readout_mode(readout_mode_)) {
      prof = tpfcore::build_tpf_gravity_profile(s, bh_mass, derived_poisson_cfg_, eps);
      prof_ptr = &prof;
    }
    tpfcore::compute_provisional_readout_with_diagnostics(
        s, 0, bh_mass, star_star, softening, source_softening_,
        readout_mode_, readout_scale_, theta_tt_scale_, theta_tr_scale_, ax, ay, diag,
        tpfcore::is_derived_tpf_radial_readout_mode(readout_mode_) ? &derived_poisson_cfg_ : nullptr,
        prof_ptr);

    double radial_closure = diag.provisional_radial_readout;
    double tangential_closure = diag.provisional_tangential_readout;
    double a_radial = ax * rx + ay * ry;
    double a_inward = -a_radial;
    double v_radial = vx * rx + vy * ry;
    double v_tangential = x * vy - y * vx;

    const char* sign_str = "zero";
    if (a_radial * v_radial > 1e-30) { sign_str = "same"; ++n_same_sign; }
    else if (a_radial * v_radial < -1e-30) { sign_str = "opposite"; ++n_opposite_sign; }

    if (radial_closure < 0.0) ++n_inward_radial;
    sum_abs_tangential += std::abs(tangential_closure);
    if (v_radial > 0.0) {
      ++n_outward_drift;
      if (a_inward > 0.0) ++n_outward_drift_with_inward_pull;
    }

    csv << std::scientific << t << "," << r << "," << diag.theta_rr << "," << diag.theta_tt << "," << diag.theta_tr
        << "," << radial_closure << "," << tangential_closure << "," << a_inward << "," << v_radial << "," << v_tangential
        << "," << sign_str << "\n";
    ++n_rows;
  }

  if (n_rows == 0) return;

  std::ofstream txt(params.output_dir + "/tpf_closure_diagnostics.txt");
  if (!txt) return;

  txt << "TPFCore closure-term decomposition (single-body)\n";
  txt << "Mode: " << readout_mode_ << ". Diagnostics only; no change to formulas or behavior.\n\n";

  txt << "--- Per-step CSV: tpf_closure_diagnostics.csv ---\n";
  txt << "Columns: time, r, theta_rr, theta_tt, theta_tr, radial_closure, tangential_closure,\n";
  txt << "  a_inward, v_radial, v_tangential, sign_radial_acc_vs_radial_vel (same/opposite/zero)\n";
  txt << "radial_closure = readout radial contribution (provisional_radial_readout); negative = inward.\n";
  txt << "tangential_closure = readout tangential contribution (provisional_tangential_readout).\n\n";

  double frac_inward = 100.0 * n_inward_radial / n_rows;
  double mean_abs_tangential = sum_abs_tangential / n_rows;
  double frac_outward_drift = 100.0 * n_outward_drift / n_rows;
  double frac_outward_with_inward_pull = (n_outward_drift > 0) ? (100.0 * n_outward_drift_with_inward_pull / n_outward_drift) : 0.0;

  txt << "--- Summary statistics ---\n";
  txt << "  Steps with radial_closure < 0 (inward): " << n_inward_radial << " / " << n_rows << " (" << std::fixed << std::setprecision(1) << frac_inward << "%)\n";
  txt << "  Mean |tangential_closure|: " << std::scientific << mean_abs_tangential << "\n";
  txt << "  Steps with v_radial > 0 (outward drift): " << n_outward_drift << " / " << n_rows << " (" << std::fixed << std::setprecision(1) << frac_outward_drift << "%)\n";
  txt << "  Of those, steps with a_inward > 0 (inward pull): " << n_outward_drift_with_inward_pull << " (" << std::fixed << std::setprecision(1) << frac_outward_with_inward_pull << "% of outward-drift steps)\n";
  txt << "  sign_radial_acc_vs_radial_vel: same=" << n_same_sign << ", opposite=" << n_opposite_sign << "\n\n";

  txt << "--- Conservative diagnostic answers ---\n";
  txt << "  Is the radial term mostly inward? ";
  if (frac_inward >= 80.0) txt << "Yes (" << std::fixed << std::setprecision(0) << frac_inward << "% of steps have radial_closure < 0).\n";
  else if (frac_inward >= 50.0) txt << "Partially (radial inward in " << frac_inward << "% of steps).\n";
  else txt << "No (radial inward in only " << frac_inward << "% of steps).\n";

  txt << "  Is the tangential/coherence term large enough to bend the trajectory? ";
  txt << "Mean |tangential_closure| = " << std::scientific << mean_abs_tangential << "; compare to |radial_closure| in CSV for relative size.\n";

  txt << "  Does the tangential term correlate with continued outward drift? ";
  txt << "Check CSV: where v_radial > 0, is tangential_closure positive or negative? (Diagnostic only; no causal claim.)\n";

  txt << "  Is the closure producing mostly outward kinematics despite inward radial pull? ";
  if (n_outward_drift_with_inward_pull > 0 && n_outward_drift > 0 && frac_outward_with_inward_pull > 50.0)
    txt << "In " << std::fixed << std::setprecision(0) << frac_outward_with_inward_pull << "% of outward-drift steps the radial acceleration is inward (a_inward > 0); trajectory continues outward.\n";
  else if (n_outward_drift_with_inward_pull > 0)
    txt << "Some steps show outward drift with inward pull (" << n_outward_drift_with_inward_pull << " steps). See CSV.\n";
  else
    txt << "When v_radial > 0, a_inward is typically not positive; see CSV for details.\n";
}

void TPFCorePackage::write_live_orbit_force_audit(const std::vector<Snapshot>& snapshots,
                                                  const Config& config,
                                                  const std::string& output_dir) const {
  if (!provisional_readout_ || snapshots.empty()) return;
  if (config.simulation_mode != SimulationMode::bh_orbit_validation) return;
  const int n_part = snapshots[0].state.n();
  if (n_part != 1) return;

  PhysicsPackage* newton = get_physics_package("Newtonian");
  if (!newton) return;

  tpfcore::TPFCoreParams params = build_params(config, output_dir);
  double softening = params.softening;
  double bh_mass = params.bh_mass;
  bool star_star = params.enable_star_star_gravity;
  double eps2 = softening * softening;

  std::ofstream csv(params.output_dir + "/tpf_live_orbit_force_audit.csv");
  std::ofstream txt(params.output_dir + "/tpf_live_orbit_force_audit.txt");
  if (!csv || !txt) return;

  csv << "step,time,x,y,vx,vy,radius,v_radial,v_tangential,"
      << "ax_tpf,ay_tpf,a_rad_tpf,a_tan_tpf,"
      << "ax_newt,ay_newt,a_rad_newt,a_tan_newt,"
      << "diff_x,diff_y,diff_rad,diff_tan\n";

  txt << "TPFCore live orbit force audit — bh_orbit_validation (Newtonian vs TPF, same evolving state)\n";
  txt << "Diagnostics only; does not alter integrator or motion law.\n";
  txt << "Positions, velocities, and accelerations are taken from the same snapshots used by the integrator (every validation_snapshot_every steps). For steps 0,1,2,5,10,20,50,100 set validation_snapshot_every=1 for a short run.\n\n";

  bool agree_step0 = true;
  int first_diverge_step = -1;
  double first_diverge_time = 0.0;
  double first_diff_rad = 0.0, first_diff_tan = 0.0;

  const double REL_TOL = 0.05;

  for (const auto& snap : snapshots) {
    const State& s = snap.state;
    double t = snap.time;
    int step = snap.step;
    double x = s.x[0], y = s.y[0], vx = s.vx[0], vy = s.vy[0];

    double r2 = x * x + y * y + eps2;
    double r = std::sqrt(r2);
    double rx = 0.0, ry = 0.0, tx = 0.0, ty = 0.0;
    if (r > 1e-30) {
      rx = x / r;
      ry = y / r;
      tx = -ry;
      ty = rx;
    }
    double v_rad = vx * rx + vy * ry;
    double v_tan = x * vy - y * vx;

    std::vector<double> ax_t(1), ay_t(1), ax_n(1), ay_n(1);
    compute_accelerations(s, bh_mass, softening, star_star, ax_t, ay_t);
    newton->compute_accelerations(s, bh_mass, softening, star_star, ax_n, ay_n);

    double a_rad_tpf = ax_t[0] * rx + ay_t[0] * ry;
    double a_tan_tpf = ax_t[0] * tx + ay_t[0] * ty;
    double a_rad_newt = ax_n[0] * rx + ay_n[0] * ry;
    double a_tan_newt = ax_n[0] * tx + ay_n[0] * ty;

    double diff_x = ax_t[0] - ax_n[0];
    double diff_y = ay_t[0] - ay_n[0];
    double diff_rad = a_rad_tpf - a_rad_newt;
    double diff_tan = a_tan_tpf - a_tan_newt;

    csv << step << "," << std::scientific << t << ","
        << x << "," << y << "," << vx << "," << vy << ","
        << r << "," << v_rad << "," << v_tan << ","
        << ax_t[0] << "," << ay_t[0] << "," << a_rad_tpf << "," << a_tan_tpf << ","
        << ax_n[0] << "," << ay_n[0] << "," << a_rad_newt << "," << a_tan_newt << ","
        << diff_x << "," << diff_y << "," << diff_rad << "," << diff_tan << "\n";

    double scale_rad = std::max(std::abs(a_rad_newt), 1e-12);
    double scale_tan = std::max(std::abs(a_tan_newt), 1e-12);
    bool rad_diff = std::abs(diff_rad) > REL_TOL * scale_rad;
    bool tan_diff = std::abs(diff_tan) > REL_TOL * scale_tan;

    if (step == 0 && (rad_diff || tan_diff))
      agree_step0 = false;

    if (first_diverge_step < 0 && (rad_diff || tan_diff)) {
      first_diverge_step = step;
      first_diverge_time = t;
      first_diff_rad = diff_rad;
      first_diff_tan = diff_tan;
    }
  }

  txt << "--- Explicit diagnostic answers ---\n";
  txt << "  Do TPF and Newtonian agree at step 0 for the actual orbit state?\n";
  txt << "    -> " << (agree_step0 ? "Yes (within 5% in both radial and tangential components).\n"
                                : "No (radial and/or tangential components differ by more than 5%).\n");

  txt << "\n  If yes, at what step do they begin to diverge materially?\n";
  if (first_diverge_step >= 0) {
    txt << "    -> First material divergence at step " << first_diverge_step
        << " (t = " << std::scientific << first_diverge_time << ").\n";
  } else {
    txt << "    -> No material divergence detected within sampled snapshots (<=5% threshold).\n";
  }

  txt << "\n  Is the divergence primarily radial, tangential, or both (at first divergence)?\n";
  if (first_diverge_step >= 0) {
    double abs_rad = std::abs(first_diff_rad);
    double abs_tan = std::abs(first_diff_tan);
    if (abs_rad > 2.0 * abs_tan)
      txt << "    -> Mostly radial.\n";
    else if (abs_tan > 2.0 * abs_rad)
      txt << "    -> Mostly tangential.\n";
    else
      txt << "    -> Both radial and tangential are comparable.\n";
  } else {
    txt << "    -> Not applicable (no material divergence detected at 5% level).\n";
  }

  txt << "\n  Is there evidence that the live orbit run is not using the same effective force path as the static force audit?\n";
  txt << "    Static audits and this live audit both use:\n";
  txt << "      - TPFCorePackage::compute_accelerations for TPF.\n";
  txt << "      - NewtonianPackage::compute_accelerations for the Newtonian benchmark.\n";
  txt << "      - Same softening and radial/tangential decomposition as the simulator.\n";
  txt << "    Any mismatch is therefore due to closure behavior on the evolving state, not a separate force path.\n\n";

  txt << "--- DECISIVE CONCLUSION ---\n";
  if (!agree_step0)
    txt << "  live run already differs at step 0 (beyond 5% in at least one component).\n";
  else if (first_diverge_step >= 0)
    txt << "  live run matches static audit initially; divergence develops later along the trajectory.\n";
  else
    txt << "  live run and static audit are consistent within the 5% threshold over sampled snapshots.\n";
}

void TPFCorePackage::write_step0_orbit_audit(const std::vector<Snapshot>& snapshots,
                                             const Config& config,
                                             const std::string& output_dir) const {
  if (snapshots.empty()) return;
  if (config.tpf_dynamics_mode != "direct_tpf") return;

  PhysicsPackage* newton = get_physics_package("Newtonian");
  if (!newton && config.simulation_mode == SimulationMode::bh_orbit_validation) return;

  tpfcore::TPFCoreParams params = build_params(config, output_dir);
  double softening = params.softening;
  double bh_mass = params.bh_mass;
  bool star_star = params.enable_star_star_gravity;
  const double eps = (source_softening_ > 0.0) ? source_softening_ : softening;
  const double eps2 = softening * softening;
  const double lambda = tpfcore::LAMBDA_4D;
  const double kappa = kappa_;

  std::ofstream txt(params.output_dir + "/tpf_step0_orbit_audit.txt");
  std::ofstream raw_csv(params.output_dir + "/direct_tpf_step0_raw_accel_audit.csv");
  std::ofstream summary(params.output_dir + "/direct_tpf_step0_raw_accel_summary.txt");
  const bool bh_only_direct_tpf_step0_decomp =
      (!star_star && vdsg_coupling_ == 0.0);
  std::ofstream decomp_csv;
  std::ofstream decomp_summary;
  if (bh_only_direct_tpf_step0_decomp) {
    decomp_csv.open(params.output_dir + "/direct_tpf_step0_decomposition_audit.csv");
    decomp_summary.open(params.output_dir + "/direct_tpf_step0_decomposition_summary.txt");
  }
  if (!raw_csv || !summary) return;

  raw_csv << "particle_index,x,y,mass,xi_x,xi_y,xi_norm,theta_xx,theta_xy,theta_yy,theta_trace,"
          << "invariant_I,b_xx,b_xy,b_yy,ax_raw,ay_raw,kappa,ax,ay,ax_newton,ay_newton\n";
  raw_csv << std::scientific << std::setprecision(16);
  summary << std::scientific << std::setprecision(16);
  if (decomp_csv) {
    decomp_csv << "particle_index,x,y,r,xi_x,xi_y,xi_mag,u_x,u_y,theta_xx,theta_xy,theta_yy,theta_zz,"
               << "theta_trace,invariant_I,c_xx,c_xy,c_yy,c_zz,proj_x,proj_y,a_raw_x,a_raw_y,a_raw_mag,"
               << "a_raw_radial_dot,a_newton_x,a_newton_y,a_newton_mag,a_newton_radial_dot,ratio_mag\n";
    decomp_csv << std::scientific << std::setprecision(16);
  }
  if (decomp_summary) {
    decomp_summary << std::scientific << std::setprecision(16);
  }

  const State& s0 = snapshots[0].state;
  std::vector<double> ax_t, ay_t, ax_n, ay_n;
  compute_accelerations(s0, bh_mass, softening, star_star, ax_t, ay_t);
  if (newton) {
    newton->compute_accelerations(s0, bh_mass, softening, star_star, ax_n, ay_n);
  } else {
    ax_n.assign(static_cast<size_t>(s0.n()), 0.0);
    ay_n.assign(static_cast<size_t>(s0.n()), 0.0);
  }

  bool strong_kappa_mismatch = false;
  double ratio_min = std::numeric_limits<double>::infinity();
  double ratio_max = -std::numeric_limits<double>::infinity();
  double ratio_sum = 0.0;
  size_t ratio_count = 0;
  size_t raw_inward_count = 0;
  size_t raw_outward_count = 0;
  double theta_trace_min = std::numeric_limits<double>::infinity();
  double theta_trace_max = -std::numeric_limits<double>::infinity();
  double theta_trace_sum = 0.0;
  double invariant_min = std::numeric_limits<double>::infinity();
  double invariant_max = -std::numeric_limits<double>::infinity();
  double invariant_sum = 0.0;
  for (int i = 0; i < s0.n(); ++i) {
    const tpfcore::FieldAtPoint field =
        tpfcore::evaluate_provisional_field_multi_source(s0, i, bh_mass, star_star, eps);
    const tpfcore::DirectTpfBaselineArtifacts artifacts =
        tpfcore::compute_direct_tpf_baseline_artifacts(field, kappa, lambda);
    const tpfcore::XiDirectedReadoutResult baseline_readout =
        tpfcore::compute_xi_directed_tensor_readout(artifacts.xi, artifacts.principal_cij);
    const tpfcore::Theta3D& theta = field.theta;
    const double theta_trace = field.theta_trace;
    const double invariant_I = field.invariant_I;
    const double b_xx = (theta.xx * theta.xx + theta.xy * theta.xy + theta.xz * theta.xz -
                         lambda * theta_trace * theta.xx - 0.5 * invariant_I);
    const double b_xy = (theta.xx * theta.xy + theta.xy * theta.yy + theta.xz * theta.yz -
                         lambda * theta_trace * theta.xy);
    const double b_yy = (theta.xy * theta.xy + theta.yy * theta.yy + theta.yz * theta.yz -
                         lambda * theta_trace * theta.yy - 0.5 * invariant_I);
    const double xi_x = field.xi.x;
    const double xi_y = field.xi.y;
    const double xi_norm = std::sqrt(xi_x * xi_x + xi_y * xi_y);
    double ax_raw = 0.0;
    double ay_raw = 0.0;
    if (xi_norm > 1e-300) {
      const double u_x = xi_x / xi_norm;
      const double u_y = xi_y / xi_norm;
      ax_raw = -(b_xx * u_x + b_xy * u_y);
      ay_raw = -(b_xy * u_x + b_yy * u_y);
    }

    const double ax_from_raw = kappa * ax_raw;
    const double ay_from_raw = kappa * ay_raw;
    const double ax_newton = ax_n[static_cast<size_t>(i)];
    const double ay_newton = ay_n[static_cast<size_t>(i)];
    raw_csv << i << ',' << s0.x[i] << ',' << s0.y[i] << ',' << s0.mass[i] << ','
            << xi_x << ',' << xi_y << ',' << xi_norm << ','
            << theta.xx << ',' << theta.xy << ',' << theta.yy << ',' << theta_trace << ','
            << invariant_I << ',' << b_xx << ',' << b_xy << ',' << b_yy << ','
            << ax_raw << ',' << ay_raw << ',' << kappa << ','
            << ax_from_raw << ',' << ay_from_raw << ','
            << ax_newton << ',' << ay_newton << '\n';

    if (decomp_csv) {
      const double x = s0.x[i];
      const double y = s0.y[i];
      const double r = std::hypot(x, y);
      const double xi_x = artifacts.xi.xi.x;
      const double xi_y = artifacts.xi.xi.y;
      const double xi_mag = std::hypot(xi_x, xi_y);
      double u_x = 0.0;
      double u_y = 0.0;
      if (xi_mag > 1e-300) {
        u_x = xi_x / xi_mag;
        u_y = xi_y / xi_mag;
      }
      const double proj_x = artifacts.principal_cij.c_xx * u_x + artifacts.principal_cij.c_xy * u_y;
      const double proj_y = artifacts.principal_cij.c_xy * u_x + artifacts.principal_cij.c_yy * u_y;
      const double a_raw_x = baseline_readout.ax;
      const double a_raw_y = baseline_readout.ay;
      const double a_raw_mag = std::hypot(a_raw_x, a_raw_y);
      double r_hat_x = 0.0;
      double r_hat_y = 0.0;
      if (r > 1e-300) {
        r_hat_x = x / r;
        r_hat_y = y / r;
      }
      const double a_raw_radial_dot = a_raw_x * r_hat_x + a_raw_y * r_hat_y;
      const double a_newton_mag = std::hypot(ax_newton, ay_newton);
      const double a_newton_radial_dot = ax_newton * r_hat_x + ay_newton * r_hat_y;
      const double ratio_mag =
          (std::isfinite(a_raw_mag) && std::isfinite(a_newton_mag) && (a_newton_mag > 1e-300))
              ? (a_raw_mag / a_newton_mag)
              : std::numeric_limits<double>::quiet_NaN();
      const bool has_ratio = std::isfinite(ratio_mag);
      decomp_csv << i << ',' << x << ',' << y << ',' << r << ','
                 << xi_x << ',' << xi_y << ',' << xi_mag << ',' << u_x << ',' << u_y << ','
                 << artifacts.theta.theta.xx << ',' << artifacts.theta.theta.xy << ','
                 << artifacts.theta.theta.yy << ',' << artifacts.theta.theta.zz << ','
                 << artifacts.theta_trace.value << ',' << artifacts.invariant_I.value << ','
                 << artifacts.principal_cij.c_xx << ',' << artifacts.principal_cij.c_xy << ','
                 << artifacts.principal_cij.c_yy << ',' << artifacts.principal_cij.c_zz << ','
                 << proj_x << ',' << proj_y << ','
                 << a_raw_x << ',' << a_raw_y << ',' << a_raw_mag << ','
                 << a_raw_radial_dot << ','
                 << ax_newton << ',' << ay_newton << ',' << a_newton_mag << ','
                 << a_newton_radial_dot << ',' << ratio_mag << '\n';

      if (has_ratio) {
        ratio_min = std::min(ratio_min, ratio_mag);
        ratio_max = std::max(ratio_max, ratio_mag);
        ratio_sum += ratio_mag;
        ++ratio_count;
      }
      if (a_raw_radial_dot < 0.0) ++raw_inward_count;
      if (a_raw_radial_dot > 0.0) ++raw_outward_count;
      theta_trace_min = std::min(theta_trace_min, artifacts.theta_trace.value);
      theta_trace_max = std::max(theta_trace_max, artifacts.theta_trace.value);
      theta_trace_sum += artifacts.theta_trace.value;
      invariant_min = std::min(invariant_min, artifacts.invariant_I.value);
      invariant_max = std::max(invariant_max, artifacts.invariant_I.value);
      invariant_sum += artifacts.invariant_I.value;
    }

    const double mag_raw = std::hypot(ax_raw, ay_raw);
    const double mag_a = std::hypot(ax_from_raw, ay_from_raw);
    const double mag_newton = std::hypot(ax_newton, ay_newton);
    const bool has_kx = std::abs(ax_raw) > 1e-300;
    const bool has_ky = std::abs(ay_raw) > 1e-300;
    const bool has_km = mag_raw > 1e-300;
    const double kappa_x = has_kx ? (ax_newton / ax_raw) : 0.0;
    const double kappa_y = has_ky ? (ay_newton / ay_raw) : 0.0;
    const double kappa_mag = has_km ? (mag_newton / mag_raw) : 0.0;
    summary << "particle " << i << ":\n";
    summary << "  |a_raw| = " << mag_raw << "\n";
    summary << "  |a| = " << mag_a << "\n";
    summary << "  |a_newton| = " << mag_newton << "\n";
    summary << "  implied_kappa_from_x = " << (has_kx ? std::to_string(kappa_x) : "undefined (ax_raw≈0)") << "\n";
    summary << "  implied_kappa_from_y = " << (has_ky ? std::to_string(kappa_y) : "undefined (ay_raw≈0)") << "\n";
    summary << "  implied_kappa_from_mag = " << (has_km ? std::to_string(kappa_mag) : "undefined (|a_raw|≈0)") << "\n";
    if (has_kx && has_ky && std::abs(kappa_x - kappa_y) > 0.1 * std::max({std::abs(kappa_x), std::abs(kappa_y), 1.0})) {
      summary << "  WARNING: implied kappa differs strongly between x and y components.\n";
      strong_kappa_mismatch = true;
    }
    if (has_km && ((has_kx && std::abs(kappa_x - kappa_mag) > 0.1 * std::max({std::abs(kappa_x), std::abs(kappa_mag), 1.0})) ||
                   (has_ky && std::abs(kappa_y - kappa_mag) > 0.1 * std::max({std::abs(kappa_y), std::abs(kappa_mag), 1.0})))) {
      summary << "  WARNING: implied kappa magnitude ratio differs strongly from component ratios.\n";
      strong_kappa_mismatch = true;
    }
    summary << "\n";
  }
  if (strong_kappa_mismatch) {
    summary << "Overall: strong implied-kappa mismatch detected across particles/components.\n";
  } else {
    summary << "Overall: no strong implied-kappa mismatch detected by the 10% threshold.\n";
  }
  if (decomp_summary) {
    const double ratio_mean = (ratio_count > 0) ? (ratio_sum / static_cast<double>(ratio_count))
                                                : std::numeric_limits<double>::quiet_NaN();
    const double ratio_min_out = (ratio_count > 0) ? ratio_min : std::numeric_limits<double>::quiet_NaN();
    const double ratio_max_out = (ratio_count > 0) ? ratio_max : std::numeric_limits<double>::quiet_NaN();
    const size_t outward_or_inward_total = raw_outward_count + raw_inward_count;
    const size_t zero_radial_dot_count = static_cast<size_t>(std::max(0, s0.n())) - outward_or_inward_total;
    const double theta_trace_mean =
        (s0.n() > 0) ? (theta_trace_sum / static_cast<double>(s0.n())) : std::numeric_limits<double>::quiet_NaN();
    const double theta_trace_min_out = (s0.n() > 0) ? theta_trace_min : std::numeric_limits<double>::quiet_NaN();
    const double theta_trace_max_out = (s0.n() > 0) ? theta_trace_max : std::numeric_limits<double>::quiet_NaN();
    const double invariant_mean =
        (s0.n() > 0) ? (invariant_sum / static_cast<double>(s0.n())) : std::numeric_limits<double>::quiet_NaN();
    const double invariant_min_out = (s0.n() > 0) ? invariant_min : std::numeric_limits<double>::quiet_NaN();
    const double invariant_max_out = (s0.n() > 0) ? invariant_max : std::numeric_limits<double>::quiet_NaN();
    decomp_summary << "|a_raw|/|a_newton| (finite only): min=" << ratio_min_out
                   << ", mean=" << ratio_mean
                   << ", max=" << ratio_max_out
                   << ", count=" << ratio_count << "\n";
    decomp_summary << "a_raw_radial_dot sign counts: inward=" << raw_inward_count
                   << ", outward=" << raw_outward_count
                   << ", zero=" << zero_radial_dot_count << "\n";
    decomp_summary << "theta_trace: min=" << theta_trace_min_out
                   << ", mean=" << theta_trace_mean
                   << ", max=" << theta_trace_max_out << "\n";
    decomp_summary << "invariant_I: min=" << invariant_min_out
                   << ", mean=" << invariant_mean
                   << ", max=" << invariant_max_out << "\n";
  }

  if (!txt) return;
  if (config.simulation_mode != SimulationMode::bh_orbit_validation) return;
  if (s0.n() != 1) return;

  const State& s = s0;
  double x = s.x[0], y = s.y[0], vx = s.vx[0], vy = s.vy[0];

  double r2 = x * x + y * y + eps2;
  double r = std::sqrt(r2);
  double rx = 0.0, ry = 0.0, tx = 0.0, ty = 0.0;
  if (r > 1e-30) {
    rx = x / r;
    ry = y / r;
    tx = -ry;
    ty = rx;
  }
  double v_radial = vx * rx + vy * ry;
  double v_tangential = x * vy - y * vx;

  compute_accelerations(s, bh_mass, softening, star_star, ax_t, ay_t);
  newton->compute_accelerations(s, bh_mass, softening, star_star, ax_n, ay_n);

  double a_rad_tpf = ax_t[0] * rx + ay_t[0] * ry;
  double a_tan_tpf = ax_t[0] * tx + ay_t[0] * ty;
  double a_rad_newton = ax_n[0] * rx + ay_n[0] * ry;
  double a_tan_newton = ax_n[0] * tx + ay_n[0] * ty;

  double diff_radial = a_rad_tpf - a_rad_newton;
  double diff_tangential = a_tan_tpf - a_tan_newton;

  double ratio_radial = 0.0;
  if (std::abs(a_rad_newton) > 1e-300)
    ratio_radial = std::abs(a_rad_tpf) / std::abs(a_rad_newton);

  tpfcore::ReadoutDiagnostics diag;
  double ax_d = 0.0, ay_d = 0.0;
  tpfcore::TpfRadialGravityProfile step0_prof;
  const tpfcore::TpfRadialGravityProfile* step0_prof_ptr = nullptr;
  if (tpfcore::is_derived_tpf_radial_readout_mode(readout_mode_)) {
    step0_prof = tpfcore::build_tpf_gravity_profile(s, bh_mass, derived_poisson_cfg_, eps);
    step0_prof_ptr = &step0_prof;
  }
  tpfcore::compute_provisional_readout_with_diagnostics(
      s, 0, bh_mass, star_star, softening, source_softening_,
      readout_mode_, readout_scale_, theta_tt_scale_, theta_tr_scale_,
      ax_d, ay_d, diag,
      tpfcore::is_derived_tpf_radial_readout_mode(readout_mode_) ? &derived_poisson_cfg_ : nullptr,
      step0_prof_ptr);
  (void)ax_d;
  (void)ay_d;

  txt << std::scientific << std::setprecision(16);
  txt << "TPFCore step-0 orbit audit (exact initial state of bh_orbit_validation)\n";
  txt << "Diagnostics only; no averaging or summaries. All values exact for step 0.\n\n";

  txt << "x = " << x << "\n";
  txt << "y = " << y << "\n";
  txt << "r = " << r << "\n\n";

  txt << "vx = " << vx << "\n";
  txt << "vy = " << vy << "\n";
  txt << "v_radial = " << v_radial << "\n";
  txt << "v_tangential = " << v_tangential << "\n\n";

  txt << "ax_tpf = " << ax_t[0] << "\n";
  txt << "ay_tpf = " << ay_t[0] << "\n";
  txt << "a_rad_tpf = " << a_rad_tpf << "\n";
  txt << "a_tan_tpf = " << a_tan_tpf << "\n\n";

  txt << "ax_newton = " << ax_n[0] << "\n";
  txt << "ay_newton = " << ay_n[0] << "\n";
  txt << "a_rad_newton = " << a_rad_newton << "\n";
  txt << "a_tan_newton = " << a_tan_newton << "\n\n";

  txt << "diff_radial = " << diff_radial << "\n";
  txt << "diff_tangential = " << diff_tangential << "\n";
  txt << "ratio_radial = |a_rad_tpf| / |a_rad_newton| = " << ratio_radial << "\n\n";

  txt << "theta_rr = " << diag.theta_rr << "\n";
  txt << "theta_tt = " << diag.theta_tt << "\n";
  txt << "theta_tr = " << diag.theta_tr << "\n";
  txt << "provisional_radial_readout = " << diag.provisional_radial_readout << "\n";
  txt << "provisional_tangential_readout = " << diag.provisional_tangential_readout << "\n\n";

  double abs_a_rad_newton = std::abs(a_rad_newton);
  double radial_mismatch_pct = 0.0;
  if (abs_a_rad_newton > 1e-300)
    radial_mismatch_pct = 100.0 * std::abs(diff_radial) / abs_a_rad_newton;

  double abs_a_tan_newton = std::abs(a_tan_newton);
  double tangential_mismatch_pct = 0.0;
  if (abs_a_tan_newton > 1e-300)
    tangential_mismatch_pct = 100.0 * std::abs(diff_tangential) / abs_a_tan_newton;

  txt << "--- Conclusion ---\n";
  txt << "Exact step-0 radial mismatch percentage = " << std::fixed << std::setprecision(4) << radial_mismatch_pct << "%\n";
  txt << "Exact step-0 tangential mismatch percentage = " << tangential_mismatch_pct << "%\n";
  if (tangential_mismatch_pct < 1.0)
    txt << "The mismatch is purely radial (tangential mismatch negligible, < 1%).\n";
  else if (tangential_mismatch_pct < 5.0)
    txt << "The mismatch is mostly radial; tangential mismatch is small but non-negligible.\n";
  else
    txt << "The mismatch is both radial and tangential (tangential non-negligible).\n";

  double r_min = params.tpfcore_probe_radius_min;
  double r_max = params.tpfcore_probe_radius_max;
  bool inside_probe_range = (r >= r_min && r <= r_max);
  bool ratio_near_one = (ratio_radial >= 0.9 && ratio_radial <= 1.1);
  if (inside_probe_range && ratio_near_one)
    txt << "This initial state sits inside the correspondence-calibration sweet spot (r in [" << r_min << ", " << r_max << "], ratio_radial near 1).\n";
  else if (inside_probe_range)
    txt << "This initial state is within the correspondence-calibration probe range r in [" << r_min << ", " << r_max << "] but ratio_radial = " << std::scientific << ratio_radial << " is not near 1 (outside sweet spot).\n";
  else
    txt << "This initial state is outside the correspondence-calibration probe range r in [" << r_min << ", " << r_max << "] (r = " << r << ").\n";
}

void TPFCorePackage::write_accel_pipeline_diagnostics(const std::vector<Snapshot>& snapshots,
                                                      const Config& config,
                                                      const std::string& output_dir) const {
  if (!pipeline_diagnostics_csv_ || !provisional_readout_) return;
  std::ostringstream path;
  path << output_dir << "/tpf_accel_pipeline_diagnostics.csv";
  std::ofstream f(path.str());
  if (!f) return;
  f << "step,time,mean_baseline_accel_mag,mean_vdsg_accel_mag,vdsg_over_baseline_ratio,shunt_events,frac_capped,tpf_global_accel_shunt_enabled,vdsg_pairs_evaluated,vdsg_min_beta_rad,vdsg_max_beta_rad,vdsg_mean_abs_beta_rad\n";
  f << std::scientific << std::setprecision(17);
  for (const auto& sn : snapshots) {
    std::vector<double> ax, ay;
    AccelPipelineStats st;
    eval_accel_pipeline(sn.state, config.bh_mass, config.softening, config.enable_star_star_gravity, ax, ay, &st);
    f << sn.step << ',' << sn.time << ',' << st.mean_baseline_mag << ',' << st.mean_vdsg_mag << ','
      << st.vdsg_over_baseline_ratio << ',' << st.shunt_events_last_step << ',' << st.frac_capped_last_step << ','
      << (shunt_enable_ ? 1 : 0) << '\n';
  }
}

void tpf_test_reset_global_accel_shunt_events() { tpfcore::reset_global_accel_shunt_events(); }

unsigned tpf_test_global_accel_shunt_events() { return tpfcore::global_accel_shunt_events(); }

}  // namespace galaxy
