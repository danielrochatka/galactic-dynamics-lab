/**
 * TPFCore package implementation.
 * Honest primitive TPF: Xi, Theta, I. Hessian-based provisional ansatz. No Newtonian substitution.
 */

#include "tpf_core_package.hpp"
#include "../../config.hpp"
#include "field_evaluation.hpp"
#include "provisional_readout.hpp"
#include "regime_diagnostics.hpp"
#include "source_ansatz.hpp"
#include "tpf_core_params.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

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
  p.tpfcore_isotropic_correction_c = config.tpfcore_isotropic_correction_c;
  p.tpfcore_probe_radius_min = config.tpfcore_probe_radius_min;
  p.tpfcore_probe_radius_max = config.tpfcore_probe_radius_max;
  p.tpfcore_probe_samples = config.tpfcore_probe_samples;
  p.tpfcore_dump_theta_profile = config.tpfcore_dump_theta_profile;
  p.tpfcore_dump_invariant_profile = config.tpfcore_dump_invariant_profile;
  p.tpfcore_dump_readout_debug = config.tpfcore_dump_readout_debug;
  p.tpfcore_c_sweep_min = config.tpfcore_c_sweep_min;
  p.tpfcore_c_sweep_max = config.tpfcore_c_sweep_max;
  p.tpfcore_c_sweep_steps = config.tpfcore_c_sweep_steps;
  p.tpfcore_c_objective = config.tpfcore_c_objective;
  p.validation_symmetric_separation = config.validation_symmetric_separation;
  return p;
}

TPFCorePackage::TPFCorePackage()
    : provisional_readout_(false),
      readout_mode_("tensor_radial_projection"),
      readout_scale_(1.0),
      theta_tt_scale_(1.0),
      theta_tr_scale_(1.0),
      isotropic_c_(0.0),
      source_softening_(0.0) {}

void TPFCorePackage::init_from_config(const Config& config) {
  provisional_readout_ = config.tpfcore_enable_provisional_readout;
  readout_mode_ = config.tpfcore_readout_mode;
  readout_scale_ = config.tpfcore_readout_scale;
  theta_tt_scale_ = config.tpfcore_theta_tt_scale;
  theta_tr_scale_ = config.tpfcore_theta_tr_scale;
  isotropic_c_ = config.tpfcore_isotropic_correction_c;
  source_softening_ = config.tpfcore_source_softening;  /* 0 => use global softening at runtime */
}

void TPFCorePackage::compute_accelerations(const State& state,
                                            double bh_mass,
                                            double softening,
                                            bool star_star,
                                            std::vector<double>& ax,
                                            std::vector<double>& ay) const {
  if (!provisional_readout_) {
    throw std::runtime_error(
        "TPFCore does not support acceleration readout unless provisional readout is enabled. "
        "Set tpfcore_enable_provisional_readout = true in config, or use Newtonian for dynamics, "
        "or run inspection modes (tpf_single_source_inspect, tpf_symmetric_pair_inspect, tpf_single_source_optimize_c).");
  }

  const int n = state.n();
  ax.assign(n, 0.0);
  ay.assign(n, 0.0);

  for (int i = 0; i < n; ++i) {
    tpfcore::compute_provisional_readout_acceleration(
        state, i, bh_mass, star_star, softening, source_softening_,
        isotropic_c_, readout_mode_, readout_scale_,
        theta_tt_scale_, theta_tr_scale_, ax[i], ay[i]);
  }
}

void TPFCorePackage::write_readout_debug(const std::vector<Snapshot>& snapshots,
                                         const Config& config,
                                         const std::string& output_dir) const {
  tpfcore::TPFCoreParams params = build_params(config, output_dir);
  if (!provisional_readout_ || !params.tpfcore_dump_readout_debug || snapshots.empty())
    return;
  tpfcore::write_readout_debug_csv(snapshots, params.output_dir,
                                   params.softening, params.bh_mass, params.enable_star_star_gravity,
                                   source_softening_, isotropic_c_,
                                   readout_mode_, readout_scale_, theta_tt_scale_, theta_tr_scale_);
}

void TPFCorePackage::run_single_source_inspect(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  TPFCoreParams params = build_params(config, output_dir);
  double r_min = params.tpfcore_probe_radius_min;
  double r_max = params.tpfcore_probe_radius_max;
  int n_samples = params.tpfcore_probe_samples;
  double eps = params.effective_source_softening;
  double m = params.bh_mass;
  double c = params.tpfcore_isotropic_correction_c;

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

    FieldAtPoint field = evaluate_provisional_field_single_source(0, 0, m, x, y, eps, c);

    radii.push_back(r);
    x_v.push_back(x);
    y_v.push_back(y);
    xi_x_v.push_back(field.xi.x);
    xi_y_v.push_back(field.xi.y);
    theta_xx_v.push_back(field.theta.xx);
    theta_xy_v.push_back(field.theta.xy);
    theta_yy_v.push_back(field.theta.yy);
    theta_trace_v.push_back(field.theta.trace());
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
      f << "  exploratory_ansatz:  isotropic correction c=" << c << " (NOT a fundamental constant)\n";
      f << "  inspection:          probe r in [" << r_min << ", " << r_max << "], n=" << n_samples << "\n";
      f << "---\n";
      f << "Ansatz: Phi=-M/sqrt(r^2+eps^2), Xi=grad Phi, Theta=Hess(Phi)+B(r)*delta, B(r)=c*M/(r^2+eps^2)^(3/2)\n";
      f << "Source: (0,0), mass=" << m << "\n";
      f << "Isotropic correction coefficient c=" << std::scientific << c << " (exploratory; not theory)\n";
      f << "Probe: +x axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << "\n";
      f << "Source softening eps=" << eps << " (numerical regularization)\n";
      f << "Provisional weak-field point-source ansatz: yes (exploratory correction)\n";
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
      f << "\nField-equation residual (R_nu = partial_i(Theta_i_nu - lambda*delta_i_nu*Theta)):\n";
      f << "  max|residual_x|=" << std::scientific << max_residual_x_abs << "\n";
      f << "  max|residual_y|=" << max_residual_y_abs << "\n";
      f << "  max residual_norm=" << max_residual_norm << "\n";
      f << "  residual_y near zero on +x axis (y=0 symmetry): " << (residual_y_near_zero ? "yes [OK]\n" : "no [check]\n");
      f << "\nSymmetry expectations on +x axis:\n";
      f << "  theta_xy should be zero (y=0 symmetry): max|theta_xy|=" << max_theta_xy_abs;
      f << (max_theta_xy_abs < 1e-10 ? " [OK]\n" : " [check]\n");
      f << "  theta_xx and theta_yy should differ (radial vs transverse): " << (theta_xx_vs_yy_differ ? "yes [OK]\n" : "no [check]\n");
      f << "  invariant_I should decay smoothly with radius\n";
    }
  }
}

void TPFCorePackage::run_single_source_optimize_c(const Config& config, const std::string& output_dir) {
  /*
   * Exploratory ansatz-tuning: numerically fit c against field-equation residual.
   * Fitted c is NOT a final paper-derived constant. This sweeps c over a range,
   * runs the single-source inspection geometry for each, and selects the c that
   * minimizes the chosen objective (max/mean/l2 residual norm).
   */
  using namespace tpfcore;

  TPFCoreParams params = build_params(config, output_dir);
  double r_min = params.tpfcore_probe_radius_min;
  double r_max = params.tpfcore_probe_radius_max;
  int n_samples = params.tpfcore_probe_samples;
  double eps = params.effective_source_softening;
  double m = params.bh_mass;
  double c_min = params.tpfcore_c_sweep_min;
  double c_max = params.tpfcore_c_sweep_max;
  int n_steps = params.tpfcore_c_sweep_steps;
  const std::string& objective_name = params.tpfcore_c_objective;

  if (r_min >= r_max || n_samples < 2 || n_steps < 2) return;

  /* Resolve objective: we minimize, so lower is better */
  int obj_code = 0;  /* 0=max, 1=mean, 2=l2 */
  if (objective_name == "mean_residual_norm") obj_code = 1;
  else if (objective_name == "l2_residual_norm") obj_code = 2;
  /* else max_residual_norm (default) */

  std::vector<double> c_vals, max_norm_vals, mean_norm_vals, l2_norm_vals, obj_vals;

  for (int step = 0; step < n_steps; ++step) {
    double frac = (n_steps > 1) ? static_cast<double>(step) / (n_steps - 1) : 0.0;
    double c = c_min + frac * (c_max - c_min);

    double sum_norm = 0.0;
    double sum_norm_sq = 0.0;
    double max_norm = 0.0;

    for (int i = 0; i < n_samples; ++i) {
      double r_frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
      double r = r_min + r_frac * (r_max - r_min);
      double x = r, y = 0.0;

      FieldAtPoint field = evaluate_provisional_field_single_source(0, 0, m, x, y, eps, c);
      double rn = field.residual.norm();
      sum_norm += rn;
      sum_norm_sq += rn * rn;
      if (rn > max_norm) max_norm = rn;
    }

    double mean_norm = sum_norm / n_samples;
    double l2_norm = std::sqrt(sum_norm_sq);

    double obj = max_norm;
    if (obj_code == 1) obj = mean_norm;
    else if (obj_code == 2) obj = l2_norm;

    c_vals.push_back(c);
    max_norm_vals.push_back(max_norm);
    mean_norm_vals.push_back(mean_norm);
    l2_norm_vals.push_back(l2_norm);
    obj_vals.push_back(obj);
  }

  /* Find best c (minimum objective) */
  size_t best_idx = 0;
  for (size_t i = 1; i < obj_vals.size(); ++i) {
    if (obj_vals[i] < obj_vals[best_idx]) best_idx = i;
  }
  double best_c = c_vals[best_idx];
  double best_obj = obj_vals[best_idx];

  /* Write c_sweep.csv */
  {
    std::ofstream f(params.output_dir + "/c_sweep.csv");
    if (f) {
      f << "c,max_residual_norm,mean_residual_norm,l2_residual_norm\n";
      for (size_t i = 0; i < c_vals.size(); ++i) {
        f << std::scientific << c_vals[i] << "," << max_norm_vals[i] << ","
          << mean_norm_vals[i] << "," << l2_norm_vals[i] << "\n";
      }
    }
  }

  /* Write c_sweep_summary.txt */
  {
    std::ofstream f(params.output_dir + "/c_sweep_summary.txt");
    if (f) {
      f << "TPFCore c-sweep (exploratory ansatz-tuning)\n";
      f << "Parameter role: c = exploratory ansatz correction (NOT a fundamental constant; NOT fixed theory).\n";
      f << "Numerically fitting c against field-equation residual. Fitted c is NOT a final paper-derived constant.\n";
      f << "Sweep range: c in [" << std::scientific << c_min << ", " << c_max << "]\n";
      f << "Number of steps: " << n_steps << "\n";
      f << "Chosen objective (minimize): " << objective_name << "\n";
      f << "Best c: " << best_c << "\n";
      f << "Best objective value: " << best_obj << "\n";
      f << "Probe geometry: single-source at origin, +x axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << "\n";
    }
  }

  std::cout << "Best c: " << std::scientific << best_c << "\n";
  std::cout << "Best objective value (" << objective_name << "): " << best_obj << "\n";
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
  double c = params.tpfcore_isotropic_correction_c;

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

      FieldAtPoint f1 = evaluate_provisional_field_single_source(d, 0, m, x, y, eps, c);
      FieldAtPoint f2 = evaluate_provisional_field_single_source(-d, 0, m, x, y, eps, c);
      FieldAtPoint field = add_provisional_fields(f1, f2);

      axis_v.push_back(ax);
      x_v.push_back(x);
      y_v.push_back(y);
      xi_x_v.push_back(field.xi.x);
      xi_y_v.push_back(field.xi.y);
      theta_xx_v.push_back(field.theta.xx);
      theta_xy_v.push_back(field.theta.xy);
      theta_yy_v.push_back(field.theta.yy);
      theta_trace_v.push_back(field.theta.trace());
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
      f << "  exploratory_ansatz:  isotropic correction c=" << c << " (NOT a fundamental constant)\n";
      f << "  inspection:          probe r in [" << r_min << ", " << r_max << "], n=" << n_samples << " per axis\n";
      f << "---\n";
      f << "Ansatz: Phi=-M/sqrt(r^2+eps^2), Xi=grad Phi, Theta=Hess(Phi)+B(r)*delta, B(r)=c*M/(r^2+eps^2)^(3/2)\n";
      f << "Source positions: (" << d << ",0) and (-" << d << ",0)\n";
      f << "Source masses: " << m << " each\n";
      f << "Isotropic correction coefficient c=" << std::scientific << c << " (exploratory; not theory)\n";
      f << "Probe geometry: +x axis and +y axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << " per axis\n";
      f << "Source softening eps=" << eps << " (numerical regularization)\n";
      f << "Provisional weak-field point-source ansatz: yes\n";
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
      f << "\nField-equation residual (R_nu = partial_i(Theta_i_nu - lambda*delta_i_nu*Theta)):\n";
      f << "  max|residual_x|=" << std::scientific << max_residual_x_abs << "\n";
      f << "  max|residual_y|=" << max_residual_y_abs << "\n";
      f << "  max residual_norm=" << max_residual_norm << "\n";
      f << "  residual_y near zero on +x axis (y=0 symmetry): " << (residual_y_near_zero_x_axis ? "yes [OK]\n" : "no [check]\n");
      f << "  residual_x near zero on +y axis (x=0 symmetry): " << (residual_x_near_zero_y_axis ? "yes [OK]\n" : "no [check]\n");
      f << "Output files: theta_profile.csv, invariant_profile.csv (axis column: x or y)\n";
    }
  }
}

void TPFCorePackage::write_regime_diagnostics(const std::vector<Snapshot>& snapshots,
                                              const Config& config,
                                              const std::string& output_dir) const {
  using namespace tpfcore;
  if (snapshots.empty()) return;

  TPFCoreParams params = build_params(config, output_dir);
  double eps = params.effective_source_softening;
  double c = params.tpfcore_isotropic_correction_c;
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
      FieldAtPoint field = evaluate_provisional_field_multi_source(s, i, bh_mass, star_star, eps, c);
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

  TPFCoreParams params = build_params(config, output_dir);
  double eps = params.effective_source_softening;
  double c = params.tpfcore_isotropic_correction_c;
  double bh_mass = params.bh_mass;
  bool star_star = params.enable_star_star_gravity;

  double sum_theta_norm = 0.0;
  double min_theta_norm = 1e300, max_theta_norm = -1e300;
  size_t n_samples = 0;

  for (const auto& snap : snapshots) {
    const State& s = snap.state;
    for (int i = 0; i < s.n(); ++i) {
      FieldAtPoint field = evaluate_provisional_field_multi_source(s, i, bh_mass, star_star, eps, c);
      double tn = theta_frobenius_norm(field.theta);
      sum_theta_norm += tn;
      if (tn < min_theta_norm) min_theta_norm = tn;
      if (tn > max_theta_norm) max_theta_norm = tn;
      ++n_samples;
    }
  }

  if (n_samples == 0) return out;
  out.valid = true;
  out.mean_theta_norm = sum_theta_norm / n_samples;
  out.max_theta_norm = max_theta_norm;
  out.min_theta_norm = min_theta_norm;
  out.n_samples = n_samples;
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
    f << "\nTrajectory classification is only computed for single-body dynamical runs (e.g. two_body_orbit).\n";
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

}  // namespace galaxy
