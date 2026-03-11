/**
 * TPFCore package implementation.
 * Honest primitive TPF: Xi, Theta, I. Hessian-based provisional ansatz. No Newtonian substitution.
 */

#include "tpf_core_package.hpp"
#include "../../config.hpp"
#include "../physics_package.hpp"  /* get_physics_package for Newtonian benchmark and live audits */
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
    txt << "TPFCore weak-field calibration (diagnostic / provisional)\n";
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
    txt << "\nSee tpf_weak_field_calibration_comparison.txt for both closure candidates.\n";
  }

  std::string comp_path = params.output_dir + "/tpf_weak_field_calibration_comparison.txt";
  std::ofstream comp(comp_path);
  if (comp) {
    comp << "TPFCore weak-field calibration: comparison of closure candidates\n";
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
    comp << "--- Shape match verdict ---\n";
    bool tr_ok = results[0].one_constant_sufficient;
    bool exp_ok = results[1].one_constant_sufficient;
    double spread_tr = results[0].ratio_spread;
    double spread_exp = results[1].ratio_spread;
    if (tr_ok && !exp_ok)
      comp << "  tr_coherence_readout has better weak-field shape match (one constant sufficient; experimental_radial_r_scaling ratio varies more with r).\n";
    else if (!tr_ok && exp_ok)
      comp << "  experimental_radial_r_scaling has better weak-field shape match (one constant sufficient; tr_coherence_readout ratio varies more with r).\n";
    else if (tr_ok && exp_ok)
      comp << "  Both closures have acceptable shape (one constant sufficient). Lower spread wins; spread tr=" << std::scientific << spread_tr << " exp=" << spread_exp << ".\n";
    else {
      comp << "  Neither closure has one constant sufficient over the probe range; lower ratio spread is better shape match.\n";
      if (spread_tr < spread_exp)
        comp << "  tr_coherence_readout is the better weak-field shape match (spread " << std::scientific << spread_tr << " < " << spread_exp << ").\n";
      else if (spread_exp < spread_tr)
        comp << "  experimental_radial_r_scaling is the better weak-field shape match (spread " << std::scientific << spread_exp << " < " << spread_tr << ").\n";
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

void TPFCorePackage::write_closure_diagnostics(const std::vector<Snapshot>& snapshots,
                                               const Config& config,
                                               const std::string& output_dir) const {
  const bool closure_diag_mode = (readout_mode_ == "tr_coherence_readout" || readout_mode_ == "experimental_radial_r_scaling");
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
    tpfcore::compute_provisional_readout_with_diagnostics(
        s, 0, bh_mass, star_star, softening, source_softening_, isotropic_c_,
        readout_mode_, readout_scale_, theta_tt_scale_, theta_tr_scale_, ax, ay, diag);

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
  if (config.simulation_mode != SimulationMode::two_body_orbit) return;
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

  txt << "TPFCore live two_body_orbit force audit (Newtonian vs TPF, same evolving state)\n";
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

}  // namespace galaxy
