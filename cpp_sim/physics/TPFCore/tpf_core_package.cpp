/**
 * TPFCore package implementation.
 * Honest primitive TPF: Xi, Theta, I. Hessian-based provisional ansatz. No Newtonian substitution.
 */

#include "tpf_core_package.hpp"
#include "../../config.hpp"
#include "provisional_readout.hpp"
#include "source_ansatz.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace galaxy {

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
  if (!provisional_readout_ || !config.tpfcore_dump_readout_debug || snapshots.empty())
    return;
  tpfcore::write_readout_debug_csv(snapshots, output_dir,
                                   config.softening, config.bh_mass, config.enable_star_star_gravity,
                                   source_softening_, isotropic_c_,
                                   readout_mode_, readout_scale_, theta_tt_scale_, theta_tr_scale_);
}

static double effective_source_softening(const Config& config) {
  return (config.tpfcore_source_softening > 0.0) ? config.tpfcore_source_softening : config.softening;
}

void TPFCorePackage::run_single_source_inspect(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  double r_min = config.tpfcore_probe_radius_min;
  double r_max = config.tpfcore_probe_radius_max;
  int n_samples = config.tpfcore_probe_samples;
  double eps = effective_source_softening(config);
  double m = config.bh_mass;
  double c = config.tpfcore_isotropic_correction_c;

  if (r_min >= r_max || n_samples < 2) return;

  std::vector<double> radii, x_v, y_v, xi_x_v, xi_y_v;
  std::vector<double> theta_xx_v, theta_xy_v, theta_yy_v, theta_trace_v, invariant_I_v;
  std::vector<double> residual_x_v, residual_y_v, residual_norm_v;

  double max_theta_xy_abs = 0.0;
  bool theta_xx_vs_yy_differ = false;
  double max_residual_x_abs = 0.0, max_residual_y_abs = 0.0, max_residual_norm = 0.0;

  for (int i = 0; i < n_samples; ++i) {
    double frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
    double r = r_min + frac * (r_max - r_min);
    double x = r, y = 0.0;

    PointSourceField f = provisional_point_source_field(0, 0, m, x, y, eps, c);
    double I = compute_invariant_I(f.theta);
    Residual2D res = provisional_point_source_residual(0, 0, m, x, y, eps, c);

    radii.push_back(r);
    x_v.push_back(x);
    y_v.push_back(y);
    xi_x_v.push_back(f.xi.x);
    xi_y_v.push_back(f.xi.y);
    theta_xx_v.push_back(f.theta.xx);
    theta_xy_v.push_back(f.theta.xy);
    theta_yy_v.push_back(f.theta.yy);
    theta_trace_v.push_back(f.theta.trace());
    invariant_I_v.push_back(I);
    residual_x_v.push_back(res.x);
    residual_y_v.push_back(res.y);
    residual_norm_v.push_back(res.norm());

    double txy = std::abs(f.theta.xy);
    if (txy > max_theta_xy_abs) max_theta_xy_abs = txy;
    if (std::abs(f.theta.xx - f.theta.yy) > 1e-14) theta_xx_vs_yy_differ = true;
    double rax = std::abs(res.x), ray = std::abs(res.y), rn = res.norm();
    if (rax > max_residual_x_abs) max_residual_x_abs = rax;
    if (ray > max_residual_y_abs) max_residual_y_abs = ray;
    if (rn > max_residual_norm) max_residual_norm = rn;
  }

  bool residual_y_near_zero = (max_residual_y_abs < 1e-10);

  if (config.tpfcore_dump_theta_profile) {
    std::ofstream f(output_dir + "/theta_profile.csv");
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

  if (config.tpfcore_dump_invariant_profile) {
    std::ofstream f(output_dir + "/invariant_profile.csv");
    if (f) {
      f << "radius,invariant_I,theta_trace,residual_norm\n";
      for (size_t i = 0; i < radii.size(); ++i) {
        f << std::scientific << radii[i] << "," << invariant_I_v[i] << "," << theta_trace_v[i]
          << "," << residual_norm_v[i] << "\n";
      }
    }
  }

  {
    std::ofstream f(output_dir + "/field_summary.txt");
    if (f) {
      f << "TPFCore single-source inspection\n";
      f << "Ansatz: Phi=-M/sqrt(r^2+eps^2), Xi=grad Phi, Theta=Hess(Phi)+B(r)*delta, B(r)=c*M/(r^2+eps^2)^(3/2)\n";
      f << "Source: (0,0), mass=" << m << "\n";
      f << "Isotropic correction coefficient c=" << std::scientific << c << "\n";
      f << "Probe: +x axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << "\n";
      f << "Source softening eps=" << eps << "\n";
      f << "Provisional weak-field point-source ansatz: yes (exploratory correction)\n";
      f << "Lambda: " << LAMBDA_4D << " (fixed, 4D)\n";
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

  double r_min = config.tpfcore_probe_radius_min;
  double r_max = config.tpfcore_probe_radius_max;
  int n_samples = config.tpfcore_probe_samples;
  double eps = effective_source_softening(config);
  double m = config.bh_mass;
  double c_min = config.tpfcore_c_sweep_min;
  double c_max = config.tpfcore_c_sweep_max;
  int n_steps = config.tpfcore_c_sweep_steps;
  const std::string& objective_name = config.tpfcore_c_objective;

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

      Residual2D res = provisional_point_source_residual(0, 0, m, x, y, eps, c);
      double rn = res.norm();
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
    std::ofstream f(output_dir + "/c_sweep.csv");
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
    std::ofstream f(output_dir + "/c_sweep_summary.txt");
    if (f) {
      f << "TPFCore c-sweep (exploratory ansatz-tuning)\n";
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

  double d = config.validation_symmetric_separation;
  double m = config.star_mass;
  double r_min = config.tpfcore_probe_radius_min;
  double r_max = config.tpfcore_probe_radius_max;
  int n_samples = config.tpfcore_probe_samples;
  double eps = effective_source_softening(config);
  double c = config.tpfcore_isotropic_correction_c;

  if (r_min >= r_max || n_samples < 2) return;

  std::vector<double> x_v, y_v, xi_x_v, xi_y_v;
  std::vector<double> theta_xx_v, theta_xy_v, theta_yy_v, theta_trace_v, invariant_I_v;
  std::vector<double> residual_x_v, residual_y_v, residual_norm_v;
  std::vector<std::string> axis_v;

  double max_residual_x_abs = 0.0, max_residual_y_abs = 0.0, max_residual_norm = 0.0;
  double max_residual_x_abs_x_axis = 0.0, max_residual_y_abs_x_axis = 0.0;
  double max_residual_x_abs_y_axis = 0.0, max_residual_y_abs_y_axis = 0.0;

  auto append = [&](const char* ax, double (*px)(double), double (*py)(double)) {
    for (int i = 0; i < n_samples; ++i) {
      double frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
      double r = r_min + frac * (r_max - r_min);
      double x = px(r), y = py(r);

      PointSourceField f1 = provisional_point_source_field(d, 0, m, x, y, eps, c);
      PointSourceField f2 = provisional_point_source_field(-d, 0, m, x, y, eps, c);
      Residual2D res1 = provisional_point_source_residual(d, 0, m, x, y, eps, c);
      Residual2D res2 = provisional_point_source_residual(-d, 0, m, x, y, eps, c);
      Residual2D res = { res1.x + res2.x, res1.y + res2.y };

      Xi2D xi = { f1.xi.x + f2.xi.x, f1.xi.y + f2.xi.y };
      Theta2D theta = {
        f1.theta.xx + f2.theta.xx,
        f1.theta.xy + f2.theta.xy,
        f1.theta.yy + f2.theta.yy
      };
      double I = compute_invariant_I(theta);

      axis_v.push_back(ax);
      x_v.push_back(x);
      y_v.push_back(y);
      xi_x_v.push_back(xi.x);
      xi_y_v.push_back(xi.y);
      theta_xx_v.push_back(theta.xx);
      theta_xy_v.push_back(theta.xy);
      theta_yy_v.push_back(theta.yy);
      theta_trace_v.push_back(theta.trace());
      invariant_I_v.push_back(I);
      residual_x_v.push_back(res.x);
      residual_y_v.push_back(res.y);
      residual_norm_v.push_back(res.norm());

      double rax = std::abs(res.x), ray = std::abs(res.y), rn = res.norm();
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
    }
  };

  append("x", [](double r) { return r; }, [](double) { return 0.0; });
  append("y", [](double) { return 0.0; }, [](double r) { return r; });

  bool residual_y_near_zero_x_axis = (max_residual_y_abs_x_axis < 1e-10);
  bool residual_x_near_zero_y_axis = (max_residual_x_abs_y_axis < 1e-10);

  if (config.tpfcore_dump_theta_profile) {
    std::ofstream f(output_dir + "/theta_profile.csv");
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

  if (config.tpfcore_dump_invariant_profile) {
    std::ofstream f(output_dir + "/invariant_profile.csv");
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
    std::ofstream f(output_dir + "/field_summary.txt");
    if (f) {
      f << "TPFCore symmetric-pair inspection\n";
      f << "Ansatz: Phi=-M/sqrt(r^2+eps^2), Xi=grad Phi, Theta=Hess(Phi)+B(r)*delta, B(r)=c*M/(r^2+eps^2)^(3/2)\n";
      f << "Source positions: (" << d << ",0) and (-" << d << ",0)\n";
      f << "Source masses: " << m << " each\n";
      f << "Isotropic correction coefficient c=" << std::scientific << c << "\n";
      f << "Probe geometry: +x axis and +y axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << " per axis\n";
      f << "Source softening eps=" << eps << "\n";
      f << "Provisional weak-field point-source ansatz: yes\n";
      f << "Lambda: " << LAMBDA_4D << " (fixed, 4D)\n";
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

}  // namespace galaxy
