/**
 * TPFCore package implementation.
 * Honest primitive TPF: Xi, Theta, I. Hessian-based provisional ansatz. No Newtonian substitution.
 */

#include "tpf_core_package.hpp"
#include "../../config.hpp"
#include "source_ansatz.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace galaxy {

TPFCorePackage::TPFCorePackage() : provisional_readout_(false) {}

void TPFCorePackage::init_from_config(const Config& config) {
  provisional_readout_ = config.tpfcore_enable_provisional_readout;
}

void TPFCorePackage::compute_accelerations(const State&,
                                            double,
                                            double,
                                            bool,
                                            std::vector<double>& ax,
                                            std::vector<double>& ay) const {
  (void)ax;
  (void)ay;
  throw std::runtime_error(
      "TPFCore does not support acceleration readout. Use Newtonian for dynamics, "
      "or run inspection modes (tpf_single_source_inspect, tpf_symmetric_pair_inspect). "
      "Provisional readout is not yet implemented.");
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

  if (r_min >= r_max || n_samples < 2) return;

  std::vector<double> radii, x_v, y_v, xi_x_v, xi_y_v;
  std::vector<double> theta_xx_v, theta_xy_v, theta_yy_v, theta_trace_v, invariant_I_v;

  double max_theta_xy_abs = 0.0;
  bool theta_xx_vs_yy_differ = false;

  for (int i = 0; i < n_samples; ++i) {
    double frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
    double r = r_min + frac * (r_max - r_min);
    double x = r, y = 0.0;

    PointSourceField f = provisional_point_source_field(0, 0, m, x, y, eps);
    double I = compute_invariant_I(f.theta);

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

    double txy = std::abs(f.theta.xy);
    if (txy > max_theta_xy_abs) max_theta_xy_abs = txy;
    if (std::abs(f.theta.xx - f.theta.yy) > 1e-14) theta_xx_vs_yy_differ = true;
  }

  if (config.tpfcore_dump_theta_profile) {
    std::ofstream f(output_dir + "/theta_profile.csv");
    if (f) {
      f << "radius,x,y,xi_x,xi_y,theta_xx,theta_xy,theta_yy,theta_trace,invariant_I\n";
      for (size_t i = 0; i < radii.size(); ++i) {
        f << std::scientific << radii[i] << "," << x_v[i] << "," << y_v[i] << ","
          << xi_x_v[i] << "," << xi_y_v[i] << ","
          << theta_xx_v[i] << "," << theta_xy_v[i] << "," << theta_yy_v[i] << ","
          << theta_trace_v[i] << "," << invariant_I_v[i] << "\n";
      }
    }
  }

  if (config.tpfcore_dump_invariant_profile) {
    std::ofstream f(output_dir + "/invariant_profile.csv");
    if (f) {
      f << "radius,invariant_I,theta_trace\n";
      for (size_t i = 0; i < radii.size(); ++i) {
        f << std::scientific << radii[i] << "," << invariant_I_v[i] << "," << theta_trace_v[i] << "\n";
      }
    }
  }

  {
    std::ofstream f(output_dir + "/field_summary.txt");
    if (f) {
      f << "TPFCore single-source inspection\n";
      f << "Ansatz: Hessian-based Phi=-M/sqrt(r^2+eps^2), Xi=grad Phi, Theta=Hess Phi\n";
      f << "Source: (0,0), mass=" << m << "\n";
      f << "Probe: +x axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << "\n";
      f << "Source softening eps=" << eps << "\n";
      f << "Provisional weak-field point-source ansatz: yes\n";
      f << "Lambda: " << LAMBDA_4D << " (fixed, 4D)\n";
      f << "\nSymmetry expectations on +x axis:\n";
      f << "  theta_xy should be zero (y=0 symmetry): max|theta_xy|=" << max_theta_xy_abs;
      f << (max_theta_xy_abs < 1e-10 ? " [OK]\n" : " [check]\n");
      f << "  theta_xx and theta_yy should differ (radial vs transverse): " << (theta_xx_vs_yy_differ ? "yes [OK]\n" : "no [check]\n");
      f << "  invariant_I should decay smoothly with radius\n";
    }
  }
}

void TPFCorePackage::run_symmetric_pair_inspect(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  double d = config.validation_symmetric_separation;
  double m = config.star_mass;
  double r_min = config.tpfcore_probe_radius_min;
  double r_max = config.tpfcore_probe_radius_max;
  int n_samples = config.tpfcore_probe_samples;
  double eps = effective_source_softening(config);

  if (r_min >= r_max || n_samples < 2) return;

  std::vector<double> x_v, y_v, xi_x_v, xi_y_v;
  std::vector<double> theta_xx_v, theta_xy_v, theta_yy_v, theta_trace_v, invariant_I_v;
  std::vector<std::string> axis_v;

  auto append = [&](const char* ax, double (*px)(double), double (*py)(double)) {
    for (int i = 0; i < n_samples; ++i) {
      double frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
      double r = r_min + frac * (r_max - r_min);
      double x = px(r), y = py(r);

      PointSourceField f1 = provisional_point_source_field(d, 0, m, x, y, eps);
      PointSourceField f2 = provisional_point_source_field(-d, 0, m, x, y, eps);

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
    }
  };

  append("x", [](double r) { return r; }, [](double) { return 0.0; });
  append("y", [](double) { return 0.0; }, [](double r) { return r; });

  if (config.tpfcore_dump_theta_profile) {
    std::ofstream f(output_dir + "/theta_profile.csv");
    if (f) {
      f << "axis,radius,x,y,xi_x,xi_y,theta_xx,theta_xy,theta_yy,theta_trace,invariant_I\n";
      for (size_t i = 0; i < axis_v.size(); ++i) {
        double r = (axis_v[i] == "x") ? x_v[i] : y_v[i];
        f << axis_v[i] << "," << std::scientific << r << "," << x_v[i] << "," << y_v[i] << ","
          << xi_x_v[i] << "," << xi_y_v[i] << ","
          << theta_xx_v[i] << "," << theta_xy_v[i] << "," << theta_yy_v[i] << ","
          << theta_trace_v[i] << "," << invariant_I_v[i] << "\n";
      }
    }
  }

  if (config.tpfcore_dump_invariant_profile) {
    std::ofstream f(output_dir + "/invariant_profile.csv");
    if (f) {
      f << "axis,radius,invariant_I,theta_trace\n";
      for (size_t i = 0; i < axis_v.size(); ++i) {
        double r = (axis_v[i] == "x") ? x_v[i] : y_v[i];
        f << axis_v[i] << "," << std::scientific << r << "," << invariant_I_v[i] << "," << theta_trace_v[i] << "\n";
      }
    }
  }

  {
    std::ofstream f(output_dir + "/field_summary.txt");
    if (f) {
      f << "TPFCore symmetric-pair inspection\n";
      f << "Ansatz: Hessian-based Phi=-M/sqrt(r^2+eps^2), Xi=grad Phi, Theta=Hess Phi\n";
      f << "Source positions: (" << d << ",0) and (-" << d << ",0)\n";
      f << "Source masses: " << m << " each\n";
      f << "Probe geometry: +x axis and +y axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << " per axis\n";
      f << "Source softening eps=" << eps << "\n";
      f << "Provisional weak-field point-source ansatz: yes\n";
      f << "Lambda: " << LAMBDA_4D << " (fixed, 4D)\n";
      f << "Output files: theta_profile.csv, invariant_profile.csv (axis column: x or y)\n";
    }
  }
}

}  // namespace galaxy
