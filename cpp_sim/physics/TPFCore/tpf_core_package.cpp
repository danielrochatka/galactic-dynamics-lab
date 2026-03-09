/**
 * TPFCore package implementation.
 * Honest primitive TPF: Theta, I. No silent Newtonian substitution.
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

void TPFCorePackage::run_single_source_inspect(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  double r_min = config.tpfcore_probe_radius_min;
  double r_max = config.tpfcore_probe_radius_max;
  int n_samples = config.tpfcore_probe_samples;
  double eps = config.softening;
  double m = config.bh_mass;  // single source at origin

  if (r_min >= r_max || n_samples < 2) return;

  std::vector<double> radii, theta_xx_v, theta_xy_v, theta_yy_v, theta_trace_v, invariant_I_v;

  for (int i = 0; i < n_samples; ++i) {
    double frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
    double r = r_min + frac * (r_max - r_min);
    double x = r, y = 0.0;  // probe along +x axis

    Theta2D theta = provisional_point_source_theta(0, 0, m, x, y, eps);
    double I = compute_invariant_I(theta);

    radii.push_back(r);
    theta_xx_v.push_back(theta.xx);
    theta_xy_v.push_back(theta.xy);
    theta_yy_v.push_back(theta.yy);
    theta_trace_v.push_back(theta.trace());
    invariant_I_v.push_back(I);
  }

  if (config.tpfcore_dump_theta_profile) {
    std::ofstream f(output_dir + "/theta_profile.csv");
    if (f) {
      f << "radius,theta_xx,theta_xy,theta_yy,theta_trace,invariant_I\n";
      for (size_t i = 0; i < radii.size(); ++i) {
        f << std::scientific << radii[i] << "," << theta_xx_v[i] << "," << theta_xy_v[i] << ","
          << theta_yy_v[i] << "," << theta_trace_v[i] << "," << invariant_I_v[i] << "\n";
      }
    }
  }

  if (config.tpfcore_dump_invariant_profile) {
    std::ofstream f(output_dir + "/invariant_profile.csv");
    if (f) {
      f << "radius,invariant_I\n";
      for (size_t i = 0; i < radii.size(); ++i) {
        f << std::scientific << radii[i] << "," << invariant_I_v[i] << "\n";
      }
    }
  }

  {
    std::ofstream f(output_dir + "/field_summary.txt");
    if (f) {
      f << "TPFCore single-source inspection\n";
      f << "Source: (0,0), mass=" << m << "\n";
      f << "Probe: +x axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << "\n";
      f << "Provisional source ansatz: yes (see source_ansatz.hpp)\n";
      f << "Lambda: " << tpfcore::LAMBDA_4D << " (fixed, 4D)\n";
    }
  }
}

void TPFCorePackage::run_symmetric_pair_inspect(const Config& config, const std::string& output_dir) {
  using namespace tpfcore;

  double d = config.validation_symmetric_separation;
  double m = config.star_mass;  // mass per source
  double r_min = config.tpfcore_probe_radius_min;
  double r_max = config.tpfcore_probe_radius_max;
  int n_samples = config.tpfcore_probe_samples;
  double eps = config.softening;

  if (r_min >= r_max || n_samples < 2) return;

  std::vector<double> radii, theta_xx_v, theta_xy_v, theta_yy_v, theta_trace_v, invariant_I_v;

  for (int i = 0; i < n_samples; ++i) {
    double frac = (n_samples > 1) ? static_cast<double>(i) / (n_samples - 1) : 0.0;
    double r = r_min + frac * (r_max - r_min);
    double x = r, y = 0.0;  // probe along +x axis

    Theta2D t1 = provisional_point_source_theta(d, 0, m, x, y, eps);
    Theta2D t2 = provisional_point_source_theta(-d, 0, m, x, y, eps);
    Theta2D theta;
    theta.xx = t1.xx + t2.xx;
    theta.xy = t1.xy + t2.xy;
    theta.yy = t1.yy + t2.yy;

    double I = compute_invariant_I(theta);

    radii.push_back(r);
    theta_xx_v.push_back(theta.xx);
    theta_xy_v.push_back(theta.xy);
    theta_yy_v.push_back(theta.yy);
    theta_trace_v.push_back(theta.trace());
    invariant_I_v.push_back(I);
  }

  if (config.tpfcore_dump_theta_profile) {
    std::ofstream f(output_dir + "/theta_profile.csv");
    if (f) {
      f << "radius,theta_xx,theta_xy,theta_yy,theta_trace,invariant_I\n";
      for (size_t i = 0; i < radii.size(); ++i) {
        f << std::scientific << radii[i] << "," << theta_xx_v[i] << "," << theta_xy_v[i] << ","
          << theta_yy_v[i] << "," << theta_trace_v[i] << "," << invariant_I_v[i] << "\n";
      }
    }
  }

  if (config.tpfcore_dump_invariant_profile) {
    std::ofstream f(output_dir + "/invariant_profile.csv");
    if (f) {
      f << "radius,invariant_I\n";
      for (size_t i = 0; i < radii.size(); ++i) {
        f << std::scientific << radii[i] << "," << invariant_I_v[i] << "\n";
      }
    }
  }

  {
    std::ofstream f(output_dir + "/field_summary.txt");
    if (f) {
      f << "TPFCore symmetric-pair inspection\n";
      f << "Sources: (" << d << ",0) and (-" << d << ",0), mass each=" << m << "\n";
      f << "Probe: +x axis, r in [" << r_min << ", " << r_max << "], n=" << n_samples << "\n";
      f << "Provisional source ansatz: yes (see source_ansatz.hpp)\n";
      f << "Lambda: " << tpfcore::LAMBDA_4D << " (fixed, 4D)\n";
    }
  }
}

}  // namespace galaxy
