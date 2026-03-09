/**
 * PROVISIONAL motion/readout layer for TPFCore.
 *
 * EXPLORATORY: This is NOT the full derived TPF dynamics.
 *
 * What is exploratory: The mapping Theta_ij → acceleration via Theta·r_hat.
 * The paper gives Xi, Theta, I; the motion law is underspecified. This readout
 * tests tensor-driven motion without reverting to Newtonian −grad(Phi).
 *
 * What we use from the ansatz: Theta_ij only (no Phi for acceleration).
 * Easy to swap out when proper TPF motion laws are derived.
 */

#include "provisional_readout.hpp"
#include "source_ansatz.hpp"
#include "../../types.hpp"
#include <cmath>

namespace galaxy {
namespace tpfcore {

static double effective_eps(double source_softening, double global_softening) {
  return (source_softening > 0.0) ? source_softening : global_softening;
}

static bool is_negated_mode(const std::string& mode) {
  return mode == "tensor_radial_projection_negated";
}

static void compute_raw_readout(const State& state,
                                 int i,
                                 double bh_mass,
                                 bool star_star,
                                 double eps,
                                 double c,
                                 double readout_scale,
                                 double& ax,
                                 double& ay,
                                 Theta2D* theta_sum,
                                 bool* has_theta_sum) {
  ax = 0.0;
  ay = 0.0;
  if (theta_sum) {
    theta_sum->xx = theta_sum->xy = theta_sum->yy = 0.0;
    *has_theta_sum = false;
  }

  const double x = state.x[i];
  const double y = state.y[i];
  const int n = state.n();

  auto add_contribution = [&](double xs, double ys, double m) {
    if (m <= 0.0) return;
    double dx = x - xs;
    double dy = y - ys;
    double r2 = dx * dx + dy * dy + eps * eps;
    double r = std::sqrt(r2);
    if (r < 1e-30) return;
    double rx = dx / r;
    double ry = dy / r;

    Theta2D theta = provisional_point_source_theta(xs, ys, m, x, y, eps, c);
    double ax_contrib = theta.xx * rx + theta.xy * ry;
    double ay_contrib = theta.xy * rx + theta.yy * ry;

    ax += readout_scale * ax_contrib;
    ay += readout_scale * ay_contrib;

    if (theta_sum) {
      theta_sum->xx += theta.xx;
      theta_sum->xy += theta.xy;
      theta_sum->yy += theta.yy;
      *has_theta_sum = true;
    }
  };

  if (bh_mass > 0.0) add_contribution(0.0, 0.0, bh_mass);
  if (star_star) {
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      add_contribution(state.x[j], state.y[j], state.mass[j]);
    }
  }
}

void compute_provisional_readout_acceleration(const State& state,
                                               int i,
                                               double bh_mass,
                                               bool star_star,
                                               double softening,
                                               double source_softening,
                                               double c,
                                               const std::string& readout_mode,
                                               double readout_scale,
                                               double& ax,
                                               double& ay) {
  if (readout_mode != "tensor_radial_projection" && readout_mode != "tensor_radial_projection_negated")
    return;

  const double eps = effective_eps(source_softening, softening);
  compute_raw_readout(state, i, bh_mass, star_star, eps, c, readout_scale, ax, ay, nullptr, nullptr);

  if (is_negated_mode(readout_mode)) {
    ax = -ax;
    ay = -ay;
  }
}

void compute_provisional_readout_with_diagnostics(const State& state,
                                                   int i,
                                                   double bh_mass,
                                                   bool star_star,
                                                   double softening,
                                                   double source_softening,
                                                   double c,
                                                   const std::string& readout_mode,
                                                   double readout_scale,
                                                   double& ax,
                                                   double& ay,
                                                   ReadoutDiagnostics& diag) {
  if (readout_mode != "tensor_radial_projection" && readout_mode != "tensor_radial_projection_negated") {
    ax = ay = 0.0;
    diag.theta_xx = diag.theta_xy = diag.theta_yy = diag.theta_trace = diag.invariant_I = 0.0;
    return;
  }

  const double eps = effective_eps(source_softening, softening);
  Theta2D theta_sum;
  bool has_theta;
  compute_raw_readout(state, i, bh_mass, star_star, eps, c, readout_scale, ax, ay, &theta_sum, &has_theta);

  if (is_negated_mode(readout_mode)) {
    ax = -ax;
    ay = -ay;
  }

  diag.ax = ax;
  diag.ay = ay;
  if (has_theta) {
    diag.theta_xx = theta_sum.xx;
    diag.theta_xy = theta_sum.xy;
    diag.theta_yy = theta_sum.yy;
    diag.theta_trace = theta_sum.trace();
    diag.invariant_I = compute_invariant_I(theta_sum);
  } else {
    diag.theta_xx = diag.theta_xy = diag.theta_yy = diag.theta_trace = diag.invariant_I = 0.0;
  }
}

}  // namespace tpfcore
}  // namespace galaxy
