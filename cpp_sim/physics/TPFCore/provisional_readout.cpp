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
  ax = 0.0;
  ay = 0.0;

  const double eps = effective_eps(source_softening, softening);
  const double x = state.x[i];
  const double y = state.y[i];
  const int n = state.n();

  /* tensor_radial_projection: a += scale * (Theta_s · r_hat) per source.
   * r_hat = unit vector from source to particle.
   * Theta · r_hat = (Theta_xx*rx + Theta_xy*ry, Theta_xy*rx + Theta_yy*ry) / r
   * where (dx,dy) = particle - source, r = |(dx,dy)|_soft. */
  if (readout_mode != "tensor_radial_projection") return;

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
    /* Theta · r_hat */
    double ax_contrib = theta.xx * rx + theta.xy * ry;
    double ay_contrib = theta.xy * rx + theta.yy * ry;

    ax += readout_scale * ax_contrib;
    ay += readout_scale * ay_contrib;
  };

  /* BH at origin */
  if (bh_mass > 0.0) add_contribution(0.0, 0.0, bh_mass);

  /* Other particles */
  if (star_star) {
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      add_contribution(state.x[j], state.y[j], state.mass[j]);
    }
  }
}

}  // namespace tpfcore
}  // namespace galaxy
