/**
 * PROVISIONAL source ansatz implementation.
 * See source_ansatz.hpp for warnings.
 */

#include "source_ansatz.hpp"
#include <cmath>

namespace galaxy {
namespace tpfcore {

Theta2D provisional_point_source_theta(double xs, double ys, double m,
                                       double x, double y, double eps) {
  double dx = x - xs;
  double dy = y - ys;
  double r_sq = dx * dx + dy * dy + eps * eps;
  double r_soft = std::sqrt(r_sq);

  // PROVISIONAL: A=1. Paper does not specify point-source coupling.
  double coeff = m / (r_soft > 1e-20 ? r_soft : 1e-20);

  Theta2D t;
  t.xx = coeff;
  t.yy = coeff;
  t.xy = 0.0;
  return t;
}

double compute_invariant_I(const Theta2D& theta) {
  double Theta_mn_Theta_mn = theta.xx * theta.xx + theta.yy * theta.yy + 2.0 * theta.xy * theta.xy;
  double Theta_trace = theta.trace();
  return Theta_mn_Theta_mn - LAMBDA_4D * Theta_trace * Theta_trace;
}

}  // namespace tpfcore
}  // namespace galaxy
