/**
 * PROVISIONAL weak-field point-source ansatz implementation.
 * Phi = -M / sqrt(r^2 + eps^2)
 * Xi_i = partial_i Phi, Theta_ij = partial_i partial_j Phi
 *
 * Closed-form derivatives:
 * R = sqrt(dx^2 + dy^2 + eps^2)
 * Phi = -M/R
 * Xi_x = M*dx/R^3, Xi_y = M*dy/R^3
 * Theta_xx = M*(1/R^3 - 3*dx^2/R^5)
 * Theta_xy = -3*M*dx*dy/R^5
 * Theta_yy = M*(1/R^3 - 3*dy^2/R^5)
 */

#include "source_ansatz.hpp"
#include <cmath>

namespace galaxy {
namespace tpfcore {

PointSourceField provisional_point_source_field(double xs, double ys, double m,
                                                double x, double y, double eps) {
  double dx = x - xs;
  double dy = y - ys;
  double r2 = dx * dx + dy * dy + eps * eps;
  double R = std::sqrt(r2);
  if (R < 1e-30) R = 1e-30;

  double R3 = R * R * R;
  double R5 = R3 * R * R;

  PointSourceField out;

  out.xi.x = m * dx / R3;
  out.xi.y = m * dy / R3;

  out.theta.xx = m * (1.0 / R3 - 3.0 * dx * dx / R5);
  out.theta.xy = -3.0 * m * dx * dy / R5;
  out.theta.yy = m * (1.0 / R3 - 3.0 * dy * dy / R5);

  return out;
}

Theta2D provisional_point_source_theta(double xs, double ys, double m,
                                       double x, double y, double eps) {
  return provisional_point_source_field(xs, ys, m, x, y, eps).theta;
}

Residual2D provisional_point_source_residual(double xs, double ys, double m,
                                             double x, double y, double eps) {
  double dx = x - xs;
  double dy = y - ys;
  double r2 = dx * dx + dy * dy + eps * eps;
  double R = std::sqrt(r2);
  if (R < 1e-30) R = 1e-30;

  double R5 = R * R * R * R * R;
  double R7 = R5 * R * R;
  const double lam = LAMBDA_4D;

  /* Third derivatives of Phi: Theta_ij = Phi_ij, so Theta_ij,k = Phi_ijk.
   * Phi = -M/R => closed-form derivatives. */
  double dTh_xx_dx = m * (-9.0 * dx / R5 + 15.0 * dx * dx * dx / R7);
  double dTh_xy_dx = -3.0 * m * (dy / R5 - 5.0 * dx * dx * dy / R7);
  double dTh_xy_dy = -3.0 * m * (dx / R5 - 5.0 * dx * dy * dy / R7);
  double dTh_yy_dy = m * (-9.0 * dy / R5 + 15.0 * dy * dy * dy / R7);
  double dTheta_dx = m * (-12.0 * dx / R5 + 15.0 * r2 * dx / R7);
  double dTheta_dy = m * (-12.0 * dy / R5 + 15.0 * r2 * dy / R7);

  /* R_x = partial_x(Theta_xx - lam*Theta) + partial_y(Theta_xy) */
  /* R_y = partial_x(Theta_xy) + partial_y(Theta_yy - lam*Theta) */
  Residual2D out;
  out.x = dTh_xx_dx - lam * dTheta_dx + dTh_xy_dy;
  out.y = dTh_xy_dx + dTh_yy_dy - lam * dTheta_dy;
  return out;
}

double compute_invariant_I(const Theta2D& theta) {
  double Theta_mn_Theta_mn = theta.xx * theta.xx + theta.yy * theta.yy + 2.0 * theta.xy * theta.xy;
  double Theta_trace = theta.trace();
  return Theta_mn_Theta_mn - LAMBDA_4D * Theta_trace * Theta_trace;
}

}  // namespace tpfcore
}  // namespace galaxy
