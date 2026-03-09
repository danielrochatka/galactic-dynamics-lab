/**
 * PROVISIONAL weak-field point-source ansatz implementation.
 * Phi = -M / sqrt(r^2 + eps^2)
 * Xi_i = partial_i Phi (unchanged)
 * Theta_ij = Hess_ij(Phi) + B(r) delta_ij
 * B(r) = c*M/(r^2+eps^2)^(3/2) = c*M/R^3 (exploratory isotropic residual-reduction term)
 *
 * Closed-form: R = sqrt(dx^2+dy^2+eps^2), Hess_xx = M*(1/R^3 - 3*dx^2/R^5), etc.
 */

#include "source_ansatz.hpp"
#include <cmath>

namespace galaxy {
namespace tpfcore {

PointSourceField provisional_point_source_field(double xs, double ys, double m,
                                                double x, double y, double eps, double c) {
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

  /* Hessian of Phi */
  double Hess_xx = m * (1.0 / R3 - 3.0 * dx * dx / R5);
  double Hess_xy = -3.0 * m * dx * dy / R5;
  double Hess_yy = m * (1.0 / R3 - 3.0 * dy * dy / R5);

  /* Isotropic correction B(r) = c*M/R^3 on diagonal only */
  double B = (c != 0.0) ? (c * m / R3) : 0.0;
  out.theta.xx = Hess_xx + B;
  out.theta.xy = Hess_xy;
  out.theta.yy = Hess_yy + B;

  return out;
}

Theta2D provisional_point_source_theta(double xs, double ys, double m,
                                       double x, double y, double eps, double c) {
  return provisional_point_source_field(xs, ys, m, x, y, eps, c).theta;
}

Residual2D provisional_point_source_residual(double xs, double ys, double m,
                                             double x, double y, double eps, double c) {
  double dx = x - xs;
  double dy = y - ys;
  double r2 = dx * dx + dy * dy + eps * eps;
  double R = std::sqrt(r2);
  if (R < 1e-30) R = 1e-30;

  double R5 = R * R * R * R * R;
  double R7 = R5 * R * R;
  const double lam = LAMBDA_4D;

  /* Third derivatives of Phi: Hess_ij,k = Phi_ijk. Theta_ij = Hess_ij + B*delta_ij. */
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

  /* Correction from B(r): Theta gets +B on diagonal, so Theta += 2B, and
   * Theta_xx-lam*Theta gains B*(1-2*lam). With lam=1/4: (1-2*lam)=0.5.
   * partial_x(B) = -3*c*m*dx/R^5, partial_y(B) = -3*c*m*dy/R^5. */
  if (c != 0.0) {
    const double fac = (1.0 - 2.0 * lam) * (-3.0 * c * m / R5);
    out.x += fac * dx;
    out.y += fac * dy;
  }
  return out;
}

double compute_invariant_I(const Theta2D& theta) {
  double Theta_mn_Theta_mn = theta.xx * theta.xx + theta.yy * theta.yy + 2.0 * theta.xy * theta.xy;
  double Theta_trace = theta.trace();
  return Theta_mn_Theta_mn - LAMBDA_4D * Theta_trace * Theta_trace;
}

}  // namespace tpfcore
}  // namespace galaxy
