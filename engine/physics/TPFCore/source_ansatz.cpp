/**
 * PROVISIONAL weak-field point-source: 3D Hessian of Phi = -m/R on z = 0,
 * R^2 = dx^2 + dy^2 + eps^2.
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

  /* 3D Hessian of Phi = -m/R at z = 0 (dz = 0): Theta_ij = m (R^2 delta_ij - 3 s_i s_j) / R^5 */
  out.theta.xx = m * (r2 - 3.0 * dx * dx) / R5;
  out.theta.xy = -3.0 * m * dx * dy / R5;
  out.theta.xz = 0.0;
  out.theta.yy = m * (r2 - 3.0 * dy * dy) / R5;
  out.theta.yz = 0.0;
  out.theta.zz = m * r2 / R5; /* = m / R^3 */

  return out;
}

Theta3D provisional_point_source_theta(double xs, double ys, double m,
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

  double dTh_xx_dx = m * (-9.0 * r2 * dx + 15.0 * dx * dx * dx) / R7;
  double dTh_yy_dy = m * (-9.0 * r2 * dy + 15.0 * dy * dy * dy) / R7;
  double dTh_xy_dx = -3.0 * m * dy * (r2 - 5.0 * dx * dx) / R7;
  double dTh_xy_dy = -3.0 * m * dx * (r2 - 5.0 * dy * dy) / R7;

  double dTheta_dx = -15.0 * m * dx * eps * eps / R7;
  double dTheta_dy = -15.0 * m * dy * eps * eps / R7;

  /* partial_z Theta_xz and Theta_yz at z=0: Theta_{xz} = -3 m dx z / R^5, etc. */
  double dTh_xz_dz = -3.0 * m * dx * r2 / R7;
  double dTh_yz_dz = -3.0 * m * dy * r2 / R7;

  Residual2D out;
  out.x = dTh_xx_dx - lam * dTheta_dx + dTh_xy_dy + dTh_xz_dz;
  out.y = dTh_xy_dx + dTh_yy_dy - lam * dTheta_dy + dTh_yz_dz;
  return out;
}

double compute_invariant_I(const Theta3D& theta) {
  double Theta_mn_Theta_mn =
      theta.xx * theta.xx + theta.yy * theta.yy + theta.zz * theta.zz +
      2.0 * (theta.xy * theta.xy + theta.xz * theta.xz + theta.yz * theta.yz);
  double Theta_trace = theta.trace();
  return Theta_mn_Theta_mn - LAMBDA_4D * Theta_trace * Theta_trace;
}

}  // namespace tpfcore
}  // namespace galaxy
