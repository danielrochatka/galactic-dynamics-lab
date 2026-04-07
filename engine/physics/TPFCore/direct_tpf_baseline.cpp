#include "direct_tpf_baseline.hpp"

#include <cmath>

namespace galaxy {
namespace tpfcore {

Eq10ThetaTrace compute_theta_trace(const Eq10ThetaTensor& theta_tensor) {
  Eq10ThetaTrace theta_trace{};
  theta_trace.value = theta_tensor.theta.trace();
  return theta_trace;
}

Eq10InvariantI compute_invariant_I(const Eq10ThetaTensor& theta_tensor) {
  Eq10InvariantI invariant_I{};
  invariant_I.value = tpfcore::compute_invariant_I(theta_tensor.theta);
  return invariant_I;
}

Eq10DeltaCPlaceholder compute_deltaC_placeholder_zero() {
  Eq10DeltaCPlaceholder delta_c{};
  delta_c.xx = 0.0;
  delta_c.xy = 0.0;
  delta_c.yy = 0.0;
  delta_c.zz = 0.0;
  return delta_c;
}

Eq10PrincipalCij compute_principal_Cij_from_eq10_baseline(const Eq10ThetaTensor& theta_tensor,
                                                          const Eq10ThetaTrace& theta_trace,
                                                          const Eq10InvariantI& invariant_I,
                                                          const Eq10DeltaCPlaceholder& delta_c,
                                                          double kappa,
                                                          double lambda) {
  const Theta3D& theta = theta_tensor.theta;
  Eq10PrincipalCij principal_cij{};
  principal_cij.c_xx =
      kappa * (theta.xx * theta.xx + theta.xy * theta.xy + theta.xz * theta.xz -
               lambda * theta_trace.value * theta.xx - 0.5 * invariant_I.value + delta_c.xx);
  principal_cij.c_xy =
      kappa * (theta.xx * theta.xy + theta.xy * theta.yy + theta.xz * theta.yz -
               lambda * theta_trace.value * theta.xy + delta_c.xy);
  principal_cij.c_yy =
      kappa * (theta.xy * theta.xy + theta.yy * theta.yy + theta.yz * theta.yz -
               lambda * theta_trace.value * theta.yy - 0.5 * invariant_I.value + delta_c.yy);
  principal_cij.c_zz =
      kappa * (theta.xz * theta.xz + theta.yz * theta.yz + theta.zz * theta.zz -
               lambda * theta_trace.value * theta.zz - 0.5 * invariant_I.value + delta_c.zz);
  return principal_cij;
}

XiDirectedReadoutResult compute_xi_directed_tensor_readout(const Eq10XiDisplacement& xi,
                                                           const Eq10PrincipalCij& principal_cij) {
  XiDirectedReadoutResult readout{};
  const double xi_norm = std::sqrt(xi.xi.x * xi.xi.x + xi.xi.y * xi.xi.y);
  if (xi_norm <= 1e-300) {
    readout.ax = 0.0;
    readout.ay = 0.0;
    return readout;
  }
  const double u_x = xi.xi.x / xi_norm;
  const double u_y = xi.xi.y / xi_norm;
  readout.ax = -(principal_cij.c_xx * u_x + principal_cij.c_xy * u_y);
  readout.ay = -(principal_cij.c_xy * u_x + principal_cij.c_yy * u_y);
  return readout;
}

DirectTpfBaselineArtifacts compute_direct_tpf_baseline_artifacts(const FieldAtPoint& field,
                                                                 double kappa,
                                                                 double lambda) {
  DirectTpfBaselineArtifacts artifacts{};
  artifacts.xi.xi = field.xi;
  artifacts.theta.theta = field.theta;
  artifacts.theta_trace = compute_theta_trace(artifacts.theta);
  artifacts.invariant_I.value = field.invariant_I;
  artifacts.delta_c = compute_deltaC_placeholder_zero();
  artifacts.principal_cij = compute_principal_Cij_from_eq10_baseline(
      artifacts.theta, artifacts.theta_trace, artifacts.invariant_I, artifacts.delta_c, kappa, lambda);
  return artifacts;
}

DirectTpfBaselineAccelerationResult compute_direct_tpf_baseline_acceleration(const FieldAtPoint& field,
                                                                             double kappa,
                                                                             double lambda) {
  const DirectTpfBaselineArtifacts artifacts = compute_direct_tpf_baseline_artifacts(field, kappa, lambda);
  const XiDirectedReadoutResult readout = compute_xi_directed_tensor_readout(artifacts.xi, artifacts.principal_cij);
  DirectTpfBaselineAccelerationResult result{};
  result.ax = readout.ax;
  result.ay = readout.ay;
  return result;
}

}  // namespace tpfcore
}  // namespace galaxy
