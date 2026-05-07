#include "readout_model_families.hpp"
#include "provisional_readout.hpp"

#include "field_evaluation.hpp"

#include <cmath>

namespace galaxy {
namespace tpfcore {
namespace readout_models {

void apply_experimental_radial_scaling_readout(const State& state,
                                               int i,
                                               double bh_mass,
                                               bool star_star,
                                               double eps,
                                               double readout_scale,
                                               double& ax,
                                               double& ay,
                                               ReadoutDiagnostics* diag) {
  const double x = state.x[i];
  const double y = state.y[i];
  const double r2 = x * x + y * y + eps * eps;
  const double r = std::sqrt(r2);
  if (r < 1e-30) {
    ax = ay = 0.0;
    if (diag) {
      diag->ax = ax;
      diag->ay = ay;
      diag->theta_rr = diag->theta_tt = diag->theta_tr = diag->theta_rr_plus_theta_tt = 0.0;
      diag->provisional_radial_readout = diag->provisional_tangential_readout = 0.0;
    }
    return;
  }
  const double rx = x / r;
  const double ry = y / r;
  const double tx = -ry;
  const double ty = rx;

  const FieldAtPoint field = evaluate_provisional_field_multi_source(state, i, bh_mass, star_star, eps);
  const CanonicalFieldObjects canonical = as_canonical_field_objects(field);
  const Theta3D& theta_sum = canonical.theta;

  const double theta_rr = rx * rx * theta_sum.xx + 2.0 * rx * ry * theta_sum.xy + ry * ry * theta_sum.yy;
  const double theta_tr = tx * (theta_sum.xx * rx + theta_sum.xy * ry) + ty * (theta_sum.xy * rx + theta_sum.yy * ry);

  // Inward radial: magnitude = readout_scale * (-theta_rr) * r (so effective ~ 1/r^2 when theta_rr ~ 1/r^3).
  // Apply as acceleration = -magnitude * r_hat (inward). r_hat = (rx, ry) points outward; Newtonian a is inward.
  const double provisional_radial = readout_scale * (-theta_rr) * r;
  const double provisional_tangential = 0.0;

  ax = -provisional_radial * rx;
  ay = -provisional_radial * ry;

  if (!diag) return;

  diag->ax = ax;
  diag->ay = ay;
  diag->theta_rr = theta_rr;
  diag->theta_tt = 0.0;
  diag->theta_tr = theta_tr;
  diag->theta_rr_plus_theta_tt = theta_rr;
  diag->provisional_radial_readout = provisional_radial;
  diag->provisional_tangential_readout = provisional_tangential;
  diag->theta_xx = theta_sum.xx;
  diag->theta_xy = theta_sum.xy;
  diag->theta_yy = theta_sum.yy;
  diag->theta_trace = canonical.theta_trace;
  diag->invariant_I = canonical.invariant_I;
  diag->theta_norm = theta_frobenius_norm(theta_sum);
}

}  // namespace readout_models
}  // namespace tpfcore
}  // namespace galaxy
