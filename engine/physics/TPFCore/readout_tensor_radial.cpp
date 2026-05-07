#include "readout_model_families.hpp"
#include "provisional_readout.hpp"

#include "field_evaluation.hpp"
#include "source_iteration.hpp"

#include <cmath>

namespace galaxy {
namespace tpfcore {
namespace readout_models {

void apply_tensor_radial_projection_readout(const State& state,
                                            int i,
                                            double bh_mass,
                                            bool star_star,
                                            double eps,
                                            double readout_scale,
                                            bool negated,
                                            double& ax,
                                            double& ay,
                                            ReadoutDiagnostics* diag) {
  ax = 0.0;
  ay = 0.0;

  Theta3D theta_sum{};
  bool has_theta_sum = false;

  const double x = state.x[i];
  const double y = state.y[i];

  auto add_contribution = [&](const GravitationalSource& source) {
    const double dx = x - source.x;
    const double dy = y - source.y;
    const double r2 = dx * dx + dy * dy + eps * eps;
    const double r = std::sqrt(r2);
    if (r < 1e-30) return;
    const double rx = dx / r;
    const double ry = dy / r;

    const FieldAtPoint field = evaluate_provisional_field_single_source(source.x, source.y, source.mass, x, y, eps);
    const Theta3D& theta = field.theta;
    const double ax_contrib = theta.xx * rx + theta.xy * ry;
    const double ay_contrib = theta.xy * rx + theta.yy * ry;

    ax += readout_scale * ax_contrib;
    ay += readout_scale * ay_contrib;

    theta_sum.xx += theta.xx;
    theta_sum.xy += theta.xy;
    theta_sum.xz += theta.xz;
    theta_sum.yy += theta.yy;
    theta_sum.yz += theta.yz;
    theta_sum.zz += theta.zz;
    has_theta_sum = true;
  };

  for_each_gravitational_source(state, i, bh_mass, star_star, add_contribution);

  if (negated) {
    ax = -ax;
    ay = -ay;
  }

  if (!diag) return;

  diag->ax = ax;
  diag->ay = ay;
  if (has_theta_sum) {
    diag->theta_xx = theta_sum.xx;
    diag->theta_xy = theta_sum.xy;
    diag->theta_yy = theta_sum.yy;
    const CanonicalFieldObjects canonical = build_canonical_field_objects(Xi2D{}, theta_sum);
    diag->theta_trace = canonical.theta_trace;
    diag->invariant_I = canonical.invariant_I;
    diag->theta_norm = theta_frobenius_norm(theta_sum);
  } else {
    diag->theta_xx = diag->theta_xy = diag->theta_yy = 0.0;
    diag->theta_trace = diag->invariant_I = diag->theta_norm = 0.0;
  }
}

}  // namespace readout_models
}  // namespace tpfcore
}  // namespace galaxy
