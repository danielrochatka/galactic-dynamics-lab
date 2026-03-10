/**
 * Field-evaluation layer: packaged Xi, Theta, I, residual from source_ansatz.
 * No formulas here; delegates to source_ansatz and bundles results.
 */

#include "field_evaluation.hpp"
#include "source_ansatz.hpp"

namespace galaxy {
namespace tpfcore {

FieldAtPoint evaluate_provisional_field_single_source(double xs, double ys, double m,
                                                      double x, double y, double eps, double c) {
  PointSourceField pf = provisional_point_source_field(xs, ys, m, x, y, eps, c);
  FieldAtPoint out;
  out.xi = pf.xi;
  out.theta = pf.theta;
  out.invariant_I = compute_invariant_I(pf.theta);
  out.has_residual = true;
  out.residual = provisional_point_source_residual(xs, ys, m, x, y, eps, c);
  return out;
}

FieldAtPoint evaluate_provisional_field_multi_source(const State& state, int i,
                                                     double bh_mass, bool star_star,
                                                     double eps, double c) {
  const double x = state.x[i];
  const double y = state.y[i];
  const int n = state.n();

  FieldAtPoint out;
  out.xi.x = out.xi.y = 0.0;
  out.theta.xx = out.theta.xy = out.theta.yy = 0.0;
  out.has_residual = false;
  out.residual.x = out.residual.y = 0.0;

  auto add = [&](double xs, double ys, double m) {
    if (m <= 0.0) return;
    PointSourceField pf = provisional_point_source_field(xs, ys, m, x, y, eps, c);
    out.xi.x += pf.xi.x;
    out.xi.y += pf.xi.y;
    out.theta.xx += pf.theta.xx;
    out.theta.xy += pf.theta.xy;
    out.theta.yy += pf.theta.yy;
  };

  if (bh_mass > 0.0) add(0.0, 0.0, bh_mass);
  if (star_star) {
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      add(state.x[j], state.y[j], state.mass[j]);
    }
  }

  out.invariant_I = compute_invariant_I(out.theta);
  return out;
}

FieldAtPoint add_provisional_fields(const FieldAtPoint& a, const FieldAtPoint& b) {
  FieldAtPoint out;
  out.xi.x = a.xi.x + b.xi.x;
  out.xi.y = a.xi.y + b.xi.y;
  out.theta.xx = a.theta.xx + b.theta.xx;
  out.theta.xy = a.theta.xy + b.theta.xy;
  out.theta.yy = a.theta.yy + b.theta.yy;
  out.invariant_I = compute_invariant_I(out.theta);
  out.has_residual = a.has_residual && b.has_residual;
  out.residual.x = a.has_residual && b.has_residual ? (a.residual.x + b.residual.x) : 0.0;
  out.residual.y = a.has_residual && b.has_residual ? (a.residual.y + b.residual.y) : 0.0;
  return out;
}

}  // namespace tpfcore
}  // namespace galaxy
