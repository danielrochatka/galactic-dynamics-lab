/**
 * Field-evaluation layer: packaged Xi, Theta, I, residual from source_ansatz.
 * No formulas here; delegates to source_ansatz and bundles results.
 */

#include "field_evaluation.hpp"
#include "source_ansatz.hpp"
#include "source_iteration.hpp"

namespace galaxy {
namespace tpfcore {

CanonicalFieldObjects package_canonical_field_objects(const Xi2D& xi, const Theta3D& theta) {
  CanonicalFieldObjects out;
  out.xi = xi;
  out.theta = theta;
  out.theta_trace = theta.trace();
  out.invariant_I = compute_invariant_I(theta);
  return out;
}

CanonicalFieldObjects package_canonical_field_objects(const FieldAtPoint& field) {
  CanonicalFieldObjects out;
  out.xi = field.xi;
  out.theta = field.theta;
  out.theta_trace = field.theta.trace();
  out.invariant_I = field.invariant_I;
  return out;
}

FieldAtPoint evaluate_provisional_field_single_source(double xs, double ys, double m,
                                                      double x, double y, double eps) {
  PointSourceField pf = provisional_point_source_field(xs, ys, m, x, y, eps);
  FieldAtPoint out;
  out.xi = pf.xi;
  out.theta = pf.theta;
  out.invariant_I = compute_invariant_I(pf.theta);
  out.has_residual = true;
  out.residual = provisional_point_source_residual(xs, ys, m, x, y, eps);
  return out;
}

FieldAtPoint evaluate_provisional_field_multi_source(const State& state, int i,
                                                     double bh_mass, bool star_star,
                                                     double eps) {
  const double x = state.x[i];
  const double y = state.y[i];

  FieldAtPoint out;
  out.xi.x = out.xi.y = 0.0;
  out.theta.xx = out.theta.xy = out.theta.xz = out.theta.yy = out.theta.yz = out.theta.zz = 0.0;
  out.has_residual = false;
  out.residual.x = out.residual.y = 0.0;

  for_each_gravitational_source(state, i, bh_mass, star_star, [&](const GravitationalSource& source) {
    const PointSourceField pf = provisional_point_source_field(source.x, source.y, source.mass, x, y, eps);
    out.xi.x += pf.xi.x;
    out.xi.y += pf.xi.y;
    out.theta.xx += pf.theta.xx;
    out.theta.xy += pf.theta.xy;
    out.theta.xz += pf.theta.xz;
    out.theta.yy += pf.theta.yy;
    out.theta.yz += pf.theta.yz;
    out.theta.zz += pf.theta.zz;
  });

  out.invariant_I = compute_invariant_I(out.theta);
  return out;
}

CanonicalFieldObjects evaluate_canonical_field_single_source(double xs, double ys, double m,
                                                             double x, double y, double eps) {
  return package_canonical_field_objects(
      evaluate_provisional_field_single_source(xs, ys, m, x, y, eps));
}

CanonicalFieldObjects evaluate_canonical_field_multi_source(const State& state, int i,
                                                            double bh_mass, bool star_star,
                                                            double eps) {
  return package_canonical_field_objects(
      evaluate_provisional_field_multi_source(state, i, bh_mass, star_star, eps));
}

FieldAtPoint add_provisional_fields(const FieldAtPoint& a, const FieldAtPoint& b) {
  FieldAtPoint out;
  out.xi.x = a.xi.x + b.xi.x;
  out.xi.y = a.xi.y + b.xi.y;
  out.theta.xx = a.theta.xx + b.theta.xx;
  out.theta.xy = a.theta.xy + b.theta.xy;
  out.theta.xz = a.theta.xz + b.theta.xz;
  out.theta.yy = a.theta.yy + b.theta.yy;
  out.theta.yz = a.theta.yz + b.theta.yz;
  out.theta.zz = a.theta.zz + b.theta.zz;
  out.invariant_I = compute_invariant_I(out.theta);
  out.has_residual = a.has_residual && b.has_residual;
  out.residual.x = a.has_residual && b.has_residual ? (a.residual.x + b.residual.x) : 0.0;
  out.residual.y = a.has_residual && b.has_residual ? (a.residual.y + b.residual.y) : 0.0;
  return out;
}

}  // namespace tpfcore
}  // namespace galaxy
