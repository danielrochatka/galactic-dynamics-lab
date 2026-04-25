#include "physics/TPFCore/tpf_4d_static_field.hpp"

#include <cmath>

namespace galaxy {
namespace tpfcore {

namespace {

Field4DAtPoint make_zero_field_4d() {
  Field4DAtPoint out{};
  out.xi = Xi4{0.0, 0.0, 0.0, 0.0};
  out.theta = Theta4{0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0,
                     0.0};
  out.theta_trace_4d = 0.0;
  out.invariant_I_4d = 0.0;
  return out;
}

}  // namespace

Field4DAtPoint evaluate_static_source_field_4d(const StaticSourcePoint4D& source,
                                               const StaticProbePoint4D& probe,
                                               double eps) {
  const double dx = probe.x - source.x;
  const double dy = probe.y - source.y;
  const double dz = probe.z - source.z;

  double r2 = dx * dx + dy * dy + dz * dz + eps * eps;
  double r = std::sqrt(r2);
  if (r < 1e-30) {
    r = 1e-30;
    r2 = r * r;
  }

  const double r3 = r2 * r;
  const double r5 = r3 * r2;

  Field4DAtPoint out{};

  out.xi.t = 0.0;
  out.xi.x = source.mass * dx / r3;
  out.xi.y = source.mass * dy / r3;
  out.xi.z = source.mass * dz / r3;

  out.theta.tt = 0.0;
  out.theta.tx = 0.0;
  out.theta.ty = 0.0;
  out.theta.tz = 0.0;
  out.theta.xt = 0.0;
  out.theta.yt = 0.0;
  out.theta.zt = 0.0;

  const double inv_r3 = 1.0 / r3;
  const double common = source.mass;

  out.theta.xx = common * (inv_r3 - 3.0 * dx * dx / r5);
  out.theta.xy = common * (-3.0 * dx * dy / r5);
  out.theta.xz = common * (-3.0 * dx * dz / r5);

  out.theta.yx = common * (-3.0 * dy * dx / r5);
  out.theta.yy = common * (inv_r3 - 3.0 * dy * dy / r5);
  out.theta.yz = common * (-3.0 * dy * dz / r5);

  out.theta.zx = common * (-3.0 * dz * dx / r5);
  out.theta.zy = common * (-3.0 * dz * dy / r5);
  out.theta.zz = common * (inv_r3 - 3.0 * dz * dz / r5);

  out.theta_trace_4d = trace_contraction_4d(out.theta);
  out.invariant_I_4d = invariant_I_4d(out.theta);
  return out;
}

Field4DAtPoint evaluate_static_sources_field_4d(const std::vector<StaticSourcePoint4D>& sources,
                                                const StaticProbePoint4D& probe,
                                                double eps) {
  Field4DAtPoint total = make_zero_field_4d();

  for (std::vector<StaticSourcePoint4D>::const_iterator it = sources.begin(); it != sources.end(); ++it) {
    const Field4DAtPoint one = evaluate_static_source_field_4d(*it, probe, eps);
    total.xi.t += one.xi.t;
    total.xi.x += one.xi.x;
    total.xi.y += one.xi.y;
    total.xi.z += one.xi.z;

    total.theta.tt += one.theta.tt;
    total.theta.tx += one.theta.tx;
    total.theta.ty += one.theta.ty;
    total.theta.tz += one.theta.tz;

    total.theta.xt += one.theta.xt;
    total.theta.xx += one.theta.xx;
    total.theta.xy += one.theta.xy;
    total.theta.xz += one.theta.xz;

    total.theta.yt += one.theta.yt;
    total.theta.yx += one.theta.yx;
    total.theta.yy += one.theta.yy;
    total.theta.yz += one.theta.yz;

    total.theta.zt += one.theta.zt;
    total.theta.zx += one.theta.zx;
    total.theta.zy += one.theta.zy;
    total.theta.zz += one.theta.zz;
  }

  total.theta_trace_4d = trace_contraction_4d(total.theta);
  total.invariant_I_4d = invariant_I_4d(total.theta);
  return total;
}

}  // namespace tpfcore
}  // namespace galaxy
