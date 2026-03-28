/**
 * Derived TPF radial: bounce closure (Sec. IX) and Hessian diagnostics.
 */

#include "derived_tpf_radial.hpp"
#include "../../config.hpp"
#include <algorithm>
#include <cmath>

namespace galaxy {
namespace tpfcore {

namespace {

Theta3D add_theta(const Theta3D& a, const Theta3D& b) {
  Theta3D o;
  o.xx = a.xx + b.xx;
  o.xy = a.xy + b.xy;
  o.xz = a.xz + b.xz;
  o.yy = a.yy + b.yy;
  o.yz = a.yz + b.yz;
  o.zz = a.zz + b.zz;
  return o;
}

}  // namespace

Theta3D evaluate_derived_theta(double mass_kg, double dx, double dy, double dz, double eps) {
  Theta3D t;
  double d_sq = dx * dx + dy * dy + dz * dz + eps * eps;
  double d = std::sqrt(d_sq);
  if (d < 1e-30) d = 1e-30;
  double d5 = d * d * d * d * d;
  double gm = TPF_G_SI * mass_kg;
  t.xx = gm * (3.0 * dx * dx - d_sq) / d5;
  t.yy = gm * (3.0 * dy * dy - d_sq) / d5;
  t.zz = gm * (3.0 * dz * dz - d_sq) / d5;
  t.xy = gm * (3.0 * dx * dy) / d5;
  t.xz = gm * (3.0 * dx * dz) / d5;
  t.yz = gm * (3.0 * dy * dz) / d5;
  return t;
}

double derived_invariant_I_contracted(const Theta3D& theta) {
  return theta.xx * theta.xx + theta.yy * theta.yy + theta.zz * theta.zz +
         2.0 * (theta.xy * theta.xy + theta.xz * theta.xz + theta.yz * theta.yz);
}

Theta3D sum_derived_theta_at_point(const State& state, double bh_mass, double px, double py, double pz,
                                   double eps) {
  Theta3D sum;
  sum.xx = sum.xy = sum.xz = sum.yy = sum.yz = sum.zz = 0.0;
  if (bh_mass > 0.0) {
    Theta3D bh = evaluate_derived_theta(bh_mass, px - 0.0, py - 0.0, pz - 0.0, eps);
    sum = add_theta(sum, bh);
  }
  const int n = state.n();
  for (int j = 0; j < n; ++j) {
    double mj = state.mass[j];
    if (mj <= 0.0) continue;
    Theta3D tj =
        evaluate_derived_theta(mj, px - state.x[j], py - state.y[j], pz - 0.0, eps);
    sum = add_theta(sum, tj);
  }
  return sum;
}

double enclosed_stellar_mass_cyl(const State& state, double r_cyl) {
  double M = 0.0;
  const int n = state.n();
  for (int j = 0; j < n; ++j) {
    double rj = std::hypot(state.x[j], state.y[j]);
    if (rj <= r_cyl) M += state.mass[j];
  }
  return M;
}

double get_tpf_mass_at_r(double M_total, double r) {
  if (M_total <= 0.0) return 0.0;
  double rr = std::abs(r);
  if (rr <= 0.0) return 0.0;

  double r_s = (2.0 * TPF_G_SI * M_total) / (c * c);
  double r3 = rr * rr * rr;
  double rs3 = r_s * r_s * r_s;

  return M_total * (r3 / (r3 + rs3));
}

TpfRadialGravityProfile build_tpf_gravity_profile(const State& state, double bh_mass, double max_radius,
                                                  int bins, double eps) {
  (void)state;
  (void)bh_mass;
  (void)eps;
  TpfRadialGravityProfile p;
  p.bins = std::max(1, bins);
  p.max_radius = std::max(max_radius, 1e-30);
  p.delta_r = p.max_radius / static_cast<double>(p.bins);
  p.r_outer.resize(static_cast<size_t>(p.bins));
  p.M_eff_enc.resize(static_cast<size_t>(p.bins));

  for (int b = 0; b < p.bins; ++b) {
    p.r_outer[static_cast<size_t>(b)] = (static_cast<double>(b) + 1.0) * p.delta_r;
    p.M_eff_enc[static_cast<size_t>(b)] = 0.0;
  }
  return p;
}

double TpfRadialGravityProfile::M_eff_at_cylindrical_r(double r_cyl) const {
  if (r_cyl <= 0.0 || bins <= 0 || r_outer.empty()) return 0.0;
  if (r_cyl >= r_outer.back()) return M_eff_enc.back();
  for (size_t k = 0; k < r_outer.size(); ++k) {
    if (r_cyl <= r_outer[k]) {
      double r_lo = (k == 0) ? 0.0 : r_outer[k - 1];
      double M_lo = (k == 0) ? 0.0 : M_eff_enc[k - 1];
      double r_hi = r_outer[k];
      double M_hi = M_eff_enc[k];
      double t = (r_cyl - r_lo) / std::max(r_hi - r_lo, 1e-30);
      return M_lo + t * (M_hi - M_lo);
    }
  }
  return M_eff_enc.back();
}

double radial_acceleration_scalar_derived(const State& state, double bh_mass, double r_cyl, double eps) {
  double r_soft_sq = r_cyl * r_cyl + eps * eps;
  if (r_soft_sq < 1e-60) r_soft_sq = 1e-60;
  double M_enc = bh_mass + enclosed_stellar_mass_cyl(state, r_cyl);
  double m_bounce = get_tpf_mass_at_r(M_enc, r_cyl);
  return -TPF_G_SI * m_bounce / r_soft_sq;
}

}  // namespace tpfcore
}  // namespace galaxy
