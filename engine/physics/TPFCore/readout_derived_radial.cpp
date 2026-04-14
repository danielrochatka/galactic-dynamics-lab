#include "readout_model_families.hpp"

#include "field_evaluation.hpp"

#include <cmath>

namespace galaxy {
namespace tpfcore {
namespace readout_models {

void apply_derived_radial_readout(const State& state,
                                  int i,
                                  double bh_mass,
                                  double eps,
                                  const DerivedTpfPoissonConfig& dcfg,
                                  const TpfRadialGravityProfile* profile_in,
                                  double readout_scale,
                                  double theta_tt_scale,
                                  double theta_tr_scale,
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
      diag->regime.clear();
    }
    return;
  }

  TpfRadialGravityProfile profile_storage;
  const TpfRadialGravityProfile* profile = profile_in;
  if (!profile) {
    profile_storage = build_tpf_gravity_profile(state, bh_mass, dcfg, eps);
    profile = &profile_storage;
  }

  const double r_cyl = std::hypot(x, y);
  const double a_s = radial_acceleration_scalar_derived(state, bh_mass, *profile, r_cyl, eps);
  ax = a_s * (x / r);
  ay = a_s * (y / r);

  const Theta3D theta_sum = sum_derived_theta_at_point(state, bh_mass, x, y, 0.0, eps);
  const double rx = x / r;
  const double ry = y / r;
  const double tx = -ry;
  const double ty = rx;
  const double theta_rr = rx * rx * theta_sum.xx + 2.0 * rx * ry * theta_sum.xy + ry * ry * theta_sum.yy;
  const double theta_tt = theta_tt_scale * (-theta_rr);
  const double theta_rr_plus_theta_tt = theta_rr + theta_tt;
  const double theta_tr =
      tx * (theta_sum.xx * rx + theta_sum.xy * ry) + ty * (theta_sum.xy * rx + theta_sum.yy * ry);
  const double provisional_tangential = readout_scale * theta_tr_scale * theta_tr;

  if (!diag) return;

  diag->ax = ax;
  diag->ay = ay;
  diag->theta_rr = theta_rr;
  diag->theta_tt = theta_tt;
  diag->theta_tr = theta_tr;
  diag->theta_rr_plus_theta_tt = theta_rr_plus_theta_tt;
  diag->provisional_radial_readout = a_s;
  diag->provisional_tangential_readout = provisional_tangential;
  diag->theta_xx = theta_sum.xx;
  diag->theta_xy = theta_sum.xy;
  diag->theta_yy = theta_sum.yy;
  const CanonicalFieldObjects canonical = build_canonical_field_objects(Xi2D{}, theta_sum);
  diag->theta_trace = canonical.theta_trace;
  diag->invariant_I = canonical.invariant_I;
  diag->theta_norm = theta_frobenius_norm(theta_sum);
  diag->regime = "derived-tpf-radial";
}

}  // namespace readout_models
}  // namespace tpfcore
}  // namespace galaxy
