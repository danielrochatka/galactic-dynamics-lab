#include "runtime_package_helpers.hpp"

#include <cmath>

namespace galaxy {
namespace tpfcore {

namespace {
unsigned g_tpf_shunt_events = 0;
}

XiWakeKinematics compute_xi_wake_kinematics(double dx,
                                            double dy,
                                            double dz,
                                            double vx_rel,
                                            double vy_rel,
                                            double vz_rel,
                                            double c_light,
                                            bool post_pass_gate) {
  XiWakeKinematics k;
  const double r_norm = std::sqrt(dx * dx + dy * dy + dz * dz);
  const double inv_r = (r_norm > 1.0e-30) ? (1.0 / r_norm) : 0.0;
  const double r_hat_x = dx * inv_r;
  const double r_hat_y = dy * inv_r;
  const double r_hat_z = dz * inv_r;
  k.v_radial = vx_rel * r_hat_x + vy_rel * r_hat_y + vz_rel * r_hat_z;
  const double vtx = vx_rel - k.v_radial * r_hat_x;
  const double vty = vy_rel - k.v_radial * r_hat_y;
  const double vtz = vz_rel - k.v_radial * r_hat_z;
  k.v_transverse = std::sqrt(vtx * vtx + vty * vty + vtz * vtz);
  k.beta_pass = k.v_transverse / c_light;
  const double vt_eps = std::max(k.v_transverse, 1.0e-30);
  if (post_pass_gate) {
    // Post-pass wake mode constants:
    // threshold = 0.10, width = 0.05 (current implementation settings, not physical constants).
    // Circular/orbiting motion with radial_ratio ~= 0 is intentionally near-null here.
    // Activation requires separating/post-pass geometry with radial_ratio above threshold.
    constexpr double kWakeGateThreshold = 0.10;
    constexpr double kWakeGateWidth = 0.05;
    const double radial_ratio = k.v_radial / vt_eps;
    k.wake_gate = 0.5 * (1.0 + std::tanh((radial_ratio - kWakeGateThreshold) / kWakeGateWidth));
  } else {
    // Continuous transverse mode preserves older behavior:
    // v_radial ~= 0 gives wake_gate ~= 0.5 (orbit-active, including near closest-pass geometry).
    k.wake_gate = 0.5 * (1.0 + std::tanh(k.v_radial / vt_eps));
  }
  k.beta_effective = k.beta_pass * k.wake_gate;
  if (k.v_transverse > 1.0e-30) {
    const double inv_vt = 1.0 / k.v_transverse;
    k.axis_x = vtx * inv_vt;
    k.axis_y = vty * inv_vt;
    k.axis_z = vtz * inv_vt;
    k.has_axis = true;
  }
  return k;
}

unsigned apply_global_accel_magnitude_shunt(const State& state,
                                            double dt,
                                            bool enable,
                                            double fraction,
                                            std::vector<double>& ax,
                                            std::vector<double>& ay) {
  g_tpf_shunt_events = 0;
  if (!enable || !(dt > 0.0) || !std::isfinite(dt) || !(fraction > 0.0)) return 0;
  const int n = state.n();
  for (int i = 0; i < n; ++i) {
    double vx = state.vx[i];
    double vy = state.vy[i];
    double v_mag = std::sqrt(vx * vx + vy * vy + 1e-30);
    double a_cap = (fraction * v_mag) / dt;
    double a_mag = std::sqrt(ax[i] * ax[i] + ay[i] * ay[i]);
    if (a_mag > a_cap && a_mag > 0.0 && std::isfinite(a_cap)) {
      double s = a_cap / a_mag;
      ax[i] *= s;
      ay[i] *= s;
      ++g_tpf_shunt_events;
    }
  }
  return g_tpf_shunt_events;
}

void reset_global_accel_shunt_events() { g_tpf_shunt_events = 0; }

unsigned global_accel_shunt_events() { return g_tpf_shunt_events; }

}  // namespace tpfcore
}  // namespace galaxy
