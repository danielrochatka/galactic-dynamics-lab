#include "runtime_package_helpers.hpp"

#include <cmath>

namespace galaxy {
namespace tpfcore {

namespace {
unsigned g_tpf_shunt_events = 0;
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
