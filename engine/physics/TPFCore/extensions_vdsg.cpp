#include "extensions_vdsg.hpp"

#include "derived_tpf_radial.hpp"
#include "source_iteration.hpp"
#include <algorithm>
#include <cmath>

namespace galaxy {
namespace tpfcore {

namespace {
const double C_SI_LIGHT = 299792458.0;
}

double vdsg_effective_coupling(double lambda0, double source_mass_kg, double baseline_mass_kg) {
  if (!(lambda0 != 0.0) || !std::isfinite(lambda0)) return lambda0;
  if (source_mass_kg <= 0.0 || baseline_mass_kg <= 0.0) return lambda0;
  constexpr double kFloorKg = 1e-99;
  const double m_src = std::max(source_mass_kg, kFloorKg);
  const double m_base = std::max(baseline_mass_kg, kFloorKg);
  const double log_s = std::log10(m_src);
  const double log_b = std::log10(m_base);
  if (!std::isfinite(log_s) || !std::isfinite(log_b)) return lambda0;
  if (std::abs(log_s) < 1e-300) return lambda0;
  const double g = lambda0 * (log_b / log_s);
  return std::isfinite(g) ? g : lambda0;
}

void accumulate_vdsg_velocity_modifier(const State& state,
                                       double bh_mass,
                                       double softening,
                                       bool star_star,
                                       double vdsg_coupling,
                                       double vdsg_mass_baseline_kg,
                                       std::vector<double>& ax,
                                       std::vector<double>& ay) {
  const int n = state.n();
  const double G = tpfcore::TPF_G_SI;
  const double eps_sq = softening * softening;

  ax.assign(n, 0.0);
  ay.assign(n, 0.0);

  if (!(vdsg_coupling != 0.0) || !std::isfinite(vdsg_coupling)) return;

  for (int i = 0; i < n; ++i) {
    for_each_gravitational_source(state, i, bh_mass, star_star, [&](const GravitationalSource& source) {
      const double dx = state.x[i] - source.x;
      const double dy = state.y[i] - source.y;
      const double r_sq = dx * dx + dy * dy;
      const double denom = r_sq + eps_sq;
      const double r_mag = std::sqrt(denom);
      if (r_mag < 1e-300) return;

      const double ux = dx / r_mag;
      const double uy = dy / r_mag;
      double dvx = state.vx[i];
      double dvy = state.vy[i];
      if (!source.is_black_hole && source.star_index >= 0) {
        dvx -= state.vx[source.star_index];
        dvy -= state.vy[source.star_index];
      }
      const double vrel_mag = std::sqrt(dvx * dvx + dvy * dvy + 1e-300);
      const double beta = vrel_mag / C_SI_LIGHT;
      const double lambda_eff = vdsg_effective_coupling(vdsg_coupling, source.mass, vdsg_mass_baseline_kg);
      const double doppler_scale = 1.0 + lambda_eff * beta;
      const double a_newt = G * source.mass / denom;
      const double accel_mag = a_newt * (doppler_scale - 1.0);

      ax[i] -= ux * accel_mag;
      ay[i] -= uy * accel_mag;
    });
  }
}

}  // namespace tpfcore
}  // namespace galaxy
