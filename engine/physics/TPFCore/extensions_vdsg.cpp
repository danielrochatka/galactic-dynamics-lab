#include "extensions_vdsg.hpp"

#include "../../softening_policy.hpp"
#include "derived_tpf_radial.hpp"
#include "source_iteration.hpp"
#include <algorithm>
#include <cmath>

namespace galaxy {
namespace tpfcore {

namespace {
const double C_SI_LIGHT = 299792458.0;

double clampd(double x, double lo, double hi) { return std::max(lo, std::min(hi, x)); }

void update_minmax(double v, double& mn, double& mx, bool first) {
  if (first) { mn = v; mx = v; return; }
  mn = std::min(mn, v);
  mx = std::max(mx, v);
}
}

bool is_valid_vdsg_mode(const std::string& mode) {
  return mode == "legacy_speed" || mode == "radial_doppler_rational" || mode == "radial_doppler_exp" ||
         mode == "radial_doppler_bounded";
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

void accumulate_vdsg_velocity_modifier(const State& state, double bh_mass, double softening, bool star_star,
                                       double vdsg_coupling, double vdsg_mass_baseline_kg, const std::string& vdsg_mode,
                                       double mass_gate_m0_kg, double mass_gate_alpha, double x_clamp,
                                       bool weak_field_gate_enable, double weak_field_a0, double weak_field_power,
                                       double bounded_amplitude, std::vector<double>& ax, std::vector<double>& ay,
                                       VdsgDiagnosticsSummary* diagnostics_out) {
  const int n = state.n();
  const double G = tpfcore::TPF_G_SI;
  ax.assign(n, 0.0);
  ay.assign(n, 0.0);
  if (!(vdsg_coupling != 0.0) || !std::isfinite(vdsg_coupling)) return;

  const bool radial_mode = (vdsg_mode != "legacy_speed");
  if (diagnostics_out) *diagnostics_out = VdsgDiagnosticsSummary{};
  if (diagnostics_out) diagnostics_out->mode_is_radial = radial_mode;

  for (int i = 0; i < n; ++i) {
    for_each_gravitational_source(state, i, bh_mass, star_star, [&](const GravitationalSource& source) {
      const double dx = state.x[i] - source.x;
      const double dy = state.y[i] - source.y;
      const double r_sq = dx * dx + dy * dy;
      if (!(r_sq > 0.0) || !std::isfinite(r_sq)) return;
      const double r_unsoft = std::sqrt(r_sq);
      if (!(r_unsoft > 0.0) || !std::isfinite(r_unsoft)) return;
      const double ux = dx / r_unsoft;
      const double uy = dy / r_unsoft;

      double dvx = state.vx[i];
      double dvy = state.vy[i];
      double source_vx = 0.0, source_vy = 0.0;
      if (!source.is_black_hole && source.star_index >= 0) {
        source_vx = state.vx[source.star_index];
        source_vy = state.vy[source.star_index];
      }
      dvx -= source_vx;
      dvy -= source_vy;

      const double denom = r_sq + softening * softening;
      if (!(denom > 0.0) || !std::isfinite(denom)) return;
      const double softened_force_mag = G * source.mass * r_unsoft / std::pow(denom, 1.5);
      if (!std::isfinite(softened_force_mag)) return;
      double factor = 1.0;

      if (!radial_mode) {
        const double vrel_mag = std::sqrt(dvx * dvx + dvy * dvy + 1e-300);
        const double beta = vrel_mag / C_SI_LIGHT;
        const double lambda_eff = vdsg_effective_coupling(vdsg_coupling, source.mass, vdsg_mass_baseline_kg);
        factor = 1.0 + lambda_eff * beta;
      } else {
        const double v_rad = dvx * ux + dvy * uy;
        const double beta_rad = v_rad / C_SI_LIGHT;
        if (!std::isfinite(beta_rad)) return;
        if (diagnostics_out) {
          const bool first = diagnostics_out->pairs_evaluated == 0;
          diagnostics_out->pairs_evaluated++;
          update_minmax(beta_rad, diagnostics_out->min_beta_rad, diagnostics_out->max_beta_rad, first);
          diagnostics_out->sum_abs_beta_rad += std::abs(beta_rad);
          if (std::abs(beta_rad) < 1e-15) diagnostics_out->beta_near_zero++;
          else if (beta_rad > 0.0) diagnostics_out->beta_positive++;
          else diagnostics_out->beta_negative++;
        }
        if (beta_rad == 0.0) return;
        const double mi = (i < static_cast<int>(state.mass.size())) ? state.mass[i] : 0.0;
        const double mj = source.mass;
        if (!(mi > 0.0) || !(mj > 0.0) || !std::isfinite(mi) || !std::isfinite(mj)) return;
        const double mu = (mi * mj) / (mi + mj);
        if (!(mu > 0.0) || !std::isfinite(mu)) return;
        double M_gate = 0.0;
        if (mass_gate_m0_kg > 0.0 && std::isfinite(mass_gate_alpha)) {
          const double q = std::pow(mu / mass_gate_m0_kg, mass_gate_alpha);
          if (std::isfinite(q) && q > 0.0) M_gate = q / (1.0 + q);
        }
        if (!(M_gate > 0.0)) return;
        double S_gate = 1.0;
        if (weak_field_gate_enable) {
          if (weak_field_a0 > 0.0 && weak_field_power > 0.0 && std::isfinite(weak_field_a0) && std::isfinite(weak_field_power)) {
            const double a_pair = G * (mi + mj) / denom;
            const double z = std::pow(a_pair / weak_field_a0, weak_field_power);
            if (std::isfinite(z) && z >= 0.0) S_gate = 1.0 / (1.0 + z);
          } else {
            S_gate = 0.0;
          }
        }
        const double x_raw = vdsg_coupling * M_gate * S_gate * beta_rad;
        double x_cap = (x_clamp > 0.0 && std::isfinite(x_clamp)) ? x_clamp : 0.25;
        if (vdsg_mode == "radial_doppler_rational" && x_cap >= 1.0) x_cap = 0.95;
        const double x = clampd(x_raw, -x_cap, x_cap);
        if (!std::isfinite(x_raw) || !std::isfinite(x)) return;
        if (diagnostics_out) {
          const bool first = diagnostics_out->pairs_evaluated == 1;
          update_minmax(x_raw, diagnostics_out->min_x_raw, diagnostics_out->max_x_raw, first);
          update_minmax(x, diagnostics_out->min_x_clamped, diagnostics_out->max_x_clamped, first);
          if (x != x_raw) diagnostics_out->x_clamped++;
        }
        if (vdsg_mode == "radial_doppler_rational") {
          factor = 1.0 / (1.0 - x);
        } else if (vdsg_mode == "radial_doppler_exp") {
          factor = std::exp(x);
        } else {
          const double A = (bounded_amplitude > 0.0 && std::isfinite(bounded_amplitude)) ? bounded_amplitude : 0.25;
          factor = 1.0 + A * std::tanh(x / A);
        }
      }
      if (!std::isfinite(factor)) return;
      const double delta_a_mag = softened_force_mag * (factor - 1.0);
      if (diagnostics_out && radial_mode) {
        const bool first = diagnostics_out->pairs_evaluated == 1;
        update_minmax(factor, diagnostics_out->min_factor, diagnostics_out->max_factor, first);
        diagnostics_out->sum_factor += factor;
        const double abs_da = std::abs(delta_a_mag);
        diagnostics_out->sum_abs_delta_accel += abs_da;
        diagnostics_out->max_abs_delta_accel = std::max(diagnostics_out->max_abs_delta_accel, abs_da);
      }
      double dax = -ux * delta_a_mag;
      double day = -uy * delta_a_mag;
      if (!std::isfinite(dax) || !std::isfinite(day)) return;
      ax[i] += dax;
      ay[i] += day;
    });
  }
}

}  // namespace tpfcore
}  // namespace galaxy
