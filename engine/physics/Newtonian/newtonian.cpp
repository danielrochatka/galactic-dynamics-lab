#include "newtonian.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace galaxy {
namespace {
constexpr double G_SI = 6.6743e-11;

struct NewtonianRegistrar {
  NewtonianRegistrar() {
    register_physics_package_factory(
        "Newtonian",
        []() -> std::unique_ptr<PhysicsPackage> {
          return std::unique_ptr<PhysicsPackage>(new NewtonianPackage());
        });
  }
};

NewtonianRegistrar s_newtonian_registrar;
}

void NewtonianPackage::compute_accelerations(const State& state,
                                             double bh_mass,
                                             double softening,
                                             bool star_star,
                                             std::vector<double>& ax,
                                             std::vector<double>& ay) const {
  const int n = state.n();
  const double tol = 1e-12;
  softening_audit_last_call_stats_ = SofteningAuditCallStats{};
  softening_audit_last_call_stats_.eps_used = softening;
  ax.assign(n, 0.0);
  ay.assign(n, 0.0);

  const double eps2 = softening * softening;

  // Black hole at origin
  for (int i = 0; i < n; ++i) {
    double rx = state.x[i], ry = state.y[i];
    double r_sq = rx * rx + ry * ry + eps2;
    if (softening_audit_enable_) softening_audit_pair(softening_audit_last_call_stats_, std::sqrt(rx * rx + ry * ry), softening, true);
    double r_mag = std::sqrt(r_sq);
    double acc_mag = G_SI * bh_mass / (r_sq * r_mag);
    const double ax_soft = -acc_mag * rx;
    const double ay_soft = -acc_mag * ry;
    ax[i] += ax_soft;
    ay[i] += ay_soft;
    if (softening_audit_enable_ && softening > 0.0) {
      const double r2_raw = rx * rx + ry * ry;
      const double r = std::sqrt(r2_raw);
      if (r > 0.0) {
        const double unsoft_mag = G_SI * bh_mass / (r2_raw * r);
        const double ax_unsoft = -unsoft_mag * rx;
        const double ay_unsoft = -unsoft_mag * ry;
        const double soft_mag = std::hypot(ax_soft, ay_soft), un_mag = std::hypot(ax_unsoft, ay_unsoft);
        const double dot = ax_soft * ax_unsoft + ay_soft * ay_unsoft;
        const double inward_dot = ax_soft * rx + ay_soft * ry;
        const bool finite = std::isfinite(soft_mag) && std::isfinite(un_mag) && std::isfinite(dot) && std::isfinite(inward_dot);
        if (!finite) ++softening_audit_last_call_stats_.force_vector_pair_nan_inf_count;
        if (soft_mag > un_mag * (1.0 + tol)) ++softening_audit_last_call_stats_.force_vector_pair_hardening_count;
        if (dot < 0.0) ++softening_audit_last_call_stats_.force_vector_pair_direction_flip_count;
        if (inward_dot > tol) ++softening_audit_last_call_stats_.force_vector_pair_inward_bh_violation_count;
        const double ratio = (un_mag > 0.0) ? (soft_mag / un_mag) : 1.0;
        if (ratio > softening_audit_last_call_stats_.force_vector_max_soft_over_unsoft_ratio) {
          softening_audit_last_call_stats_.force_vector_max_soft_over_unsoft_ratio = ratio;
          softening_audit_last_call_stats_.force_vector_max_ratio_distance = r;
        }
      }
    }
  }

  if (!star_star) {
    if (softening_audit_enable_) {
      softening_audit_finalize(softening_audit_last_call_stats_);
      softening_audit_merge_run(softening_audit_run_stats_, softening_audit_last_call_stats_);
    }
    return;
  }

  // Pairwise star-star (O(n^2))
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      double dx = state.x[j] - state.x[i];
      double dy = state.y[j] - state.y[i];
      double r_sq = dx * dx + dy * dy + eps2;
      if (softening_audit_enable_) softening_audit_pair(softening_audit_last_call_stats_, std::sqrt(dx * dx + dy * dy), softening, false);
      double r_mag = std::sqrt(r_sq);
      double acc_mag = G_SI * state.mass[j] / (r_sq * r_mag);
      const double ax_soft = acc_mag * dx;
      const double ay_soft = acc_mag * dy;
      ax[i] += ax_soft;
      ay[i] += ay_soft;
      if (softening_audit_enable_ && softening > 0.0) {
        const double r2_raw = dx * dx + dy * dy;
        const double r = std::sqrt(r2_raw);
        if (r > 0.0) {
          const double unsoft_mag = G_SI * state.mass[j] / (r2_raw * r);
          const double ax_unsoft = unsoft_mag * dx;
          const double ay_unsoft = unsoft_mag * dy;
          const double soft_mag = std::hypot(ax_soft, ay_soft), un_mag = std::hypot(ax_unsoft, ay_unsoft);
          const double dot = ax_soft * ax_unsoft + ay_soft * ay_unsoft;
          const bool finite = std::isfinite(soft_mag) && std::isfinite(un_mag) && std::isfinite(dot);
          if (!finite) ++softening_audit_last_call_stats_.force_vector_pair_nan_inf_count;
          if (soft_mag > un_mag * (1.0 + tol)) ++softening_audit_last_call_stats_.force_vector_pair_hardening_count;
          if (dot < 0.0) ++softening_audit_last_call_stats_.force_vector_pair_direction_flip_count;
          const double ratio = (un_mag > 0.0) ? (soft_mag / un_mag) : 1.0;
          if (ratio > softening_audit_last_call_stats_.force_vector_max_soft_over_unsoft_ratio) {
            softening_audit_last_call_stats_.force_vector_max_soft_over_unsoft_ratio = ratio;
            softening_audit_last_call_stats_.force_vector_max_ratio_distance = r;
          }
        }
      }
    }
  }
  if (softening_audit_enable_ && softening > 0.0) {
    std::vector<double> ax_u(n, 0.0), ay_u(n, 0.0);
    for (int i = 0; i < n; ++i) {
      const double rx = state.x[i], ry = state.y[i];
      const double r2 = rx * rx + ry * ry;
      if (r2 > 0.0) {
        const double invr3 = 1.0 / (r2 * std::sqrt(r2));
        ax_u[i] -= G_SI * bh_mass * rx * invr3;
        ay_u[i] -= G_SI * bh_mass * ry * invr3;
      }
    }
    if (star_star) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          if (i == j) continue;
          const double dx = state.x[j] - state.x[i];
          const double dy = state.y[j] - state.y[i];
          const double r2 = dx * dx + dy * dy;
          if (r2 <= 0.0) continue;
          const double invr3 = 1.0 / (r2 * std::sqrt(r2));
          ax_u[i] += G_SI * state.mass[j] * dx * invr3;
          ay_u[i] += G_SI * state.mass[j] * dy * invr3;
        }
      }
    }
    softening_audit_last_call_stats_.net_particle_count = static_cast<std::uint64_t>(n);
    std::vector<double> ratios;
    ratios.reserve(n);
    for (int i = 0; i < n; ++i) {
      const double s = std::hypot(ax[i], ay[i]);
      const double u = std::hypot(ax_u[i], ay_u[i]);
      const double ratio = (u > 0.0) ? (s / u) : ((s == 0.0) ? 1.0 : std::numeric_limits<double>::infinity());
      ratios.push_back(ratio);
      if (ratio > (1.0 + tol)) ++softening_audit_last_call_stats_.net_hardening_count;
      if (ax[i] * ax_u[i] + ay[i] * ay_u[i] < 0.0) ++softening_audit_last_call_stats_.net_direction_flip_count;
    }
    std::sort(ratios.begin(), ratios.end());
    softening_audit_last_call_stats_.net_ratio_median = ratios[ratios.size() / 2];
    softening_audit_last_call_stats_.net_ratio_p95 = ratios[static_cast<std::size_t>(0.95 * (ratios.size() - 1))];
    softening_audit_last_call_stats_.net_ratio_max = ratios.back();
  }
  if (softening_audit_enable_) {
    softening_audit_finalize(softening_audit_last_call_stats_);
    softening_audit_merge_run(softening_audit_run_stats_, softening_audit_last_call_stats_);
  }
}

double NewtonianPackage::compute_potential_energy(const State& state,
                                                  double bh_mass,
                                                  double softening,
                                                  bool star_star) const {
  const int n = state.n();
  double pe = 0.0;
  const double eps2 = softening * softening;

  for (int i = 0; i < n; ++i) {
    double r = std::sqrt(state.x[i] * state.x[i] + state.y[i] * state.y[i] + eps2);
    pe -= G_SI * bh_mass * state.mass[i] / r;
  }
  if (!star_star) return pe;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dx = state.x[j] - state.x[i];
      double dy = state.y[j] - state.y[i];
      double r = std::sqrt(dx * dx + dy * dy + eps2);
      pe -= G_SI * state.mass[i] * state.mass[j] / r;
    }
  }
  return pe;
}

}  // namespace galaxy
