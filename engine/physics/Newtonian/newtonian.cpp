#include "newtonian.hpp"
#include <cmath>

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
  softening_audit_stats_ = SofteningAuditStats{};
  softening_audit_stats_.eps_used = softening;
  ax.assign(n, 0.0);
  ay.assign(n, 0.0);

  const double eps2 = softening * softening;

  // Black hole at origin
  for (int i = 0; i < n; ++i) {
    double rx = state.x[i], ry = state.y[i];
    double r_sq = rx * rx + ry * ry + eps2;
    if (softening_audit_enable_) softening_audit_pair(softening_audit_stats_, std::sqrt(rx * rx + ry * ry), softening, true);
    double r_mag = std::sqrt(r_sq);
    double acc_mag = G_SI * bh_mass / (r_sq * r_mag);
    ax[i] -= acc_mag * rx;
    ay[i] -= acc_mag * ry;
  }

  if (!star_star) return;

  // Pairwise star-star (O(n^2))
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      double dx = state.x[j] - state.x[i];
      double dy = state.y[j] - state.y[i];
      double r_sq = dx * dx + dy * dy + eps2;
      if (softening_audit_enable_) softening_audit_pair(softening_audit_stats_, std::sqrt(dx * dx + dy * dy), softening, false);
      double r_mag = std::sqrt(r_sq);
      double acc_mag = G_SI * state.mass[j] / (r_sq * r_mag);
      ax[i] += acc_mag * dx;
      ay[i] += acc_mag * dy;
    }
  }
  if (softening_audit_enable_) softening_audit_finalize(softening_audit_stats_);
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
