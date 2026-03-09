#include "newtonian.hpp"
#include <cmath>

namespace galaxy {

void NewtonianPackage::compute_accelerations(const State& state,
                                             double bh_mass,
                                             double softening,
                                             bool star_star,
                                             std::vector<double>& ax,
                                             std::vector<double>& ay) const {
  const int n = state.n();
  ax.assign(n, 0.0);
  ay.assign(n, 0.0);

  const double eps2 = softening * softening;

  // Black hole at origin
  for (int i = 0; i < n; ++i) {
    double rx = state.x[i], ry = state.y[i];
    double r_sq = rx * rx + ry * ry + eps2;
    double r_mag = std::sqrt(r_sq);
    double acc_mag = bh_mass / (r_sq * r_mag);
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
      double r_mag = std::sqrt(r_sq);
      double acc_mag = state.mass[j] / (r_sq * r_mag);
      ax[i] += acc_mag * dx;
      ay[i] += acc_mag * dy;
    }
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
    pe -= bh_mass * state.mass[i] / r;
  }
  if (!star_star) return pe;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dx = state.x[j] - state.x[i];
      double dy = state.y[j] - state.y[i];
      double r = std::sqrt(dx * dx + dy * dy + eps2);
      pe -= state.mass[i] * state.mass[j] / r;
    }
  }
  return pe;
}

}  // namespace galaxy
