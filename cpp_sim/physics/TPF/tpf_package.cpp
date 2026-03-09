/**
 * TPF weak-field correspondence package.
 *
 * Implements the linearized weak-field sector only:
 *   curl-free displacement xi = grad(phi)
 *   Theta_ij = partial_i partial_j phi, Theta = nabla^2 phi
 *   Linearized trace equation: nabla^2 phi = alpha * rho(x)
 *
 * We solve via superposition of Green-function contributions (no grid).
 * Acceleration is an OPERATIONAL READOUT from -grad(phi), not the primitive ontology.
 */

#include "tpf_package.hpp"
#include "../../config.hpp"
#include <cmath>

namespace galaxy {

namespace {

const double PI = 3.14159265358979323846;
const double FOUR_PI = 4.0 * PI;
const double INV_FOUR_PI = 1.0 / FOUR_PI;

}  // namespace

void TPFPackage::init_from_config(const Config& config) {
  if (config.tpf_match_newtonian_scale) {
    alpha_ = FOUR_PI;
    match_newtonian_ = true;
  } else if (config.tpf_alpha > 0.0) {
    alpha_ = config.tpf_alpha;
    match_newtonian_ = false;
  } else {
    alpha_ = FOUR_PI;
    match_newtonian_ = true;
  }
  softening_override_ = config.tpf_softening;
}

void TPFPackage::compute_accelerations(const State& state,
                                       double bh_mass,
                                       double softening,
                                       bool star_star,
                                       std::vector<double>& ax,
                                       std::vector<double>& ay) const {
  const int n = state.n();
  ax.assign(n, 0.0);
  ay.assign(n, 0.0);

  const double eps = (softening_override_ > 0.0) ? softening_override_ : softening;
  const double eps2 = eps * eps;
  const double coeff = alpha_ * INV_FOUR_PI;

  // Operational readout: a = -grad(phi). For phi = -coeff * sum m_j/r_j (weak-field Green),
  // grad(phi) = coeff * sum m_j * (r-r_j) / |r-r_j|^3, so a = -grad(phi) = -coeff * sum m_j * (r-r_j)/r^3.
  // Unit vector from source to field point: (r_i - r_j) / r. Acceleration points toward source (attractive).
  // a_i += coeff * m_j * (r_j - r_i) / r^3 = coeff * m_j * r_hat / r^2 (magnitude).

  // BH at origin
  for (int i = 0; i < n; ++i) {
    double rx = state.x[i], ry = state.y[i];
    double r_sq = rx * rx + ry * ry + eps2;
    double r_cubed = r_sq * std::sqrt(r_sq);
    double acc_mag = coeff * bh_mass / r_cubed;
    ax[i] -= acc_mag * rx;
    ay[i] -= acc_mag * ry;
  }

  if (!star_star) return;

  // Star-star: superposition from all other particles
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      double dx = state.x[j] - state.x[i];
      double dy = state.y[j] - state.y[i];
      double r_sq = dx * dx + dy * dy + eps2;
      double r_cubed = r_sq * std::sqrt(r_sq);
      double acc_mag = coeff * state.mass[j] / r_cubed;
      ax[i] += acc_mag * dx;
      ay[i] += acc_mag * dy;
    }
  }
}

double TPFPackage::compute_potential_energy(const State& state,
                                            double bh_mass,
                                            double softening,
                                            bool star_star) const {
  const int n = state.n();
  const double eps = (softening_override_ > 0.0) ? softening_override_ : softening;
  const double eps2 = eps * eps;
  const double coeff = alpha_ * INV_FOUR_PI;

  // Weak-field scalar: phi = -coeff * sum_{j!=i} m_j / r_ij.
  // PE = (1/2) sum_i m_i * phi_i (factor 1/2 to avoid double counting pairs).
  double pe = 0.0;

  for (int i = 0; i < n; ++i) {
    double r = std::sqrt(state.x[i] * state.x[i] + state.y[i] * state.y[i] + eps2);
    pe -= coeff * bh_mass * state.mass[i] / r;
  }
  if (!star_star) return pe;

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dx = state.x[j] - state.x[i];
      double dy = state.y[j] - state.y[i];
      double r = std::sqrt(dx * dx + dy * dy + eps2);
      pe -= coeff * state.mass[i] * state.mass[j] / r;
    }
  }
  return pe;
}

}  // namespace galaxy
