#ifndef GALAXY_PHYSICS_TPFCORE_SOURCE_ANSATZ_HPP
#define GALAXY_PHYSICS_TPFCORE_SOURCE_ANSATZ_HPP

#include <cmath>

/**
 * PROVISIONAL weak-field point-source ansatz for TPFCore.
 *
 * Uses the paper's weak-field point-source construction:
 *   Xi_i = partial_i Phi
 *   Theta_ij = Hess_ij(Phi) + B(r) delta_ij
 * with a softened point-source scalar Phi for numerical regularization.
 *
 * Phi = -M / sqrt(r^2 + eps^2)  (softened 1/r)
 * B(r) = c * M / (r^2 + eps^2)^(3/2)  (exploratory isotropic correction; c configurable, default 0)
 * r^2 = dx^2 + dy^2, dx = x - xs, dy = y - ys
 *
 * WARNING: This is a provisional weak-field point-source ansatz, NOT the full
 * nonlinear TPF source law. Kept for inspection-only use. The B(r) term is an
 * exploratory residual-reduction correction, not a derived TPF source law.
 */

namespace galaxy {
namespace tpfcore {

/** Displacement field Xi in 2D plane. */
struct Xi2D {
  double x, y;
};

/** Theta tensor components in 2D plane (xx, xy, yy). */
struct Theta2D {
  double xx, xy, yy;
  double trace() const { return xx + yy; }
};

/** Combined Xi and Theta from a single evaluation. */
struct PointSourceField {
  Xi2D xi;
  Theta2D theta;
};

/**
 * PROVISIONAL weak-field point-source: Phi = -M/sqrt(r^2+eps^2).
 * Xi_i = partial_i Phi (unchanged).
 * Theta_ij = Hess_ij(Phi) + B(r) delta_ij, with B(r) = c*M/(r^2+eps^2)^(3/2).
 *
 * Source at (xs, ys) with mass m; field point (x, y); eps = softening; c = isotropic correction coefficient.
 */
PointSourceField provisional_point_source_field(double xs, double ys, double m,
                                                double x, double y, double eps, double c = 0.0);

/** Legacy: Theta only (superposition by caller). */
Theta2D provisional_point_source_theta(double xs, double ys, double m,
                                       double x, double y, double eps, double c = 0.0);

/** Lambda = 1/4 in 4D (paper-fixed, not tunable). */
constexpr double LAMBDA_4D = 0.25;

/** Residual vector R_nu = partial_i (Theta_i_nu - lambda delta_i_nu Theta). */
struct Residual2D {
  double x, y;
  double norm() const { return std::sqrt(x * x + y * y); }
};

/**
 * Analytic residual for the configuration equation nabla_mu(Theta^mu_nu - lambda delta^mu_nu Theta).
 * Spatial weak-field form: R_nu = partial_x(...) + partial_y(...).
 * Computed from closed-form derivatives of Phi and B(r). Single source only.
 */
Residual2D provisional_point_source_residual(double xs, double ys, double m,
                                             double x, double y, double eps, double c = 0.0);

/**
 * Invariant I = Theta_{mu nu} Theta^{mu nu} - lambda * Theta^2.
 * For 2D Euclidean: Theta_mn Theta^mn = Theta_xx^2 + Theta_yy^2 + 2*Theta_xy^2.
 */
double compute_invariant_I(const Theta2D& theta);

}  // namespace tpfcore
}  // namespace galaxy

#endif
