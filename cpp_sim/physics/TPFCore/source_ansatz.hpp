#ifndef GALAXY_PHYSICS_TPFCORE_SOURCE_ANSATZ_HPP
#define GALAXY_PHYSICS_TPFCORE_SOURCE_ANSATZ_HPP

/**
 * PROVISIONAL source ansatz for TPFCore.
 *
 * The paper defines Theta_{mu nu} = nabla_mu Xi_nu but does not fully specify
 * the point-source primitive field. This module isolates a clearly labeled
 * PROVISIONAL ansatz until the paper provides a complete specification.
 *
 * WARNING: This is NOT derived from the full variational TPF equations.
 * It is a placeholder to enable inspection of the invariant I.
 */

namespace galaxy {
namespace tpfcore {

/** Theta tensor components in 2D plane (xx, xy, yy). */
struct Theta2D {
  double xx, xy, yy;
  double trace() const { return xx + yy; }
};

/**
 * PROVISIONAL: Point-source contribution to Theta at a field point.
 * Source at (xs, ys) with mass m; field point (x, y); softening eps.
 *
 * Ansatz: Radially isotropic Theta_xx = Theta_yy = A*m/(r_soft), Theta_xy = 0,
 * with r_soft = sqrt(dx^2 + dy^2 + eps^2). A=1 provisional (dimensionally G-like).
 * This is NOT the paper's full point-source specification.
 */
Theta2D provisional_point_source_theta(double xs, double ys, double m,
                                       double x, double y, double eps);

/** Lambda = 1/4 in 4D (paper-fixed, not tunable). */
constexpr double LAMBDA_4D = 0.25;

/**
 * Invariant I = Theta_{mu nu} Theta^{mu nu} - lambda * Theta^2.
 * For 2D Euclidean: Theta_mn Theta^mn = Theta_xx^2 + Theta_yy^2 + 2*Theta_xy^2.
 */
double compute_invariant_I(const Theta2D& theta);

}  // namespace tpfcore
}  // namespace galaxy

#endif
