#ifndef GALAXY_PHYSICS_TPFCORE_SOURCE_ANSATZ_HPP
#define GALAXY_PHYSICS_TPFCORE_SOURCE_ANSATZ_HPP

#include <cmath>

/**
 * PROVISIONAL weak-field point-source ansatz for TPFCore.
 *
 * True 3D scalar potential evaluated on the z = 0 simulation plane (source and
 * field point both at z = 0):
 *   Phi = -M / R,   R = sqrt(dx^2 + dy^2 + eps^2)
 * where dx = x - xs, dy = y - ys, and eps > 0 is isotropic softening (same as
 * embedding a perpendicular offset in R). Equivalently R is the 3D distance with
 * z = z_s = 0 and a constant eps^2 added under the radical.
 *
 *   Xi_i = partial_i Phi   (i = x, y; Xi_z = 0 on the plane)
 *   Theta_ij = Hess_ij(Phi)   (full symmetric 3x3; Theta_xz = Theta_yz = 0 on z = 0)
 *
 * The field-equation residual uses the full spatial divergence in 3D, including
 * partial_z Theta_{xz} and partial_z Theta_{yz} (nonzero as limits even when
 * Theta_xz = Theta_yz = 0 on the plane). For Phi = -m/R with R^2 = r^2 + eps^2,
 * the residual vanishes identically when eps = 0 away from the source; with
 * softening it is O(eps^2).
 */

namespace galaxy {
namespace tpfcore {

/** Displacement field Xi in 2D plane (simulation). */
struct Xi2D {
  double x, y;
};

/** Symmetric Theta tensor on the z = 0 slice (3D Hessian of Phi). */
struct Theta3D {
  double xx, xy, xz, yy, yz, zz;
  double trace() const { return xx + yy + zz; }
};

/** Combined Xi and Theta from a single evaluation. */
struct PointSourceField {
  Xi2D xi;
  Theta3D theta;
};

/**
 * Single point source at (xs, ys); field at (x, y); eps = softening.
 * Theta is the 3D Hessian of Phi = -m/R with R^2 = dx^2 + dy^2 + eps^2.
 */
PointSourceField provisional_point_source_field(double xs, double ys, double m,
                                                double x, double y, double eps);

/** Legacy: Theta only (superposition by caller). */
Theta3D provisional_point_source_theta(double xs, double ys, double m,
                                         double x, double y, double eps);

/** Lambda = 1/4 in 4D. Fixed by manuscript structure; must not become tunable. */
constexpr double LAMBDA_4D = 0.25;

/** Residual vector for the x,y components of the configuration equation on the plane. */
struct Residual2D {
  double x, y;
  double norm() const { return std::sqrt(x * x + y * y); }
};

/**
 * R_nu = partial_i (Theta_i_nu - lambda delta_i_nu Theta) for nu = x, y, with
 * i summed over x, y, z. On z = 0, Theta_xz = Theta_yz = 0 but partial_z of
 * those components contributes. Single source only.
 */
Residual2D provisional_point_source_residual(double xs, double ys, double m,
                                             double x, double y, double eps);

/**
 * Invariant I = Theta_{mu nu} Theta^{mu nu} - lambda * Theta^2.
 * Spatial Euclidean: sum of squares of all nine symmetric entries minus lambda trace^2.
 */
double compute_invariant_I(const Theta3D& theta);

/**
 * Frobenius norm sqrt(sum_ij Theta_ij^2) for the symmetric 3x3 tensor.
 */
inline double theta_frobenius_norm(const Theta3D& theta) {
  return std::sqrt(theta.xx * theta.xx + theta.yy * theta.yy + theta.zz * theta.zz +
                   2.0 * (theta.xy * theta.xy + theta.xz * theta.xz + theta.yz * theta.yz));
}

}  // namespace tpfcore
}  // namespace galaxy

#endif
