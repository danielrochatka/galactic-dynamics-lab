#ifndef GALAXY_PHYSICS_TPFCORE_FIELD_EVALUATION_HPP
#define GALAXY_PHYSICS_TPFCORE_FIELD_EVALUATION_HPP

/**
 * Field-evaluation layer for TPFCore.
 *
 * Packages the provisional field at a sample point: Xi, Theta, invariant I,
 * and residual (when available). Built from source_ansatz only; no motion/readout.
 * Downstream (inspection, provisional readout) consume this layer.
 */

#include "source_ansatz.hpp"
#include "../../types.hpp"

namespace galaxy {
namespace tpfcore {

/** Evaluated provisional field at one point (single- or multi-source). */
struct FieldAtPoint {
  Xi2D xi;
  Theta3D theta;
  double invariant_I;
  bool has_residual;
  Residual2D residual;
};

/**
 * Evaluate provisional field at (x, y) from a single point source.
 * Residual is set (has_residual = true); uses analytic formula.
 */
FieldAtPoint evaluate_provisional_field_single_source(double xs, double ys, double m,
                                                      double x, double y, double eps);

/**
 * Evaluate provisional field at particle i from all sources (BH + optional star_star).
 * Superposed Xi and Theta; I from combined Theta. Residual not set (has_residual = false).
 */
FieldAtPoint evaluate_provisional_field_multi_source(const State& state, int i,
                                                     double bh_mass, bool star_star,
                                                     double eps);

/**
 * Combine two single-source field evaluations (e.g. symmetric pair).
 * Xi and Theta summed; I from combined Theta; residual summed when both have_residual.
 */
FieldAtPoint add_provisional_fields(const FieldAtPoint& a, const FieldAtPoint& b);

}  // namespace tpfcore
}  // namespace galaxy

#endif
