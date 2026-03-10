#ifndef GALAXY_PHYSICS_TPFCORE_REGIME_DIAGNOSTICS_HPP
#define GALAXY_PHYSICS_TPFCORE_REGIME_DIAGNOSTICS_HPP

/**
 * Regime / field-strength diagnostics for TPFCore.
 *
 * Reporting only: no equation switching, no strong-field model.
 * Heuristic thresholds so researchers can see low-intensity (weak-like),
 * transitional, or high-intensity (provisional ansatz caution) regimes.
 */

namespace galaxy {
namespace tpfcore {

/** Theta Frobenius below this: low-intensity (weak-like) regime. */
constexpr double THETA_NORM_LOW_MAX = 1.0;
/** Theta Frobenius in [LOW_MAX, TRANSITIONAL_MAX): transitional. >= TRANSITIONAL_MAX: high-intensity. */
constexpr double THETA_NORM_TRANSITIONAL_MAX = 10.0;

/** Regime label from Theta Frobenius norm. Conservative wording. */
inline const char* regime_label_from_theta_norm(double theta_norm) {
  if (theta_norm < THETA_NORM_LOW_MAX) return "low-intensity";
  if (theta_norm < THETA_NORM_TRANSITIONAL_MAX) return "transitional";
  return "high-intensity (provisional ansatz caution)";
}

}  // namespace tpfcore
}  // namespace galaxy

#endif
