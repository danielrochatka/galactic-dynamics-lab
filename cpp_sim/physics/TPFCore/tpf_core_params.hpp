#ifndef GALAXY_PHYSICS_TPFCORE_TPF_CORE_PARAMS_HPP
#define GALAXY_PHYSICS_TPFCORE_TPF_CORE_PARAMS_HPP

/**
 * Package-local parameter bundle for TPFCore.
 *
 * Parameters are grouped and commented by scientific role so the code is
 * honest about what is theory, regularization, exploratory, or provisional.
 *
 * Role classification:
 *   - Fixed theory:     lambda = 1/4 (LAMBDA_4D in source_ansatz.hpp; not tunable).
 *   - Numerical reg.:   source softening (eps) for Phi; not part of the theory.
 *   - Provisional:      readout closure (mode, scale, theta_tt/tr); VDSG params live on Config (separate path).
 *   - Inspection:       probe geometry and dump flags; one-time/calibration use.
 *
 * The simulator passes Config; the package (tpf_core_package.cpp) converts
 * Config -> TPFCoreParams. All other TPFCore code consumes TPFCoreParams or
 * primitive args, not Config.
 */

#include <string>

#include "../../config.hpp"

namespace galaxy {
namespace tpfcore {

struct TPFCoreParams {
  /* --- Simulator context (passed through for dynamics/output) --- */
  std::string output_dir;
  double softening = kDefaultSofteningM;
  double bh_mass = kDefaultBhMassKg;
  double star_mass = kSolarMassKg;
  bool enable_star_star_gravity = true;

  /* --- Numerical regularization (eps for Phi; not theory) --- */
  double tpfcore_source_softening = 0.0;   /**< If > 0, use this as eps; else use global softening. */
  double effective_source_softening = kDefaultSofteningM; /**< Resolved: tpfcore_source_softening > 0 ? that : softening. */

  /* --- Inspection / calibration (probe geometry and dump flags) --- */
  /** Defaults keep v11_weak_field_correspondence earth_moon_line_of_centers valid (max < D ~ 3.84e8 m). */
  double tpfcore_probe_radius_min = 1.0;
  double tpfcore_probe_radius_max = 50.0;
  int tpfcore_probe_samples = 100;
  bool tpfcore_dump_theta_profile = true;
  bool tpfcore_dump_invariant_profile = true;

  /* --- Provisional readout (experimental closure; not source-theory) --- */
  bool tpfcore_dump_readout_debug = true;

  /* --- Validation geometry --- */
  double validation_symmetric_separation = 7.48e10;
};

}  // namespace tpfcore
}  // namespace galaxy

#endif
