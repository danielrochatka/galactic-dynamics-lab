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
 *   - Exploratory:      isotropic correction c; NOT a fundamental constant.
 *   - Provisional:      readout closure (mode, scale, theta_tt/tr); experimental.
 *   - Inspection:       probe geometry and dump flags; one-time/calibration use.
 *
 * The simulator passes Config; the package (tpf_core_package.cpp) converts
 * Config -> TPFCoreParams. All other TPFCore code consumes TPFCoreParams or
 * primitive args, not Config.
 */

#include <string>

namespace galaxy {
namespace tpfcore {

struct TPFCoreParams {
  /* --- Simulator context (passed through for dynamics/output) --- */
  std::string output_dir;
  double softening = 1.0;
  double bh_mass = 1000.0;
  double star_mass = 0.05;
  bool enable_star_star_gravity = true;

  /* --- Numerical regularization (eps for Phi; not theory) --- */
  double tpfcore_source_softening = 0.0;   /**< If > 0, use this as eps; else use global softening. */
  double effective_source_softening = 1.0; /**< Resolved: tpfcore_source_softening > 0 ? that : softening. */

  /* --- Exploratory ansatz correction (c in B(r)=c*M/(r^2+eps^2)^(3/2); NOT a fundamental constant) --- */
  double tpfcore_isotropic_correction_c = 0.0;

  /* --- Inspection / calibration (probe geometry and dump flags) --- */
  double tpfcore_probe_radius_min = 1.0;
  double tpfcore_probe_radius_max = 50.0;
  int tpfcore_probe_samples = 100;
  bool tpfcore_dump_theta_profile = true;
  bool tpfcore_dump_invariant_profile = true;

  /* --- Provisional readout (experimental closure; not source-theory) --- */
  bool tpfcore_dump_readout_debug = true;

  /* --- C-sweep: exploratory calibration tool (fitted c is NOT a paper constant) --- */
  double tpfcore_c_sweep_min = -0.5;
  double tpfcore_c_sweep_max = 1.0;
  int tpfcore_c_sweep_steps = 31;
  std::string tpfcore_c_objective = "max_residual_norm";

  /* --- Validation geometry --- */
  double validation_symmetric_separation = 30.0;
};

}  // namespace tpfcore
}  // namespace galaxy

#endif
