#ifndef GALAXY_PHYSICS_TPFCORE_READOUT_DIAGNOSTICS_HPP
#define GALAXY_PHYSICS_TPFCORE_READOUT_DIAGNOSTICS_HPP

#include <string>

namespace galaxy {
namespace tpfcore {

/** Per-particle readout diagnostics: a, Theta, I, and derived quantities. */
struct ReadoutDiagnostics {
  double ax, ay;
  double theta_xx, theta_xy, theta_yy, theta_trace, invariant_I;
  /** Theta Frobenius norm (configuration intensity); for regime diagnostics. */
  double theta_norm = 0.0;
  /* Derived-radial closure: diagnostic theta components (not added to ax, ay on that path). */
  double theta_rr = 0.0;
  double theta_tt = 0.0;
  double theta_tr = 0.0;
  double theta_rr_plus_theta_tt = 0.0;
  double provisional_radial_readout = 0.0;
  double provisional_tangential_readout = 0.0;
  /** Optional regime label when populated (often empty for non–derived-radial modes). */
  std::string regime;
};

}  // namespace tpfcore
}  // namespace galaxy

#endif
