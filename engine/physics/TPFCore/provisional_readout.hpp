#ifndef GALAXY_PHYSICS_TPFCORE_PROVISIONAL_READOUT_HPP
#define GALAXY_PHYSICS_TPFCORE_PROVISIONAL_READOUT_HPP

/**
 * PROVISIONAL motion/readout layer for TPFCore.
 *
 * EXPLORATORY closures downstream of the source ansatz (see readout_closure.hpp).
 * Not the full derived TPF dynamics.
 *
 * Modes (dispatch in compute_provisional_readout_acceleration):
 * - tensor_radial_projection / _negated: per-source Theta projected along r_hat, superposed.
 * - tr_coherence_readout, derived_tpf_radial_readout: same hybrid radial closure (a_s r̂ from
 *   radial_acceleration_scalar_derived; κ–I ledger). Extra theta_tt/tr terms are diagnostics only.
 * - experimental_radial_r_scaling: separate radial closure from theta_rr.
 *
 * Integrator note: TPFCorePackage::compute_accelerations may fill ax, ay from VDSG instead;
 * when tpf_vdsg_coupling != 0, readout closures here are not used for ax, ay on that path.
 */

#include "../../types.hpp"
#include "derived_tpf_radial.hpp"
#include <string>
#include <vector>

namespace galaxy {
struct Config;
namespace tpfcore {

/**
 * Provisional readout acceleration for one particle (readout path only; not used for ax, ay when VDSG active).
 *
 * Derived radial modes: pass derived_poisson; optional derived_profile avoids rebuilding the radial
 * profile each particle (batch from TPFCorePackage).
 */
void compute_provisional_readout_acceleration(const State& state,
                                               int i,
                                               double bh_mass,
                                               bool star_star,
                                               double softening,
                                               double source_softening,
                                               const std::string& readout_mode,
                                               double readout_scale,
                                               double theta_tt_scale,
                                               double theta_tr_scale,
                                               double& ax,
                                               double& ay,
                                               const DerivedTpfPoissonConfig* derived_poisson = nullptr,
                                               const TpfRadialGravityProfile* derived_profile = nullptr);

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

void compute_provisional_readout_with_diagnostics(const State& state,
                                                   int i,
                                                   double bh_mass,
                                                   bool star_star,
                                                   double softening,
                                                   double source_softening,
                                                   const std::string& readout_mode,
                                                   double readout_scale,
                                                   double theta_tt_scale,
                                                   double theta_tr_scale,
                                                   double& ax,
                                                   double& ay,
                                                   ReadoutDiagnostics& diag,
                                                   const DerivedTpfPoissonConfig* derived_poisson = nullptr,
                                                   const TpfRadialGravityProfile* derived_profile = nullptr);

void write_readout_debug_csv(const std::vector<Snapshot>& snapshots,
                             const std::string& output_dir,
                             double softening,
                             double bh_mass,
                             bool star_star,
                             double source_softening,
                             const std::string& readout_mode,
                             double readout_scale,
                             double theta_tt_scale,
                             double theta_tr_scale,
                             const DerivedTpfPoissonConfig& derived_poisson = DerivedTpfPoissonConfig());

}  // namespace tpfcore
}  // namespace galaxy

#endif
