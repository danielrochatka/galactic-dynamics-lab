#ifndef GALAXY_PHYSICS_TPFCORE_PROVISIONAL_READOUT_HPP
#define GALAXY_PHYSICS_TPFCORE_PROVISIONAL_READOUT_HPP

/**
 * PROVISIONAL readout/diagnostic helper layer for TPFCore.
 *
 * Branch status: the active canonical runtime route is tpf_xi_theta_v1.
 * Active v1 motion is Xi_total-driven (a = -K_xi * Xi_total_spatial), with Theta retained
 * as diagnostic-only configuration-gradient output.
 *
 * The provisional functions declared here are retained temporarily for readout/diagnostic
 * tooling cleanup and are not active v1 motion mechanics.
 *
 * Modes (dispatch in compute_provisional_readout_acceleration):
 * - tensor_radial_projection / _negated: per-source Theta projected along r_hat, superposed.
 * - tr_coherence_readout, derived_tpf_radial_readout: same hybrid radial closure (a_s r̂ from
 *   radial_acceleration_scalar_derived; κ–I ledger). Extra theta_tt/tr terms are diagnostics only.
 * - experimental_radial_r_scaling: separate radial closure from theta_rr.
 */

#include "../../types.hpp"
#include "derived_tpf_radial.hpp"
#include "readout_diagnostics.hpp"
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
