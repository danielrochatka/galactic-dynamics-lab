#ifndef GALAXY_PHYSICS_TPFCORE_PROVISIONAL_READOUT_HPP
#define GALAXY_PHYSICS_TPFCORE_PROVISIONAL_READOUT_HPP

/**
 * PROVISIONAL motion/readout layer for TPFCore.
 *
 * EXPLORATORY: This is NOT the full derived TPF dynamics. It maps the local
 * tensor field (Theta_ij) into an acceleration-like vector for the simulator,
 * without reverting to Newtonian scalar-potential dynamics.
 *
 * Design:
 * - Motion derived from local TPF tensor structure (Theta), NOT from Phi or -grad(Phi).
 * - Closures are downstream of the ansatz; see readout_closure.hpp for the boundary.
 * - Supported modes: tensor_radial_projection, tensor_radial_projection_negated,
 *   tr_coherence_readout / derived_tpf_radial_readout (manuscript Hessian + 1D radial Poisson profile),
 *   experimental_radial_r_scaling (experimental; provisional ansatz Theta, radial-only with r-scaling).
 */

#include "../../types.hpp"
#include "derived_tpf_radial.hpp"
#include <string>
#include <vector>

namespace galaxy {
struct Config;
namespace tpfcore {

/**
 * Compute provisional readout acceleration for one particle from all sources.
 * EXPLORATORY: not the full TPF dynamics. Tensor-driven, no Phi gradient.
 *
 * Supported modes: tensor_radial_projection, tensor_radial_projection_negated,
 * tr_coherence_readout / derived_tpf_radial_readout (1D radial Poisson + invariant I; uses derived_poisson),
 * experimental_radial_r_scaling (experimental radial-only closure with r-scaling).
 *
 * For derived radial modes: pass derived_poisson; optional derived_profile avoids rebuilding the radial profile
 * (required for efficient N-body timesteps when nullptr, profile is built each call — use batch API from TPFCorePackage).
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
  /* tr_coherence_readout only (0 otherwise): paper-aligned t-r provisional terms */
  double theta_rr = 0.0;
  double theta_tt = 0.0;
  double theta_tr = 0.0;
  double theta_rr_plus_theta_tt = 0.0;
  double provisional_radial_readout = 0.0;
  double provisional_tangential_readout = 0.0;
  /** tr_coherence_readout: TPF vs Newtonian regime label (empty for other modes). */
  std::string regime;
};

/**
 * Compute readout acceleration and Theta for diagnostics.
 * Same logic as compute_provisional_readout_acceleration but also returns Theta, I, and (for tr_coherence) t-r terms.
 */
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

/**
 * Write tpf_readout_debug.csv for dynamical runs.
 * Column layout and mode-specific fields are owned by the readout module.
 * Caller passes resolved run params (no Config dependency in readout).
 */
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
