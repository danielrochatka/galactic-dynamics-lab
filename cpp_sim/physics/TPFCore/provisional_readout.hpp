#ifndef GALAXY_PHYSICS_TPFCORE_PROVISIONAL_READOUT_HPP
#define GALAXY_PHYSICS_TPFCORE_PROVISIONAL_READOUT_HPP

/**
 * PROVISIONAL motion/readout layer for TPFCore.
 *
 * This is EXPLORATORY and NOT the full derived TPF dynamics. It maps the local
 * tensor field (Theta_ij) into an acceleration-like vector for the simulator,
 * without reverting to Newtonian scalar-potential dynamics.
 *
 * Design:
 * - Motion derived from local TPF tensor structure (Theta), NOT from Phi or -grad(Phi).
 * - Tensor-driven readout: Theta_ij contracted with a local directional vector.
 * - Explicitly isolated from source_ansatz and superposition logic.
 * - Easy to swap out later when proper TPF motion laws are derived.
 *
 * Supported readout_mode: tensor_radial_projection
 * - Per source: Theta_s from that source at particle location.
 * - r_hat = unit vector from source to particle.
 * - Contribution: scale * (Theta_s · r_hat)  [tensor-vector contraction]
 * - Superpose over all sources (BH + other particles when star_star).
 */

#include "../../types.hpp"
#include <string>

namespace galaxy {
struct Config;
namespace tpfcore {

/**
 * Compute provisional readout acceleration for one particle from all sources.
 * EXPLORATORY: not the full TPF dynamics. Tensor-driven, no Phi gradient.
 *
 * @param state     Full state (particle positions, masses)
 * @param i         Index of particle to compute acceleration for
 * @param bh_mass   Fixed BH mass at origin (use 0 if no BH)
 * @param star_star Include contributions from other particles
 * @param softening Global softening
 * @param c         Isotropic correction coefficient (from ansatz)
 * @param readout_mode "tensor_radial_projection" supported
 * @param readout_scale Scale factor for readout magnitude
 * @param[out] ax   x-component of acceleration
 * @param[out] ay   y-component of acceleration
 */
void compute_provisional_readout_acceleration(const State& state,
                                               int i,
                                               double bh_mass,
                                               bool star_star,
                                               double softening,
                                               double source_softening,
                                               double c,
                                               const std::string& readout_mode,
                                               double readout_scale,
                                               double& ax,
                                               double& ay);

}  // namespace tpfcore
}  // namespace galaxy

#endif
