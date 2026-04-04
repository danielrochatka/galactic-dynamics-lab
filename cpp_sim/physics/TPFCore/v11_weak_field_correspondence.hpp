#ifndef GALAXY_PHYSICS_TPFCORE_V11_WEAK_FIELD_CORRESPONDENCE_HPP
#define GALAXY_PHYSICS_TPFCORE_V11_WEAK_FIELD_CORRESPONDENCE_HPP

#include "../../config.hpp"
#include "../../types.hpp"
#include <string>
#include <vector>

namespace galaxy {

/**
 * Standalone audit: manuscript v11 static weak-field correspondence construction only.
 * - axis_monopole (default): Ξ, Θ, I, and the principal part of C_{μν} from Eq. (10) with ΔC_{μν} omitted.
 * - earth_moon_line_of_centers: Sec. XI line-of-centers φ (Eq.44) / aTPF (Eq.45) vs Newtonian (Eq.46) CSV.
 * Does not integrate particles or supply accelerations from TPF field equations.
 */
void run_v11_weak_field_correspondence_audit(const Config& config, const std::string& output_dir);

/**
 * Paper-backed weak-field correspondence dynamics helper (Eq. 42-44 scalar superposition).
 * Limited to static/quasi-static correspondence use; not full Eq. (9)/(10)+DeltaC dynamics.
 */
void compute_v11_weak_field_correspondence_accelerations(const State& state,
                                                         double bh_mass,
                                                         double softening,
                                                         bool star_star,
                                                         double alpha_si,
                                                         std::vector<double>& ax,
                                                         std::vector<double>& ay);

}  // namespace galaxy

#endif
