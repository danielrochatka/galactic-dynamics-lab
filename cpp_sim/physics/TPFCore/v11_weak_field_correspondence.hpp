#ifndef GALAXY_PHYSICS_TPFCORE_V11_WEAK_FIELD_CORRESPONDENCE_HPP
#define GALAXY_PHYSICS_TPFCORE_V11_WEAK_FIELD_CORRESPONDENCE_HPP

#include "../../config.hpp"
#include <string>

namespace galaxy {

/**
 * Standalone audit: manuscript v11 static weak-field correspondence construction only.
 * Computes Ξ, Θ, I, and the principal part of C_{μν} from Eq. (10) with ΔC_{μν} omitted.
 * Does not integrate particles or supply accelerations from TPF field equations.
 */
void run_v11_weak_field_correspondence_audit(const Config& config, const std::string& output_dir);

}  // namespace galaxy

#endif
