#ifndef GALAXY_FORCE_COMPARE_HPP
#define GALAXY_FORCE_COMPARE_HPP

#include "config.hpp"
#include <string>

namespace galaxy {

/**
 * Newtonian-vs-TPF acceleration comparison diagnostic (diagnostics only; no physics change).
 * Requires physics_package=TPFCore. Writes tpf_newtonian_force_compare.csv
 * and tpf_newtonian_force_compare.txt to output_dir.
 */
void run_tpf_newtonian_force_compare(const Config& config, const std::string& output_dir);

/**
 * Quarantined on tpf_xi_theta_v1 because it depends on provisional readout
 * diagnostics that are not part of the v1 Xi/Theta runtime path.
 */
bool run_tpf_diagnostic_consistency_audit(const Config& config, const std::string& output_dir);

}  // namespace galaxy

#endif
