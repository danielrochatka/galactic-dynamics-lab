#ifndef GALAXY_FORCE_COMPARE_HPP
#define GALAXY_FORCE_COMPARE_HPP

#include "config.hpp"
#include <string>

namespace galaxy {

/**
 * Newtonian-vs-TPF acceleration comparison diagnostic (diagnostics only; no physics change).
 * Requires physics_package=TPFCore and tpfcore_enable_provisional_readout=true so TPF
 * uses the current readout mode and calibrated scale. Writes tpf_newtonian_force_compare.csv
 * and tpf_newtonian_force_compare.txt to output_dir.
 */
void run_tpf_newtonian_force_compare(const Config& config, const std::string& output_dir);

/**
 * Consistency audit between weak_field_calibration and force_compare diagnostics.
 * Same points, side-by-side intermediates; decisive conclusion. Diagnostics only.
 */
void run_tpf_diagnostic_consistency_audit(const Config& config, const std::string& output_dir);

}  // namespace galaxy

#endif
