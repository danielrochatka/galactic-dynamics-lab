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

}  // namespace galaxy

#endif
