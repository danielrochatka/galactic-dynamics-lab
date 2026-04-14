#ifndef GALAXY_PHYSICS_TPFCORE_EXTENSIONS_VDSG_HPP
#define GALAXY_PHYSICS_TPFCORE_EXTENSIONS_VDSG_HPP

#include "../../types.hpp"
#include <vector>

namespace galaxy {
namespace tpfcore {

double vdsg_effective_coupling(double lambda0, double source_mass_kg, double baseline_mass_kg);

void accumulate_vdsg_velocity_modifier(const State& state,
                                       double bh_mass,
                                       double softening,
                                       bool star_star,
                                       double vdsg_coupling,
                                       double vdsg_mass_baseline_kg,
                                       std::vector<double>& ax,
                                       std::vector<double>& ay);

}  // namespace tpfcore
}  // namespace galaxy

#endif
