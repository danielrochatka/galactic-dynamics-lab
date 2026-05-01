#ifndef GALAXY_PHYSICS_TPFCORE_EXTENSIONS_VDSG_HPP
#define GALAXY_PHYSICS_TPFCORE_EXTENSIONS_VDSG_HPP

#include "../../types.hpp"
#include <string>
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
                                       const std::string& vdsg_mode,
                                       double mass_gate_m0_kg,
                                       double mass_gate_alpha,
                                       double x_clamp,
                                       bool weak_field_gate_enable,
                                       double weak_field_a0,
                                       double weak_field_power,
                                       double bounded_amplitude,
                                       std::vector<double>& ax,
                                       std::vector<double>& ay);

}  // namespace tpfcore
}  // namespace galaxy

#endif
