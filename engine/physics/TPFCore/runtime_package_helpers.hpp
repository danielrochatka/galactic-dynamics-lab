#ifndef GALAXY_PHYSICS_TPFCORE_RUNTIME_PACKAGE_HELPERS_HPP
#define GALAXY_PHYSICS_TPFCORE_RUNTIME_PACKAGE_HELPERS_HPP

#include "../../types.hpp"
#include <vector>

namespace galaxy {
namespace tpfcore {

unsigned apply_global_accel_magnitude_shunt(const State& state,
                                            double dt,
                                            bool enable,
                                            double fraction,
                                            std::vector<double>& ax,
                                            std::vector<double>& ay);

void reset_global_accel_shunt_events();
unsigned global_accel_shunt_events();

}  // namespace tpfcore
}  // namespace galaxy

#endif
