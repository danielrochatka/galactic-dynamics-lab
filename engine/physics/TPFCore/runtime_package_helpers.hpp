#ifndef GALAXY_PHYSICS_TPFCORE_RUNTIME_PACKAGE_HELPERS_HPP
#define GALAXY_PHYSICS_TPFCORE_RUNTIME_PACKAGE_HELPERS_HPP

#include "../../types.hpp"
#include <vector>

namespace galaxy {
namespace tpfcore {

struct XiWakeKinematics {
  double v_radial = 0.0;
  double v_transverse = 0.0;
  double beta_pass = 0.0;
  double wake_gate = 0.0;
  double beta_effective = 0.0;
  double axis_x = 0.0;
  double axis_y = 0.0;
  double axis_z = 0.0;
  bool has_axis = false;
};

XiWakeKinematics compute_xi_wake_kinematics(double dx,
                                            double dy,
                                            double dz,
                                            double vx_rel,
                                            double vy_rel,
                                            double vz_rel,
                                            double c_light,
                                            bool post_pass_gate);

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
