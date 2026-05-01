#ifndef GALAXY_PHYSICS_TPFCORE_RUNTIME_PACKAGE_HELPERS_HPP
#define GALAXY_PHYSICS_TPFCORE_RUNTIME_PACKAGE_HELPERS_HPP

#include "../../types.hpp"
#include <vector>

namespace galaxy {
namespace tpfcore {

struct XiWakeKinematics {
  // v_radial = dot(v_rel, r_hat) (positive for separating motion along source->target radial direction).
  double v_radial = 0.0;
  // v_transverse = ||v_rel - v_radial * r_hat||.
  double v_transverse = 0.0;
  // beta_pass = v_transverse / c.
  double beta_pass = 0.0;
  // wake_gate is mode-dependent: post-pass threshold gate (wake mode) or continuous gate.
  double wake_gate = 0.0;
  // beta_effective = beta_pass * wake_gate.
  double beta_effective = 0.0;
  // axis_* is a transverse-passing direction diagnostic; not necessarily the final metric deformation axis.
  double axis_x = 0.0;
  double axis_y = 0.0;
  double axis_z = 0.0;
  bool has_axis = false;
};

// compute_xi_wake_kinematics(..., post_pass_gate):
// - post_pass_gate=true  => metric_transverse_wake (post-pass/separation threshold gate).
// - post_pass_gate=false => metric_transverse_continuous (legacy continuous transverse gate).
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
