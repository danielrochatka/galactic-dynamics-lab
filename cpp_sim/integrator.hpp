#ifndef GALAXY_INTEGRATOR_HPP
#define GALAXY_INTEGRATOR_HPP

#include "types.hpp"

namespace galaxy {

// Single velocity Verlet step (same as Python):
// x_new = x + v*dt + 0.5*a*dt^2
// v_new = v + 0.5*(a + a_new)*dt
// Modifies state in place; ax, ay are scratch buffers (size n).
void velocity_verlet_step(State& state,
                          double bh_mass,
                          double softening,
                          bool star_star,
                          double dt,
                          std::vector<double>& ax,
                          std::vector<double>& ay);

}  // namespace galaxy

#endif
