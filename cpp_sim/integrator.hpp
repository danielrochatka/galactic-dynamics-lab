#ifndef GALAXY_INTEGRATOR_HPP
#define GALAXY_INTEGRATOR_HPP

#include "types.hpp"
#include "physics/physics_package.hpp"

namespace galaxy {

// Single velocity Verlet step (same as Python):
// x_new = x + v*dt + 0.5*a*dt^2
// v_new = v + 0.5*(a + a_new)*dt
// Uses the given physics package for acceleration. Modifies state in place; ax, ay are scratch buffers.
void velocity_verlet_step(State& state,
                          const PhysicsPackage* physics,
                          double bh_mass,
                          double softening,
                          bool star_star,
                          double dt,
                          std::vector<double>& ax,
                          std::vector<double>& ay);

}  // namespace galaxy

#endif
