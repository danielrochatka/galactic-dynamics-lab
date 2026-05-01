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

// Semi-implicit Euler (symplectic Euler, kick-drift):
// a_n = a(x_n, v_n), v_{n+1}=v_n+a_n*dt, x_{n+1}=x_n+v_{n+1}*dt.
// Preferred when acceleration depends on velocity.
void semi_implicit_euler_step(State& state,
                              const PhysicsPackage* physics,
                              double bh_mass,
                              double softening,
                              bool star_star,
                              double dt,
                              std::vector<double>& ax,
                              std::vector<double>& ay);

}  // namespace galaxy

#endif
