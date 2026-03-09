#ifndef GALAXY_PHYSICS_PACKAGE_HPP
#define GALAXY_PHYSICS_PACKAGE_HPP

/**
 * Physics package interface.
 * Every physics package (Newtonian, custom, etc.) must implement this interface.
 * The simulator and integrator depend only on this interface, not on any specific implementation.
 */

#include "../types.hpp"
#include <memory>
#include <string>
#include <vector>

namespace galaxy {

class PhysicsPackage {
 public:
  virtual ~PhysicsPackage() = default;

  /** Package identifier (e.g. "Newtonian"). Must match config physics_package name. */
  virtual const char* name() const = 0;

  /**
   * Compute accelerations for all particles. Required.
   * Fills ax[i], ay[i] for each particle i. Same signature used by the integrator.
   */
  virtual void compute_accelerations(const State& state,
                                     double bh_mass,
                                     double softening,
                                     bool star_star,
                                     std::vector<double>& ax,
                                     std::vector<double>& ay) const = 0;

  /**
   * Optional: total gravitational potential energy for diagnostics/validation.
   * Default returns 0. Override in packages that define a potential.
   */
  virtual double compute_potential_energy(const State& state,
                                         double bh_mass,
                                         double softening) const {
    (void)state;
    (void)bh_mass;
    (void)softening;
    return 0.0;
  }

  /** Optional: called once before simulation (e.g. load tables). Default no-op. */
  virtual void init() {}

  /** Optional: validation/reporting name or hook. Default returns name(). */
  virtual const char* validation_name() const { return name(); }
};

/** Generic kinetic energy (same for all packages): 0.5 * sum(m_i * v_i^2). */
double compute_kinetic_energy(const State& state);

/**
 * Registry: get a physics package by name.
 * Returns nullptr if unknown. Caller must not delete the pointer (static/registry-owned).
 */
PhysicsPackage* get_physics_package(const std::string& name);

/** Returns true if the named package is available. */
bool has_physics_package(const std::string& name);

}  // namespace galaxy

#endif
