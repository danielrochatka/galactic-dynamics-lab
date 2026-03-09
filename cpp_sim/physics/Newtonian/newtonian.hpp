#ifndef GALAXY_PHYSICS_NEWTONIAN_HPP
#define GALAXY_PHYSICS_NEWTONIAN_HPP

#include "../physics_package.hpp"

namespace galaxy {

/** Newtonian gravity package: BH at origin + optional pairwise star-star with softening. */
class NewtonianPackage : public PhysicsPackage {
 public:
  const char* name() const override { return "Newtonian"; }

  void compute_accelerations(const State& state,
                             double bh_mass,
                             double softening,
                             bool star_star,
                             std::vector<double>& ax,
                             std::vector<double>& ay) const override;

  double compute_potential_energy(const State& state,
                                 double bh_mass,
                                 double softening) const override;
};

}  // namespace galaxy

#endif
