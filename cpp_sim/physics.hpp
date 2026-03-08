#ifndef GALAXY_PHYSICS_HPP
#define GALAXY_PHYSICS_HPP

#include "types.hpp"

namespace galaxy {

// Newtonian acceleration: BH at origin + optional pairwise star-star.
// Same formula as Python: a_i -= bh_mass * r_i / (r_sq * r_mag), r_sq = r^2 + eps^2.
// Star-star: a_i += sum_{j!=i} m_j * dr_ij / (r_sq * r_mag).
void compute_accelerations(const State& state,
                           double bh_mass,
                           double softening,
                           bool star_star,
                           std::vector<double>& ax,
                           std::vector<double>& ay);

double compute_kinetic_energy(const State& state);
double compute_potential_energy(const State& state, double bh_mass, double softening);

}  // namespace galaxy

#endif
