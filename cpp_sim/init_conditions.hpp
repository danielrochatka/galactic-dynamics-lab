#ifndef GALAXY_INIT_CONDITIONS_HPP
#define GALAXY_INIT_CONDITIONS_HPP

#include "config.hpp"
#include "types.hpp"

namespace galaxy {

// Galaxy: annulus disk, template-driven ICs (see galaxy_init.hpp); v_circ from enclosed mass or TPF derived radial.
void init_galaxy_disk(const Config& config, State& state);

// Earth–Moon SI benchmark for two_body_orbit: n=2, clears and resizes state arrays.
void init_two_body(const Config& config, State& state);

// Legacy validation: one particle at (r0,0) with v from speed_ratio * v_circ(BH); BH mass from config.
void init_two_body_star_around_bh(const Config& config, State& state);

// Symmetric pair: two stars at (-a,0) and (a,0), v = (0, ±v).
void init_symmetric_pair(const Config& config, State& state);

// Small-N: ring with seed 42 (match Python), n in [3,10].
void init_small_n(const Config& config, State& state);

}  // namespace galaxy

#endif
