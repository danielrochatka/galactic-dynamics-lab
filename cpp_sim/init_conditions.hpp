#ifndef GALAXY_INIT_CONDITIONS_HPP
#define GALAXY_INIT_CONDITIONS_HPP

#include "config.hpp"
#include "types.hpp"

namespace galaxy {

// Galaxy: rotating disk, uniform in area, v_circ from enclosed mass (BH + stars inside r).
void init_galaxy_disk(const Config& config, State& state, unsigned seed = 12345);

// Two-body: one star at (r0, 0), v = (0, speed_ratio * v_circ).
void init_two_body(const Config& config, State& state);

// Symmetric pair: two stars at (-a,0) and (a,0), v = (0, ±v).
void init_symmetric_pair(const Config& config, State& state);

// Small-N: ring with seed 42 (match Python), n in [3,10].
void init_small_n(const Config& config, State& state);

}  // namespace galaxy

#endif
