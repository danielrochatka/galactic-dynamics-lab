#ifndef GALAXY_SIMULATION_HPP
#define GALAXY_SIMULATION_HPP

#include "config.hpp"
#include "types.hpp"
#include <vector>

namespace galaxy {

// Run simulation and collect snapshots (step 0, then every snapshot_every).
// Uses config.dt, config.n_steps, config.snapshot_every, config.bh_mass,
// config.softening, config.enable_star_star_gravity.
std::vector<Snapshot> run_simulation(const Config& config,
                                     State state,
                                     int n_steps,
                                     int snapshot_every);

}  // namespace galaxy

#endif
