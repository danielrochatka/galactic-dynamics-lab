#ifndef GALAXY_SIMULATION_HPP
#define GALAXY_SIMULATION_HPP

#include "config.hpp"
#include "types.hpp"
#include <functional>
#include <vector>

namespace galaxy {

// Progress callback: (step, n_steps, sim_time). Called every progress_interval steps when set.
using ProgressCallback = std::function<void(int step, int n_steps, double sim_time)>;

// Run simulation and collect snapshots (step 0, then every snapshot_every).
// When progress_interval > 0 and progress_callback is callable, invokes it every progress_interval steps.
std::vector<Snapshot> run_simulation(const Config& config,
                                     State state,
                                     int n_steps,
                                     int snapshot_every,
                                     ProgressCallback progress_callback = nullptr,
                                     int progress_interval = 0);

}  // namespace galaxy

#endif
