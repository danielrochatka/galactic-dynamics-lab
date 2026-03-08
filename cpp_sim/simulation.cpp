#include "simulation.hpp"
#include "integrator.hpp"
#include <algorithm>

namespace galaxy {

std::vector<Snapshot> run_simulation(const Config& config,
                                     State state,
                                     int n_steps,
                                     int snapshot_every,
                                     ProgressCallback progress_callback,
                                     int progress_interval) {
  std::vector<Snapshot> snapshots;
  Snapshot initial;
  initial.step = 0;
  initial.time = 0.0;
  initial.state = state;
  snapshots.push_back(std::move(initial));

  std::vector<double> ax, ay;
  const double dt = config.dt;
  const double bh_mass = config.bh_mass;
  const double softening = config.softening;
  const bool star_star = config.enable_star_star_gravity;
  const bool use_progress = (progress_interval > 0 && progress_callback);

  for (int step = 1; step <= n_steps; ++step) {
    velocity_verlet_step(state, bh_mass, softening, star_star, dt, ax, ay);

    if (step % snapshot_every == 0) {
      Snapshot snap;
      snap.step = step;
      snap.time = step * dt;
      snap.state = state;
      snapshots.push_back(std::move(snap));
    }

    if (use_progress && (step % progress_interval == 0 || step == n_steps)) {
      progress_callback(step, n_steps, step * dt);
    }
  }

  return snapshots;
}

}  // namespace galaxy
