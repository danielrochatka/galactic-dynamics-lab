#include "simulation.hpp"
#include "integrator.hpp"
#include <algorithm>

namespace galaxy {

std::vector<Snapshot> run_simulation(const Config& config,
                                     State state,
                                     const PhysicsPackage* physics,
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

  const bool cooling_on = physics->cooling_active(config);
  const bool xi_kernel_deformation_active =
      (config.physics_package == "TPFCore" && config.tpf_dynamics_mode == "tpf_xi_theta_v1" &&
       config.tpf_4d_xi_kernel_mode != "off" && config.tpf_4d_xi_kernel_coupling != 0.0);
  const int cooling_steps = cooling_on ? physics->cooling_steps(config, n_steps) : 0;

  for (int step = 1; step <= n_steps; ++step) {
    if (xi_kernel_deformation_active) {
      semi_implicit_euler_step(state, physics, bh_mass, softening, star_star, dt, ax, ay);
    } else {
      velocity_verlet_step(state, physics, bh_mass, softening, star_star, dt, ax, ay);
    }

    if (cooling_on) {
      physics->apply_cooling_step(state, config, step, cooling_steps);
    }

    /* Suppress snapshot collection during cooling (saves memory and disk when snapshots are written). */
    const bool skip_snapshot_for_cooling =
        cooling_on && physics->suppress_snapshot_for_cooling(config, step, cooling_steps);
    if (!skip_snapshot_for_cooling && step % snapshot_every == 0) {
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
