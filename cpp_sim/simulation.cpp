#include "simulation.hpp"
#include "integrator.hpp"
#include <algorithm>
#include <cmath>

namespace galaxy {

namespace {

/** 1% radial velocity damping per step (TPF artificial cooling). */
constexpr double kTpfRadialCoolingDampingFactor = 0.99;

void apply_radial_cooling_damping(State& state, double damping_factor) {
  const int n = state.n();
  for (int i = 0; i < n; ++i) {
    double x = state.x[i];
    double y = state.y[i];
    double r = std::sqrt(x * x + y * y);
    if (r > 0.0) {
      double vx = state.vx[i];
      double vy = state.vy[i];
      double v_rad = (vx * x + vy * y) / r;
      double damped_v_rad = v_rad * damping_factor;
      double v_rad_diff = v_rad - damped_v_rad;
      state.vx[i] -= v_rad_diff * (x / r);
      state.vy[i] -= v_rad_diff * (y / r);
    }
  }
}

}  // namespace

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

  const bool tpf_cooling_on =
      (config.physics_package == "TPFCore" && config.tpf_cooling_fraction > 0.0);
  const int cooling_steps = tpf_cooling_on
      ? std::min(n_steps, std::max(0, static_cast<int>(n_steps * config.tpf_cooling_fraction)))
      : 0;

  for (int step = 1; step <= n_steps; ++step) {
    velocity_verlet_step(state, physics, bh_mass, softening, star_star, dt, ax, ay);

    if (tpf_cooling_on && step < cooling_steps) {
      apply_radial_cooling_damping(state, kTpfRadialCoolingDampingFactor);
    }

    /* Suppress snapshot collection during cooling (saves memory and disk when snapshots are written). */
    const bool skip_snapshot_for_cooling = tpf_cooling_on && step < cooling_steps;
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
