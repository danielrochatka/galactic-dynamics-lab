#include "integrator.hpp"
#include "physics.hpp"
#include <cstddef>

namespace galaxy {

void velocity_verlet_step(State& state,
                         double bh_mass,
                         double softening,
                         bool star_star,
                         double dt,
                         std::vector<double>& ax,
                         std::vector<double>& ay) {
  const int n = state.n();
  if (ax.size() != static_cast<size_t>(n)) ax.resize(n);
  if (ay.size() != static_cast<size_t>(n)) ay.resize(n);

  compute_accelerations(state, bh_mass, softening, star_star, ax, ay);

  const double dt2 = dt * dt;
  for (int i = 0; i < n; ++i) {
    state.x[i] += state.vx[i] * dt + 0.5 * ax[i] * dt2;
    state.y[i] += state.vy[i] * dt + 0.5 * ay[i] * dt2;
  }

  std::vector<double> ax_new(n), ay_new(n);
  compute_accelerations(state, bh_mass, softening, star_star, ax_new, ay_new);

  for (int i = 0; i < n; ++i) {
    state.vx[i] += 0.5 * (ax[i] + ax_new[i]) * dt;
    state.vy[i] += 0.5 * (ay[i] + ay_new[i]) * dt;
  }
}

}  // namespace galaxy
