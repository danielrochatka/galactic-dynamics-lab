#ifndef GALAXY_PHYSICS_TPFCORE_SOURCE_ITERATION_HPP
#define GALAXY_PHYSICS_TPFCORE_SOURCE_ITERATION_HPP

#include "../../types.hpp"

namespace galaxy {
namespace tpfcore {

struct GravitationalSource {
  double x = 0.0;
  double y = 0.0;
  double mass = 0.0;
  bool is_black_hole = false;
  int star_index = -1;
};

template <typename Fn>
inline void for_each_gravitational_source(const State& state,
                                          int target_index,
                                          double bh_mass,
                                          bool star_star,
                                          Fn&& fn) {
  if (bh_mass > 0.0) {
    GravitationalSource source;
    source.mass = bh_mass;
    source.is_black_hole = true;
    source.star_index = -1;
    fn(source);
  }

  if (!star_star) {
    return;
  }

  const int n = state.n();
  for (int j = 0; j < n; ++j) {
    if (j == target_index) {
      continue;
    }
    const double mj = state.mass[j];
    if (mj <= 0.0) {
      continue;
    }
    GravitationalSource source;
    source.x = state.x[j];
    source.y = state.y[j];
    source.mass = mj;
    source.is_black_hole = false;
    source.star_index = j;
    fn(source);
  }
}

}  // namespace tpfcore
}  // namespace galaxy

#endif
