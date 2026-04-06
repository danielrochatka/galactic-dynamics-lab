#ifndef GALAXY_TYPES_HPP
#define GALAXY_TYPES_HPP

#include <vector>

namespace galaxy {

// Struct-of-arrays for cache-friendly layout; easy to parallelize over i.
struct State {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> mass;

  int n() const { return static_cast<int>(x.size()); }
  void resize(int n) {
    x.resize(n);
    y.resize(n);
    vx.resize(n);
    vy.resize(n);
    mass.resize(n);
  }
};

// One snapshot: step, time, and state copy (for output).
struct Snapshot {
  int step = 0;
  double time = 0.0;
  State state;
};

}  // namespace galaxy

#endif
