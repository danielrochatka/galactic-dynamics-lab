#ifndef GALAXY_SOFTENING_POLICY_HPP
#define GALAXY_SOFTENING_POLICY_HPP
#include "config.hpp"
#include "types.hpp"
#include <string>
namespace galaxy {
struct ResolvedSoftening {
  double effective_softening = 0.0;
  std::string mode;
  std::string profile;
  double mean_separation = 0.0;
  double radius_inner_used = 0.0;
  double radius_outer_used = 0.0;
  int dimension = 0;
  double factor = 0.0;
  std::string source;
  bool max_capped = false;
  bool min_floored = false;
  std::string max_cap_source;
};
ResolvedSoftening resolve_softening(const Config& cfg, const State& state);

double plummer_softening_scale(double r_sq, double softening);
void apply_plummer_softening(double dx,
                             double dy,
                             double softening,
                             double& ax,
                             double& ay);
}
#endif
