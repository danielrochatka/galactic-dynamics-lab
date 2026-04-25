#ifndef GALAXY_PHYSICS_TPFCORE_TPF_4D_STATIC_FIELD_HPP
#define GALAXY_PHYSICS_TPFCORE_TPF_4D_STATIC_FIELD_HPP

#include <vector>

#include "physics/TPFCore/tpf_4d_field.hpp"

namespace galaxy {
namespace tpfcore {

struct StaticSourcePoint4D {
  double mass;
  double x;
  double y;
  double z;
};

struct StaticProbePoint4D {
  double x;
  double y;
  double z;
};

struct Field4DAtPoint {
  Xi4 xi;
  Theta4 theta;
  double theta_trace_4d;
  double invariant_I_4d;
};

Field4DAtPoint evaluate_static_source_field_4d(const StaticSourcePoint4D& source,
                                               const StaticProbePoint4D& probe,
                                               double eps);

Field4DAtPoint evaluate_static_sources_field_4d(const std::vector<StaticSourcePoint4D>& sources,
                                                const StaticProbePoint4D& probe,
                                                double eps);

}  // namespace tpfcore
}  // namespace galaxy

#endif
