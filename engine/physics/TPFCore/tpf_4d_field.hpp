#ifndef GALAXY_PHYSICS_TPFCORE_TPF_4D_FIELD_HPP
#define GALAXY_PHYSICS_TPFCORE_TPF_4D_FIELD_HPP

#include <cstddef>

#include "physics/TPFCore/source_ansatz.hpp"

namespace galaxy {
namespace tpfcore {

struct Xi4 {
  double t, x, y, z;

  double& component(std::size_t mu);
  double component(std::size_t mu) const;
};

struct Theta4 {
  double tt, tx, ty, tz;
  double xt, xx, xy, xz;
  double yt, yx, yy, yz;
  double zt, zx, zy;
  double zz;

  double& component(std::size_t mu, std::size_t nu);
  double component(std::size_t mu, std::size_t nu) const;
};

double metric_signature_sign(std::size_t mu);
double trace_contraction_4d(const Theta4& theta);
double theta_double_contraction_4d(const Theta4& theta);
double invariant_I_4d(const Theta4& theta);

}  // namespace tpfcore
}  // namespace galaxy

#endif
