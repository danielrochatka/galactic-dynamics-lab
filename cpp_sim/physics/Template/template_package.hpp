#ifndef GALAXY_PHYSICS_TEMPLATE_PACKAGE_HPP
#define GALAXY_PHYSICS_TEMPLATE_PACKAGE_HPP

/**
 * Template / stub physics package.
 * Copy this folder, rename to your package (e.g. MyCustomPhysics), and implement
 * at least name() and compute_accelerations(). See this directory's README.md.
 */

#include "../physics_package.hpp"

namespace galaxy {

class TemplatePackage : public PhysicsPackage {
 public:
  const char* name() const override { return "Template"; }

  void compute_accelerations(const State& state,
                             double bh_mass,
                             double softening,
                             bool star_star,
                             std::vector<double>& ax,
                             std::vector<double>& ay) const override {
    (void)bh_mass;
    (void)softening;
    (void)star_star;
    // TODO: implement your acceleration law here
    const int n = state.n();
    ax.assign(n, 0.0);
    ay.assign(n, 0.0);
  }
};

}  // namespace galaxy

#endif
