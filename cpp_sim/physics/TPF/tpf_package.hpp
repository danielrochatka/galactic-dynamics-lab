#ifndef GALAXY_PHYSICS_TPF_PACKAGE_HPP
#define GALAXY_PHYSICS_TPF_PACKAGE_HPP

/**
 * TPF weak-field correspondence package.
 *
 * Implements ONLY the paper's linearized weak-field / quasi-static sector.
 * NOT a full variational TPF solver. See README.md.
 */

#include "../physics_package.hpp"

namespace galaxy {

class TPFPackage : public PhysicsPackage {
 public:
  TPFPackage() : alpha_(4.0 * 3.14159265358979323846), softening_override_(0.0), match_newtonian_(true) {}

  const char* name() const override { return "TPF"; }

  void init_from_config(const Config& config) override;

  void compute_accelerations(const State& state,
                            double bh_mass,
                            double softening,
                            bool star_star,
                            std::vector<double>& ax,
                            std::vector<double>& ay) const override;

  double compute_potential_energy(const State& state,
                                  double bh_mass,
                                  double softening,
                                  bool star_star = true) const override;

  /** Resolved alpha used for this run (e.g. for run_info). */
  double resolved_alpha() const { return alpha_; }
  /** Whether Newtonian-scale normalization was applied. */
  bool match_newtonian_scale() const { return match_newtonian_; }

 private:
  double alpha_;
  double softening_override_;  // > 0 means use this instead of global
  bool match_newtonian_;
};

}  // namespace galaxy

#endif
