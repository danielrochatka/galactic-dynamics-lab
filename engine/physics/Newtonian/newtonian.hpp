#ifndef GALAXY_PHYSICS_NEWTONIAN_HPP
#define GALAXY_PHYSICS_NEWTONIAN_HPP

#include "../physics_package.hpp"
#include "../../softening_audit.hpp"

namespace galaxy {

/** Newtonian gravity package: BH at origin + optional pairwise star-star with softening. */
class NewtonianPackage : public PhysicsPackage {
 public:
  const char* name() const override { return "Newtonian"; }
  void init_from_config(const Config& config) override { softening_audit_enable_ = config.softening_audit_enable; }

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
  const SofteningAuditStats& softening_audit_stats() const { return softening_audit_stats_; }
 private:
  bool softening_audit_enable_ = false;
  mutable SofteningAuditStats softening_audit_stats_;
};

}  // namespace galaxy

#endif
