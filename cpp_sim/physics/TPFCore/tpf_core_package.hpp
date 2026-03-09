#ifndef GALAXY_PHYSICS_TPFCORE_PACKAGE_HPP
#define GALAXY_PHYSICS_TPFCORE_PACKAGE_HPP

/**
 * TPFCore: Honest primitive TPF structure package.
 *
 * Implements Xi^mu, Theta_{mu nu}, I = Theta_mn Theta^mn - lambda Theta^2.
 * Lambda fixed at 1/4 in 4D.
 *
 * NOT the removed weak-field Newtonian-like package.
 * NOT a full nonlinear/dynamic TPF solver.
 * Inspection-first; dynamics require provisional readout (not yet implemented).
 */

#include "../physics_package.hpp"
#include "source_ansatz.hpp"
#include <string>
#include <vector>

namespace galaxy {

class TPFCorePackage : public PhysicsPackage {
 public:
  TPFCorePackage();

  const char* name() const override { return "TPFCore"; }

  void init_from_config(const Config& config) override;

  /** TPFCore does NOT implement acceleration readout unless provisional (not yet impl). Always throws. */
  void compute_accelerations(const State& state,
                            double bh_mass,
                            double softening,
                            bool star_star,
                            std::vector<double>& ax,
                            std::vector<double>& ay) const override;

  double compute_potential_energy(const State&, double, double, bool) const override { return 0.0; }

  bool provisional_readout_enabled() const { return provisional_readout_; }
  bool provisional_source_ansatz_in_use() const { return true; }  // source_ansatz is always provisional

  /** Run single-source inspection: one source at origin, probe along +x. */
  void run_single_source_inspect(const Config& config, const std::string& output_dir);

  /** Run symmetric-pair inspection: sources at (+d,0) and (-d,0), probe along axes. */
  void run_symmetric_pair_inspect(const Config& config, const std::string& output_dir);

  /**
   * Run c-sweep optimization: numerically fit c against field-equation residual.
   * Exploratory ansatz-tuning tool. Fitted c is NOT a final paper-derived constant.
   */
  void run_single_source_optimize_c(const Config& config, const std::string& output_dir);

 private:
  bool provisional_readout_;
};

}  // namespace galaxy

#endif
