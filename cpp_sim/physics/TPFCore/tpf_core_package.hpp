#ifndef GALAXY_PHYSICS_TPFCORE_PACKAGE_HPP
#define GALAXY_PHYSICS_TPFCORE_PACKAGE_HPP

/**
 * TPFCore: Honest primitive TPF structure package.
 *
 * Implements Xi^mu, Theta_{mu nu}, I = Theta_mn Theta^mn - lambda Theta^2.
 * Lambda = 1/4 in 4D (fixed theory; not tunable).
 *
 * Parameter roles: fixed theory (lambda); numerical regularization (source eps);
 * exploratory ansatz (isotropic correction c, NOT a fundamental constant);
 * provisional experimental (readout mode/scale/theta_tt/theta_tr).
 *
 * When tpfcore_enable_provisional_readout=true, a PROVISIONAL motion/readout layer
 * maps Theta into acceleration. EXPLORATORY—not the full derived TPF dynamics.
 * Inspection-first; dynamics require provisional readout to be enabled.
 */

#include "../../types.hpp"
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

  /** When provisional readout enabled: tensor-driven acceleration. Otherwise throws. */
  void compute_accelerations(const State& state,
                            double bh_mass,
                            double softening,
                            bool star_star,
                            std::vector<double>& ax,
                            std::vector<double>& ay) const override;

  double compute_potential_energy(const State&, double, double, bool) const override { return 0.0; }

  bool provisional_readout_enabled() const { return provisional_readout_; }
  const std::string& readout_mode() const { return readout_mode_; }
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

  /** Write tpf_readout_debug.csv for dynamical runs when tpfcore_dump_readout_debug. */
  void write_readout_debug(const std::vector<Snapshot>& snapshots,
                           const Config& config,
                           const std::string& output_dir) const;

  /** Write tpf_regime_diagnostics.txt for dynamical runs (field intensity, I, regime distribution). */
  void write_regime_diagnostics(const std::vector<Snapshot>& snapshots,
                                const Config& config,
                                const std::string& output_dir) const;

  /** Write tpf_trajectory_diagnostics.txt for dynamical runs (geometry/time-series classification). */
  void write_trajectory_diagnostics(const std::vector<Snapshot>& snapshots,
                                    const Config& config,
                                    const std::string& output_dir) const;

 private:
  bool provisional_readout_;
  std::string readout_mode_;
  double readout_scale_;
  double theta_tt_scale_;
  double theta_tr_scale_;
  double isotropic_c_;
  double source_softening_;
};

}  // namespace galaxy

#endif
