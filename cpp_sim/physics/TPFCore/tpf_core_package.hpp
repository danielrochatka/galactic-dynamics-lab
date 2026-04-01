#ifndef GALAXY_PHYSICS_TPFCORE_PACKAGE_HPP
#define GALAXY_PHYSICS_TPFCORE_PACKAGE_HPP

/**
 * TPFCore: Honest primitive TPF structure package.
 *
 * Implements Xi^mu, Theta_{mu nu}, I = Theta_mn Theta^mn - lambda Theta^2.
 * Lambda = 1/4 in 4D (fixed theory; not tunable).
 *
 * Parameter roles: fixed theory (lambda); numerical regularization (source eps);
 * provisional readout knobs (mode/scale/theta_tt/theta_tr); VDSG coupling (exploratory SI path).
 *
 * compute_accelerations requires tpfcore_enable_provisional_readout (gate to this API).
 * When tpf_vdsg_coupling != 0, ax, ay are filled by VDSG (provisional readout closures do not supply them).
 * When VDSG is off, ax, ay come from provisional_readout closures for tpfcore_readout_mode.
 * Active branch identity: run_info / render_manifest (active_dynamics_branch, acceleration_code_path).
 */

#include "../../accel_pipeline_stats.hpp"
#include "../../types.hpp"
#include "../physics_package.hpp"
#include "derived_tpf_radial.hpp"
#include "source_ansatz.hpp"
#include <string>
#include <vector>

namespace galaxy {

class TPFCorePackage : public PhysicsPackage {
 public:
  TPFCorePackage();

  const char* name() const override { return "TPFCore"; }

  void init_from_config(const Config& config) override;

  /**
   * Particle accelerations. Requires provisional readout enabled (else throws).
   * If tpf_vdsg_coupling != 0: ax, ay from VDSG SI centripetal path.
   * Else: ax, ay from compute_provisional_readout_acceleration for configured mode.
   */
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

  /** Run weak-field calibration: compare TPF provisional radial acceleration to Newtonian benchmark. */
  void run_weak_field_calibration(const Config& config, const std::string& output_dir);

  /** Run symmetric-pair inspection: sources at (+d,0) and (-d,0), probe along axes. */
  void run_symmetric_pair_inspect(const Config& config, const std::string& output_dir);

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

  /** Write tpf_closure_diagnostics (tr_coherence_readout, single-body only): closure-term decomposition. */
  void write_closure_diagnostics(const std::vector<Snapshot>& snapshots,
                                 const Config& config,
                                 const std::string& output_dir) const;

  /** Per-snapshot pipeline metrics CSV (readout baseline, VDSG modifier, shunt); optional end-of-run pass. */
  void write_accel_pipeline_diagnostics(const std::vector<Snapshot>& snapshots,
                                          const Config& config,
                                          const std::string& output_dir) const;

  /** Last integrator-step pipeline stats (updated every compute_accelerations). */
  const AccelPipelineStats& last_accel_pipeline_stats() const { return last_pipeline_; }

  /** Live orbit force audit for bh_orbit_validation (Newtonian vs TPF for the actual evolving state). */
  void write_live_orbit_force_audit(const std::vector<Snapshot>& snapshots,
                                    const Config& config,
                                    const std::string& output_dir) const;

  /** Exact step-0 orbit audit: raw numbers for bh_orbit_validation initial state only. */
  void write_step0_orbit_audit(const std::vector<Snapshot>& snapshots,
                               const Config& config,
                               const std::string& output_dir) const;

  /** Trajectory metrics + class from snapshots (single-body only). For sweep harness. */
  struct TrajectorySummary {
    bool valid = false;
    double r_initial = 0.0, r_final = 0.0, r_min = 0.0, r_max = 0.0;
    double radial_drift = 0.0, revolutions = 0.0;
    std::string trajectory_class;
  };
  TrajectorySummary compute_trajectory_summary(const std::vector<Snapshot>& snapshots) const;

  /** Regime stats from snapshots (mean/max theta norm, regime fractions). For sweep harness. */
  struct RegimeSummary {
    bool valid = false;
    double mean_theta_norm = 0.0, max_theta_norm = 0.0, min_theta_norm = 0.0;
    size_t n_samples = 0;
    double frac_low = 0.0, frac_transitional = 0.0, frac_high = 0.0;
  };
  RegimeSummary compute_regime_summary(const std::vector<Snapshot>& snapshots,
                                       const Config& config,
                                       const std::string& output_dir) const;

 private:
  bool provisional_readout_;
  std::string readout_mode_;
  double readout_scale_;
  double theta_tt_scale_;
  double theta_tr_scale_;
  double source_softening_;
  double vdsg_coupling_;
  /** Resolved M_ref (kg): explicit tpf_vdsg_mass_baseline_kg or star_mass when baseline key <= 0. */
  double vdsg_mass_baseline_resolved_kg_;
  double simulation_dt_;
  bool shunt_enable_;
  double shunt_fraction_;
  bool pipeline_diagnostics_csv_;
  mutable AccelPipelineStats last_pipeline_;
  tpfcore::DerivedTpfPoissonConfig derived_poisson_cfg_;

  void eval_accel_pipeline(const State& state,
                           double bh_mass,
                           double softening,
                           bool star_star,
                           std::vector<double>& ax,
                           std::vector<double>& ay,
                           AccelPipelineStats* stats_out) const;
};

/** Test-only: reset before compute_accelerations; counts per-particle caps in last apply_global_accel_magnitude_shunt. */
void tpf_test_reset_global_accel_shunt_events();
unsigned tpf_test_global_accel_shunt_events();

}  // namespace galaxy

#endif
