#ifndef GALAXY_PHYSICS_TPFCORE_PACKAGE_HPP
#define GALAXY_PHYSICS_TPFCORE_PACKAGE_HPP

/**
 * TPFCore: Honest primitive TPF structure package.
 *
 * Active runtime route on this branch: `tpf_xi_theta_v1`.
 * - Xi_total is accumulated from per-source Xi contributions.
 * - Theta is the full unsymmetrized spatial Jacobian Theta_ij = d_j Xi_i,total.
 * - Motion update uses Xi_total only; Theta is diagnostic-only.
 * Isolated 4D math/static residual benchmark helpers live in tpf_4d_field.* / tpf_4d_static_* and are not wired into dynamics.
 */

#include "../../types.hpp"
#include "../physics_package.hpp"
#include "derived_tpf_radial.hpp"
#include "source_ansatz.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace galaxy {

class TPFCorePackage : public PhysicsPackage {
 public:
  TPFCorePackage();

  const char* name() const override { return "TPFCore"; }

  void init_from_config(const Config& config) override;
  void sync_config_from_package_options(Config& config, const std::string& changed_key) const override;
  PackageModeInfo resolve_mode_token(const std::string& mode_token) const override;
  bool supports_utility_mode(SimulationMode mode) const override;
  bool run_utility_mode(const Config& config, const std::string& output_dir) override;
  bool cooling_active(const Config& config) const override;
  int cooling_steps(const Config& config, int n_steps) const override;
  void apply_cooling_step(State& state, const Config& config, int step, int cooling_steps) const override;
  bool suppress_snapshot_for_cooling(const Config& config, int step, int cooling_steps) const override;
  bool write_post_run_diagnostics(const std::vector<Snapshot>& snapshots,
                                  const Config& config,
                                  const std::string& output_dir) override;
  std::vector<PackageMetadataEntry> run_info_metadata(const Config& config,
                                                      const std::string& package_defaults_path) const override;
  std::vector<PackageMetadataEntry> run_info_supplement_metadata(const Config& config) const override;
  std::vector<PackageMetadataEntry> run_info_runtime_metadata(const Config& config,
                                                              const RunInfoContext& context) const override;
  void write_run_info_section(std::ostream& out,
                              const Config& config,
                              const RunInfoContext& context,
                              PackageRunInfoSection section) const override;
  std::vector<PackageMetadataEntry> render_metadata(const Config& config) const override;
  std::string render_active_dynamics_branch(const Config& config) const override;
  std::string render_active_metrics_branch(const Config& config) const override;
  std::string render_acceleration_code_path(const Config& config) const override;
  std::vector<PackageMetadataEntry> resolved_runtime_metadata(const Config& config,
                                                              SimulationMode mode) const override;

  /** Particle accelerations for `tpf_xi_theta_v1`: Xi_total-driven motion update; Theta is diagnostic-only. */
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

  /** Run correspondence calibration: compare provisional readout radial acceleration to Newtonian benchmark. */
  void run_weak_field_calibration(const Config& config, const std::string& output_dir);

  /** Run symmetric-pair inspection: sources at (+d,0) and (-d,0), probe along axes. */
  void run_symmetric_pair_inspect(const Config& config, const std::string& output_dir);
  /** Run source-first field benchmark and dump plane probe grid CSV. */
  void run_source_field_benchmark(const Config& config, const std::string& output_dir);
  /** Run static 4D residual benchmark harness and dump summary + view-plane slice artifacts. */
  void run_4d_static_residual_benchmark(const Config& config, const std::string& output_dir);
  /** Run static 4D field -> principal tensor -> probe-motion readout benchmark harness. */
  void run_4d_static_motion_readout_benchmark(const Config& config, const std::string& output_dir);
  /** Run dynamic probe-motion benchmark using Xi-direct acceleration from fixed-source static 4D field evaluation. */
  void run_4d_xi_motion_probe_benchmark(const Config& config, const std::string& output_dir);

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

  struct XiRuntimeCounters {
    std::uint64_t theta_evaluations = 0;
    std::uint64_t invariant_I_evaluations = 0;
    std::uint64_t xi_last_call_pair_evaluations = 0;
    std::uint64_t xi_total_pair_evaluations = 0;
  };
  XiRuntimeCounters xi_runtime_counters() const { return xi_runtime_counters_; }
  struct XiThetaV1Sample {
    double xi_x = 0.0;
    double xi_y = 0.0;
    double xi_z = 0.0;
    double theta[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  };
  XiThetaV1Sample last_xi_theta_v1_sample() const { return last_xi_theta_v1_sample_; }

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
  std::string tpf_dynamics_mode_;
  bool provisional_readout_;
  std::string readout_mode_;
  double readout_scale_;
  double theta_tt_scale_;
  double theta_tr_scale_;
  double source_softening_;
  /** Internal paper-baseline coupling used by direct_tpf tensor principal-part route. */
  double kappa_;
  double weak_field_correspondence_alpha_si_;
  double vdsg_coupling_;
  /** Resolved M_ref (kg): explicit tpf_vdsg_mass_baseline_kg or star_mass when baseline key <= 0. */
  double vdsg_mass_baseline_resolved_kg_;
  std::string vdsg_mode_;
  double vdsg_mass_gate_m0_kg_;
  double vdsg_mass_gate_alpha_;
  double vdsg_x_clamp_;
  bool vdsg_weak_field_gate_enable_;
  double vdsg_weak_field_a0_;
  double vdsg_weak_field_power_;
  double vdsg_bounded_amplitude_;
  double simulation_dt_;
  double cooling_fraction_;
  bool shunt_enable_;
  double shunt_fraction_;
  /** Legacy derived-radial closure ledger configuration (still fed from flat config keys). */
  tpfcore::DerivedTpfPoissonConfig derived_poisson_cfg_;

  void compute_xi_kernel_deformed_accelerations(const State& state,
                                                double bh_mass,
                                                double softening,
                                                bool star_star,
                                                std::vector<double>& ax,
                                                std::vector<double>& ay) const;
  void validate_xi_kernel_runtime_config() const;
  bool xi_kernel_deformation_active() const;

  double xi_motion_readout_scale_;
  std::string xi_kernel_mode_;
  double xi_kernel_coupling_;
  double xi_kernel_beta_power_;
  std::string xi_kernel_factor_mode_;
  double xi_kernel_metric_min_;
  double xi_kernel_metric_max_;
  std::string xi_temporal_mode_;
  double xi_temporal_coupling_;
  double xi_source_speed_x_;
  double xi_source_speed_y_;
  double xi_source_speed_z_;
  mutable XiRuntimeCounters xi_runtime_counters_;
  mutable XiThetaV1Sample last_xi_theta_v1_sample_;
};

/** Test-only: reset before compute_accelerations; counts per-particle caps in last apply_global_accel_magnitude_shunt. */
void tpf_test_reset_global_accel_shunt_events();
unsigned tpf_test_global_accel_shunt_events();

}  // namespace galaxy

#endif
