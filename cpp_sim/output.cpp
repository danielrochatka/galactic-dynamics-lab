#include "output.hpp"
#include "accel_pipeline_stats.hpp"
#include "galaxy_init.hpp"
#include "git_provenance.hpp"
#include "render_audit.hpp"
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace galaxy {

void write_run_info(const std::string& output_dir,
                    const Config& config,
                    int n_steps_done,
                    int n_snapshots,
                    int n_particles,
                    const std::string& run_config_path,
                    const std::string& package_defaults_path,
                    const GalaxyInitAudit* galaxy_init_audit,
                    const CoolingAuditInfo* cooling_audit,
                    const AccelPipelineStats* tpf_pipeline) {
  std::ostringstream path;
  path << output_dir << "/run_info.txt";
  std::ofstream f(path.str());
  if (!f) return;

  const GitProvenance gp = resolve_git_provenance();

  int n_star = (n_particles >= 0) ? n_particles : config.n_stars;

  /* Resolved config banner at top so it is obvious what was used */
  f << "=== Resolved config (built-in < package defaults < run config < CLI) ===\n";
  f << "run_config\t" << (run_config_path.empty() ? "(none)" : run_config_path) << "\n";
  f << "package_defaults\t" << (package_defaults_path.empty() ? "(none)" : package_defaults_path) << "\n";
  f << "physics_package\t" << config.physics_package << "\n";
  f << "physics_package_compare\t" << config.physics_package_compare << "\n";
  f << "simulation_mode\t" << mode_to_string(config.simulation_mode) << "\n";
  f << "tpf_analysis_mode\t" << config.tpf_analysis_mode << "\n";
  f << "artifact_naming_convention\t<mode>__<physics>__<scope>__<quantity>__<stage>.<ext>\n";
  f << "artifact_label_scope_primary\tprimary=main interpretation outputs (mode-dependent)\n";
  f << "artifact_label_scope_secondary\tsecondary=lab/origin radial diagnostics where applicable\n";
  f << "render_overlay_mode\t" << config.render_overlay_mode << "\n";
  f << "display_distance_unit\t" << config.display_distance_unit << "\n";
  f << "display_time_unit\t" << config.display_time_unit << "\n";
  f << "display_velocity_unit\t" << config.display_velocity_unit << "\n";
  f << "display_units_in_overlay\t" << (config.display_units_in_overlay ? 1 : 0) << "\n";
  f << "display_show_unit_reference\t" << (config.display_show_unit_reference ? 1 : 0) << "\n";
  f << "active_dynamics_branch\t" << compute_active_dynamics_branch(config) << "\n";
  f << "active_metrics_branch\t" << compute_active_metrics_branch(config) << "\n";
  f << "acceleration_code_path\t" << compute_acceleration_code_path(config) << "\n";
  f << "output_dir\t" << output_dir << "\n";
  f << "n_stars\t" << n_star << "\n";
  f << "bh_mass\t" << config.bh_mass << "\n";
  if (config.simulation_mode == galaxy::SimulationMode::galaxy) {
    f << "galaxy_init_template\t" << config.galaxy_init_template << "\n";
    f << "galaxy_init_seed\t" << config.galaxy_init_seed << "\n";
    f << "galaxy_init_master_chaos\t" << config.galaxy_init_master_chaos << "\n";
    f << "galaxy_init_position_noise\t" << config.galaxy_init_position_noise << "\n";
    f << "galaxy_init_velocity_angle_noise\t" << config.galaxy_init_velocity_angle_noise << "\n";
    f << "galaxy_init_velocity_magnitude_noise\t" << config.galaxy_init_velocity_magnitude_noise << "\n";
    f << "galaxy_init_clumpiness\t" << config.galaxy_init_clumpiness << "\n";
    f << "galaxy_init_num_clumps\t" << config.galaxy_init_num_clumps << "\n";
    f << "galaxy_init_clump_radius_fraction\t" << config.galaxy_init_clump_radius_fraction << "\n";
    f << "galaxy_init_m2_amplitude\t" << config.galaxy_init_m2_amplitude << "\n";
    f << "galaxy_init_m3_amplitude\t" << config.galaxy_init_m3_amplitude << "\n";
    f << "galaxy_init_bar_amplitude\t" << config.galaxy_init_bar_amplitude << "\n";
    f << "galaxy_init_bar_axis_ratio\t" << config.galaxy_init_bar_axis_ratio << "\n";
    f << "galaxy_init_spiral_amplitude\t" << config.galaxy_init_spiral_amplitude << "\n";
    f << "galaxy_init_spiral_winding\t" << config.galaxy_init_spiral_winding << "\n";
    f << "galaxy_init_spiral_phase\t" << config.galaxy_init_spiral_phase << "\n";
    f << "velocity_noise\t" << config.velocity_noise << "\n";
  }
  if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
    f << "\n=== v11 weak-field correspondence audit (manuscript v11; not integrator dynamics) ===\n";
    f << "v11_weak_field_correspondence_benchmark\t" << config.v11_weak_field_correspondence_benchmark << "\n";
    if (config.v11_weak_field_correspondence_benchmark == "earth_moon_line_of_centers") {
      f << "v11_audit_claim\tearth_moon_weak_field_line_of_centers_correspondence_benchmark_only_not_full_tpf_solver\n";
      f << "v11_extensions_note\tVDSG_off_DeltaC_omitted_no_direct_tpf_no_orbit_integrator_paper_extensions_off\n";
      f << "v11_phi_input\tweak_field_Poisson_correspondence_phi_Eq44_point_mass_line_exterior_only\n";
      f << "v11_eq9_audited\t0\n";
      f << "v11_eq9_note\tnot_applicable_earth_moon_benchmark_is_1D_phi_Eq44_to_acceleration_Eq45_vs_Newtonian_Eq46\n";
      f << "v11_earth_moon_benchmark_png_note\tauto-generated_when_plot_v11_earth_moon_line_benchmark_py_succeeds_requires_python3_numpy_matplotlib\n";
      f << "v11_earth_moon_benchmark_png_files\ttpf_v11_earth_moon_line_correspondence_benchmark_compare.png;"
            "tpf_v11_earth_moon_line_correspondence_benchmark_difference.png;"
            "tpf_v11_earth_moon_line_correspondence_benchmark_earth_zoom.png;"
            "tpf_v11_earth_moon_line_correspondence_benchmark_normalized_shape.png\n";
    } else {
      f << "v11_audit_claim\tstatic_weak_field_tensor_correspondence_on_positive_z_axis_only\n";
      f << "v11_eq9_audited\t1\n";
      f << "v11_eq9_note\tflat_static_axis_residual_Rj_equals_1minus_lambda_dj_Theta_see_CSV_eq9_columns\n";
      f << "v11_phi_input\tcorrespondence_benchmark_-GM_over_sqrt_z2_plus_eps2_softening_is_numerical_only\n";
      f << "v11_eq10_C_principal_scaling\tkappa_multiplies_entire_principal_bracket_including_minus_half_gI_DeltaC_omitted\n";
      f << "v11_eq10_weak_field_C00_identity\tC00_principal_SI_equals_kappa_times_I_over_2_when_Theta0mu_zero_g00_minus_1\n";
    }
    f << "v11_delta_c_computed\t0\n";
    f << "v11_delta_c_note\tDeltaC_mu_nu omitted per manuscript v11 (connection-variation terms deferred)\n";
    f << "=== End v11 weak-field audit header ===\n\n";
  }
  if (config.physics_package == "TPFCore") {
    if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
      f << "v11_audit_tpfcore_dynamics_note\tno particle integration; TPFCore acceleration API (legacy_readout / "
           "direct_tpf) not used in this mode\n";
      f << "tpf_dynamics_mode_configured_in_layered_config\t" << config.tpf_dynamics_mode << "\n";
      f << "tpf_dynamics_mode_effective_for_this_run\tnone_audit_only_no_integrator\n";
      f << "tpfcore_enable_provisional_readout_configured\t" << (config.tpfcore_enable_provisional_readout ? 1 : 0)
         << "\n";
      f << "tpfcore_enable_provisional_readout_effective_for_this_run\t0_unused_in_v11_audit_mode\n";
      f << "tpfcore_readout_mode_configured\t" << config.tpfcore_readout_mode << "\n";
      f << "v11_audit_readout_scalars_note\ttpfcore_readout_scale/kappa/etc. below are inherited layered-config "
           "artifacts; not used for particle accelerations in this mode\n";
    } else {
      f << "tpf_dynamics_mode\t" << config.tpf_dynamics_mode << "\n";
      f << "tpfcore_enable_provisional_readout\t" << (config.tpfcore_enable_provisional_readout ? 1 : 0) << "\n";
      f << "tpfcore_readout_mode\t" << config.tpfcore_readout_mode << "\n";
    }
    f << "tpfcore_readout_scale\t" << config.tpfcore_readout_scale << "\n";
    f << "tpfcore_theta_tt_scale\t" << config.tpfcore_theta_tt_scale << "\n";
    f << "tpfcore_theta_tr_scale\t" << config.tpfcore_theta_tr_scale << "\n";
    f << "tpf_kappa\t" << config.tpf_kappa << "\n";
    f << "tpf_vdsg_coupling\t" << config.tpf_vdsg_coupling << "\n";
    f << "tpf_vdsg_mass_baseline_kg\t" << config.tpf_vdsg_mass_baseline_kg << "\n";
    f << "tpf_poisson_bins\t" << config.tpf_poisson_bins << "\n";
    f << "tpf_poisson_max_radius\t" << config.tpf_poisson_max_radius << "\n";
    f << "tpf_cooling_fraction\t" << config.tpf_cooling_fraction << "\n";
    f << "tpf_global_accel_shunt_enable\t" << (config.tpf_global_accel_shunt_enable ? 1 : 0) << "\n";
    f << "tpf_global_accel_shunt_fraction\t" << config.tpf_global_accel_shunt_fraction << "\n";
    f << "tpf_accel_pipeline_diagnostics_csv\t" << (config.tpf_accel_pipeline_diagnostics_csv ? 1 : 0) << "\n";
  }
  f << "=== End resolved config ===\n\n";

  f << "=== Code provenance (git) ===\n";
  f << "git_commit_full\t" << gp.git_commit_full << "\n";
  f << "git_commit_short\t" << gp.git_commit_short << "\n";
  f << "git_branch\t" << gp.git_branch << "\n";
  f << "git_tag\t" << gp.git_tag << "\n";
  f << "git_dirty\t" << (gp.git_dirty ? 1 : 0) << "\n";
  f << "code_version_label\t" << gp.code_version_label << "\n";
  f << "=== End code provenance ===\n\n";

  if (config.physics_package == "TPFCore") {
    if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
      f << "=== TPFCore (v11 correspondence audit only: dynamics routing not operative) ===\n";
      f << "note\tlayered config may still set tpf_dynamics_mode, tpfcore_enable_provisional_readout, readout_mode; "
           "they are inherited only and not used for particle accelerations in this audit-only mode\n";
      f << "fixed_theory\tlambda=1/4 (LAMBDA_4D; fixed in code; not tunable)\n";
      f << "numerical_regularization\ttpfcore_source_softening, effective_source_softening (eps for Phi) — not used "
           "by all v11 benchmarks\n";
      f << "inspection\ttpfcore_probe_radius_min/max, probe_samples (sampling range for correspondence audit)\n";
      f << "=== End TPFCore v11 audit note ===\n\n";
    } else {
      f << "=== TPFCore parameter roles (theory vs regularization vs exploratory vs provisional) ===\n";
      f << "fixed_theory\tlambda=1/4 (LAMBDA_4D; fixed in code; not tunable)\n";
      f << "numerical_regularization\ttpfcore_source_softening, effective_source_softening (eps for Phi)\n";
      f << "dynamics_routing\ttpf_dynamics_mode (legacy_readout vs direct_tpf); legacy_readout uses "
           "tpfcore_enable_provisional_readout as gate; direct_tpf does not use that gate (stub throws until implemented)\n";
      f << "provisional_readout\ttpfcore_enable_provisional_readout (gate to legacy_readout accelerations), readout_mode "
           "(configured label; may differ from integrator ax,ay path when tpf_vdsg_coupling != 0 on legacy_readout), "
           "readout_scale, theta_tt_scale, theta_tr_scale, dump_readout_debug (experimental readout closures; "
           "diagnostics)\n";
      f << "inspection\ttpfcore_probe_radius_min/max, probe_samples, dump_theta_profile, dump_invariant_profile\n";
      f << "=== End TPFCore parameter roles ===\n\n";
    }
  }

  f << "dt\t" << config.dt << "\n";
  f << "n_steps\t" << n_steps_done << "\n";
  f << "snapshot_every\t" << config.snapshot_every << "\n";
  f << "softening\t" << config.softening << "\n";
  f << "star_mass\t" << config.star_mass << "\n";
  f << "bh_mass\t" << config.bh_mass << "\n";
  f << "enable_star_star_gravity\t" << (config.enable_star_star_gravity ? 1 : 0) << "\n";
  f << "galaxy_radius\t" << config.galaxy_radius << "\n";
  f << "plot_animation_dynamic_zoom\t" << (config.plot_animation_dynamic_zoom ? 1 : 0) << "\n";
  f << "plot_skip_initial_steps\t" << config.plot_skip_initial_steps << "\n";
  f << "plot_skip_initial_snapshots\t" << config.plot_skip_initial_snapshots << "\n";
  f << "diagnostic_cutoff_radius\t" << config.diagnostic_cutoff_radius << "\n";
  f << "render_overlay_mode\t" << config.render_overlay_mode << "\n";
  f << "display_distance_unit\t" << config.display_distance_unit << "\n";
  f << "display_time_unit\t" << config.display_time_unit << "\n";
  f << "display_velocity_unit\t" << config.display_velocity_unit << "\n";
  f << "display_units_in_overlay\t" << (config.display_units_in_overlay ? 1 : 0) << "\n";
  f << "display_show_unit_reference\t" << (config.display_show_unit_reference ? 1 : 0) << "\n";
  f << "active_dynamics_branch\t" << compute_active_dynamics_branch(config) << "\n";
  f << "active_metrics_branch\t" << compute_active_metrics_branch(config) << "\n";
  f << "acceleration_code_path\t" << compute_acceleration_code_path(config) << "\n";
  f << "total_simulated_time\t" << (n_steps_done * config.dt) << "\n";
  f << "number_of_snapshots\t" << n_snapshots << "\n";
  if (config.simulation_mode == SimulationMode::earth_moon_benchmark ||
      config.simulation_mode == SimulationMode::bh_orbit_validation ||
      config.simulation_mode == SimulationMode::two_body_orbit) {
    f << "\n=== Two-body run (postprocess diagnostics) ===\n";
    f << "diagnostics_primary_pair_frame\tplot_cpp_run.py: diagnostic_pair_*.png, diagnostic_com_radius.png, "
         "diagnostic_relative_angular_momentum_z.png, diagnostic_relative_energy.png (Newtonian equivalent), "
         "diagnostic_two_body_timeseries.csv, two_body_diagnostics_README.txt\n";
    f << "diagnostics_primary_pair_frame_mode_aware_alias\tplot_cpp_run.py also writes "
         "<mode>__<physics>__primary__<quantity>__<stage> aliases for these files\n";
    f << "diagnostics_secondary_lab_frame\tOrigin-radial plots (diagnostic_median_radius.png, etc.) are secondary; "
         "titles prefixed when regenerated for these modes\n";
    if (config.simulation_mode == SimulationMode::bh_orbit_validation) {
      f << "bh_orbit_validation_postprocess\tplot_cpp_run.py: same primary pair PNGs + footnotes; "
           "bh_orbit_trajectory_xy.png, bh_orbit_trajectory_xy_zoom.png, bh_orbit_separation_extrema.png "
           "(experimental; not paper correspondence; not direct_tpf; compare Newtonian vs TPFCore legacy_readout with VDSG off)\n";
    }
    f << "=== End two-body diagnostics note ===\n\n";
  }
  if (cooling_audit) {
    f << "cooling_active\t" << (cooling_audit->cooling_active ? 1 : 0) << "\n";
    f << "cooling_steps\t" << cooling_audit->cooling_steps << "\n";
    f << "cooling_end_step\t" << cooling_audit->cooling_end_step << "\n";
    f << "first_saved_snapshot_step\t" << cooling_audit->first_saved_snapshot_step << "\n";
    f << "first_saved_snapshot_time\t" << cooling_audit->first_saved_snapshot_time << "\n";
  }
  if (config.physics_package == "TPFCore" && tpf_pipeline && tpf_pipeline->valid) {
    f << "\n=== TPF acceleration pipeline (last integrator step) ===\n";
    f << "tpf_last_mean_baseline_accel_mag\t" << tpf_pipeline->mean_baseline_mag << "\n";
    f << "tpf_last_mean_vdsg_accel_mag\t" << tpf_pipeline->mean_vdsg_mag << "\n";
    f << "tpf_last_vdsg_over_baseline_ratio\t" << tpf_pipeline->vdsg_over_baseline_ratio << "\n";
    f << "tpf_last_shunt_events\t" << tpf_pipeline->shunt_events_last_step << "\n";
    f << "tpf_last_frac_capped\t" << tpf_pipeline->frac_capped_last_step << "\n";
    f << "tpf_last_global_accel_shunt_enabled\t" << (tpf_pipeline->shunt_enabled ? 1 : 0) << "\n";
    f << "tpf_last_global_accel_shunt_fraction\t" << tpf_pipeline->shunt_fraction << "\n";
    f << "=== End TPF acceleration pipeline ===\n";
  }
  f << "n_stars\t" << n_star << "\n";
  f << "simulation_mode\t" << static_cast<int>(config.simulation_mode) << "\n";
  f << "physics_package\t" << config.physics_package << "\n";
  f << "physics_package_compare\t" << config.physics_package_compare << "\n";
  if (config.physics_package == "TPFCore") {
    double src_eps = (config.tpfcore_source_softening > 0.0) ? config.tpfcore_source_softening : config.softening;
    if (config.simulation_mode != SimulationMode::tpf_v11_weak_field_correspondence) {
      f << "tpf_dynamics_mode\t" << config.tpf_dynamics_mode << "\n";
      f << "tpfcore_enable_provisional_readout\t" << (config.tpfcore_enable_provisional_readout ? 1 : 0) << "\n";
    } else {
      f << "v11_audit_repeat_note\ttpf_dynamics_mode / provisional readout: see resolved config section above "
           "(configured vs effective; not operative)\n";
    }
    f << "tpfcore_provisional_source_ansatz\t1\n";
    f << "tpfcore_source_softening\t" << src_eps << "\n";
    f << "tpfcore_probe_radius_min\t" << config.tpfcore_probe_radius_min << "\n";
    f << "tpfcore_probe_radius_max\t" << config.tpfcore_probe_radius_max << "\n";
    f << "tpfcore_probe_samples\t" << config.tpfcore_probe_samples << "\n";
    f << "tpfcore_residual_method\tanalytic\n";
    f << "tpfcore_residual_step\t" << config.tpfcore_residual_step << "\n";
    if (config.simulation_mode != SimulationMode::tpf_v11_weak_field_correspondence) {
      f << "tpfcore_readout_mode\t" << config.tpfcore_readout_mode << "\n";
    } else {
      f << "v11_audit_repeat_readout_mode\tsee tpfcore_readout_mode_configured in resolved config section above\n";
    }
    f << "tpfcore_readout_scale\t" << config.tpfcore_readout_scale << "\n";
    f << "tpfcore_readout_scale_note\tweak-field calibrated effective scale (K_eff); not proof of final TPF dynamics\n";
    f << "tpfcore_theta_tt_scale\t" << config.tpfcore_theta_tt_scale << "\n";
    f << "tpfcore_theta_tr_scale\t" << config.tpfcore_theta_tr_scale << "\n";
    f << "tpf_kappa\t" << config.tpf_kappa << "\n";
    f << "tpf_vdsg_coupling\t" << config.tpf_vdsg_coupling << "\n";
    f << "tpf_vdsg_mass_baseline_kg\t" << config.tpf_vdsg_mass_baseline_kg << "\n";
    f << "tpf_poisson_bins\t" << config.tpf_poisson_bins << "\n";
    f << "tpf_poisson_max_radius\t" << config.tpf_poisson_max_radius << "\n";
    f << "tpf_cooling_fraction\t" << config.tpf_cooling_fraction << "\n";
    f << "tpf_global_accel_shunt_enable\t" << (config.tpf_global_accel_shunt_enable ? 1 : 0) << "\n";
    f << "tpf_global_accel_shunt_fraction\t" << config.tpf_global_accel_shunt_fraction << "\n";
    f << "tpf_accel_pipeline_diagnostics_csv\t" << (config.tpf_accel_pipeline_diagnostics_csv ? 1 : 0) << "\n";
    f << "tpfcore_dump_readout_debug\t" << (config.tpfcore_dump_readout_debug ? 1 : 0) << "\n";
    f << "tpf_regime_diagnostics\tsee tpf_regime_diagnostics.txt (dynamical runs with provisional readout)\n";
    f << "tpf_trajectory_diagnostics\tsee tpf_trajectory_diagnostics.txt (dynamical runs; single-body only)\n";
    f << "tpf_closure_diagnostics\tsee tpf_closure_diagnostics.csv, tpf_closure_diagnostics.txt (tr_coherence_readout, single-body dynamical runs)\n";
    if (config.simulation_mode == galaxy::SimulationMode::tpf_two_body_sweep)
      f << "tpf_sweep_summary\tsee tpf_sweep_summary.csv, tpf_sweep_summary.txt\n";
    if (config.simulation_mode == galaxy::SimulationMode::bh_orbit_validation)
      f << "detailed_follow_up\tTPFCore bh_orbit_validation: full snapshots + regime + trajectory diagnostics (use after tpf_two_body_sweep to inspect a chosen speed_ratio)\n";
  }
  if (!run_config_path.empty())
    f << "config_loaded_run\t" << run_config_path << "\n";
  if (!package_defaults_path.empty())
    f << "config_loaded_package_defaults\t" << package_defaults_path << "\n";

  if (galaxy_init_audit && galaxy_init_audit->valid &&
      config.simulation_mode == galaxy::SimulationMode::galaxy) {
    f << "\n=== Galaxy initialization (resolved; audit) ===\n";
    f << "galaxy_init_template\t" << galaxy_init_audit->template_name << "\n";
    f << "galaxy_init_template_defaults_used\t" << (galaxy_init_audit->template_defaults_used ? 1 : 0)
      << "\n";
    f << "galaxy_init_seed\t" << galaxy_init_audit->seed << "\n";
    f << "galaxy_init_master_chaos\t" << galaxy_init_audit->master_chaos << "\n";
    f << "master_chaos_scales_position_noise\t" << (galaxy_init_audit->master_scales_position_noise ? 1 : 0)
      << "\n";
    f << "master_chaos_scales_velocity_angle_noise\t"
      << (galaxy_init_audit->master_scales_velocity_angle_noise ? 1 : 0) << "\n";
    f << "master_chaos_scales_velocity_magnitude_noise\t"
      << (galaxy_init_audit->master_scales_velocity_magnitude_noise ? 1 : 0) << "\n";
    f << "galaxy_init_position_noise_raw\t" << galaxy_init_audit->raw_position_noise << "\n";
    f << "galaxy_init_velocity_angle_noise_raw\t" << galaxy_init_audit->raw_velocity_angle_noise << "\n";
    f << "galaxy_init_velocity_magnitude_noise_raw\t" << galaxy_init_audit->raw_velocity_magnitude_noise
      << "\n";
    f << "eff_position_noise\t" << galaxy_init_audit->eff_position_noise << "\n";
    f << "eff_velocity_angle_noise_rad\t" << galaxy_init_audit->eff_velocity_angle_noise_rad << "\n";
    f << "eff_velocity_magnitude_noise\t" << galaxy_init_audit->eff_velocity_magnitude_noise << "\n";
    f << "used_new_state_noise\t" << (galaxy_init_audit->used_new_state_noise ? 1 : 0) << "\n";
    f << "used_legacy_velocity_noise\t" << (galaxy_init_audit->used_legacy_velocity_noise ? 1 : 0) << "\n";
    f << "velocity_noise_config\t" << config.velocity_noise
      << "  (legacy Cartesian vx,vy if new noises zero)\n";
    f << "galaxy_init_m2_amplitude_raw\t" << galaxy_init_audit->galaxy_init_m2_amplitude_raw << "\n";
    f << "galaxy_init_m3_amplitude_raw\t" << galaxy_init_audit->galaxy_init_m3_amplitude_raw << "\n";
    f << "galaxy_init_bar_amplitude_raw\t" << galaxy_init_audit->galaxy_init_bar_amplitude_raw << "\n";
    f << "galaxy_init_bar_axis_ratio_raw\t" << galaxy_init_audit->galaxy_init_bar_axis_ratio_raw << "\n";
    f << "galaxy_init_spiral_amplitude_raw\t" << galaxy_init_audit->galaxy_init_spiral_amplitude_raw << "\n";
    f << "galaxy_init_spiral_winding_raw\t" << galaxy_init_audit->galaxy_init_spiral_winding_raw << "\n";
    f << "galaxy_init_spiral_phase_raw\t" << galaxy_init_audit->galaxy_init_spiral_phase_raw << "\n";
    f << "galaxy_init_clumpiness_raw\t" << galaxy_init_audit->galaxy_init_clumpiness_raw << "\n";
    f << "galaxy_init_num_clumps_raw\t" << galaxy_init_audit->galaxy_init_num_clumps_raw << "\n";
    f << "galaxy_init_clump_radius_fraction_raw\t" << galaxy_init_audit->galaxy_init_clump_radius_fraction_raw
      << "\n";
    f << "velocity_noise_raw\t" << galaxy_init_audit->velocity_noise_raw << "\n";
    f << "structured_m2\t" << (galaxy_init_audit->structured_m2 ? 1 : 0) << "\n";
    f << "structured_m3\t" << (galaxy_init_audit->structured_m3 ? 1 : 0) << "\n";
    f << "structured_bar\t" << (galaxy_init_audit->structured_bar ? 1 : 0) << "\n";
    f << "structured_spiral\t" << (galaxy_init_audit->structured_spiral ? 1 : 0) << "\n";
    f << "structured_clumps\t" << (galaxy_init_audit->structured_clumps ? 1 : 0) << "\n";
    f << "galaxy_init_clumpiness\t" << galaxy_init_audit->galaxy_init_clumpiness << "\n";
    f << "galaxy_init_num_clumps\t" << galaxy_init_audit->galaxy_init_num_clumps << "\n";
    f << "galaxy_init_clump_radius_fraction\t" << galaxy_init_audit->galaxy_init_clump_radius_fraction
      << "\n";
    f << "galaxy_init_m2_amplitude\t" << galaxy_init_audit->galaxy_init_m2_amplitude << "\n";
    f << "galaxy_init_m3_amplitude\t" << galaxy_init_audit->galaxy_init_m3_amplitude << "\n";
    f << "galaxy_init_bar_amplitude\t" << galaxy_init_audit->galaxy_init_bar_amplitude << "\n";
    f << "galaxy_init_bar_axis_ratio\t" << galaxy_init_audit->galaxy_init_bar_axis_ratio << "\n";
    f << "galaxy_init_spiral_amplitude\t" << galaxy_init_audit->galaxy_init_spiral_amplitude << "\n";
    f << "galaxy_init_spiral_winding\t" << galaxy_init_audit->galaxy_init_spiral_winding << "\n";
    f << "galaxy_init_spiral_phase\t" << galaxy_init_audit->galaxy_init_spiral_phase << "\n";
    f << "structured_weight_cap_used\t" << galaxy_init_audit->weight_w_max << "\n";
    f << "rejection_sampling_fallbacks\t" << galaxy_init_audit->rejection_fallbacks << "\n";
    f << "clump_centers_count\t" << galaxy_init_audit->clump_centers_xy.size() << "\n";
    for (size_t i = 0; i < galaxy_init_audit->clump_centers_xy.size(); ++i) {
      f << "clump_center_" << i << "_x\t" << galaxy_init_audit->clump_centers_xy[i].first << "\n";
      f << "clump_center_" << i << "_y\t" << galaxy_init_audit->clump_centers_xy[i].second << "\n";
    }
    for (size_t i = 0; i < galaxy_init_audit->template_defaults_log.applied.size(); ++i) {
      f << "template_default_applied_" << i << "\t" << galaxy_init_audit->template_defaults_log.applied[i]
        << "\n";
    }
    for (size_t i = 0; i < galaxy_init_audit->template_defaults_log.warnings.size(); ++i) {
      f << "template_warning_" << i << "\t" << galaxy_init_audit->template_defaults_log.warnings[i] << "\n";
    }
    f << "galaxy_init_diagnostics_file\tgalaxy_init_diagnostics.txt\n";
    f << "=== End galaxy initialization ===\n";
  }
}

void write_snapshots(const std::string& output_dir,
                     const std::vector<Snapshot>& snapshots) {
  /* Decimal precision: default ostringstream uses 6 significant digits in scientific format, which
   * quantizes Earth–Moon-scale coordinates (~1e8 m) to ~10–100 m steps. Consecutive snapshots then
   * repeat identical x or y text while vx,vy still tick, producing sawtooth-like artifacts in
   * post-processed separation, energy, and Lz. Use max_digits10 for IEEE-754 double round-trip. */
  const int csv_prec = std::numeric_limits<double>::max_digits10;
  for (const auto& snap : snapshots) {
    std::ostringstream fname;
    fname << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snap.step << ".csv";
    std::ofstream f(fname.str());
    if (!f) continue;

    f << std::scientific << std::setprecision(csv_prec);
    f << "# step," << snap.step << ",time," << snap.time << "\n";
    f << "i,x,y,vx,vy,mass\n";
    const State& s = snap.state;
    for (int i = 0; i < s.n(); ++i) {
      f << i << "," << s.x[i] << "," << s.y[i] << "," << s.vx[i] << "," << s.vy[i] << "," << s.mass[i]
        << "\n";
    }
  }
}

}  // namespace galaxy
