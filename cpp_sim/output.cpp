#include "output.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace galaxy {

void write_run_info(const std::string& output_dir,
                    const Config& config,
                    int n_steps_done,
                    int n_snapshots,
                    int n_particles,
                    const std::string& run_config_path,
                    const std::string& package_defaults_path) {
  std::ostringstream path;
  path << output_dir << "/run_info.txt";
  std::ofstream f(path.str());
  if (!f) return;

  int n_star = (n_particles >= 0) ? n_particles : config.n_stars;

  /* Resolved config banner at top so it is obvious what was used */
  f << "=== Resolved config (built-in < package defaults < run config < CLI) ===\n";
  f << "run_config\t" << (run_config_path.empty() ? "(none)" : run_config_path) << "\n";
  f << "package_defaults\t" << (package_defaults_path.empty() ? "(none)" : package_defaults_path) << "\n";
  f << "physics_package\t" << config.physics_package << "\n";
  f << "simulation_mode\t" << mode_to_string(config.simulation_mode) << "\n";
  f << "output_dir\t" << output_dir << "\n";
  f << "n_stars\t" << n_star << "\n";
  f << "bh_mass\t" << config.bh_mass << "\n";
  if (config.physics_package == "TPFCore") {
    f << "tpfcore_enable_provisional_readout\t" << (config.tpfcore_enable_provisional_readout ? 1 : 0) << "\n";
    f << "tpfcore_readout_mode\t" << config.tpfcore_readout_mode << "\n";
    f << "tpfcore_readout_scale\t" << config.tpfcore_readout_scale << "\n";
    f << "tpfcore_theta_tt_scale\t" << config.tpfcore_theta_tt_scale << "\n";
    f << "tpfcore_theta_tr_scale\t" << config.tpfcore_theta_tr_scale << "\n";
    f << "tpf_kappa\t" << config.tpf_kappa << "\n";
    f << "tpf_vdsg_coupling\t" << config.tpf_vdsg_coupling << "\n";
    f << "tpf_vdsg_mass_baseline_kg\t" << config.tpf_vdsg_mass_baseline_kg << "\n";
    f << "tpf_poisson_bins\t" << config.tpf_poisson_bins << "\n";
    f << "tpf_poisson_max_radius\t" << config.tpf_poisson_max_radius << "\n";
    f << "tpf_cooling_fraction\t" << config.tpf_cooling_fraction << "\n";
  }
  f << "=== End resolved config ===\n\n";

  if (config.physics_package == "TPFCore") {
    f << "=== TPFCore parameter roles (theory vs regularization vs exploratory vs provisional) ===\n";
    f << "fixed_theory\tlambda=1/4 (LAMBDA_4D; manuscript structure; not tunable)\n";
    f << "numerical_regularization\ttpfcore_source_softening, effective_source_softening (eps for Phi)\n";
    f << "provisional_readout\ttpfcore_enable_provisional_readout, readout_mode, readout_scale (weak-field calibrated effective scale), theta_tt_scale, theta_tr_scale, dump_readout_debug (experimental closure)\n";
    f << "inspection\ttpfcore_probe_radius_min/max, probe_samples, dump_theta_profile, dump_invariant_profile\n";
    f << "=== End TPFCore parameter roles ===\n\n";
  }

  f << "dt\t" << config.dt << "\n";
  f << "n_steps\t" << n_steps_done << "\n";
  f << "snapshot_every\t" << config.snapshot_every << "\n";
  f << "softening\t" << config.softening << "\n";
  f << "star_mass\t" << config.star_mass << "\n";
  f << "bh_mass\t" << config.bh_mass << "\n";
  f << "enable_star_star_gravity\t" << (config.enable_star_star_gravity ? 1 : 0) << "\n";
  f << "galaxy_radius\t" << config.galaxy_radius << "\n";
  f << "total_simulated_time\t" << (n_steps_done * config.dt) << "\n";
  f << "number_of_snapshots\t" << n_snapshots << "\n";
  f << "n_stars\t" << n_star << "\n";
  f << "simulation_mode\t" << static_cast<int>(config.simulation_mode) << "\n";
  f << "physics_package\t" << config.physics_package << "\n";
  if (config.physics_package == "TPFCore") {
    double src_eps = (config.tpfcore_source_softening > 0.0) ? config.tpfcore_source_softening : config.softening;
    f << "tpfcore_enable_provisional_readout\t" << (config.tpfcore_enable_provisional_readout ? 1 : 0) << "\n";
    f << "tpfcore_provisional_source_ansatz\t1\n";
    f << "tpfcore_source_softening\t" << src_eps << "\n";
    f << "tpfcore_probe_radius_min\t" << config.tpfcore_probe_radius_min << "\n";
    f << "tpfcore_probe_radius_max\t" << config.tpfcore_probe_radius_max << "\n";
    f << "tpfcore_probe_samples\t" << config.tpfcore_probe_samples << "\n";
    f << "tpfcore_residual_method\tanalytic\n";
    f << "tpfcore_residual_step\t" << config.tpfcore_residual_step << "\n";
    f << "tpfcore_readout_mode\t" << config.tpfcore_readout_mode << "\n";
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
    f << "tpfcore_dump_readout_debug\t" << (config.tpfcore_dump_readout_debug ? 1 : 0) << "\n";
    f << "tpf_regime_diagnostics\tsee tpf_regime_diagnostics.txt (dynamical runs with provisional readout)\n";
    f << "tpf_trajectory_diagnostics\tsee tpf_trajectory_diagnostics.txt (dynamical runs; single-body only)\n";
    f << "tpf_closure_diagnostics\tsee tpf_closure_diagnostics.csv, tpf_closure_diagnostics.txt (tr_coherence_readout, single-body dynamical runs)\n";
    if (config.simulation_mode == galaxy::SimulationMode::tpf_two_body_sweep)
      f << "tpf_sweep_summary\tsee tpf_sweep_summary.csv, tpf_sweep_summary.txt\n";
    if (config.simulation_mode == galaxy::SimulationMode::two_body_orbit)
      f << "detailed_follow_up\tTPFCore two_body_orbit: full snapshots + regime + trajectory diagnostics (use after tpf_two_body_sweep to inspect a chosen speed_ratio)\n";
  }
  if (!run_config_path.empty())
    f << "config_loaded_run\t" << run_config_path << "\n";
  if (!package_defaults_path.empty())
    f << "config_loaded_package_defaults\t" << package_defaults_path << "\n";
}

void write_snapshots(const std::string& output_dir,
                     const std::vector<Snapshot>& snapshots) {
  for (const auto& snap : snapshots) {
    std::ostringstream fname;
    fname << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snap.step << ".csv";
    std::ofstream f(fname.str());
    if (!f) continue;

    f << "# step," << snap.step << ",time," << std::scientific << snap.time << "\n";
    f << "i,x,y,vx,vy,mass\n";
    const State& s = snap.state;
    for (int i = 0; i < s.n(); ++i) {
      f << i << ","
        << std::scientific << s.x[i] << "," << s.y[i] << ","
        << s.vx[i] << "," << s.vy[i] << "," << s.mass[i] << "\n";
    }
  }
}

}  // namespace galaxy
