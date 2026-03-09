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
    f << "tpfcore_isotropic_correction_c\t" << config.tpfcore_isotropic_correction_c << "\n";
  }
  f << "=== End resolved config ===\n\n";

  f << "dt\t" << config.dt << "\n";
  f << "n_steps\t" << n_steps_done << "\n";
  f << "snapshot_every\t" << config.snapshot_every << "\n";
  f << "softening\t" << config.softening << "\n";
  f << "star_mass\t" << config.star_mass << "\n";
  f << "bh_mass\t" << config.bh_mass << "\n";
  f << "enable_star_star_gravity\t" << (config.enable_star_star_gravity ? 1 : 0) << "\n";
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
    f << "tpfcore_isotropic_correction_c\t" << config.tpfcore_isotropic_correction_c << "\n";
    f << "tpfcore_readout_mode\t" << config.tpfcore_readout_mode << "\n";
    f << "tpfcore_readout_scale\t" << config.tpfcore_readout_scale << "\n";
    f << "tpfcore_dump_readout_debug\t" << (config.tpfcore_dump_readout_debug ? 1 : 0) << "\n";
    f << "tpfcore_c_sweep_min\t" << config.tpfcore_c_sweep_min << "\n";
    f << "tpfcore_c_sweep_max\t" << config.tpfcore_c_sweep_max << "\n";
    f << "tpfcore_c_sweep_steps\t" << config.tpfcore_c_sweep_steps << "\n";
    f << "tpfcore_c_objective\t" << config.tpfcore_c_objective << "\n";
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
