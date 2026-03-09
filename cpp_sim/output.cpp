#include "output.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace galaxy {

namespace {
const double FOUR_PI = 4.0 * 3.14159265358979323846;
}

void write_run_info(const std::string& output_dir,
                    const Config& config,
                    int n_steps_done,
                    int n_snapshots,
                    int n_particles) {
  std::ostringstream path;
  path << output_dir << "/run_info.txt";
  std::ofstream f(path.str());
  if (!f) return;

  int n_star = (n_particles >= 0) ? n_particles : config.n_stars;

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
  if (config.physics_package == "TPF") {
    double alpha = config.tpf_match_newtonian_scale ? FOUR_PI
        : (config.tpf_alpha > 0.0 ? config.tpf_alpha : FOUR_PI);
    f << "tpf_alpha\t" << std::scientific << alpha << "\n";
    f << "tpf_match_newtonian_scale\t" << (config.tpf_match_newtonian_scale ? 1 : 0) << "\n";
  }
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
