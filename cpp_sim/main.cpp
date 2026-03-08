#include "config.hpp"
#include "init_conditions.hpp"
#include "output.hpp"
#include "physics.hpp"
#include "simulation.hpp"
#include "types.hpp"

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path, mode) _mkdir(path)
#else
#include <sys/stat.h>
#define MKDIR(path, mode) mkdir(path, mode)
#endif

namespace {

std::string run_id_from_time() {
  auto now = std::chrono::system_clock::now();
  auto t = std::chrono::system_clock::to_time_t(now);
  std::tm* bt = std::localtime(&t);
  std::ostringstream os;
  os << std::put_time(bt, "%Y%m%d_%H%M%S");
  return os.str();
}

bool ensure_dir(const std::string& path) {
  return MKDIR(path.c_str(), 0755) == 0 || errno == EEXIST;
}

double L_z_total(const galaxy::State& s) {
  double L = 0;
  for (int i = 0; i < s.n(); ++i)
    L += s.mass[i] * (s.x[i] * s.vy[i] - s.y[i] * s.vx[i]);
  return L;
}

}  // namespace

int main(int argc, char** argv) {
  galaxy::Config config;
  config.run_id = run_id_from_time();
  config.output_dir = "outputs/" + config.run_id;

  if (argc >= 2) {
    try {
      config.simulation_mode = galaxy::parse_mode(argv[1]);
    } catch (const std::exception& e) {
      std::cerr << e.what() << "\nAllowed: galaxy, two_body_orbit, symmetric_pair, small_n_conservation, timestep_convergence\n";
      return 1;
    }
  }

  if (!ensure_dir("outputs") || !ensure_dir(config.output_dir)) {
    std::cerr << "Failed to create output dir " << config.output_dir << "\n";
    return 1;
  }

  std::cout << "Galaxy N-body (C++)\n";
  std::cout << "Output directory: " << config.output_dir << "\n";
  std::cout << "Mode: " << static_cast<int>(config.simulation_mode) << "\n";

  galaxy::State state;
  int n_steps = config.n_steps;
  int snapshot_every = config.snapshot_every;
  double bh_mass = config.bh_mass;

  switch (config.simulation_mode) {
    case galaxy::SimulationMode::galaxy: {
      galaxy::init_galaxy_disk(config, state);
      std::cout << "Stars: " << config.n_stars << ", steps: " << n_steps << "\n";
      break;
    }
    case galaxy::SimulationMode::two_body_orbit: {
      galaxy::init_two_body(config, state);
      config.enable_star_star_gravity = false;
      n_steps = config.validation_n_steps;
      snapshot_every = config.validation_snapshot_every;
      std::cout << "Two-body orbit: r0=" << config.validation_two_body_radius
                << " speed_ratio=" << config.validation_two_body_speed_ratio << "\n";
      break;
    }
    case galaxy::SimulationMode::symmetric_pair: {
      galaxy::init_symmetric_pair(config, state);
      if (!config.validation_symmetric_include_bh) bh_mass = 0.0;
      config.bh_mass = bh_mass;  // so run_simulation uses it
      n_steps = config.validation_n_steps;
      snapshot_every = config.validation_snapshot_every;
      std::cout << "Symmetric pair: a=" << config.validation_symmetric_separation
                << " include_bh=" << config.validation_symmetric_include_bh << "\n";
      break;
    }
    case galaxy::SimulationMode::small_n_conservation: {
      galaxy::init_small_n(config, state);
      config.n_stars = state.n();
      n_steps = config.validation_n_steps;
      snapshot_every = config.validation_snapshot_every;
      std::cout << "Small-N: n=" << state.n() << "\n";
      break;
    }
    case galaxy::SimulationMode::timestep_convergence: {
      // Run two_body at dt, dt/2, dt/4 (same total time)
      double total_time = config.validation_n_steps * config.dt;
      std::vector<double> dts = {config.dt, config.dt / 2, config.dt / 4};
      config.enable_star_star_gravity = false;

      std::cout << "Timestep convergence (two_body), total_time=" << total_time << "\n";
      std::ofstream summary(config.output_dir + "/validation_timestep_convergence.txt");
      summary << "Timestep convergence (two_body_orbit)\n";
      summary << "dt\tfinal_x\tfinal_y\tfinal_r\tL_z\tE_drift\n";

      double E0 = 0;
      galaxy::State state0;
      galaxy::init_two_body(config, state0);
      E0 = galaxy::compute_kinetic_energy(state0) + galaxy::compute_potential_energy(state0, config.bh_mass, config.softening);

      for (double dt : dts) {
        int steps = static_cast<int>(std::round(total_time / dt));
        galaxy::Config c2 = config;
        c2.dt = dt;
        galaxy::State s0;
        galaxy::init_two_body(c2, s0);
        auto snaps = galaxy::run_simulation(c2, s0, steps, std::max(1, config.validation_snapshot_every));
        const auto& last = snaps.back().state;
        double r_final = std::sqrt(last.x[0] * last.x[0] + last.y[0] * last.y[0]);
        double Lz = L_z_total(last);
        double E_final = galaxy::compute_kinetic_energy(last) + galaxy::compute_potential_energy(last, config.bh_mass, config.softening);
        double E_drift = std::abs(E_final - E0);
        summary << std::scientific << dt << "\t" << last.x[0] << "\t" << last.y[0] << "\t"
                << r_final << "\t" << Lz << "\t" << E_drift << "\n";
        std::cout << "  dt=" << dt << " steps=" << steps << " final_r=" << r_final << " L_z=" << Lz << " E_drift=" << E_drift << "\n";
      }
      summary.close();
      std::cout << "Wrote " << config.output_dir << "/validation_timestep_convergence.txt\n";
      std::cout << "Output directory: " << config.output_dir << "\n";
      return 0;
    }
  }

  auto snapshots = galaxy::run_simulation(config, state, n_steps, snapshot_every);

  galaxy::write_run_info(config.output_dir, config, n_steps, static_cast<int>(snapshots.size()), state.n());
  galaxy::write_snapshots(config.output_dir, snapshots);

  std::cout << "Snapshots: " << snapshots.size() << "\n";
  std::cout << "Wrote " << config.output_dir << "/run_info.txt and snapshot_*.csv\n";
  std::cout << "Output directory: " << config.output_dir << "\n";

  return 0;
}
