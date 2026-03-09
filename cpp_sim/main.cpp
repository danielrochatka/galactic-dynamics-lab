#include "config.hpp"
#include "init_conditions.hpp"
#include "output.hpp"
#include "physics/physics_package.hpp"
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
#include <io.h>
#define MKDIR(path, mode) _mkdir(path)
#define IS_STDOUT_TERMINAL() (_isatty(_fileno(stdout)) != 0)
#else
#include <sys/stat.h>
#include <unistd.h>
#define MKDIR(path, mode) mkdir(path, mode)
#define IS_STDOUT_TERMINAL() (isatty(STDOUT_FILENO) != 0)
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

// Format seconds as HH:MM:SS or MM:SS
std::string format_elapsed(double sec) {
  int s = static_cast<int>(sec + 0.5);
  if (s < 0) s = 0;
  int m = s / 60;
  s %= 60;
  int h = m / 60;
  m %= 60;
  std::ostringstream os;
  os << std::setfill('0');
  if (h > 0)
    os << std::setw(2) << h << ":" << std::setw(2) << m << ":" << std::setw(2) << s;
  else
    os << std::setw(2) << m << ":" << std::setw(2) << s;
  return os.str();
}

}  // namespace

int main(int argc, char** argv) {
  // 1. Find run config path (for probe and later full load)
  std::string run_config_path = galaxy::find_run_config_path();

  // 2. Probe run config for physics_package (needed to load package defaults)
  std::string physics_pkg = "Newtonian";
  if (!run_config_path.empty()) {
    std::string probed = galaxy::probe_config_key(run_config_path, "physics_package");
    if (!probed.empty()) physics_pkg = probed;
  }

  // 3. Layered load: built-in -> package defaults -> run config
  galaxy::Config config;

  std::string package_defaults_path = galaxy::find_package_defaults_path(physics_pkg);
  if (!package_defaults_path.empty()) {
    galaxy::load_config_file(package_defaults_path, config);
    std::cout << "Loaded package defaults: " << package_defaults_path << "\n";
  }

  if (!run_config_path.empty()) {
    galaxy::load_config_file(run_config_path, config);
    std::cout << "Loaded run config: " << run_config_path << "\n";
  }

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

  if (config.physics_package.empty())
    config.physics_package = "Newtonian";
  galaxy::PhysicsPackage* physics = galaxy::get_physics_package(config.physics_package);
  if (!physics) {
    std::cerr << "Unknown physics_package: '" << config.physics_package << "'. Available: Newtonian, TPF (add more in physics/registry.cpp).\n";
    return 1;
  }
  physics->init_from_config(config);

  if (config.physics_package == "TPF") {
    std::cout << "Physics: TPF weak-field correspondence package\n";
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
      std::cout << "Running galaxy: n_stars=" << config.n_stars
                << ", n_steps=" << n_steps
                << ", dt=" << config.dt
                << ", sim_time=" << (n_steps * config.dt) << "\n";
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
      E0 = galaxy::compute_kinetic_energy(state0) + physics->compute_potential_energy(state0, config.bh_mass, config.softening, config.enable_star_star_gravity);

      for (double dt : dts) {
        int steps = static_cast<int>(std::round(total_time / dt));
        galaxy::Config c2 = config;
        c2.dt = dt;
        galaxy::State s0;
        galaxy::init_two_body(c2, s0);
        auto snaps = galaxy::run_simulation(c2, s0, physics, steps, std::max(1, config.validation_snapshot_every));
        const auto& last = snaps.back().state;
        double r_final = std::sqrt(last.x[0] * last.x[0] + last.y[0] * last.y[0]);
        double Lz = L_z_total(last);
        double E_final = galaxy::compute_kinetic_energy(last) + physics->compute_potential_energy(last, config.bh_mass, config.softening, config.enable_star_star_gravity);
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

  // Progress reporting: only for galaxy (long runs); in-place line when stdout is a terminal
  int progress_interval = 0;
  galaxy::ProgressCallback progress_callback;
  bool progress_to_terminal = false;
  if (config.simulation_mode == galaxy::SimulationMode::galaxy && n_steps > 0) {
    progress_interval = std::max(1, std::min(1000, n_steps / 100));
    auto start_wall = std::chrono::steady_clock::now();
    progress_to_terminal = IS_STDOUT_TERMINAL();
    progress_callback = [start_wall, &config, progress_to_terminal](int step, int n_steps, double sim_time) {
      auto now = std::chrono::steady_clock::now();
      double elapsed_sec = 1e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_wall).count();
      double eta_sec = (step > 0 && step < n_steps)
          ? (elapsed_sec / step) * (n_steps - step)
          : 0.0;
      double pct = 100.0 * step / n_steps;
      if (progress_to_terminal) {
        std::cout << "\r[ " << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%] "
                  << "step " << step << "/" << n_steps
                  << ", sim t=" << std::setprecision(2) << sim_time
                  << ", elapsed=" << format_elapsed(elapsed_sec)
                  << ", eta=" << format_elapsed(eta_sec) << "    " << std::flush;
      } else {
        std::cout << "[ " << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%] "
                  << "step " << step << "/" << n_steps
                  << ", sim t=" << std::setprecision(2) << sim_time
                  << ", elapsed=" << format_elapsed(elapsed_sec)
                  << ", eta=" << format_elapsed(eta_sec) << "\n"
                  << std::flush;
      }
    };
  }

  auto run_start = std::chrono::steady_clock::now();
  auto snapshots = galaxy::run_simulation(config, state, physics, n_steps, snapshot_every,
                                          progress_callback, progress_interval);
  double run_elapsed_sec = 1e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now() - run_start).count();

  if (config.simulation_mode == galaxy::SimulationMode::galaxy) {
    if (progress_to_terminal)
      std::cout << "\n";
    std::cout << "Completed in " << format_elapsed(run_elapsed_sec)
              << ". Snapshots: " << snapshots.size()
              << ". Output: " << config.output_dir << "\n";
  }

  if (config.save_run_info) {
    galaxy::write_run_info(config.output_dir, config, n_steps, static_cast<int>(snapshots.size()), state.n(),
                           run_config_path, package_defaults_path);
    std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
  }
  if (config.save_snapshots) {
    galaxy::write_snapshots(config.output_dir, snapshots);
    std::cout << "Wrote " << config.output_dir << "/snapshot_*.csv\n";
  }

  std::cout << "Snapshots: " << snapshots.size() << "\n";
  std::cout << "Output directory: " << config.output_dir << "\n";

  return 0;
}
