#include "config.hpp"
#include "init_conditions.hpp"
#include "output.hpp"
#include "physics/physics_package.hpp"
#include "physics/TPFCore/tpf_core_package.hpp"
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
  // 1. Find run config path (root configs/ only; cpp_sim/configs/ is not used)
  std::string run_config_path = galaxy::find_run_config_path();
  if (!galaxy::check_run_config_canonical(run_config_path))
    return 1;
  if (!run_config_path.empty()) {
    std::cout << "Run config selected: " << run_config_path << "\n";
  } else {
    std::cout << "Run config selected: (none)\n";
  }

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
  }
  if (!run_config_path.empty()) {
    galaxy::load_config_file(run_config_path, config);
  }

  config.run_id = run_id_from_time();
  config.output_dir = "outputs/" + config.run_id;

  bool cli_override_mode = false;
  if (argc >= 2) {
    try {
      config.simulation_mode = galaxy::parse_mode(argv[1]);
      cli_override_mode = true;
      std::cout << "CLI override applied: simulation_mode=" << galaxy::mode_to_string(config.simulation_mode) << "\n";
    } catch (const std::exception& e) {
      std::cerr << e.what() << "\nAllowed: galaxy, two_body_orbit, symmetric_pair, small_n_conservation, timestep_convergence, tpf_single_source_inspect, tpf_symmetric_pair_inspect, tpf_single_source_optimize_c, tpf_two_body_sweep, tpf_weak_field_calibration\n";
      return 1;
    }
  }

  /* Probe run config for consistency check (what did the file say?) */
  std::string run_cfg_physics = run_config_path.empty() ? "" : galaxy::probe_config_key(run_config_path, "physics_package");
  std::string run_cfg_mode_str = run_config_path.empty() ? "" : galaxy::probe_config_key(run_config_path, "simulation_mode");

  /* Hard failure: run config says TPFCore + two_body_orbit but resolved is galaxy with no CLI override */
  if (!run_config_path.empty() && cli_override_mode == false) {
    if (run_cfg_physics == "TPFCore" && run_cfg_mode_str == "two_body_orbit" &&
        config.simulation_mode == galaxy::SimulationMode::galaxy) {
      std::cerr << "Config mismatch: " << run_config_path << " specifies physics_package=TPFCore and simulation_mode=two_body_orbit, "
                << "but resolved simulation_mode is galaxy. Refusing to run. Fix config precedence or run: ./galaxy_sim two_body_orbit\n";
      return 1;
    }
  }

  /* Startup banner: resolved config */
  std::cout << "--- Resolved config ---\n";
  std::cout << "RUN CONFIG: " << (run_config_path.empty() ? "(none)" : run_config_path) << "\n";
  std::cout << "PACKAGE DEFAULTS: " << (package_defaults_path.empty() ? "(none)" : package_defaults_path) << "\n";
  std::cout << "PHYSICS PACKAGE: " << config.physics_package << "\n";
  std::cout << "SIMULATION MODE: " << galaxy::mode_to_string(config.simulation_mode) << "\n";
  std::cout << "OUTPUT DIR: " << config.output_dir << "\n";
  std::cout << "n_stars: " << config.n_stars << "  bh_mass: " << config.bh_mass << "\n";
  if (config.physics_package == "TPFCore") {
    std::cout << "tpfcore_enable_provisional_readout: " << (config.tpfcore_enable_provisional_readout ? "true" : "false")
              << "  tpfcore_readout_mode: " << config.tpfcore_readout_mode
              << "  tpfcore_isotropic_correction_c: " << config.tpfcore_isotropic_correction_c;
    if (config.tpfcore_readout_mode == "tr_coherence_readout")
      std::cout << "  theta_tt_scale: " << config.tpfcore_theta_tt_scale << "  theta_tr_scale: " << config.tpfcore_theta_tr_scale;
    std::cout << "\n";
  }
  std::cout << "------------------------\n";

  if (config.physics_package.empty())
    config.physics_package = "Newtonian";
  galaxy::PhysicsPackage* physics = galaxy::get_physics_package(config.physics_package);
  if (!physics) {
    std::cerr << "Unknown physics_package: '" << config.physics_package << "'. Available: Newtonian, TPFCore (add more in physics/registry.cpp).\n";
    return 1;
  }
  physics->init_from_config(config);

  if (config.physics_package == "TPFCore") {
    galaxy::TPFCorePackage* tpf = dynamic_cast<galaxy::TPFCorePackage*>(physics);
    std::cout << "Physics: TPFCore (primitive TPF structure)\n";
    std::cout << "  Hessian-based provisional ansatz: Phi=-M/sqrt(r^2+eps^2), Theta=Hess(Phi)+B(r)*delta\n";
    std::cout << "  Parameter roles: lambda=1/4 (fixed theory) | eps (numerical regularization) | c (exploratory ansatz) | readout (provisional experimental)\n";
    std::cout << "  Provisional readout: " << (tpf && tpf->provisional_readout_enabled() ? "enabled" : "disabled");
    if (tpf && tpf->provisional_readout_enabled()) {
      std::cout << " (readout mode: " << tpf->readout_mode() << ", scale=" << config.tpfcore_readout_scale << " [weak-field calibrated])";
      if (config.tpfcore_readout_mode == "tr_coherence_readout")
        std::cout << " [exploratory t-r; theta_tt_scale=" << config.tpfcore_theta_tt_scale << ", theta_tr_scale=" << config.tpfcore_theta_tr_scale << "]";
      if (config.tpfcore_dump_readout_debug)
        std::cout << " [readout debug CSV enabled]";
    }
    std::cout << "\n";

    if (!tpf->provisional_readout_enabled()) {
      bool is_dynamical = (config.simulation_mode == galaxy::SimulationMode::galaxy ||
                           config.simulation_mode == galaxy::SimulationMode::two_body_orbit ||
                           config.simulation_mode == galaxy::SimulationMode::symmetric_pair ||
                           config.simulation_mode == galaxy::SimulationMode::small_n_conservation ||
                           config.simulation_mode == galaxy::SimulationMode::timestep_convergence);
      if (is_dynamical) {
        std::cerr << "TPFCore does not support dynamical modes (galaxy, two_body_orbit, etc.) unless provisional readout is enabled.\n";
        std::cerr << "Use physics_package = Newtonian for dynamics, or run inspection modes: tpf_single_source_inspect, tpf_symmetric_pair_inspect, tpf_single_source_optimize_c.\n";
        return 1;
      }
    }
  }

  if (!ensure_dir("outputs") || !ensure_dir(config.output_dir)) {
    std::cerr << "Failed to create output dir " << config.output_dir << "\n";
    return 1;
  }

  std::cout << "Galaxy N-body (C++)\n";

  galaxy::TPFCorePackage* tpfcore = dynamic_cast<galaxy::TPFCorePackage*>(physics);

  if (config.simulation_mode == galaxy::SimulationMode::tpf_single_source_inspect) {
    if (!tpfcore) {
      std::cerr << "tpf_single_source_inspect requires physics_package = TPFCore.\n";
      return 1;
    }
    std::cout << "Inspection mode: tpf_single_source_inspect\n";
    std::cout << "Residual checking: enabled (analytic)\n";
    std::cout << "Probe: +x axis, r in [" << config.tpfcore_probe_radius_min << ", " << config.tpfcore_probe_radius_max << "], n=" << config.tpfcore_probe_samples << "\n";
    tpfcore->run_single_source_inspect(config, config.output_dir);
    std::cout << "Wrote " << config.output_dir << "/theta_profile.csv, invariant_profile.csv, field_summary.txt\n";
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    return 0;
  }

  if (config.simulation_mode == galaxy::SimulationMode::tpf_symmetric_pair_inspect) {
    if (!tpfcore) {
      std::cerr << "tpf_symmetric_pair_inspect requires physics_package = TPFCore.\n";
      return 1;
    }
    std::cout << "Inspection mode: tpf_symmetric_pair_inspect\n";
    std::cout << "Residual checking: enabled (analytic)\n";
    std::cout << "Probe: +x and +y axes, r in [" << config.tpfcore_probe_radius_min << ", " << config.tpfcore_probe_radius_max << "], n=" << config.tpfcore_probe_samples << " per axis\n";
    tpfcore->run_symmetric_pair_inspect(config, config.output_dir);
    std::cout << "Wrote " << config.output_dir << "/theta_profile.csv, invariant_profile.csv, field_summary.txt\n";
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    return 0;
  }

  if (config.simulation_mode == galaxy::SimulationMode::tpf_single_source_optimize_c) {
    if (!tpfcore) {
      std::cerr << "tpf_single_source_optimize_c requires physics_package = TPFCore.\n";
      return 1;
    }
    std::cout << "Utility mode: tpf_single_source_optimize_c (exploratory ansatz-tuning)\n";
    std::cout << "Sweep range: c in [" << config.tpfcore_c_sweep_min << ", " << config.tpfcore_c_sweep_max << "]\n";
    std::cout << "Number of steps: " << config.tpfcore_c_sweep_steps << "\n";
    std::cout << "Chosen objective: " << config.tpfcore_c_objective << "\n";
    tpfcore->run_single_source_optimize_c(config, config.output_dir);
    std::cout << "Wrote " << config.output_dir << "/c_sweep.csv, c_sweep_summary.txt\n";
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    return 0;
  }

  if (config.simulation_mode == galaxy::SimulationMode::tpf_two_body_sweep) {
    if (!tpfcore) {
      std::cerr << "tpf_two_body_sweep requires physics_package = TPFCore.\n";
      return 1;
    }
    if (!tpfcore->provisional_readout_enabled()) {
      std::cerr << "tpf_two_body_sweep requires tpfcore_enable_provisional_readout = true.\n";
      return 1;
    }
    if (!ensure_dir("outputs") || !ensure_dir(config.output_dir)) {
      std::cerr << "Failed to create output dir " << config.output_dir << "\n";
      return 1;
    }
    /* Conservative default: sweep validation_two_body_speed_ratio only. */
    const double speed_ratios[] = { 0.5, 0.8, 1.0, 1.2, 1.5 };
    const int n_speed = sizeof(speed_ratios) / sizeof(speed_ratios[0]);

    std::cout << "Sweep mode: tpf_two_body_sweep (TPFCore two-body trajectory classification)\n";
    std::cout << "Sweeping speed_ratio over " << n_speed << " values; reusing existing simulation + diagnostics.\n";

    struct SweepRow {
      double speed_ratio;
      double r_initial, r_final, r_min, r_max, radial_drift, revolutions;
      std::string trajectory_class;
      double mean_theta_norm, max_theta_norm;
    };
    std::vector<SweepRow> rows;

    for (int i = 0; i < n_speed; ++i) {
      galaxy::Config c = config;
      c.validation_two_body_speed_ratio = speed_ratios[i];
      c.enable_star_star_gravity = false;
      physics->init_from_config(c);

      galaxy::State state;
      galaxy::init_two_body(c, state);
      int n_steps = c.validation_n_steps;
      int snapshot_every = c.validation_snapshot_every;

      auto snapshots = galaxy::run_simulation(c, state, physics, n_steps, snapshot_every, nullptr, 0);

      auto traj = tpfcore->compute_trajectory_summary(snapshots);
      auto regime = tpfcore->compute_regime_summary(snapshots, c, config.output_dir);

      SweepRow row = {};
      row.speed_ratio = speed_ratios[i];
      if (traj.valid) {
        row.r_initial = traj.r_initial;
        row.r_final = traj.r_final;
        row.r_min = traj.r_min;
        row.r_max = traj.r_max;
        row.radial_drift = traj.radial_drift;
        row.revolutions = traj.revolutions;
        row.trajectory_class = traj.trajectory_class;
      }
      if (regime.valid) {
        row.mean_theta_norm = regime.mean_theta_norm;
        row.max_theta_norm = regime.max_theta_norm;
      }
      rows.push_back(row);
    }

    std::string out_dir = config.output_dir;
    std::ofstream csv(out_dir + "/tpf_sweep_summary.csv");
    if (csv) {
      csv << "speed_ratio,r_initial,r_final,r_min,r_max,radial_drift,revolutions,trajectory_class,mean_theta_norm,max_theta_norm\n";
      for (const auto& row : rows) {
        csv << std::scientific << row.speed_ratio << "," << row.r_initial << "," << row.r_final << ","
            << row.r_min << "," << row.r_max << "," << row.radial_drift << "," << row.revolutions << ","
            << row.trajectory_class << "," << row.mean_theta_norm << "," << row.max_theta_norm << "\n";
      }
    }
    std::ofstream txt(out_dir + "/tpf_sweep_summary.txt");
    if (txt) {
      txt << "TPFCore two-body parameter sweep (exploratory; downstream of physics)\n";
      txt << "Swept parameter: validation_two_body_speed_ratio\n";
      txt << "Values: ";
      for (int i = 0; i < n_speed; ++i) txt << (i ? ", " : "") << speed_ratios[i];
      txt << "\n\n";
      txt << "speed_ratio\tr_initial\tr_final\tr_min\tr_max\tr_drift\trevolutions\tclass\tmean_theta_norm\tmax_theta_norm\n";
      for (const auto& row : rows) {
        txt << std::scientific << row.speed_ratio << "\t" << row.r_initial << "\t" << row.r_final << "\t"
            << row.r_min << "\t" << row.r_max << "\t" << row.radial_drift << "\t" << row.revolutions << "\t"
            << row.trajectory_class << "\t" << row.mean_theta_norm << "\t" << row.max_theta_norm << "\n";
      }
    }
    std::cout << "Wrote " << out_dir << "/tpf_sweep_summary.csv, tpf_sweep_summary.txt\n";
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    return 0;
  }

  if (config.simulation_mode == galaxy::SimulationMode::tpf_weak_field_calibration) {
    if (!tpfcore) {
      std::cerr << "tpf_weak_field_calibration requires physics_package = TPFCore.\n";
      return 1;
    }
    if (!tpfcore->provisional_readout_enabled()) {
      std::cerr << "tpf_weak_field_calibration requires tpfcore_enable_provisional_readout = true.\n";
      return 1;
    }
    if (!ensure_dir("outputs") || !ensure_dir(config.output_dir)) {
      std::cerr << "Failed to create output dir " << config.output_dir << "\n";
      return 1;
    }
    std::cout << "Calibration mode: tpf_weak_field_calibration (TPF vs Newtonian benchmark)\n";
    tpfcore->run_weak_field_calibration(config, config.output_dir);
    std::cout << "Wrote " << config.output_dir << "/tpf_weak_field_calibration.csv, tpf_weak_field_calibration.txt\n";
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    return 0;
  }

  galaxy::State state;
  int n_steps = config.n_steps;
  int snapshot_every = config.snapshot_every;
  double bh_mass = config.bh_mass;

  switch (config.simulation_mode) {
    case galaxy::SimulationMode::tpf_single_source_inspect:
    case galaxy::SimulationMode::tpf_symmetric_pair_inspect:
    case galaxy::SimulationMode::tpf_single_source_optimize_c:
    case galaxy::SimulationMode::tpf_two_body_sweep:
    case galaxy::SimulationMode::tpf_weak_field_calibration:
      std::cerr << "Internal error: inspection/utility/sweep/calibration modes should have returned earlier.\n";
      return 1;
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

  if (config.physics_package == "TPFCore") {
    galaxy::TPFCorePackage* tpf = dynamic_cast<galaxy::TPFCorePackage*>(physics);
    if (tpf && tpf->provisional_readout_enabled()) {
      tpf->write_readout_debug(snapshots, config, config.output_dir);
      if (config.tpfcore_dump_readout_debug)
        std::cout << "Wrote " << config.output_dir << "/tpf_readout_debug.csv\n";
      tpf->write_regime_diagnostics(snapshots, config, config.output_dir);
      std::cout << "Wrote " << config.output_dir << "/tpf_regime_diagnostics.txt\n";
      tpf->write_trajectory_diagnostics(snapshots, config, config.output_dir);
      std::cout << "Wrote " << config.output_dir << "/tpf_trajectory_diagnostics.txt\n";
      tpf->write_closure_diagnostics(snapshots, config, config.output_dir);
      if (config.physics_package == "TPFCore" && config.tpfcore_readout_mode == "tr_coherence_readout" && snapshots[0].state.n() == 1)
        std::cout << "Wrote " << config.output_dir << "/tpf_closure_diagnostics.csv, tpf_closure_diagnostics.txt\n";
    }
  }

  std::cout << "Snapshots: " << snapshots.size() << "\n";
  std::cout << "Output directory: " << config.output_dir << "\n";

  return 0;
}
