#include "accel_pipeline_stats.hpp"
#include "config.hpp"
#include "compare_orchestration.hpp"
#include "force_compare.hpp"
#include "galaxy_init.hpp"
#include "git_provenance.hpp"
#include "init_conditions.hpp"
#include "render_audit.hpp"
#include "resolved_scenario.hpp"
#include "output.hpp"
#include "physics/physics_package.hpp"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "physics/TPFCore/v11_weak_field_correspondence.hpp"
#include "simulation.hpp"
#include "types.hpp"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdint>
#include <cstdio>

#ifdef _WIN32
#include <direct.h>
#include <io.h>
#define MKDIR(path, mode) _mkdir(path)
#define IS_STDOUT_TERMINAL() (_isatty(_fileno(stdout)) != 0)
#else
#include <sys/stat.h>
#include <sys/wait.h>
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

/** Galaxy step progress; stage_tag e.g. "left"/"right" in compare mode, empty for single run. */
galaxy::ProgressCallback make_galaxy_step_progress_callback(
    std::chrono::steady_clock::time_point start_wall,
    bool progress_to_terminal,
    const std::string& stage_tag) {
  return [start_wall, progress_to_terminal, stage_tag](int step, int n_steps, double sim_time) {
    auto now = std::chrono::steady_clock::now();
    double elapsed_sec =
        1e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_wall).count();
    double eta_sec =
        (step > 0 && step < n_steps) ? (elapsed_sec / step) * (n_steps - step) : 0.0;
    double pct = 100.0 * step / n_steps;
    const bool single_line_terminal_progress = progress_to_terminal && stage_tag.empty();
    if (single_line_terminal_progress) {
      // Clear the full line before redrawing to prevent remnants when text width shrinks.
      std::cout << "\r\033[2K[ ";
      if (!stage_tag.empty()) std::cout << stage_tag << " ";
      std::cout << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%] "
                << "step " << step << "/" << n_steps << ", sim t=" << std::setprecision(2) << sim_time
                << ", elapsed=" << format_elapsed(elapsed_sec) << ", eta=" << format_elapsed(eta_sec)
                << "    " << std::flush;
    } else {
      std::cout << "[ ";
      if (!stage_tag.empty()) std::cout << stage_tag << " ";
      std::cout << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%] "
                << "step " << step << "/" << n_steps << ", sim t=" << std::setprecision(2) << sim_time
                << ", elapsed=" << format_elapsed(elapsed_sec) << ", eta=" << format_elapsed(eta_sec) << "\n"
                << std::flush;
    }
  };
}

std::string sanitize_label(const std::string& in) {
  std::string out;
  out.reserve(in.size());
  for (char c : in) {
    if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9'))
      out.push_back(c);
    else
      out.push_back('_');
  }
  return out;
}

std::uint64_t fnv1a_u64_append(std::uint64_t h, const void* p, std::size_t n) {
  const unsigned char* b = static_cast<const unsigned char*>(p);
  for (std::size_t i = 0; i < n; ++i) {
    h ^= static_cast<std::uint64_t>(b[i]);
    h *= 1099511628211ull;
  }
  return h;
}

std::string state_fingerprint_hex(const galaxy::State& s) {
  std::uint64_t h = 1469598103934665603ull;
  const int n = s.n();
  h = fnv1a_u64_append(h, &n, sizeof(n));
  for (int i = 0; i < n; ++i) {
    h = fnv1a_u64_append(h, &s.x[i], sizeof(double));
    h = fnv1a_u64_append(h, &s.y[i], sizeof(double));
    h = fnv1a_u64_append(h, &s.vx[i], sizeof(double));
    h = fnv1a_u64_append(h, &s.vy[i], sizeof(double));
    h = fnv1a_u64_append(h, &s.mass[i], sizeof(double));
  }
  std::ostringstream os;
  os << std::hex << std::setfill('0') << std::setw(16) << h;
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

  bool auto_plot = false;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--plot")
      auto_plot = true;
  }

  /* First positional is the simulation mode only when it is not a long option (--key=value).
   * Otherwise mode comes from layered config (package defaults + run config). This allows
   * e.g. simulation_mode=earth_moon_benchmark in the run config with ./galaxy_sim --output_dir=... */
  bool cli_override_mode = false;
  int first_cli_config_idx = 2;
  if (argc >= 2) {
    std::string a1 = argv[1];
    if (a1.size() >= 2 && a1[0] == '-' && a1[1] == '-') {
      first_cli_config_idx = 1;
    } else {
      try {
        config.simulation_mode = galaxy::parse_mode(argv[1]);
        cli_override_mode = true;
        std::cout << "CLI override applied: simulation_mode=" << galaxy::mode_to_string(config.simulation_mode) << "\n";
      } catch (const std::exception& e) {
        std::cerr << e.what() << "\nAllowed: galaxy, earth_moon_benchmark, bh_orbit_validation, two_body_orbit (deprecated), "
                     "symmetric_pair, small_n_conservation, timestep_convergence, tpf_single_source_inspect, "
                     "tpf_symmetric_pair_inspect, tpf_two_body_sweep, tpf_weak_field_calibration, "
                     "tpf_newtonian_force_compare, tpf_diagnostic_consistency_audit, tpf_bound_orbit_sweep, "
                     "tpf_v11_weak_field_correspondence\n";
        return 1;
      }
    }
  }

  for (int i = first_cli_config_idx; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--plot") continue;
    if (a.size() < 4 || a.substr(0, 2) != "--") continue;
    std::size_t eq = a.find('=');
    if (eq == std::string::npos) {
      std::cerr << "CLI config: expected --key=value, got: " << a << "\n";
      return 1;
    }
    std::string key = a.substr(2, eq - 2);
    std::string val = a.substr(eq + 1);
    try {
      if (!galaxy::apply_config_kv(key, val, config)) {
        std::cerr << "Unknown CLI config key: " << key << "\n";
        return 1;
      }
    } catch (const std::exception& e) {
      std::cerr << "CLI config error (" << key << "): " << e.what() << "\n";
      return 1;
    }
  }

  /* Probe run config for consistency check (what did the file say?) */
  std::string run_cfg_physics = run_config_path.empty() ? "" : galaxy::probe_config_key(run_config_path, "physics_package");
  std::string run_cfg_mode_str = run_config_path.empty() ? "" : galaxy::probe_config_key(run_config_path, "simulation_mode");

  /* Hard failure: run config says TPFCore + dynamical validation mode but resolved is galaxy with no CLI override */
  if (!run_config_path.empty() && cli_override_mode == false) {
    if (run_cfg_physics == "TPFCore" &&
        (run_cfg_mode_str == "two_body_orbit" || run_cfg_mode_str == "earth_moon_benchmark" ||
         run_cfg_mode_str == "bh_orbit_validation") &&
        config.simulation_mode == galaxy::SimulationMode::galaxy) {
      std::cerr << "Config mismatch: " << run_config_path << " specifies physics_package=TPFCore and simulation_mode="
                << run_cfg_mode_str << ", but resolved simulation_mode is galaxy. Refusing to run. "
                << "Fix config precedence or run e.g.: ./galaxy_sim earth_moon_benchmark\n";
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
    if (config.simulation_mode == galaxy::SimulationMode::tpf_v11_weak_field_correspondence) {
      std::cout << "v11 correspondence audit: tpf_dynamics_mode / provisional readout are configured-inherited only "
                   "(not operative); no TPFCore particle accelerations this run.\n";
      std::cout << "  (configured: tpf_dynamics_mode=" << config.tpf_dynamics_mode
                << ", tpfcore_enable_provisional_readout=" << (config.tpfcore_enable_provisional_readout ? "true" : "false")
                << ", tpfcore_readout_mode=" << config.tpfcore_readout_mode << ")\n";
    } else {
      std::cout << "tpf_dynamics_mode: " << config.tpf_dynamics_mode << "  "
                << "tpfcore_enable_provisional_readout: " << (config.tpfcore_enable_provisional_readout ? "true" : "false")
                << "  tpfcore_readout_mode: " << config.tpfcore_readout_mode;
      if (config.tpfcore_readout_mode == "tr_coherence_readout")
        std::cout << "  theta_tt_scale: " << config.tpfcore_theta_tt_scale << "  theta_tr_scale: " << config.tpfcore_theta_tr_scale;
      std::cout << "\n";
    }
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
    std::cout << "  Provisional source ansatz: 3D Phi=-M/R on z=0 (R^2=dx^2+dy^2+eps^2), Theta=Hess_3D(Phi)\n";
    std::cout << "  Parameter roles: lambda=1/4 (fixed theory) | eps (numerical regularization) | readout (provisional experimental)\n";
    if (config.simulation_mode == galaxy::SimulationMode::tpf_v11_weak_field_correspondence) {
      std::cout << "  v11 audit: TPFCore acceleration routing not used (no legacy_readout / direct_tpf dynamics this run).\n";
    } else {
      std::cout << "  Provisional readout: " << (tpf && tpf->provisional_readout_enabled() ? "enabled" : "disabled");
      if (tpf && tpf->provisional_readout_enabled()) {
        std::cout << " (readout mode: " << tpf->readout_mode() << ", scale=" << config.tpfcore_readout_scale << " [weak-field calibrated])";
        if (config.tpfcore_readout_mode == "tr_coherence_readout")
          std::cout << " [exploratory t-r; theta_tt_scale=" << config.tpfcore_theta_tt_scale << ", theta_tr_scale=" << config.tpfcore_theta_tr_scale << "]";
        if (config.tpfcore_dump_readout_debug)
          std::cout << " [readout debug CSV enabled]";
      }
      std::cout << "\n";
    }

    if (config.tpf_dynamics_mode == "legacy_readout" && !tpf->provisional_readout_enabled()) {
      bool is_dynamical =
          (config.simulation_mode == galaxy::SimulationMode::galaxy ||
           config.simulation_mode == galaxy::SimulationMode::two_body_orbit ||
           config.simulation_mode == galaxy::SimulationMode::earth_moon_benchmark ||
           config.simulation_mode == galaxy::SimulationMode::bh_orbit_validation ||
           config.simulation_mode == galaxy::SimulationMode::symmetric_pair ||
           config.simulation_mode == galaxy::SimulationMode::small_n_conservation ||
           config.simulation_mode == galaxy::SimulationMode::timestep_convergence);
      if (is_dynamical) {
        std::cerr << "TPFCore legacy_readout does not support dynamical modes (galaxy, earth_moon_benchmark, "
                     "bh_orbit_validation, etc.) unless provisional readout is enabled.\n";
        std::cerr << "Set tpfcore_enable_provisional_readout = true, or use tpf_dynamics_mode = direct_tpf for the "
                     "future direct path, or physics_package = Newtonian, or inspection modes: "
                     "tpf_single_source_inspect, tpf_symmetric_pair_inspect.\n";
        return 1;
      }
    }
  }

  if (!ensure_dir("outputs") || !ensure_dir(config.output_dir)) {
    std::cerr << "Failed to create output dir " << config.output_dir << "\n";
    return 1;
  }
  auto write_resolved_artifacts = [](const galaxy::Config& cfg) {
    galaxy::ResolvedScenario r = galaxy::resolve_scenario(cfg);
    galaxy::write_resolved_scenario_artifacts(cfg.output_dir, r);
  };

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
    write_resolved_artifacts(config);
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
    write_resolved_artifacts(config);
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
      galaxy::init_two_body_star_around_bh(c, state);
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
    write_resolved_artifacts(config);
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
    std::cout << "Wrote " << config.output_dir << "/tpf_weak_field_calibration_comparison.csv, tpf_weak_field_calibration_comparison.txt\n";
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    write_resolved_artifacts(config);
    return 0;
  }

  if (config.simulation_mode == galaxy::SimulationMode::tpf_newtonian_force_compare) {
    if (!tpfcore) {
      std::cerr << "tpf_newtonian_force_compare requires physics_package = TPFCore.\n";
      return 1;
    }
    if (!tpfcore->provisional_readout_enabled()) {
      std::cerr << "tpf_newtonian_force_compare requires tpfcore_enable_provisional_readout = true.\n";
      return 1;
    }
    if (!ensure_dir("outputs") || !ensure_dir(config.output_dir)) {
      std::cerr << "Failed to create output dir " << config.output_dir << "\n";
      return 1;
    }
    std::cout << "Diagnostic: Newtonian vs TPF acceleration comparison (same positions/states)\n";
    tpfcore->init_from_config(config);
    galaxy::run_tpf_newtonian_force_compare(config, config.output_dir);
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    write_resolved_artifacts(config);
    return 0;
  }

  if (config.simulation_mode == galaxy::SimulationMode::tpf_diagnostic_consistency_audit) {
    if (!ensure_dir("outputs") || !ensure_dir(config.output_dir)) {
      std::cerr << "Failed to create output dir " << config.output_dir << "\n";
      return 1;
    }
    std::cout << "Diagnostic: tpf_diagnostic_consistency_audit (weak_field vs force_compare intermediates)\n";
    galaxy::run_tpf_diagnostic_consistency_audit(config, config.output_dir);
    std::cout << "Wrote " << config.output_dir << "/tpf_diagnostic_consistency_audit.csv, tpf_diagnostic_consistency_audit.txt\n";
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    write_resolved_artifacts(config);
    return 0;
  }

  if (config.simulation_mode == galaxy::SimulationMode::tpf_bound_orbit_sweep) {
    if (!tpfcore) {
      std::cerr << "tpf_bound_orbit_sweep requires physics_package = TPFCore.\n";
      return 1;
    }
    if (!tpfcore->provisional_readout_enabled()) {
      std::cerr << "tpf_bound_orbit_sweep requires tpfcore_enable_provisional_readout = true.\n";
      return 1;
    }
    if (config.tpfcore_readout_mode != "experimental_radial_r_scaling") {
      std::cerr << "tpf_bound_orbit_sweep is for experimental_radial_r_scaling only; current mode is " << config.tpfcore_readout_mode << ".\n";
      return 1;
    }
    if (!ensure_dir("outputs") || !ensure_dir(config.output_dir)) {
      std::cerr << "Failed to create output dir " << config.output_dir << "\n";
      return 1;
    }

    const double speed_ratios[] = { 0.40, 0.425, 0.45, 0.475, 0.49, 0.495, 0.499, 0.5, 0.505, 0.51, 0.525, 0.55 };
    const int n_speed = static_cast<int>(sizeof(speed_ratios) / sizeof(speed_ratios[0]));
    const double APOCENTER_TOL = 0.02;
    const double R0 = config.validation_two_body_radius;

    struct BoundSweepRow {
      double speed_ratio;
      double r_initial, r_final, r_min, r_max, radial_drift, revolutions;
      std::string trajectory_class;
      double frac_low, frac_transitional, frac_high;
      bool remained_inside_initial_apocenter;
      double returned_near_start;
      double periapsis_ratio, apocenter_ratio, eccentricity_like;
    };
    std::vector<BoundSweepRow> rows;

    std::cout << "Bound orbit sweep: experimental_radial_r_scaling, speed_ratio x " << n_speed << ", same total time\n";

    for (int i = 0; i < n_speed; ++i) {
      galaxy::Config c = config;
      c.validation_two_body_speed_ratio = speed_ratios[i];
      c.enable_star_star_gravity = false;
      physics->init_from_config(c);

      galaxy::State state;
      galaxy::init_two_body_star_around_bh(c, state);
      int n_steps = c.validation_n_steps;
      int snapshot_every = c.validation_snapshot_every;

      auto snapshots = galaxy::run_simulation(c, state, physics, n_steps, snapshot_every, nullptr, 0);

      auto traj = tpfcore->compute_trajectory_summary(snapshots);
      auto regime = tpfcore->compute_regime_summary(snapshots, c, config.output_dir);

      BoundSweepRow row = {};
      row.speed_ratio = speed_ratios[i];
      if (traj.valid) {
        row.r_initial = traj.r_initial;
        row.r_final = traj.r_final;
        row.r_min = traj.r_min;
        row.r_max = traj.r_max;
        row.radial_drift = traj.radial_drift;
        row.revolutions = traj.revolutions;
        row.trajectory_class = traj.trajectory_class;

        row.remained_inside_initial_apocenter = (traj.r_max <= traj.r_initial + APOCENTER_TOL * traj.r_initial);
        row.returned_near_start = (traj.r_initial > 1e-30) ? (std::abs(traj.r_final - traj.r_initial) / traj.r_initial) : 0.0;
        row.periapsis_ratio = (traj.r_initial > 1e-30) ? (traj.r_min / traj.r_initial) : 0.0;
        row.apocenter_ratio = (traj.r_initial > 1e-30) ? (traj.r_max / traj.r_initial) : 0.0;
        double sum_r = traj.r_max + traj.r_min;
        row.eccentricity_like = (sum_r > 1e-30) ? ((traj.r_max - traj.r_min) / sum_r) : 0.0;
      }
      if (regime.valid) {
        row.frac_low = regime.frac_low;
        row.frac_transitional = regime.frac_transitional;
        row.frac_high = regime.frac_high;
      }
      rows.push_back(row);
    }

    std::string out_dir = config.output_dir;
    std::ofstream csv(out_dir + "/tpf_bound_orbit_sweep.csv");
    if (csv) {
      csv << "speed_ratio,r_initial,r_final,r_min,r_max,radial_drift,revolutions,trajectory_class,"
          << "frac_low,frac_transitional,frac_high,remained_inside_apocenter,returned_near_start,"
          << "periapsis_ratio,apocenter_ratio,eccentricity_like\n";
      for (const auto& row : rows) {
        csv << std::scientific << row.speed_ratio << "," << row.r_initial << "," << row.r_final << ","
            << row.r_min << "," << row.r_max << "," << row.radial_drift << "," << row.revolutions << ","
            << "\"" << row.trajectory_class << "\","
            << row.frac_low << "," << row.frac_transitional << "," << row.frac_high << ","
            << (row.remained_inside_initial_apocenter ? "1" : "0") << "," << row.returned_near_start << ","
            << row.periapsis_ratio << "," << row.apocenter_ratio << "," << row.eccentricity_like << "\n";
      }
    }

    std::ofstream txt(out_dir + "/tpf_bound_orbit_sweep.txt");
    if (txt) {
      txt << "TPFCore bound-orbit sweep (experimental_radial_r_scaling, fixed closure)\n";
      txt << "Orchestration/reporting only; same formulas and total time as a single bh_orbit_validation run.\n";
      txt << "Speed list: 0.40, 0.425, 0.45, 0.475, 0.49, 0.495, 0.499, 0.5, 0.505, 0.51, 0.525, 0.55\n";
      txt << "Initial radius: " << R0 << ", n_steps: " << config.validation_n_steps << ", snapshot_every: " << config.validation_snapshot_every << "\n\n";

      txt << "Note: The legacy trajectory_class may still say \"strongly drifting\" or \"unclear\" even when the new boundness indicators (remained_inside_initial_apocenter, returned_near_start, periapsis_ratio) show a non-escaping bound-like orbit. The old bounded-band classifier is unchanged.\n\n";

      std::vector<size_t> by_inside(rows.size()), by_returned(rows.size()), by_periapsis(rows.size()), by_regime(rows.size());
      for (size_t k = 0; k < rows.size(); ++k) by_inside[k] = by_returned[k] = by_periapsis[k] = by_regime[k] = k;
      std::sort(by_inside.begin(), by_inside.end(), [&rows](size_t a, size_t b) {
        return rows[a].remained_inside_initial_apocenter != rows[b].remained_inside_initial_apocenter
          ? rows[a].remained_inside_initial_apocenter
          : (rows[a].apocenter_ratio < rows[b].apocenter_ratio);
      });
      std::sort(by_returned.begin(), by_returned.end(), [&rows](size_t a, size_t b) {
        return rows[a].returned_near_start < rows[b].returned_near_start;
      });
      std::sort(by_periapsis.begin(), by_periapsis.end(), [&rows](size_t a, size_t b) {
        return rows[a].periapsis_ratio > rows[b].periapsis_ratio;
      });
      std::sort(by_regime.begin(), by_regime.end(), [&rows](size_t a, size_t b) {
        return rows[a].frac_high < rows[b].frac_high;
      });

      txt << "--- Rank by non-escape / remained_inside_initial_apocenter (then by apocenter_ratio) ---\n";
      for (size_t j = 0; j < rows.size(); ++j) {
        size_t k = by_inside[j];
        const auto& row = rows[k];
        txt << "  " << (j + 1) << ". speed_ratio=" << std::scientific << row.speed_ratio
            << " inside_apocenter=" << (row.remained_inside_initial_apocenter ? "yes" : "no")
            << " apocenter_ratio=" << row.apocenter_ratio << " class=" << row.trajectory_class << "\n";
      }
      txt << "\n--- Rank by closeness of final radius to initial (returned_near_start, lower is better) ---\n";
      for (size_t j = 0; j < rows.size(); ++j) {
        size_t k = by_returned[j];
        const auto& row = rows[k];
        txt << "  " << (j + 1) << ". speed_ratio=" << std::scientific << row.speed_ratio
            << " returned_near_start=" << row.returned_near_start << " r_final=" << row.r_final << "\n";
      }
      txt << "\n--- Rank by shallowness of periapsis dip (periapsis_ratio, higher is better) ---\n";
      for (size_t j = 0; j < rows.size(); ++j) {
        size_t k = by_periapsis[j];
        const auto& row = rows[k];
        txt << "  " << (j + 1) << ". speed_ratio=" << std::scientific << row.speed_ratio
            << " periapsis_ratio=" << row.periapsis_ratio << " r_min=" << row.r_min << "\n";
      }
      txt << "\n--- Rank by stability of regime mix (lower frac_high preferred) ---\n";
      for (size_t j = 0; j < rows.size(); ++j) {
        size_t k = by_regime[j];
        const auto& row = rows[k];
        txt << "  " << (j + 1) << ". speed_ratio=" << std::scientific << row.speed_ratio
            << " frac_low=" << row.frac_low << " frac_trans=" << row.frac_transitional << " frac_high=" << row.frac_high << "\n";
      }
    }

    std::cout << "Wrote " << out_dir << "/tpf_bound_orbit_sweep.csv, tpf_bound_orbit_sweep.txt\n";
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    write_resolved_artifacts(config);
    return 0;
  }

  if (config.simulation_mode == galaxy::SimulationMode::tpf_v11_weak_field_correspondence) {
    if (!tpfcore) {
      std::cerr << "tpf_v11_weak_field_correspondence requires physics_package = TPFCore.\n";
      return 1;
    }
    if (config.tpf_analysis_mode != "v11_weak_field_correspondence") {
      std::cerr << "tpf_v11_weak_field_correspondence requires tpf_analysis_mode = v11_weak_field_correspondence.\n";
      return 1;
    }
    if (config.tpf_vdsg_coupling != 0.0 || !std::isfinite(config.tpf_vdsg_coupling)) {
      std::cerr << "v11 weak-field correspondence audit requires tpf_vdsg_coupling = 0 (VDSG is not in manuscript v11).\n";
      return 1;
    }
    const bool v11_em = (config.v11_weak_field_correspondence_benchmark == "earth_moon_line_of_centers");
    if (!v11_em && !(config.bh_mass > 0.0)) {
      std::cerr << "tpf_v11_weak_field_correspondence (axis_monopole) requires bh_mass > 0 (point-mass source M in kg).\n";
      return 1;
    }
    if (v11_em) {
      const double D = config.v11_em_mean_distance_m;
      if (!(config.tpfcore_probe_radius_min > 0.0) ||
          !(config.tpfcore_probe_radius_max > config.tpfcore_probe_radius_min) ||
          !(config.tpfcore_probe_radius_max < D)) {
        std::cerr << "earth_moon_line_of_centers: require 0 < tpfcore_probe_radius_min < "
                     "tpfcore_probe_radius_max < v11_em_mean_distance_m (line-of-centers interior sampling).\n";
        return 1;
      }
    } else if (!(config.tpfcore_probe_radius_min > 0.0) ||
               !(config.tpfcore_probe_radius_max > config.tpfcore_probe_radius_min)) {
      std::cerr << "tpf_v11_weak_field_correspondence requires tpfcore_probe_radius_min > 0 and "
                   "tpfcore_probe_radius_max > tpfcore_probe_radius_min (positive z axis samples).\n";
      return 1;
    }
    std::cout << "Audit mode: tpf_v11_weak_field_correspondence (manuscript v11 static weak-field correspondence only)\n";
    std::cout << "  Benchmark: " << config.v11_weak_field_correspondence_benchmark << "\n";
    std::cout << "  No particle integration; Delta C_mu_nu omitted; VDSG off.\n";
    galaxy::run_v11_weak_field_correspondence_audit(config, config.output_dir);
    if (v11_em) {
      std::cout << "Wrote " << config.output_dir << "/tpf_v11_earth_moon_line_correspondence_benchmark.csv\n";
      std::cout << "Wrote " << config.output_dir << "/tpf_v11_earth_moon_line_correspondence_benchmark_summary.txt\n";
      std::cout << "Wrote " << config.output_dir << "/tpf_v11_earth_moon_line_correspondence_benchmark.gnu\n";
    } else {
      std::cout << "Wrote " << config.output_dir << "/tpf_v11_weak_field_correspondence_profile.csv\n";
      std::cout << "Wrote " << config.output_dir << "/tpf_v11_weak_field_correspondence_summary.txt\n";
    }
    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, 0, 0, 0, run_config_path, package_defaults_path);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    galaxy::write_render_manifest(config.output_dir, config, 0, 0, 0, nullptr);
    std::cout << "Wrote " << config.output_dir << "/render_manifest.json, render_manifest.txt\n";
    write_resolved_artifacts(config);
    return 0;
  }

  switch (config.simulation_mode) {
    case galaxy::SimulationMode::tpf_single_source_inspect:
    case galaxy::SimulationMode::tpf_symmetric_pair_inspect:
    case galaxy::SimulationMode::tpf_two_body_sweep:
    case galaxy::SimulationMode::tpf_weak_field_calibration:
    case galaxy::SimulationMode::tpf_newtonian_force_compare:
    case galaxy::SimulationMode::tpf_diagnostic_consistency_audit:
    case galaxy::SimulationMode::tpf_bound_orbit_sweep:
    case galaxy::SimulationMode::tpf_v11_weak_field_correspondence:
      std::cerr << "Internal error: inspection/utility/sweep/calibration/audit modes should have returned earlier.\n";
      return 1;
    default:
      break;
  }

  const galaxy::Config configured_after_layering = config;
  galaxy::ResolvedScenario resolved = galaxy::resolve_scenario(configured_after_layering);
  config = resolved.config;
  galaxy::State state = resolved.initial_state;
  int n_steps = resolved.effective_n_steps;
  int snapshot_every = resolved.effective_snapshot_every;

  if (config.simulation_mode == galaxy::SimulationMode::galaxy) {
    galaxy::write_galaxy_init_diagnostics(config.output_dir, state, config,
                                          galaxy::last_galaxy_init_audit());
    std::cout << "Galaxy IC: template=" << galaxy::last_galaxy_init_audit().template_name
              << ", seed=" << galaxy::last_galaxy_init_audit().seed;
    if (galaxy::last_galaxy_init_audit().used_new_state_noise)
      std::cout << ", noise=new (pos/angle/mag)";
    else if (galaxy::last_galaxy_init_audit().used_legacy_velocity_noise)
      std::cout << ", noise=legacy velocity_noise";
    else
      std::cout << ", noise=none";
    std::cout << "\n";
    std::cout << "Wrote " << config.output_dir << "/galaxy_init_diagnostics.txt\n";
    if (!config.save_snapshots)
      std::cout << "Wrote " << config.output_dir << "/galaxy_init_snapshot.csv (snapshots disabled)\n";
    std::cout << "Running galaxy: n_stars=" << config.n_stars
              << ", n_steps=" << n_steps
              << ", dt=" << config.dt
              << ", sim_time=" << (n_steps * config.dt) << "\n";
  } else if (config.simulation_mode == galaxy::SimulationMode::two_body_orbit ||
             config.simulation_mode == galaxy::SimulationMode::earth_moon_benchmark) {
    std::cout << "Earth–Moon benchmark (SI), n=" << state.n()
              << " (pairwise gravity; bh_mass cleared by scenario resolver)\n";
  } else if (config.simulation_mode == galaxy::SimulationMode::bh_orbit_validation) {
    std::cout << "BH orbit validation: n=" << state.n() << " star, r0=" << config.validation_two_body_radius
              << " m, speed_ratio=" << config.validation_two_body_speed_ratio
              << " (star–star gravity off by scenario resolver)\n";
  } else if (config.simulation_mode == galaxy::SimulationMode::symmetric_pair) {
    std::cout << "Symmetric pair: a=" << config.validation_symmetric_separation
              << " include_bh=" << config.validation_symmetric_include_bh << "\n";
  } else if (config.simulation_mode == galaxy::SimulationMode::small_n_conservation) {
    std::cout << "Small-N: n=" << state.n() << "\n";
  }

  switch (config.simulation_mode) {
    case galaxy::SimulationMode::timestep_convergence: {
      // Run two_body at dt, dt/2, dt/4 (same total time)
      double total_time = config.validation_n_steps * config.dt;
      std::vector<double> dts = {config.dt, config.dt / 2, config.dt / 4};
      config.enable_star_star_gravity = false;

      std::cout << "Timestep convergence (star_around_bh IC), total_time=" << total_time << "\n";
      std::ofstream summary(config.output_dir + "/validation_timestep_convergence.txt");
      summary << "Timestep convergence (star_around_bh IC; same initializer as bh_orbit_validation)\n";
      summary << "dt\tfinal_x\tfinal_y\tfinal_r\tL_z\tE_drift\n";

      double E0 = 0;
      galaxy::State state0;
      galaxy::init_two_body_star_around_bh(config, state0);
      E0 = galaxy::compute_kinetic_energy(state0) + physics->compute_potential_energy(state0, config.bh_mass, config.softening, config.enable_star_star_gravity);

      for (double dt : dts) {
        int steps = static_cast<int>(std::round(total_time / dt));
        galaxy::Config c2 = config;
        c2.dt = dt;
        galaxy::State s0;
        galaxy::init_two_body_star_around_bh(c2, s0);
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
      write_resolved_artifacts(config);
      return 0;
    }
    default:
      break;
  }

  write_resolved_artifacts(config);

  const bool compare_mode_requested =
      (config.simulation_mode == galaxy::SimulationMode::galaxy &&
       !config.physics_package_compare.empty());
  const bool compare_same_package =
      (compare_mode_requested && config.physics_package_compare == config.physics_package);
  if (compare_same_package) {
    std::cout << "Warning: physics_package_compare equals physics_package ("
              << config.physics_package << "); falling back to single-package run.\n";
  }

  if (compare_mode_requested && !compare_same_package) {
    const std::string compare_parent_dir = config.output_dir;
    const std::string left_dir = compare_parent_dir + "/left_" + sanitize_label(config.physics_package);
    const std::string right_dir = compare_parent_dir + "/right_" + sanitize_label(config.physics_package_compare);
    if (!ensure_dir("outputs") || !ensure_dir(compare_parent_dir) ||
        !ensure_dir(left_dir) || !ensure_dir(right_dir)) {
      std::cerr << "Failed to create compare output directories under " << compare_parent_dir << "\n";
      return 1;
    }
    std::cout << "Compare run directories (created):\n  " << left_dir << "\n  " << right_dir << "\n";

    galaxy::Config left_cfg = config;
    left_cfg.output_dir = left_dir;
    left_cfg.run_id = config.run_id + "_left";
    left_cfg.physics_package = config.physics_package;

    galaxy::Config right_cfg = config;
    right_cfg.output_dir = right_dir;
    right_cfg.run_id = config.run_id + "_right";
    right_cfg.physics_package = config.physics_package_compare;
    right_cfg.physics_package_compare.clear();

    galaxy::PhysicsPackage* left_physics_probe = galaxy::get_physics_package(left_cfg.physics_package);
    galaxy::PhysicsPackage* right_physics_probe = galaxy::get_physics_package(right_cfg.physics_package);
    if (!left_physics_probe || !right_physics_probe) {
      std::cerr << "Compare mode failed: unknown package(s): left=" << left_cfg.physics_package
                << ", right=" << right_cfg.physics_package << "\n";
      return 1;
    }

    const std::string ic_hash = state_fingerprint_hex(state);

    std::cout << "Compare mode: left=" << left_cfg.physics_package
              << "  right=" << right_cfg.physics_package << "\n";
    std::cout << "Shared IC fingerprint (fnv1a64): " << ic_hash << "\n";

    auto write_side_outputs =
        [&](const galaxy::Config& side_cfg,
            galaxy::PhysicsPackage* side_physics,
            const std::vector<galaxy::Snapshot>& side_snaps,
            const std::string& side_defaults_path) {
          const bool cooling_active =
              (side_cfg.physics_package == "TPFCore" && side_cfg.tpf_cooling_fraction > 0.0);
          const int cooling_steps = cooling_active
              ? std::min(n_steps, std::max(0, static_cast<int>(n_steps * side_cfg.tpf_cooling_fraction)))
              : 0;
          galaxy::CoolingAuditInfo cooling_audit;
          cooling_audit.cooling_active = cooling_active;
          cooling_audit.cooling_steps = cooling_steps;
          cooling_audit.cooling_end_step = std::max(0, cooling_steps - 1);
          if (!side_snaps.empty()) {
            const galaxy::Snapshot* first_saved = nullptr;
            for (const auto& snap : side_snaps) {
              if (snap.step > 0) {
                first_saved = &snap;
                break;
              }
            }
            if (!first_saved) first_saved = &side_snaps.front();
            cooling_audit.first_saved_snapshot_step = first_saved->step;
            cooling_audit.first_saved_snapshot_time = first_saved->time;
          }

          const galaxy::AccelPipelineStats* tpf_pipeline_stats = nullptr;
          if (side_cfg.physics_package == "TPFCore") {
            if (auto* tpf_pkg = dynamic_cast<galaxy::TPFCorePackage*>(side_physics)) {
              if (tpf_pkg->provisional_readout_enabled() && tpf_pkg->last_accel_pipeline_stats().valid)
                tpf_pipeline_stats = &tpf_pkg->last_accel_pipeline_stats();
            }
          }

          if (side_cfg.save_run_info) {
            galaxy::write_run_info(side_cfg.output_dir, side_cfg, n_steps, static_cast<int>(side_snaps.size()),
                                   state.n(), run_config_path, side_defaults_path,
                                   nullptr, nullptr,
                                   &galaxy::last_galaxy_init_audit(), &cooling_audit, tpf_pipeline_stats);
            galaxy::write_render_manifest(side_cfg.output_dir, side_cfg, n_steps,
                                          static_cast<int>(side_snaps.size()), state.n(),
                                          &galaxy::last_galaxy_init_audit());
          }
          if (side_cfg.save_snapshots)
            galaxy::write_snapshots(side_cfg.output_dir, side_snaps);
          write_resolved_artifacts(side_cfg);
          galaxy::write_galaxy_init_diagnostics(side_cfg.output_dir, state, side_cfg,
                                                galaxy::last_galaxy_init_audit());
        };
    auto run_compare_side =
        [&](const char* side_tag, const galaxy::Config& side_cfg) -> int {
          galaxy::PhysicsPackage* side_physics = galaxy::get_physics_package(side_cfg.physics_package);
          if (!side_physics) {
            std::cerr << "Compare mode failed in " << side_tag << " side: unknown package "
                      << side_cfg.physics_package << "\n";
            return 2;
          }
          side_physics->init_from_config(side_cfg);
          galaxy::State side_state = state;
          int compare_progress_interval = 0;
          if (n_steps > 0) {
            compare_progress_interval = std::max(1, std::min(1000, n_steps / 100));
          }
          const bool progress_to_terminal = IS_STDOUT_TERMINAL();
          auto start_wall = std::chrono::steady_clock::now();
          galaxy::ProgressCallback side_progress =
              make_galaxy_step_progress_callback(start_wall, progress_to_terminal, side_tag);
          auto side_snapshots = galaxy::run_simulation(side_cfg, side_state, side_physics,
                                                       n_steps, snapshot_every, side_progress,
                                                       compare_progress_interval);
          if (progress_to_terminal && n_steps > 0) std::cout << "\n";
          write_side_outputs(side_cfg, side_physics, side_snapshots,
                             galaxy::find_package_defaults_path(side_cfg.physics_package));
          return 0;
        };

    const bool compare_parallel_enabled =
        galaxy::should_run_compare_parallel(compare_mode_requested, compare_same_package,
#ifdef _WIN32
                                            false
#else
                                            true
#endif
        );

    if (n_steps > 0) {
      std::cout << "Running compare simulations (" << n_steps << " steps each)."
                << (compare_parallel_enabled ? " [parallel process mode]\n" : " [sequential mode]\n")
                << std::flush;
    }

    if (compare_parallel_enabled) {
#ifdef _WIN32
      std::cerr << "Process-parallel compare is unavailable on this platform; using sequential mode.\n";
      if (run_compare_side("left", left_cfg) != 0) return 1;
      if (run_compare_side("right", right_cfg) != 0) return 1;
#else
      const std::string left_log = compare_parent_dir + "/left_run.log";
      const std::string right_log = compare_parent_dir + "/right_run.log";
      const bool show_live_compare_progress = IS_STDOUT_TERMINAL();
      if (show_live_compare_progress) {
        std::cout << "Parallel compare enabled; streaming child progress to terminal.\n";
      } else {
        std::cout << "Parallel compare enabled; child logs:\n  " << left_log << "\n  " << right_log << "\n";
      }

      pid_t left_pid = fork();
      if (left_pid == 0) {
        if (!show_live_compare_progress) {
          FILE* lf = std::freopen(left_log.c_str(), "w", stdout);
          FILE* le = std::freopen(left_log.c_str(), "a", stderr);
          if (!lf || !le) _exit(90);
        }
        const int rc = run_compare_side("left", left_cfg);
        _exit(rc);
      }
      if (left_pid < 0) {
        std::perror("fork(left)");
        return 1;
      }

      pid_t right_pid = fork();
      if (right_pid == 0) {
        if (!show_live_compare_progress) {
          FILE* rf = std::freopen(right_log.c_str(), "w", stdout);
          FILE* re = std::freopen(right_log.c_str(), "a", stderr);
          if (!rf || !re) _exit(91);
        }
        const int rc = run_compare_side("right", right_cfg);
        _exit(rc);
      }
      if (right_pid < 0) {
        std::perror("fork(right)");
        int left_status = 0;
        (void)waitpid(left_pid, &left_status, 0);
        return 1;
      }

      std::cout << "Started compare children: left pid=" << static_cast<long>(left_pid)
                << ", right pid=" << static_cast<long>(right_pid) << "\n";
      int left_status = 0;
      int right_status = 0;
      if (waitpid(left_pid, &left_status, 0) < 0) {
        std::perror("waitpid(left)");
        return 1;
      }
      if (waitpid(right_pid, &right_status, 0) < 0) {
        std::perror("waitpid(right)");
        return 1;
      }

      std::string left_err;
      std::string right_err;
      const bool left_ok = galaxy::child_exit_ok(left_status, "left", &left_err);
      const bool right_ok = galaxy::child_exit_ok(right_status, "right", &right_err);
      if (!left_ok || !right_ok) {
        if (!left_ok) std::cerr << "Compare parallel failure: " << left_err << "\n";
        if (!right_ok) std::cerr << "Compare parallel failure: " << right_err << "\n";
        std::cerr << "Compare run failed; manifests/plot skipped.\n";
        return 1;
      }
#endif
    } else {
      if (run_compare_side("left", left_cfg) != 0) return 1;
      if (run_compare_side("right", right_cfg) != 0) return 1;
    }

    const galaxy::GitProvenance gp = galaxy::resolve_git_provenance();
    {
      std::ofstream jf(compare_parent_dir + "/compare_manifest.json");
      if (jf) {
        jf << "{\n"
           << "  \"schema\": \"galaxy_compare_manifest_v1\",\n"
           << "  \"compare_run_id\": \"" << config.run_id << "\",\n"
           << "  \"primary_package\": \"" << left_cfg.physics_package << "\",\n"
           << "  \"compare_package\": \"" << right_cfg.physics_package << "\",\n"
           << "  \"left_dir\": \"" << left_dir << "\",\n"
           << "  \"right_dir\": \"" << right_dir << "\",\n"
           << "  \"ic_seed\": " << left_cfg.galaxy_init_seed << ",\n"
           << "  \"ic_fingerprint_fnv1a64\": \"" << ic_hash << "\",\n"
           << "  \"git_commit_full\": \"" << gp.git_commit_full << "\",\n"
           << "  \"git_commit_short\": \"" << gp.git_commit_short << "\",\n"
           << "  \"git_branch\": \"" << gp.git_branch << "\",\n"
           << "  \"git_tag\": \"" << gp.git_tag << "\",\n"
           << "  \"git_dirty\": " << (gp.git_dirty ? "true" : "false") << ",\n"
           << "  \"code_version_label\": \"" << gp.code_version_label << "\"\n"
           << "}\n";
      }
      std::ofstream tf(compare_parent_dir + "/compare_manifest.txt");
      if (tf) {
        tf << "compare_run_id\t" << config.run_id << "\n";
        tf << "primary_package\t" << left_cfg.physics_package << "\n";
        tf << "compare_package\t" << right_cfg.physics_package << "\n";
        tf << "left_dir\t" << left_dir << "\n";
        tf << "right_dir\t" << right_dir << "\n";
        tf << "ic_seed\t" << left_cfg.galaxy_init_seed << "\n";
        tf << "ic_fingerprint_fnv1a64\t" << ic_hash << "\n";
        tf << "git_commit_full\t" << gp.git_commit_full << "\n";
        tf << "git_commit_short\t" << gp.git_commit_short << "\n";
        tf << "git_branch\t" << gp.git_branch << "\n";
        tf << "git_tag\t" << gp.git_tag << "\n";
        tf << "git_dirty\t" << (gp.git_dirty ? 1 : 0) << "\n";
        tf << "code_version_label\t" << gp.code_version_label << "\n";
      }
    }

    std::cout << "Wrote compare manifests in " << compare_parent_dir << "\n";
    std::cout << "Left output: " << left_dir << "\n";
    std::cout << "Right output: " << right_dir << "\n";

    if (auto_plot) {
      std::cout << "Rendering compare figures (plot_cpp_compare.py)...\n" << std::flush;
      // Run from cpp_sim cwd: script lives at repo_root/plot_cpp_compare.py
      const std::string dev_py = "../dev/bin/python3";
      const bool dev_py_exists = static_cast<bool>(std::ifstream(dev_py).good());
      const std::string py = dev_py_exists ? dev_py : "python3";
      std::string cmd = py + " ../plot_cpp_compare.py " + compare_parent_dir;
      int ret = std::system(cmd.c_str());
      if (ret != 0) {
        std::cerr << "Warning: compare renderer returned non-zero exit code. Command was:\n  " << cmd << "\n";
      } else {
        std::cout << "Compare render finished (see galaxy_initial_compare.png, galaxy_final_compare.png in "
                  << compare_parent_dir << ").\n";
      }
    } else {
      std::cout
          << "\nCompare side-by-side PNGs/animation were not generated (run without --plot).\n"
          << "From the cpp_sim directory:\n  python3 ../plot_cpp_compare.py " << compare_parent_dir << "\n"
          << "From the repository root:\n  python3 plot_cpp_compare.py cpp_sim/" << compare_parent_dir << "\n"
          << "Or re-run with --plot to render automatically.\n";
    }
    return 0;
  }

  // Progress reporting: only for galaxy (long runs); in-place line when stdout is a terminal
  int progress_interval = 0;
  galaxy::ProgressCallback progress_callback;
  bool progress_to_terminal = false;
  if (config.simulation_mode == galaxy::SimulationMode::galaxy && n_steps > 0) {
    progress_interval = std::max(1, std::min(1000, n_steps / 100));
    auto start_wall = std::chrono::steady_clock::now();
    progress_to_terminal = IS_STDOUT_TERMINAL();
    progress_callback = make_galaxy_step_progress_callback(start_wall, progress_to_terminal, "");
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

  {
    const bool cooling_active =
        (config.physics_package == "TPFCore" && config.tpf_cooling_fraction > 0.0);
    const int cooling_steps = cooling_active
        ? std::min(n_steps, std::max(0, static_cast<int>(n_steps * config.tpf_cooling_fraction)))
        : 0;
    galaxy::CoolingAuditInfo cooling_audit;
    cooling_audit.cooling_active = cooling_active;
    cooling_audit.cooling_steps = cooling_steps;
    cooling_audit.cooling_end_step = std::max(0, cooling_steps - 1);
    if (!snapshots.empty()) {
      const galaxy::Snapshot* first_saved = nullptr;
      for (const auto& snap : snapshots) {
        if (snap.step > 0) {
          first_saved = &snap;
          break;
        }
      }
      if (!first_saved) first_saved = &snapshots.front();
      cooling_audit.first_saved_snapshot_step = first_saved->step;
      cooling_audit.first_saved_snapshot_time = first_saved->time;
    }

    const galaxy::AccelPipelineStats* tpf_pipeline_stats = nullptr;
    if (config.physics_package == "TPFCore") {
      if (auto* tpf_pkg = dynamic_cast<galaxy::TPFCorePackage*>(physics)) {
        if (tpf_pkg->provisional_readout_enabled() && tpf_pkg->last_accel_pipeline_stats().valid)
          tpf_pipeline_stats = &tpf_pkg->last_accel_pipeline_stats();
      }
    }

    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, n_steps, static_cast<int>(snapshots.size()), state.n(),
                             run_config_path, package_defaults_path,
                             &configured_after_layering,
                             &resolved,
                             config.simulation_mode == galaxy::SimulationMode::galaxy
                                 ? &galaxy::last_galaxy_init_audit()
                                 : nullptr,
                             &cooling_audit,
                             tpf_pipeline_stats);
      std::cout << "Wrote " << config.output_dir << "/run_info.txt\n";
    }
    if (config.save_run_info && config.simulation_mode == galaxy::SimulationMode::galaxy) {
      galaxy::write_render_manifest(config.output_dir, config, n_steps, static_cast<int>(snapshots.size()),
                                  state.n(), &galaxy::last_galaxy_init_audit());
      std::cout << "Wrote " << config.output_dir << "/render_manifest.json, render_manifest.txt\n";
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
        if (config.physics_package == "TPFCore" && snapshots[0].state.n() == 1 &&
            (config.tpfcore_readout_mode == "tr_coherence_readout" || config.tpfcore_readout_mode == "experimental_radial_r_scaling"))
          std::cout << "Wrote " << config.output_dir << "/tpf_closure_diagnostics.csv, tpf_closure_diagnostics.txt\n";
        if (config.simulation_mode == galaxy::SimulationMode::bh_orbit_validation && snapshots[0].state.n() == 1) {
          tpf->write_step0_orbit_audit(snapshots, config, config.output_dir);
          std::cout << "Wrote " << config.output_dir << "/tpf_step0_orbit_audit.txt\n";
        }
        if (config.tpfcore_live_orbit_force_audit) {
          tpf->write_live_orbit_force_audit(snapshots, config, config.output_dir);
          std::cout << "Wrote " << config.output_dir << "/tpf_live_orbit_force_audit.csv, tpf_live_orbit_force_audit.txt\n";
        }
        if (config.tpf_accel_pipeline_diagnostics_csv) {
          tpf->write_accel_pipeline_diagnostics(snapshots, config, config.output_dir);
          std::cout << "Wrote " << config.output_dir << "/tpf_accel_pipeline_diagnostics.csv\n";
        }
      }
      if (tpf && config.tpf_dynamics_mode == "direct_tpf" && !snapshots.empty()) {
        tpf->write_step0_orbit_audit(snapshots, config, config.output_dir);
        std::cout << "Wrote " << config.output_dir
                  << "/direct_tpf_step0_raw_accel_audit.csv, direct_tpf_step0_raw_accel_summary.txt\n";
      }
    }
  }

  std::cout << "Snapshots: " << snapshots.size() << "\n";
  std::cout << "Output directory: " << config.output_dir << "\n";

  std::string output_dir = config.output_dir;
  if (auto_plot) {
    std::cout.flush();
    std::cerr.flush();
    std::cout << "Simulation complete. Rendering animation..." << std::endl;
    std::string cmd = "cd .. && ./dev/bin/python plot_cpp_run.py cpp_sim/" + output_dir;
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
      std::cerr << "Warning: Python rendering script returned non-zero exit code." << std::endl;
    }
  }

  return 0;
}
