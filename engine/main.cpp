#include "config.hpp"
#include "physics/TPFCore/tpf_core_config.hpp"
#include "compare_orchestration.hpp"
#include "force_compare.hpp"
#include "galaxy_init.hpp"
#include "git_provenance.hpp"
#include "init_conditions.hpp"
#include "render_audit.hpp"
#include "resolved_scenario.hpp"
#include "output.hpp"
#include "preflight.hpp"
#include "progress_time.hpp"
#include "physics/physics_package.hpp"
#include "physics/Newtonian/newtonian.hpp"
#include "simulation.hpp"
#include "types.hpp"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <cstdint>
#include <cstdio>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#include <io.h>
#define MKDIR(path, mode) _mkdir(path)
#define IS_STDOUT_TERMINAL() (_isatty(_fileno(stdout)) != 0)
#define IS_STDIN_TERMINAL() (_isatty(_fileno(stdin)) != 0)
#else
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#define MKDIR(path, mode) mkdir(path, mode)
#define IS_STDOUT_TERMINAL() (isatty(STDOUT_FILENO) != 0)
#define IS_STDIN_TERMINAL() (isatty(STDIN_FILENO) != 0)
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

bool ensure_dir_recursive(const std::string& path) {
  if (path.empty()) return false;
  std::string cur;
  std::size_t i = 0;
  if (path[0] == '/') {
    cur = '/';
    i = 1;
  }
  while (i <= path.size()) {
    const std::size_t j = path.find('/', i);
    const std::string part = path.substr(i, j == std::string::npos ? std::string::npos : j - i);
    if (!part.empty() && part != ".") {
      if (!cur.empty() && cur[cur.size() - 1] != '/') cur += '/';
      cur += part;
      if (!ensure_dir(cur)) return false;
    }
    if (j == std::string::npos) break;
    i = j + 1;
  }
  return true;
}

std::string shell_single_quote(const std::string& raw) {
  std::string out = "'";
  for (std::size_t i = 0; i < raw.size(); ++i) {
    if (raw[i] == '\'') {
      out += "'\"'\"'";
    } else {
      out += raw[i];
    }
  }
  out += "'";
  return out;
}

bool file_exists(const std::string& path) {
  std::ifstream f(path.c_str());
  return static_cast<bool>(f);
}

std::vector<std::string> find_compare_side_by_side_pngs(const std::string& dir,
                                                         const std::string& stage) {
  std::vector<std::string> out;
  DIR* dp = opendir(dir.c_str());
  if (!dp) return out;
  const std::string stage_token = stage + "_side_by_side";
  const std::string ext = ".png";
  for (dirent* ent = readdir(dp); ent != nullptr; ent = readdir(dp)) {
    const std::string name(ent->d_name);
    if (name.find("galaxy_compare__") != std::string::npos &&
        name.find(stage_token) != std::string::npos &&
        name.size() >= ext.size() &&
        name.compare(name.size() - ext.size(), ext.size(), ext) == 0) {
      out.push_back(name);
    }
  }
  closedir(dp);
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<std::string> existing_tpf_4d_static_plot_pngs(const std::string& output_dir) {
  const char* candidates[] = {
      "tpf_4d_static_residual_xy_normalized_residual.png",
      "tpf_4d_static_residual_xy_residual_spatial_norm.png",
      "tpf_4d_static_residual_xy_xi_spatial_norm.png",
      "tpf_4d_static_residual_xy_invariant_I.png",
      "tpf_4d_static_residual_xz_normalized_residual.png",
      "tpf_4d_static_residual_xz_residual_spatial_norm.png",
      "tpf_4d_static_residual_xz_xi_spatial_norm.png",
      "tpf_4d_static_residual_xz_invariant_I.png",
      "tpf_4d_static_residual_yz_normalized_residual.png",
      "tpf_4d_static_residual_yz_residual_spatial_norm.png",
      "tpf_4d_static_residual_yz_xi_spatial_norm.png",
      "tpf_4d_static_residual_yz_invariant_I.png",
      "tpf_4d_static_residual_xy_xi_quiver.png",
      "tpf_4d_static_residual_xz_xi_quiver.png",
      "tpf_4d_static_residual_yz_xi_quiver.png",
  };
  std::vector<std::string> found;
  for (std::size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); ++i) {
    const std::string full = output_dir + "/" + candidates[i];
    if (file_exists(full)) found.push_back(candidates[i]);
  }
  return found;
}

void write_resolved_artifacts(const galaxy::Config& config) {
  const galaxy::ResolvedScenario resolved = galaxy::resolve_scenario(config);
  galaxy::write_resolved_scenario_artifacts(config.output_dir, resolved);
}

void finalize_utility_mode_run(const galaxy::Config& config,
                               const std::string& run_config_path,
                               const std::string& package_defaults_path) {
  const galaxy::ResolvedScenario resolved = galaxy::resolve_scenario(config);
  galaxy::write_resolved_scenario_artifacts(config.output_dir, resolved);
  if (config.save_run_info) {
    galaxy::write_run_info(config.output_dir, config, 0, 0, resolved.initial_state.n(),
                           run_config_path, package_defaults_path, &config, &resolved);
  }
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

bool compare_dual_line_progress_enabled(bool progress_to_terminal, const std::string& stage_tag) {
  return progress_to_terminal &&
         std::getenv("GALAXY_COMPARE_DUAL_LINE_PROGRESS") != nullptr &&
         (stage_tag == "left" || stage_tag == "right");
}

/** Galaxy step progress; stage_tag e.g. "left"/"right" in compare mode, empty for single run. */
galaxy::ProgressCallback make_galaxy_step_progress_callback(
    std::chrono::steady_clock::time_point start_wall,
    bool progress_to_terminal,
    const std::string& stage_tag,
    const std::string& display_time_unit) {
  return [start_wall, progress_to_terminal, stage_tag, display_time_unit](int step, int n_steps, double sim_time) {
    auto now = std::chrono::steady_clock::now();
    double elapsed_sec =
        1e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_wall).count();
    double eta_sec =
        (step > 0 && step < n_steps) ? (elapsed_sec / step) * (n_steps - step) : 0.0;
    double pct = 100.0 * step / n_steps;
    const bool compare_dual_line_mode =
        compare_dual_line_progress_enabled(progress_to_terminal, stage_tag);
    const std::string formatted_sim_time =
        galaxy::format_sim_time_for_progress(sim_time, display_time_unit);
    const bool single_line_terminal_progress = progress_to_terminal;
    if (single_line_terminal_progress) {
      // Clear the full line before redrawing to prevent remnants when text width shrinks.
      if (compare_dual_line_mode) {
        if (stage_tag == "left")
          std::cout << "\033[2A\r\033[2K[ ";
        else
          std::cout << "\033[1A\r\033[2K[ ";
      } else {
        std::cout << "\r\033[2K[ ";
      }
      if (!stage_tag.empty()) std::cout << stage_tag << " ";
      std::cout << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%] "
                << "step " << step << "/" << n_steps << ", sim t=" << formatted_sim_time
                << ", elapsed=" << format_elapsed(elapsed_sec) << ", eta=" << format_elapsed(eta_sec);
      if (compare_dual_line_mode) {
        if (stage_tag == "left")
          std::cout << "    \033[2B" << std::flush;
        else
          std::cout << "    \033[1B" << std::flush;
      } else {
        std::cout << "    " << std::flush;
      }
    } else {
      std::cout << "[ ";
      if (!stage_tag.empty()) std::cout << stage_tag << " ";
      std::cout << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%] "
                << "step " << step << "/" << n_steps << ", sim t=" << formatted_sim_time
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

void write_galaxy_step0_accel_audit(const galaxy::Config& config,
                                    galaxy::PhysicsPackage* physics,
                                    const galaxy::State& state) {
  if (config.simulation_mode != galaxy::SimulationMode::galaxy) return;
  if (!physics) return;

  const auto tpf_cfg = build_tpfcore_config(config);
  std::ofstream csv(config.output_dir + "/galaxy_step0_accel_audit.csv");
  if (!csv) return;
  csv << std::scientific << std::setprecision(17);
  csv << "particle_index,radius,vx0,vy0,"
      << "active_route,applied_acceleration,decomposition_note,"
      << "ax_reference_direct_tpf,ay_reference_direct_tpf,ax_reference_vdsg,ay_reference_vdsg,"
      << "ax_total,ay_total,ax_bh_only,ay_bh_only,ax_star_star_only,ay_star_star_only\n";

  std::vector<double> ax_total, ay_total;
  physics->compute_accelerations(state, config.bh_mass, config.softening, config.enable_star_star_gravity, ax_total, ay_total);

  std::vector<double> ax_bh_only, ay_bh_only;
  physics->compute_accelerations(state, config.bh_mass, config.softening, false, ax_bh_only, ay_bh_only);

  std::vector<double> ax_star_only, ay_star_only;
  physics->compute_accelerations(state, 0.0, config.softening, config.enable_star_star_gravity, ax_star_only, ay_star_only);

  std::vector<double> ax_reference_direct_tpf(state.n(), std::numeric_limits<double>::quiet_NaN());
  std::vector<double> ay_reference_direct_tpf(state.n(), std::numeric_limits<double>::quiet_NaN());
  std::vector<double> ax_reference_vdsg(state.n(), std::numeric_limits<double>::quiet_NaN());
  std::vector<double> ay_reference_vdsg(state.n(), std::numeric_limits<double>::quiet_NaN());
  std::string active_route = config.physics_package;
  std::string applied_accel_formula = "n/a";
  std::string decomposition_note = "reference decomposition unavailable";

  if (config.physics_package == "TPFCore") {
    active_route = tpf_cfg.tpf_dynamics_mode;
    if (tpf_cfg.tpf_dynamics_mode == "tpf_xi_theta_v1") {
      applied_accel_formula = "a=-K_xi*Xi_total_spatial";
      decomposition_note =
          "v1 route: Xi_total-driven motion; Theta=grad(Xi_total) is diagnostic-only";
    } else {
      applied_accel_formula = "unsupported non-v1 dynamics mode";
      decomposition_note = "non-v1 decomposition disabled on this branch";

      galaxy::Config cfg_direct = config;
      // Step-0 reference decomposition branch:
      //  1) force direct_tpf dynamics,
      //  2) force tensor_radial_projection readout implementation,
      //  3) zero VDSG for the direct/baseline branch,
      //  4) compute additive VDSG = (total-with-configured-VDSG) - (direct baseline).
      cfg_direct.tpf_dynamics_mode = "tpf_xi_theta_v1";
      cfg_direct.tpfcore_enable_provisional_readout = false;
      cfg_direct.tpfcore_readout_mode = "tensor_radial_projection";
      cfg_direct.tpfcore_readout_scale = 1.0;
      cfg_direct.tpfcore_theta_tt_scale = 1.0;
      cfg_direct.tpfcore_theta_tr_scale = 1.0;
      cfg_direct.tpf_global_accel_shunt_enable = false;
      cfg_direct.tpf_cooling_fraction = 0.0;
      cfg_direct.tpf_vdsg_coupling = 0.0;

      galaxy::Config cfg_total = cfg_direct;
      cfg_total.tpf_vdsg_coupling = 0.0;

      galaxy::PhysicsPackage* tpf_direct = galaxy::get_physics_package("TPFCore");
      if (tpf_direct) {
        tpf_direct->init_from_config(cfg_direct);
        tpf_direct->compute_accelerations(state, config.bh_mass, config.softening,
                                          config.enable_star_star_gravity, ax_reference_direct_tpf,
                                          ay_reference_direct_tpf);
      }

      std::vector<double> ax_total_direct, ay_total_direct;
      galaxy::PhysicsPackage* tpf_total = galaxy::get_physics_package("TPFCore");
      if (tpf_total) {
        tpf_total->init_from_config(cfg_total);
        tpf_total->compute_accelerations(state, config.bh_mass, config.softening,
                                         config.enable_star_star_gravity, ax_total_direct, ay_total_direct);
      }
      for (int i = 0; i < state.n(); ++i) {
        ax_reference_vdsg[static_cast<size_t>(i)] =
            ax_total_direct[static_cast<size_t>(i)] - ax_reference_direct_tpf[static_cast<size_t>(i)];
        ay_reference_vdsg[static_cast<size_t>(i)] =
            ay_total_direct[static_cast<size_t>(i)] - ay_reference_direct_tpf[static_cast<size_t>(i)];
      }
    }
  }

  for (int i = 0; i < state.n(); ++i) {
    const double radius = std::hypot(state.x[i], state.y[i]);
    csv << i << ',' << radius << ',' << state.vx[i] << ',' << state.vy[i] << ','
        << active_route << ',' << applied_accel_formula << ',' << decomposition_note << ','
        << ax_reference_direct_tpf[static_cast<size_t>(i)] << ',' << ay_reference_direct_tpf[static_cast<size_t>(i)]
        << ','
        << ax_reference_vdsg[static_cast<size_t>(i)] << ',' << ay_reference_vdsg[static_cast<size_t>(i)] << ','
        << ax_total[static_cast<size_t>(i)] << ',' << ay_total[static_cast<size_t>(i)] << ','
        << ax_bh_only[static_cast<size_t>(i)] << ',' << ay_bh_only[static_cast<size_t>(i)] << ','
        << ax_star_only[static_cast<size_t>(i)] << ',' << ay_star_only[static_cast<size_t>(i)] << '\n';
  }
}

}  // namespace

int main(int argc, char** argv) {
  // 1. Find run config path (root configs/ only; engine/configs/ is not used)
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
  const bool enable_softening_trace = (std::getenv("GALAXY_STARTUP_SOFTENING_TRACE") != nullptr);
  auto print_key_occurrences = [&](const std::string& path, const std::string& key) {
    if (!enable_softening_trace) return;
    if (path.empty()) {
      std::cout << "[startup_diag][softening_trace] key=" << key << " matches: (run config not found)\n";
      return;
    }
    const std::vector<galaxy::ConfigKeyOccurrence> matches =
        galaxy::scan_config_key_occurrences(path, key);
    if (matches.empty()) {
      std::cout << "[startup_diag][softening_trace] key=" << key << " matches: (none)\n";
      return;
    }
    std::cout << "[startup_diag][softening_trace] key=" << key
              << " matches (" << matches.size() << "):\n";
    for (const auto& match : matches) {
      std::cout << "  - line " << match.line_number << ": " << match.value << "\n";
    }
  };
  if (enable_softening_trace)
    std::cout << "[startup_diag][softening_trace] run_config_path="
              << (run_config_path.empty() ? "(none)" : run_config_path) << "\n";
  print_key_occurrences(run_config_path, "softening");
  print_key_occurrences(run_config_path, "tpfcore_source_softening");

  std::string package_defaults_path = galaxy::find_package_defaults_path(physics_pkg);
  if (!package_defaults_path.empty()) {
    galaxy::load_config_file(package_defaults_path, config);
  }
  if (enable_softening_trace) {
    std::cout << "[startup_diag][softening_trace] final softening after package defaults: "
              << config.softening << "\n";
    std::cout << "[startup_diag][softening_trace] final tpfcore_source_softening after package defaults: "
              << config.tpfcore_source_softening << "\n";
  }
  if (!run_config_path.empty()) {
    galaxy::load_config_file(run_config_path, config);
  }
  if (enable_softening_trace) {
    std::cout << "[startup_diag][softening_trace] final softening after run config: "
              << config.softening << "\n";
    std::cout << "[startup_diag][softening_trace] final tpfcore_source_softening after run config: "
              << config.tpfcore_source_softening << "\n";
  }

  config.run_id = run_id_from_time();
  config.output_dir = "../outputs/" + config.run_id;

  bool auto_plot = false;
  bool assume_yes = false;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--plot")
      auto_plot = true;
    if (std::string(argv[i]) == "--yes")
      assume_yes = true;
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
                     "tpf_symmetric_pair_inspect, tpf_source_field_benchmark, tpf_4d_static_residual_benchmark, "
                     "tpf_4d_static_motion_readout_benchmark, tpf_4d_xi_motion_probe_benchmark, "
                     "tpf_weak_field_calibration, tpf_newtonian_force_compare, "
                     "tpf_diagnostic_consistency_audit.\n";
        return 1;
      }
    }
  }

  for (int i = first_cli_config_idx; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--plot" || a == "--yes") continue;
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
  if (galaxy::simulation_mode_requires_output_dir(config.simulation_mode)) {
    if (!ensure_dir_recursive(config.output_dir)) {
      std::cerr << "Failed to create output directory for simulation_mode="
                << galaxy::mode_to_string(config.simulation_mode)
                << ": " << config.output_dir << "\n";
      return 1;
    }
  }

  if (galaxy::is_tpf_utility_mode(config.simulation_mode)) {
    galaxy::PhysicsPackage* utility_physics = galaxy::get_physics_package(config.physics_package);
    if (!utility_physics) {
      std::cerr << "Unknown physics package: " << config.physics_package << "\n";
      return 1;
    }
    utility_physics->init_from_config(config);
    if (!utility_physics->supports_utility_mode(config.simulation_mode)) {
      std::cerr << "simulation_mode=" << galaxy::mode_to_string(config.simulation_mode)
                << " is not supported by physics_package=" << config.physics_package << "\n";
      return 1;
    }
    if (!utility_physics->run_utility_mode(config, config.output_dir)) {
      return 1;
    }
    if (config.simulation_mode == galaxy::SimulationMode::tpf_4d_static_residual_benchmark && auto_plot) {
      const std::string dev_py = "../dev/bin/python3";
      const bool dev_py_exists = static_cast<bool>(std::ifstream(dev_py).good());
      const std::string py = dev_py_exists ? dev_py : "python3";
      const std::string cmd = py + " ../plot_tpf_4d_static_residual.py " + shell_single_quote(config.output_dir);
      const int ret = std::system(cmd.c_str());
      if (ret != 0) {
        std::cerr << "Warning: optional tpf_4d_static_residual plot script returned non-zero exit code.\n";
      }
      const std::vector<std::string> pngs = existing_tpf_4d_static_plot_pngs(config.output_dir);
      if (pngs.empty()) {
        std::cout << "plot script completed but no expected PNGs were found/generated\n";
      } else {
        std::cout << "Generated optional PNGs in " << config.output_dir << "\n";
      }
    }
    finalize_utility_mode_run(config, run_config_path, package_defaults_path);
    return 0;
  }

  galaxy::PhysicsPackage* physics = galaxy::get_physics_package(config.physics_package);
  if (!physics) {
    std::cerr << "Unknown physics package: " << config.physics_package << "\n";
    return 1;
  }
  physics->init_from_config(config);

  const galaxy::Config configured_after_layering = config;
  galaxy::ResolvedScenario resolved = galaxy::resolve_scenario(configured_after_layering);
  config = resolved.config;
  galaxy::State state = resolved.initial_state;
  int n_steps = resolved.effective_n_steps;
  int snapshot_every = resolved.effective_snapshot_every;

  galaxy::GalaxyPreflightSummary galaxy_preflight;
  if (config.simulation_mode == galaxy::SimulationMode::galaxy) {
    galaxy_preflight = galaxy::build_galaxy_preflight_summary(config, resolved);
    galaxy::print_galaxy_preflight_summary(galaxy_preflight);
    if (!galaxy_preflight.warnings.empty()) {
      std::cout << "Galaxy preflight warnings:\n";
      for (std::size_t i = 0; i < galaxy_preflight.warnings.size(); ++i)
        std::cout << "  - " << galaxy_preflight.warnings[i] << "\n";

      const bool interactive = IS_STDOUT_TERMINAL() && IS_STDIN_TERMINAL();
      if (interactive && !assume_yes) {
        std::cout << "Continue anyway? [y/N] ";
        std::string ans;
        std::getline(std::cin, ans);
        if (!(ans == "y" || ans == "Y")) {
          std::cerr << "Aborting due to preflight warnings.\n";
          return 1;
        }
      }
    }
  }

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
    write_galaxy_step0_accel_audit(config, physics, state);
    std::cout << "Wrote " << config.output_dir << "/galaxy_step0_accel_audit.csv\n";
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
    if (!ensure_dir_recursive(compare_parent_dir) ||
        !ensure_dir_recursive(left_dir) || !ensure_dir_recursive(right_dir)) {
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
              side_physics->cooling_active(side_cfg);
          const int cooling_steps = cooling_active
              ? side_physics->cooling_steps(side_cfg, n_steps)
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

          if (side_cfg.save_run_info) {
            galaxy::write_run_info(side_cfg.output_dir, side_cfg, n_steps, static_cast<int>(side_snaps.size()),
                                   state.n(), run_config_path, side_defaults_path,
                                   nullptr, nullptr,
                                   &galaxy::last_galaxy_init_audit(), &cooling_audit, nullptr);
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
              make_galaxy_step_progress_callback(start_wall, progress_to_terminal, side_tag, side_cfg.display_time_unit);
          auto side_snapshots = galaxy::run_simulation(side_cfg, side_state, side_physics,
                                                       n_steps, snapshot_every, side_progress,
                                                       compare_progress_interval);
          const bool compare_dual_line_mode_active =
              (std::getenv("GALAXY_COMPARE_DUAL_LINE_PROGRESS") != nullptr);
          if (progress_to_terminal && n_steps > 0 && !compare_dual_line_mode_active) std::cout << "\n";
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
        // Reserve two terminal lines for dual-line in-place progress updates without
        // showing placeholder labels.
        std::cout << "\n\n";
      } else {
        std::cout << "Parallel compare enabled; child logs:\n  " << left_log << "\n  " << right_log << "\n";
      }

      pid_t left_pid = fork();
      if (left_pid == 0) {
        if (!show_live_compare_progress) {
          FILE* lf = std::freopen(left_log.c_str(), "w", stdout);
          FILE* le = std::freopen(left_log.c_str(), "a", stderr);
          if (!lf || !le) _exit(90);
        } else {
          setenv("GALAXY_COMPARE_DUAL_LINE_PROGRESS", "1", 1);
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
        } else {
          setenv("GALAXY_COMPARE_DUAL_LINE_PROGRESS", "1", 1);
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
      if (show_live_compare_progress) std::cout << "\n";
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
      // Run from engine cwd: script lives at repo_root/plot_cpp_compare.py
      const std::string dev_py = "../dev/bin/python3";
      const bool dev_py_exists = static_cast<bool>(std::ifstream(dev_py).good());
      const std::string py = dev_py_exists ? dev_py : "python3";
      std::string cmd = py + " ../plot_cpp_compare.py " + compare_parent_dir;
      int ret = std::system(cmd.c_str());
      if (ret != 0) {
        std::cerr << "Warning: compare renderer returned non-zero exit code. Command was:\n  " << cmd << "\n";
      } else {
        const std::string legacy_initial = compare_parent_dir + "/galaxy_initial_compare.png";
        const std::string legacy_final = compare_parent_dir + "/galaxy_final_compare.png";
        const auto mode_initial = find_compare_side_by_side_pngs(compare_parent_dir, "initial");
        const auto mode_final = find_compare_side_by_side_pngs(compare_parent_dir, "final");
        std::vector<std::string> missing;
        if (!file_exists(legacy_initial)) missing.push_back("galaxy_initial_compare.png");
        if (!file_exists(legacy_final)) missing.push_back("galaxy_final_compare.png");
        if (mode_initial.empty()) missing.push_back("galaxy_compare__*__initial_side_by_side.png");
        if (mode_final.empty()) missing.push_back("galaxy_compare__*__final_side_by_side.png");
        if (!missing.empty()) {
          std::cerr << "Error: compare renderer reported success but expected files are missing in "
                    << compare_parent_dir << ": ";
          for (std::size_t i = 0; i < missing.size(); ++i) {
            if (i) std::cerr << ", ";
            std::cerr << missing[i];
          }
          std::cerr << "\n";
          return 1;
        }
        std::cout << "Compare render finished. Generated PNGs in " << compare_parent_dir << ":\n"
                  << "  galaxy_initial_compare.png\n"
                  << "  galaxy_final_compare.png\n";
        for (const auto& f : mode_initial) std::cout << "  " << f << "\n";
        for (const auto& f : mode_final) std::cout << "  " << f << "\n";
      }
    } else {
      std::cout
          << "\nCompare side-by-side PNGs/animation were not generated (run without --plot).\n"
          << "From the engine directory:\n  python3 ../plot_cpp_compare.py " << compare_parent_dir << "\n"
          << "From the repository root:\n  python3 plot_cpp_compare.py " << compare_parent_dir << "\n"
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
    progress_callback = make_galaxy_step_progress_callback(start_wall, progress_to_terminal, "", config.display_time_unit);
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
        physics->cooling_active(config);
    const int cooling_steps = cooling_active
        ? physics->cooling_steps(config, n_steps)
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

    if (config.save_run_info) {
      galaxy::write_run_info(config.output_dir, config, n_steps, static_cast<int>(snapshots.size()), state.n(),
                             run_config_path, package_defaults_path,
                             &configured_after_layering,
                             &resolved,
                             config.simulation_mode == galaxy::SimulationMode::galaxy
                                 ? &galaxy::last_galaxy_init_audit()
                                 : nullptr,
                             &cooling_audit,
                             nullptr);
      if (config.simulation_mode == galaxy::SimulationMode::galaxy)
        galaxy::append_galaxy_preflight_to_run_info(config.output_dir, galaxy_preflight);
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

    if (!physics->write_post_run_diagnostics(snapshots, config, config.output_dir)) {
      return 1;
    }
  }

  std::cout << "Snapshots: " << snapshots.size() << "\n";
  std::cout << "Output directory: " << config.output_dir << "\n";

  std::string output_dir = config.output_dir;
  if (auto_plot) {
    std::cout.flush();
    std::cerr.flush();
    std::cout << "Simulation complete. Rendering animation..." << std::endl;
    // Run from engine cwd (same approach as compare renderer): script lives at repo_root/plot_cpp_run.py.
    const std::string dev_py = "../dev/bin/python3";
    const bool dev_py_exists = static_cast<bool>(std::ifstream(dev_py).good());
    const std::string py = dev_py_exists ? dev_py : "python3";
    std::string cmd = py + " ../plot_cpp_run.py " + output_dir;
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
      std::cerr << "Warning: Python rendering script returned non-zero exit code." << std::endl;
    }
  }

  return 0;
}
