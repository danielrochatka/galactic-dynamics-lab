#include "config.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iomanip>

namespace galaxy {

namespace {

std::string trim(const std::string& s) {
  auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end == std::string::npos ? std::string::npos : end - start + 1);
}

bool parse_bool(const std::string& v) {
  std::string s = trim(v);
  if (s == "1" || s == "true" || s == "yes") return true;
  if (s == "0" || s == "false" || s == "no") return false;
  return false;
}

bool file_exists(const std::string& path) {
  std::ifstream f(path);
  return f.good();
}

std::string lower_copy(std::string s) {
  for (char& c : s) {
    if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
  }
  return s;
}

bool in_list(const std::string& s, const std::initializer_list<const char*>& vals) {
  for (const char* v : vals)
    if (s == v) return true;
  return false;
}

}  // namespace

SimulationMode parse_mode(const std::string& s) {
  std::string t = trim(s);
  if (t == "galaxy") return SimulationMode::galaxy;
  if (t == "earth_moon_benchmark") return SimulationMode::earth_moon_benchmark;
  if (t == "bh_orbit_validation") return SimulationMode::bh_orbit_validation;
  if (t == "two_body_orbit") {
    std::cerr << "Warning: simulation_mode \"two_body_orbit\" is deprecated; use \"earth_moon_benchmark\" "
                 "(same Earth–Moon SI benchmark).\n";
    return SimulationMode::earth_moon_benchmark;
  }
  if (t == "symmetric_pair") return SimulationMode::symmetric_pair;
  if (t == "small_n_conservation") return SimulationMode::small_n_conservation;
  if (t == "timestep_convergence") return SimulationMode::timestep_convergence;
  if (t == "tpf_single_source_inspect") return SimulationMode::tpf_single_source_inspect;
  if (t == "tpf_symmetric_pair_inspect") return SimulationMode::tpf_symmetric_pair_inspect;
  if (t == "tpf_two_body_sweep") return SimulationMode::tpf_two_body_sweep;
  if (t == "tpf_weak_field_calibration") return SimulationMode::tpf_weak_field_calibration;
  if (t == "tpf_newtonian_force_compare") return SimulationMode::tpf_newtonian_force_compare;
  if (t == "tpf_diagnostic_consistency_audit") return SimulationMode::tpf_diagnostic_consistency_audit;
  if (t == "tpf_bound_orbit_sweep") return SimulationMode::tpf_bound_orbit_sweep;
  if (t == "tpf_v11_weak_field_correspondence") return SimulationMode::tpf_v11_weak_field_correspondence;
  throw std::runtime_error("Unknown simulation_mode: " + s);
}

std::string mode_to_string(SimulationMode m) {
  switch (m) {
    case SimulationMode::galaxy: return "galaxy";
    /* Legacy enum value (int 1 in old run_info): same physics as earth_moon_benchmark. */
    case SimulationMode::two_body_orbit: return "earth_moon_benchmark";
    case SimulationMode::earth_moon_benchmark: return "earth_moon_benchmark";
    case SimulationMode::bh_orbit_validation: return "bh_orbit_validation";
    case SimulationMode::symmetric_pair: return "symmetric_pair";
    case SimulationMode::small_n_conservation: return "small_n_conservation";
    case SimulationMode::timestep_convergence: return "timestep_convergence";
    case SimulationMode::tpf_single_source_inspect: return "tpf_single_source_inspect";
    case SimulationMode::tpf_symmetric_pair_inspect: return "tpf_symmetric_pair_inspect";
    case SimulationMode::tpf_two_body_sweep: return "tpf_two_body_sweep";
    case SimulationMode::tpf_weak_field_calibration: return "tpf_weak_field_calibration";
    case SimulationMode::tpf_newtonian_force_compare: return "tpf_newtonian_force_compare";
    case SimulationMode::tpf_diagnostic_consistency_audit: return "tpf_diagnostic_consistency_audit";
    case SimulationMode::tpf_bound_orbit_sweep: return "tpf_bound_orbit_sweep";
    case SimulationMode::tpf_v11_weak_field_correspondence: return "tpf_v11_weak_field_correspondence";
  }
  return "unknown";
}

bool apply_config_kv(const std::string& key, const std::string& val, Config& config) {
  if (key == "simulation_mode") {
    config.simulation_mode = parse_mode(val);
    return true;
  }
  if (key == "n_stars") {
    config.n_stars = std::stoi(val);
    return true;
  }
  if (key == "star_mass") {
    config.star_mass = std::stod(val);
    return true;
  }
  if (key == "bh_mass") {
    config.bh_mass = std::stod(val);
    return true;
  }
  if (key == "inner_radius") {
    config.inner_radius = std::stod(val);
    return true;
  }
  if (key == "outer_radius") {
    config.outer_radius = std::stod(val);
    return true;
  }
  if (key == "galaxy_radius") {
    config.galaxy_radius = std::stod(val);
    return true;
  }
  if (key == "dt") {
    config.dt = std::stod(val);
    return true;
  }
  if (key == "n_steps") {
    config.n_steps = std::stoi(val);
    return true;
  }
  if (key == "snapshot_every") {
    config.snapshot_every = std::stoi(val);
    return true;
  }
  if (key == "softening") {
    config.softening = std::stod(val);
    return true;
  }
  if (key == "enable_star_star_gravity") {
    config.enable_star_star_gravity = parse_bool(val);
    return true;
  }
  if (key == "physics_package") {
    config.physics_package = val;
    return true;
  }
  if (key == "physics_package_compare") {
    config.physics_package_compare = trim(val);
    return true;
  }
  if (key == "compare_parallel") {
    config.compare_parallel = parse_bool(val);
    return true;
  }
  if (key == "tpfcore_enable_provisional_readout") {
    config.tpfcore_enable_provisional_readout = parse_bool(val);
    return true;
  }
  if (key == "tpfcore_readout_mode") {
    config.tpfcore_readout_mode = trim(val);
    return true;
  }
  if (key == "tpfcore_readout_scale") {
    config.tpfcore_readout_scale = std::stod(val);
    return true;
  }
  if (key == "tpfcore_theta_tt_scale") {
    config.tpfcore_theta_tt_scale = std::stod(val);
    return true;
  }
  if (key == "tpfcore_theta_tr_scale") {
    config.tpfcore_theta_tr_scale = std::stod(val);
    return true;
  }
  if (key == "tpf_dynamics_mode") {
    std::string s = trim(val);
    if (s == "weak_field_correspondence") s = "v11_weak_field_truncation";
    if (s != "legacy_readout" && s != "v11_weak_field_truncation" && s != "direct_tpf") {
      throw std::runtime_error(
          "tpf_dynamics_mode must be legacy_readout, v11_weak_field_truncation, or direct_tpf, got: " + val);
    }
    config.tpf_dynamics_mode = s;
    return true;
  }
  if (key == "tpf_weak_field_correspondence_alpha_si") {
    config.tpf_weak_field_correspondence_alpha_si = std::stod(val);
    return true;
  }
  if (key == "tpf_analysis_mode") {
    std::string s = trim(val);
    if (s != "none" && s != "v11_weak_field_correspondence") {
      throw std::runtime_error("tpf_analysis_mode must be none or v11_weak_field_correspondence, got: " + val);
    }
    config.tpf_analysis_mode = s;
    return true;
  }
  if (key == "v11_weak_field_correspondence_benchmark") {
    std::string s = trim(val);
    if (s != "axis_monopole" && s != "earth_moon_line_of_centers") {
      throw std::runtime_error(
          "v11_weak_field_correspondence_benchmark must be axis_monopole or earth_moon_line_of_centers, got: " + val);
    }
    config.v11_weak_field_correspondence_benchmark = s;
    return true;
  }
  if (key == "v11_em_mass_earth_kg") {
    config.v11_em_mass_earth_kg = std::stod(val);
    return true;
  }
  if (key == "v11_em_mass_moon_kg") {
    config.v11_em_mass_moon_kg = std::stod(val);
    return true;
  }
  if (key == "v11_em_mean_distance_m") {
    config.v11_em_mean_distance_m = std::stod(val);
    return true;
  }
  if (key == "v11_em_sidereal_period_s") {
    config.v11_em_sidereal_period_s = std::stod(val);
    return true;
  }
  if (key == "v11_em_calib_surface_radius_m") {
    config.v11_em_calib_surface_radius_m = std::stod(val);
    return true;
  }
  if (key == "v11_em_calib_surface_g_m_s2") {
    config.v11_em_calib_surface_g_m_s2 = std::stod(val);
    return true;
  }
  if (key == "tpfcore_closure_kappa" || key == "tpf_kappa") {
    config.tpf_kappa = std::stod(val);
    return true;
  }
  if (key == "tpf_vdsg_coupling") {
    config.tpf_vdsg_coupling = std::stod(val);
    return true;
  }
  if (key == "tpf_vdsg_mass_baseline_kg") {
    config.tpf_vdsg_mass_baseline_kg = std::stod(val);
    return true;
  }
  if (key == "tpf_global_accel_shunt_enable") {
    config.tpf_global_accel_shunt_enable = parse_bool(val);
    return true;
  }
  if (key == "tpf_global_accel_shunt_fraction") {
    config.tpf_global_accel_shunt_fraction = std::stod(val);
    return true;
  }
  if (key == "tpf_accel_pipeline_diagnostics_csv") {
    config.tpf_accel_pipeline_diagnostics_csv = parse_bool(val);
    return true;
  }
  if (key == "tpf_poisson_bins") {
    config.tpf_poisson_bins = std::stoi(val);
    return true;
  }
  if (key == "tpf_poisson_max_radius") {
    config.tpf_poisson_max_radius = std::stod(val);
    return true;
  }
  if (key == "tpf_cooling_fraction") {
    config.tpf_cooling_fraction = std::stod(val);
    return true;
  }
  if (key == "tpfcore_dump_readout_debug") {
    config.tpfcore_dump_readout_debug = parse_bool(val);
    return true;
  }
  if (key == "tpfcore_probe_radius_min") {
    config.tpfcore_probe_radius_min = std::stod(val);
    return true;
  }
  if (key == "tpfcore_probe_radius_max") {
    config.tpfcore_probe_radius_max = std::stod(val);
    return true;
  }
  if (key == "tpfcore_probe_samples") {
    config.tpfcore_probe_samples = std::stoi(val);
    return true;
  }
  if (key == "tpfcore_dump_invariant_profile") {
    config.tpfcore_dump_invariant_profile = parse_bool(val);
    return true;
  }
  if (key == "tpfcore_dump_theta_profile") {
    config.tpfcore_dump_theta_profile = parse_bool(val);
    return true;
  }
  if (key == "tpfcore_source_softening") {
    config.tpfcore_source_softening = std::stod(val);
    return true;
  }
  if (key == "tpfcore_residual_step") {
    config.tpfcore_residual_step = std::stod(val);
    return true;
  }
  if (key == "tpfcore_live_orbit_force_audit") {
    config.tpfcore_live_orbit_force_audit = parse_bool(val);
    return true;
  }
  if (key == "galaxy_init_template") {
    config.galaxy_init_template = trim(val);
    return true;
  }
  if (key == "galaxy_init_seed") {
    config.galaxy_init_seed = static_cast<unsigned>(std::stoul(val));
    return true;
  }
  if (key == "galaxy_init_position_noise") {
    config.galaxy_init_position_noise = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_velocity_angle_noise") {
    config.galaxy_init_velocity_angle_noise = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_velocity_magnitude_noise") {
    config.galaxy_init_velocity_magnitude_noise = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_clumpiness") {
    config.galaxy_init_clumpiness = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_num_clumps") {
    config.galaxy_init_num_clumps = std::stoi(val);
    return true;
  }
  if (key == "galaxy_init_clump_radius_fraction") {
    config.galaxy_init_clump_radius_fraction = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_m2_amplitude") {
    config.galaxy_init_m2_amplitude = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_m3_amplitude") {
    config.galaxy_init_m3_amplitude = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_bar_amplitude") {
    config.galaxy_init_bar_amplitude = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_bar_axis_ratio") {
    config.galaxy_init_bar_axis_ratio = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_spiral_amplitude") {
    config.galaxy_init_spiral_amplitude = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_spiral_winding") {
    config.galaxy_init_spiral_winding = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_spiral_phase") {
    config.galaxy_init_spiral_phase = std::stod(val);
    return true;
  }
  if (key == "galaxy_init_master_chaos") {
    config.galaxy_init_master_chaos = std::stod(val);
    return true;
  }
  if (key == "velocity_noise") {
    config.velocity_noise = std::stod(val);
    return true;
  }
  if (key == "initial_velocity_scale") {
    config.initial_velocity_scale = std::stod(val);
    return true;
  }
  if (key == "save_snapshots") {
    config.save_snapshots = parse_bool(val);
    return true;
  }
  if (key == "save_run_info") {
    config.save_run_info = parse_bool(val);
    return true;
  }
  if (key == "plot_animation_dynamic_zoom") {
    config.plot_animation_dynamic_zoom = parse_bool(val);
    return true;
  }
  if (key == "plot_skip_initial_steps") {
    config.plot_skip_initial_steps = std::stoi(val);
    return true;
  }
  if (key == "plot_skip_initial_snapshots") {
    config.plot_skip_initial_snapshots = std::stoi(val);
    return true;
  }
  if (key == "diagnostic_cutoff_radius") {
    config.diagnostic_cutoff_radius = std::stod(val);
    return true;
  }
  if (key == "render_overlay_mode") {
    std::string lo = lower_copy(trim(val));
    if (lo == "none" || lo == "minimal" || lo == "audit_full") {
      config.render_overlay_mode = lo;
      return true;
    }
    throw std::runtime_error("render_overlay_mode must be none, minimal, or audit_full");
  }
  if (key == "display_distance_unit") {
    const std::string s = trim(val);
    if (in_list(s, {"auto", "m", "km", "AU", "ly", "pc", "kpc"})) {
      config.display_distance_unit = s;
      return true;
    }
    throw std::runtime_error("display_distance_unit must be auto, m, km, AU, ly, pc, or kpc");
  }
  if (key == "display_time_unit") {
    const std::string s = trim(val);
    if (in_list(s, {"auto", "s", "min", "hr", "day", "yr", "kyr", "Myr"})) {
      config.display_time_unit = s;
      return true;
    }
    throw std::runtime_error("display_time_unit must be auto, s, min, hr, day, yr, kyr, or Myr");
  }
  if (key == "display_velocity_unit") {
    const std::string s = trim(val);
    if (in_list(s, {"auto", "m/s", "km/s"})) {
      config.display_velocity_unit = s;
      return true;
    }
    throw std::runtime_error("display_velocity_unit must be auto, m/s, or km/s");
  }
  if (key == "display_units_in_overlay") {
    config.display_units_in_overlay = parse_bool(val);
    return true;
  }
  if (key == "display_show_unit_reference") {
    config.display_show_unit_reference = parse_bool(val);
    return true;
  }
  if (key == "tpf_gdd_coupling") {
    config.tpf_vdsg_coupling = std::stod(val);
    return true;
  }
  if (key == "validation_two_body_radius") {
    config.validation_two_body_radius = std::stod(val);
    return true;
  }
  if (key == "validation_two_body_speed_ratio") {
    config.validation_two_body_speed_ratio = std::stod(val);
    return true;
  }
  if (key == "validation_earth_mass") {
    config.validation_earth_mass = std::stod(val);
    return true;
  }
  if (key == "validation_moon_mass") {
    config.validation_moon_mass = std::stod(val);
    return true;
  }
  if (key == "validation_earth_moon_distance") {
    config.validation_earth_moon_distance = std::stod(val);
    return true;
  }
  if (key == "validation_moon_tangential_speed") {
    config.validation_moon_tangential_speed = std::stod(val);
    return true;
  }
  if (key == "validation_symmetric_include_bh") {
    config.validation_symmetric_include_bh = parse_bool(val);
    return true;
  }
  if (key == "validation_symmetric_separation") {
    config.validation_symmetric_separation = std::stod(val);
    return true;
  }
  if (key == "validation_symmetric_speed") {
    config.validation_symmetric_speed = std::stod(val);
    return true;
  }
  if (key == "validation_small_n") {
    config.validation_small_n = std::stoi(val);
    return true;
  }
  if (key == "validation_n_steps") {
    config.validation_n_steps = std::stoi(val);
    return true;
  }
  if (key == "validation_snapshot_every") {
    config.validation_snapshot_every = std::stoi(val);
    return true;
  }
  if (key == "output_dir") {
    config.output_dir = trim(val);
    return true;
  }
  return false;
}

bool load_config_file(const std::string& path, Config& config) {
  std::ifstream f(path);
  if (!f) return false;

  std::string line;
  int line_num = 0;
  while (std::getline(f, line)) {
    ++line_num;
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;
    auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    std::string key = trim(line.substr(0, eq));
    std::string val = trim(line.substr(eq + 1));
    if (key.empty()) continue;

    try {
      if (!apply_config_kv(key, val, config)) {
        /* unknown key: ignore (same as before) */
      }
    } catch (const std::exception& e) {
      throw std::runtime_error("Config error in " + path + " line " + std::to_string(line_num) +
          " (key=" + key + "): " + e.what());
    }
  }
  return true;
}

std::string probe_config_key(const std::string& path, const std::string& key) {
  std::ifstream f(path);
  if (!f) return "";

  std::string line;
  while (std::getline(f, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;
    auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    std::string k = trim(line.substr(0, eq));
    if (k == key)
      return trim(line.substr(eq + 1));
  }
  return "";
}

// Root configs only: try repo-root configs/ (../configs when run from cpp_sim, configs when run from repo root).
// We never use cpp_sim/configs/; if it exists we warn or fail (see check_run_config_canonical).
std::string find_run_config_path() {
  // Optional explicit path for CI / integration tests (takes precedence over my.local.cfg).
  const char* env = std::getenv("GALAXY_RUN_CONFIG");
  if (env && env[0] != '\0') {
    std::string p = trim(std::string(env));
    if (file_exists(p)) return p;
  }
  // 1. Root configs when run from cpp_sim/
  if (file_exists("../configs/my.local.cfg")) return "../configs/my.local.cfg";
  if (file_exists("../configs/local/my.local.cfg")) return "../configs/local/my.local.cfg";
  // 2. Root configs when run from repo root (only if ../configs is not the repo root = we are at repo root)
  if (!file_exists("../configs/example.cfg")) {
    if (file_exists("configs/my.local.cfg")) return "configs/my.local.cfg";
    if (file_exists("configs/local/my.local.cfg")) return "configs/local/my.local.cfg";
  }
  return "";
}

// If cpp_sim/configs/ exists (configs/ when cwd is cpp_sim): warn when we're using root configs, fail when we'd have used it.
// Returns false if caller should exit(1).
bool check_run_config_canonical(const std::string& run_config_path) {
  bool have_local = file_exists("configs/my.local.cfg") || file_exists("configs/local/my.local.cfg");
  if (run_config_path.empty()) {
    if (have_local) {
      std::cerr << "Error: run config found under configs/ but run configs must live in repository root configs/, "
                   "not cpp_sim/configs/. Use root configs/ (e.g. ../configs/my.local.cfg when running from cpp_sim).\n";
      return false;
    }
    return true;
  }
  if (have_local && run_config_path.find("../configs/") == 0u) {
    std::cout << "Warning: configs/ (cpp_sim/configs/) exists but is ignored. Run configs must live in root configs/ only.\n";
  }
  return true;
}

// Package defaults only in cpp_sim/physics/<Package>/defaults.cfg (relative to cpp_sim when run from cpp_sim).
std::string find_package_defaults_path(const std::string& package_name) {
  const std::string path = "physics/" + package_name + "/defaults.cfg";
  if (file_exists(path)) return path;
  return "";
}

std::vector<std::pair<std::string, std::string>> serialize_config_kv(const Config& config) {
  auto b = [](bool v) { return v ? std::string("1") : std::string("0"); };
  auto i = [](int v) { return std::to_string(v); };
  auto u = [](unsigned v) { return std::to_string(v); };
  auto d = [](double v) {
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
  };
  std::vector<std::pair<std::string, std::string>> kv;
  kv.reserve(80);
  kv.emplace_back("simulation_mode", mode_to_string(config.simulation_mode));
  kv.emplace_back("n_stars", i(config.n_stars));
  kv.emplace_back("star_mass", d(config.star_mass));
  kv.emplace_back("bh_mass", d(config.bh_mass));
  kv.emplace_back("inner_radius", d(config.inner_radius));
  kv.emplace_back("outer_radius", d(config.outer_radius));
  kv.emplace_back("galaxy_radius", d(config.galaxy_radius));
  kv.emplace_back("dt", d(config.dt));
  kv.emplace_back("n_steps", i(config.n_steps));
  kv.emplace_back("snapshot_every", i(config.snapshot_every));
  kv.emplace_back("softening", d(config.softening));
  kv.emplace_back("enable_star_star_gravity", b(config.enable_star_star_gravity));
  kv.emplace_back("physics_package", config.physics_package);
  kv.emplace_back("physics_package_compare", config.physics_package_compare);
  kv.emplace_back("compare_parallel", b(config.compare_parallel));
  kv.emplace_back("tpf_dynamics_mode", config.tpf_dynamics_mode);
  kv.emplace_back("tpf_weak_field_correspondence_alpha_si", d(config.tpf_weak_field_correspondence_alpha_si));
  kv.emplace_back("tpf_analysis_mode", config.tpf_analysis_mode);
  kv.emplace_back("v11_weak_field_correspondence_benchmark", config.v11_weak_field_correspondence_benchmark);
  kv.emplace_back("v11_em_mass_earth_kg", d(config.v11_em_mass_earth_kg));
  kv.emplace_back("v11_em_mass_moon_kg", d(config.v11_em_mass_moon_kg));
  kv.emplace_back("v11_em_mean_distance_m", d(config.v11_em_mean_distance_m));
  kv.emplace_back("v11_em_sidereal_period_s", d(config.v11_em_sidereal_period_s));
  kv.emplace_back("v11_em_calib_surface_radius_m", d(config.v11_em_calib_surface_radius_m));
  kv.emplace_back("v11_em_calib_surface_g_m_s2", d(config.v11_em_calib_surface_g_m_s2));
  kv.emplace_back("tpfcore_enable_provisional_readout", b(config.tpfcore_enable_provisional_readout));
  kv.emplace_back("tpfcore_readout_mode", config.tpfcore_readout_mode);
  kv.emplace_back("tpfcore_readout_scale", d(config.tpfcore_readout_scale));
  kv.emplace_back("tpfcore_theta_tt_scale", d(config.tpfcore_theta_tt_scale));
  kv.emplace_back("tpfcore_theta_tr_scale", d(config.tpfcore_theta_tr_scale));
  kv.emplace_back("tpf_kappa", d(config.tpf_kappa));
  kv.emplace_back("tpf_vdsg_coupling", d(config.tpf_vdsg_coupling));
  kv.emplace_back("tpf_vdsg_mass_baseline_kg", d(config.tpf_vdsg_mass_baseline_kg));
  kv.emplace_back("tpf_global_accel_shunt_enable", b(config.tpf_global_accel_shunt_enable));
  kv.emplace_back("tpf_global_accel_shunt_fraction", d(config.tpf_global_accel_shunt_fraction));
  kv.emplace_back("tpf_accel_pipeline_diagnostics_csv", b(config.tpf_accel_pipeline_diagnostics_csv));
  kv.emplace_back("tpf_poisson_bins", i(config.tpf_poisson_bins));
  kv.emplace_back("tpf_poisson_max_radius", d(config.tpf_poisson_max_radius));
  kv.emplace_back("tpf_cooling_fraction", d(config.tpf_cooling_fraction));
  kv.emplace_back("tpfcore_dump_readout_debug", b(config.tpfcore_dump_readout_debug));
  kv.emplace_back("tpfcore_live_orbit_force_audit", b(config.tpfcore_live_orbit_force_audit));
  kv.emplace_back("tpfcore_probe_radius_min", d(config.tpfcore_probe_radius_min));
  kv.emplace_back("tpfcore_probe_radius_max", d(config.tpfcore_probe_radius_max));
  kv.emplace_back("tpfcore_probe_samples", i(config.tpfcore_probe_samples));
  kv.emplace_back("tpfcore_dump_invariant_profile", b(config.tpfcore_dump_invariant_profile));
  kv.emplace_back("tpfcore_dump_theta_profile", b(config.tpfcore_dump_theta_profile));
  kv.emplace_back("tpfcore_source_softening", d(config.tpfcore_source_softening));
  kv.emplace_back("tpfcore_residual_step", d(config.tpfcore_residual_step));
  kv.emplace_back("galaxy_init_template", config.galaxy_init_template);
  kv.emplace_back("galaxy_init_seed", u(config.galaxy_init_seed));
  kv.emplace_back("galaxy_init_position_noise", d(config.galaxy_init_position_noise));
  kv.emplace_back("galaxy_init_velocity_angle_noise", d(config.galaxy_init_velocity_angle_noise));
  kv.emplace_back("galaxy_init_velocity_magnitude_noise", d(config.galaxy_init_velocity_magnitude_noise));
  kv.emplace_back("galaxy_init_clumpiness", d(config.galaxy_init_clumpiness));
  kv.emplace_back("galaxy_init_num_clumps", i(config.galaxy_init_num_clumps));
  kv.emplace_back("galaxy_init_clump_radius_fraction", d(config.galaxy_init_clump_radius_fraction));
  kv.emplace_back("galaxy_init_m2_amplitude", d(config.galaxy_init_m2_amplitude));
  kv.emplace_back("galaxy_init_m3_amplitude", d(config.galaxy_init_m3_amplitude));
  kv.emplace_back("galaxy_init_bar_amplitude", d(config.galaxy_init_bar_amplitude));
  kv.emplace_back("galaxy_init_bar_axis_ratio", d(config.galaxy_init_bar_axis_ratio));
  kv.emplace_back("galaxy_init_spiral_amplitude", d(config.galaxy_init_spiral_amplitude));
  kv.emplace_back("galaxy_init_spiral_winding", d(config.galaxy_init_spiral_winding));
  kv.emplace_back("galaxy_init_spiral_phase", d(config.galaxy_init_spiral_phase));
  kv.emplace_back("galaxy_init_master_chaos", d(config.galaxy_init_master_chaos));
  kv.emplace_back("velocity_noise", d(config.velocity_noise));
  kv.emplace_back("initial_velocity_scale", d(config.initial_velocity_scale));
  kv.emplace_back("save_snapshots", b(config.save_snapshots));
  kv.emplace_back("save_run_info", b(config.save_run_info));
  kv.emplace_back("plot_animation_dynamic_zoom", b(config.plot_animation_dynamic_zoom));
  kv.emplace_back("plot_skip_initial_steps", i(config.plot_skip_initial_steps));
  kv.emplace_back("plot_skip_initial_snapshots", i(config.plot_skip_initial_snapshots));
  kv.emplace_back("diagnostic_cutoff_radius", d(config.diagnostic_cutoff_radius));
  kv.emplace_back("render_overlay_mode", config.render_overlay_mode);
  kv.emplace_back("display_distance_unit", config.display_distance_unit);
  kv.emplace_back("display_time_unit", config.display_time_unit);
  kv.emplace_back("display_velocity_unit", config.display_velocity_unit);
  kv.emplace_back("display_units_in_overlay", b(config.display_units_in_overlay));
  kv.emplace_back("display_show_unit_reference", b(config.display_show_unit_reference));
  kv.emplace_back("validation_two_body_radius", d(config.validation_two_body_radius));
  kv.emplace_back("validation_two_body_speed_ratio", d(config.validation_two_body_speed_ratio));
  kv.emplace_back("validation_earth_mass", d(config.validation_earth_mass));
  kv.emplace_back("validation_moon_mass", d(config.validation_moon_mass));
  kv.emplace_back("validation_earth_moon_distance", d(config.validation_earth_moon_distance));
  kv.emplace_back("validation_moon_tangential_speed", d(config.validation_moon_tangential_speed));
  kv.emplace_back("validation_symmetric_include_bh", b(config.validation_symmetric_include_bh));
  kv.emplace_back("validation_symmetric_separation", d(config.validation_symmetric_separation));
  kv.emplace_back("validation_symmetric_speed", d(config.validation_symmetric_speed));
  kv.emplace_back("validation_small_n", i(config.validation_small_n));
  kv.emplace_back("validation_n_steps", i(config.validation_n_steps));
  kv.emplace_back("validation_snapshot_every", i(config.validation_snapshot_every));
  kv.emplace_back("output_dir", config.output_dir);
  kv.emplace_back("run_id", config.run_id);
  return kv;
}

}  // namespace galaxy
