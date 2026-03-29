#include "config.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

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

}  // namespace

SimulationMode parse_mode(const std::string& s) {
  std::string t = trim(s);
  if (t == "galaxy") return SimulationMode::galaxy;
  if (t == "two_body_orbit") return SimulationMode::two_body_orbit;
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
  throw std::runtime_error("Unknown simulation_mode: " + s);
}

std::string mode_to_string(SimulationMode m) {
  switch (m) {
    case SimulationMode::galaxy: return "galaxy";
    case SimulationMode::two_body_orbit: return "two_body_orbit";
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
  if (key == "tpf_kappa") {
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
  if (key == "validation_two_body_radius") {
    config.validation_two_body_radius = std::stod(val);
    return true;
  }
  if (key == "validation_two_body_speed_ratio") {
    config.validation_two_body_speed_ratio = std::stod(val);
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

}  // namespace galaxy
