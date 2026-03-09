#include "config.hpp"
#include <fstream>
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

}  // namespace

SimulationMode parse_mode(const std::string& s) {
  std::string t = trim(s);
  if (t == "galaxy") return SimulationMode::galaxy;
  if (t == "two_body_orbit") return SimulationMode::two_body_orbit;
  if (t == "symmetric_pair") return SimulationMode::symmetric_pair;
  if (t == "small_n_conservation") return SimulationMode::small_n_conservation;
  if (t == "timestep_convergence") return SimulationMode::timestep_convergence;
  throw std::runtime_error("Unknown simulation_mode: " + s);
}

bool load_config_file(const std::string& path, Config& config) {
  std::ifstream f(path);
  if (!f) return false;

  std::string line;
  while (std::getline(f, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;
    auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    std::string key = trim(line.substr(0, eq));
    std::string val = trim(line.substr(eq + 1));
    if (key.empty()) continue;

    if (key == "simulation_mode") { config.simulation_mode = parse_mode(val); continue; }
    if (key == "n_stars") { config.n_stars = std::stoi(val); continue; }
    if (key == "star_mass") { config.star_mass = std::stod(val); continue; }
    if (key == "bh_mass") { config.bh_mass = std::stod(val); continue; }
    if (key == "inner_radius") { config.inner_radius = std::stod(val); continue; }
    if (key == "outer_radius") { config.outer_radius = std::stod(val); continue; }
    if (key == "dt") { config.dt = std::stod(val); continue; }
    if (key == "n_steps") { config.n_steps = std::stoi(val); continue; }
    if (key == "snapshot_every") { config.snapshot_every = std::stoi(val); continue; }
    if (key == "softening") { config.softening = std::stod(val); continue; }
    if (key == "enable_star_star_gravity") { config.enable_star_star_gravity = parse_bool(val); continue; }
    if (key == "physics_package") { config.physics_package = val; continue; }
    if (key == "velocity_noise") { config.velocity_noise = std::stod(val); continue; }
    if (key == "initial_velocity_scale") { config.initial_velocity_scale = std::stod(val); continue; }
    if (key == "save_snapshots") { config.save_snapshots = parse_bool(val); continue; }
    if (key == "save_run_info") { config.save_run_info = parse_bool(val); continue; }
    if (key == "validation_two_body_radius") { config.validation_two_body_radius = std::stod(val); continue; }
    if (key == "validation_two_body_speed_ratio") { config.validation_two_body_speed_ratio = std::stod(val); continue; }
    if (key == "validation_symmetric_include_bh") { config.validation_symmetric_include_bh = parse_bool(val); continue; }
    if (key == "validation_symmetric_separation") { config.validation_symmetric_separation = std::stod(val); continue; }
    if (key == "validation_symmetric_speed") { config.validation_symmetric_speed = std::stod(val); continue; }
    if (key == "validation_small_n") { config.validation_small_n = std::stoi(val); continue; }
    if (key == "validation_n_steps") { config.validation_n_steps = std::stoi(val); continue; }
    if (key == "validation_snapshot_every") { config.validation_snapshot_every = std::stoi(val); continue; }
  }
  return true;
}

}  // namespace galaxy
