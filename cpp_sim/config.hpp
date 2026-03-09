#ifndef GALAXY_CONFIG_HPP
#define GALAXY_CONFIG_HPP

#include <string>

namespace galaxy {

// Simulation mode: matches Python VALIDATION_MODES + "galaxy"
enum class SimulationMode {
  galaxy,
  two_body_orbit,
  symmetric_pair,
  small_n_conservation,
  timestep_convergence
};

SimulationMode parse_mode(const std::string& s);

struct Config;

// Load key=value pairs from a .cfg file into config. Returns true if file was read.
// Keys match Config member names (e.g. n_steps, simulation_mode). Unknown keys are skipped.
bool load_config_file(const std::string& path, Config& config);

struct Config {
  SimulationMode simulation_mode = SimulationMode::galaxy;

  int n_stars = 5000;
  double star_mass = 0.05;
  double bh_mass = 1000.0;

  double inner_radius = 5.0;
  double outer_radius = 50.0;

  double dt = 0.01;
  int n_steps = 50000;
  int snapshot_every = 50;

  double softening = 1.0;
  bool enable_star_star_gravity = true;

  /** Physics package name (e.g. "Newtonian"). Must match a registered package. Default: Newtonian. */
  std::string physics_package = "Newtonian";

  double velocity_noise = 0.05;
  double initial_velocity_scale = 1.0;

  bool save_snapshots = true;
  bool save_run_info = true;

  // Validation-only
  double validation_two_body_radius = 20.0;
  double validation_two_body_speed_ratio = 1.0;
  bool validation_symmetric_include_bh = true;
  double validation_symmetric_separation = 30.0;
  double validation_symmetric_speed = 4.0;
  int validation_small_n = 5;
  int validation_n_steps = 5000;
  int validation_snapshot_every = 5;

  std::string output_dir = "outputs";
  std::string run_id;
};

}  // namespace galaxy

#endif
