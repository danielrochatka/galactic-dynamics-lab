#ifndef GALAXY_CONFIG_HPP
#define GALAXY_CONFIG_HPP

#include <string>

namespace galaxy {

// Simulation mode: matches Python VALIDATION_MODES + "galaxy" + TPFCore inspection modes
enum class SimulationMode {
  galaxy,
  two_body_orbit,
  symmetric_pair,
  small_n_conservation,
  timestep_convergence,
  tpf_single_source_inspect,
  tpf_symmetric_pair_inspect,
  tpf_single_source_optimize_c,
  tpf_two_body_sweep,
  tpf_weak_field_calibration
};

SimulationMode parse_mode(const std::string& s);

/** Return canonical string for mode (for display and run_info). */
std::string mode_to_string(SimulationMode m);

struct Config;

// Load key=value pairs from a .cfg file into config. Returns true if file was read.
// Throws on malformed values (invalid number, etc.). Unknown keys are skipped.
bool load_config_file(const std::string& path, Config& config);

// Probe a config file for a single key. Returns value if found, else empty string.
// Returns empty if file does not exist or key not found.
std::string probe_config_key(const std::string& path, const std::string& key);

// Find first existing run config path (root configs/ only). Returns empty if none found.
std::string find_run_config_path();

// If cpp_sim/configs/ exists: warn when using root configs, fail when run config was only there. Returns false to exit(1).
bool check_run_config_canonical(const std::string& run_config_path);

// Find package defaults path for given package name. Returns path if exists, else empty.
std::string find_package_defaults_path(const std::string& package_name);

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

  /** Physics package name (e.g. "Newtonian", "TPFCore"). Must match a registered package. Default: Newtonian. */
  std::string physics_package = "Newtonian";

  /** TPFCore only: enable provisional acceleration readout for dynamical modes. Default false. */
  bool tpfcore_enable_provisional_readout = false;
  /** TPFCore readout: mode. tensor_radial_projection (spatial) or tr_coherence_readout (paper t-r). */
  std::string tpfcore_readout_mode = "tensor_radial_projection";
  /** TPFCore readout: scale factor for magnitude. Default: weak-field calibrated effective scale (K_eff from tpf_weak_field_calibration; not proof of final TPF dynamics). */
  double tpfcore_readout_scale = 0.2046442;
  /** TPFCore tr_coherence_readout: Theta_tt balancing companion scale. Default 1.0. */
  double tpfcore_theta_tt_scale = 1.0;
  /** TPFCore tr_coherence_readout: Theta_tr mixed coupling scale. Default 1.0. */
  double tpfcore_theta_tr_scale = 1.0;
  /** TPFCore readout: dump debug CSV (tpf_readout_debug.csv) for dynamical runs. Default true. */
  bool tpfcore_dump_readout_debug = true;
  /** TPFCore inspection: probe radius min. */
  double tpfcore_probe_radius_min = 1.0;
  /** TPFCore inspection: probe radius max. */
  double tpfcore_probe_radius_max = 50.0;
  /** TPFCore inspection: number of probe samples along axis. */
  int tpfcore_probe_samples = 100;
  /** TPFCore inspection: write invariant_profile.csv. */
  bool tpfcore_dump_invariant_profile = true;
  /** TPFCore inspection: write theta_profile.csv. */
  bool tpfcore_dump_theta_profile = true;
  /** TPFCore source ansatz: softening for Phi. If <= 0, use global softening. */
  double tpfcore_source_softening = 0.0;
  /** TPFCore inspection: step size for numerical residual (if used). Not used when analytic. Default 1e-6. */
  double tpfcore_residual_step = 1e-6;
  /** TPFCore ansatz: isotropic correction B(r)=c*M/(r^2+eps^2)^(3/2). Default 0.0 reproduces pure Hessian. */
  double tpfcore_isotropic_correction_c = 0.0;
  /** TPFCore c-sweep: min c (default -0.5). */
  double tpfcore_c_sweep_min = -0.5;
  /** TPFCore c-sweep: max c (default 1.0). */
  double tpfcore_c_sweep_max = 1.0;
  /** TPFCore c-sweep: number of steps (default 31). */
  int tpfcore_c_sweep_steps = 31;
  /** TPFCore c-sweep: objective to minimize. max_residual_norm, mean_residual_norm, or l2_residual_norm. */
  std::string tpfcore_c_objective = "max_residual_norm";

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
