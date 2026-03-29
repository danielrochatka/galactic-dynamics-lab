#ifndef GALAXY_CONFIG_HPP
#define GALAXY_CONFIG_HPP

#include <string>

namespace galaxy {

/** SI speed of light (m/s). TPF bounce closure (Sec. IX, Eq. 31). */
constexpr double c = 299792458.0;

// Simulation mode: matches Python VALIDATION_MODES + "galaxy" + TPFCore inspection modes
enum class SimulationMode {
  galaxy,
  two_body_orbit,
  symmetric_pair,
  small_n_conservation,
  timestep_convergence,
  tpf_single_source_inspect,
  tpf_symmetric_pair_inspect,
  tpf_two_body_sweep,
  tpf_weak_field_calibration,
  tpf_newtonian_force_compare,
  tpf_diagnostic_consistency_audit,
  tpf_bound_orbit_sweep
};

SimulationMode parse_mode(const std::string& s);

/** Return canonical string for mode (for display and run_info). */
std::string mode_to_string(SimulationMode m);

struct Config;

// Load key=value pairs from a .cfg file into config. Returns true if file was read.
// Throws on malformed values (invalid number, etc.). Unknown keys are skipped.
bool load_config_file(const std::string& path, Config& config);

/** Apply one key=value (same rules as .cfg lines). Returns true if key was recognized. Throws on bad value. */
bool apply_config_kv(const std::string& key, const std::string& val, Config& config);

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
  /** Galaxy disk: stars sampled uniformly in area between 0.05 * galaxy_radius and galaxy_radius (m). */
  double galaxy_radius = 50.0;

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
  /** TPFCore readout: scale factor for magnitude. Was 0.2046442 (INVALIDATED by benchmark mismatch). Pending recalibration; default 1.0. */
  double tpfcore_readout_scale = 1.0;
  /** TPFCore tr_coherence_readout: Theta_tt balancing companion scale. Default 1.0. */
  double tpfcore_theta_tt_scale = 1.0;
  /** TPFCore tr_coherence_readout: Theta_tr mixed coupling scale. Default 1.0. */
  double tpfcore_theta_tr_scale = 1.0;
  /**
   * TPFCore hybrid radial: coupling κ from geometric invariant I to effective mass density ρ_eff = κ I (SI ledger).
   * Default 1e32.
   */
  double tpf_kappa = 1.0e32;
  /**
   * TPFCore VDSG: global strength λ in doppler_scale = 1 + λ_eff (v·r̂)/c per interaction.
   * λ_eff is mass-normalized (see tpf_vdsg_mass_baseline_kg). Total acceleration still passes the
   * 0.1% |v|/dt shunt after summing neighbors. Default tuned with TPF_G_SI.
   */
  double tpf_vdsg_coupling = 1.0e-20;
  /**
   * VDSG mass baseline M_ref (kg) for log normalization: λ_eff = λ · log10(M_ref) / log10(M_source).
   * If <= 0, uses star_mass at runtime (same units as simulation). Heuristic / provisional closure.
   */
  double tpf_vdsg_mass_baseline_kg = 0.0;
  /** TPFCore derived radial profile: number of radial bins (grid for diagnostics / interpolation). Default 100. */
  int tpf_poisson_bins = 100;
  /** TPFCore derived radial profile outer radius (m); <=0 uses galaxy_radius. */
  double tpf_poisson_max_radius = 0.0;
  /**
   * TPFCore dynamical runs: fraction of n_steps in the artificial radial cooling phase (1% radial
   * damping per step). Snapshots are not recorded during this phase. Ignored when physics_package != TPFCore.
   */
  double tpf_cooling_fraction = 0.2;
  /** TPFCore readout: dump debug CSV (tpf_readout_debug.csv) for dynamical runs. Default true. */
  bool tpfcore_dump_readout_debug = true;
  /** TPFCore diagnostics: enable live two_body_orbit Newtonian-vs-TPF force audit. */
  bool tpfcore_live_orbit_force_audit = false;
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

  double velocity_noise = 0.05;
  double initial_velocity_scale = 1.0;

  bool save_snapshots = true;
  bool save_run_info = true;
  /** Post-process hint for plot_cpp_run.py: smooth animation viewport vs per-frame target (velocity-gated). */
  bool plot_animation_dynamic_zoom = true;

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
