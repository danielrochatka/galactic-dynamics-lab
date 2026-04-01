#ifndef GALAXY_CONFIG_HPP
#define GALAXY_CONFIG_HPP

#include <string>

namespace galaxy {

/** SI speed of light (m/s). TPF bounce closure (Sec. IX, Eq. 31). */
constexpr double c = 299792458.0;

// Simulation mode: see config.py VALIDATION_MODES (Python) + "galaxy" + TPFCore inspection modes
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
  tpf_bound_orbit_sweep,
  /** v11 static weak-field correspondence audit (Ξ,Θ,I,C principal); not particle dynamics. */
  tpf_v11_weak_field_correspondence,
  /** Canonical Earth–Moon SI benchmark (same IC as legacy two_body_orbit string). */
  earth_moon_benchmark,
  /** Single star + central mass; uses init_two_body_star_around_bh. */
  bh_orbit_validation
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
  /**
   * Optional secondary package for declared side-by-side compare workflow.
   * Empty => single-package run (default). Non-empty + different from physics_package => compare mode.
   */
  std::string physics_package_compare = "";

  /**
   * TPFCore only: how dynamical accelerations are produced.
   * - legacy_readout (default): provisional readout closures (+ optional VDSG); requires tpfcore_enable_provisional_readout for dynamics.
   * - direct_tpf: future direct TPF dynamics path (not implemented yet); does not use the provisional-readout gate; VDSG must be 0.
   */
  std::string tpf_dynamics_mode = "legacy_readout";

  /**
   * TPFCore audit/analysis layer (no substitute for dynamics).
   * none (default) | v11_weak_field_correspondence — manuscript v11 static weak-field tensor correspondence only.
   */
  std::string tpf_analysis_mode = "none";

  /** TPFCore only: gate for legacy_readout dynamics (required for that path; ignored for direct_tpf). VDSG may own ax, ay when tpf_vdsg_coupling != 0. Default false. */
  bool tpfcore_enable_provisional_readout = false;
  /** TPFCore configured readout mode string (see provisional_readout.cpp). May label metrics while VDSG supersedes ax, ay. */
  std::string tpfcore_readout_mode = "tensor_radial_projection";
  /** TPFCore readout: scale factor for magnitude. Was 0.2046442 (INVALIDATED by benchmark mismatch). Pending recalibration; default 1.0. */
  double tpfcore_readout_scale = 1.0;
  /** TPFCore derived-radial readout modes: theta_tt scale in diagnostics (derived closure; not added to ax, ay). Default 1.0. */
  double tpfcore_theta_tt_scale = 1.0;
  /** TPFCore derived-radial readout modes: theta_tr scale in diagnostics (not added to ax, ay). Default 1.0. */
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
  /**
   * TPFCore: optional global |a| cap after readout + VDSG (fraction of |v|/dt per particle).
   * Default **false** so λ=0 runs are a clean readout baseline without this stabilization path.
   * Enable explicitly when investigating numerical caps; not implied by tpf_vdsg_coupling.
   */
  bool tpf_global_accel_shunt_enable = false;
  /** TPFCore: cap magnitude = tpf_global_accel_shunt_fraction * |v| / dt (when shunt enabled). Default 0.001 (0.1%). */
  double tpf_global_accel_shunt_fraction = 0.001;
  /**
   * TPFCore dynamical runs: write tpf_accel_pipeline_diagnostics.csv (per-snapshot pipeline metrics).
   * Default true; set false to skip the extra pass over snapshots at end of run.
   */
  bool tpf_accel_pipeline_diagnostics_csv = true;
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
  /** TPFCore diagnostics: enable live Newtonian-vs-TPF force audit for bh_orbit_validation runs. */
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

  /**
   * Galaxy mode: named IC template (see galaxy_init.hpp). Each template selects structured density
   * seeds (m2/m3/bar/spiral/clumps). When structured or noise parameters are still at their neutral
   * built-in defaults (0.0, 1.0 axis ratio, etc.), apply_galaxy_init_template_defaults supplies
   * modest template-specific values so the preset name matches visible structure. Explicit nonzero
   * user settings are never overwritten. Random “chaos” uses galaxy_init_* noise keys and
   * galaxy_init_master_chaos (noise-only scaler).
   */
  std::string galaxy_init_template = "symmetric_disk";
  /** mt19937 seed for galaxy placement and velocity perturbations (reproducible). Default 12345 matches legacy hardcoded seed. */
  unsigned galaxy_init_seed = 12345u;
  /**
   * Position jitter (Gaussian): delta_r ~ N(0, noise * (R - r_min)), delta_theta ~ N(0, noise * pi).
   * Scaled by galaxy_init_master_chaos only (not clump/mode amplitudes).
   */
  double galaxy_init_position_noise = 0.0;
  /** RMS rotation (rad) from pure tangential: mix tangent/radial directions by delta ~ N(0, noise^2). Scaled by master_chaos. */
  double galaxy_init_velocity_angle_noise = 0.0;
  /**
   * Fractional speed spread: v *= (1 + noise * N(0,1)) after direction set; floor 0.05 on scale to avoid sign flip.
   * Scaled by master_chaos.
   */
  double galaxy_init_velocity_magnitude_noise = 0.0;
  /** Probability [0,1] a star is drawn as a clump member (clumpy_disk). Not scaled by master_chaos. */
  double galaxy_init_clumpiness = 0.0;
  /** Number of clump centers (deterministic positions from RNG stream). */
  int galaxy_init_num_clumps = 8;
  /** Clump Gaussian sigma as fraction of galaxy_radius. */
  double galaxy_init_clump_radius_fraction = 0.08;
  /** Mode-2 azimuthal density weight ~ (1 + amp * cos(2 theta)). Weak seed only. Not scaled by master_chaos. */
  double galaxy_init_m2_amplitude = 0.0;
  /** Mode-3: (1 + amp * cos(3 theta)). */
  double galaxy_init_m3_amplitude = 0.0;
  /** Bar-like weight ~ (1 + amp * cos(2 theta)) * inner concentration; plus optional axis stretch. */
  double galaxy_init_bar_amplitude = 0.0;
  /** Bar stretch: x *= sqrt(ratio), y /= sqrt(ratio), then clipped to annulus; 1.0 = no stretch. */
  double galaxy_init_bar_axis_ratio = 1.0;
  /** Two-arm spiral envelope: weight ~ (1 + amp * cos(2*theta + winding*log(r/r_min) + phase)). */
  double galaxy_init_spiral_amplitude = 0.0;
  /** Dimensionless winding on log(r/r_min) in spiral phase. */
  double galaxy_init_spiral_winding = 1.0;
  /** Added to spiral phase (rad). */
  double galaxy_init_spiral_phase = 0.0;
  /**
   * Multiplies only galaxy_init_position_noise, galaxy_init_velocity_angle_noise,
   * galaxy_init_velocity_magnitude_noise. Does not scale clumpiness, m2/m3, bar, or spiral amplitudes.
   */
  double galaxy_init_master_chaos = 1.0;

  /**
   * Legacy isotropic Cartesian velocity noise (fraction of v_circ per component) when all new galaxy_init_* noises are zero.
   * If any of position/angle/magnitude noise (after master_chaos) is > 0, new pipeline is used instead for velocities.
   */
  double velocity_noise = 0.05;
  double initial_velocity_scale = 1.0;

  bool save_snapshots = true;
  bool save_run_info = true;
  /** Post-process hint for plot_cpp_run.py: smooth animation viewport vs per-frame target (velocity-gated). */
  bool plot_animation_dynamic_zoom = false;
  /** plot_cpp_run.py: plotting-only burn-in filter: ignore snapshots with step < this value. Default 0 (no filtering). */
  int plot_skip_initial_steps = 0;
  /** plot_cpp_run.py: plotting-only burn-in filter: after step filtering, drop first N remaining snapshots. Default 0 (no filtering). */
  int plot_skip_initial_snapshots = 0;
  /** plot_cpp_run.py diagnostics cutoff radius (m). <=0 means unset and plotter falls back to galaxy_radius if present. */
  double diagnostic_cutoff_radius = 0.0;
  /**
   * plot_cpp_run.py: text overlay on galaxy_initial/final PNG and animation frames.
   * none | minimal | audit_full. Default none keeps PNG/MP4 pixel-identical to pre-overlay behavior;
   * set minimal or audit_full in run config for audit labels.
   */
  std::string render_overlay_mode = "none";

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
