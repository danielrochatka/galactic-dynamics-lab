#ifndef GALAXY_CONFIG_HPP
#define GALAXY_CONFIG_HPP

#include <string>
#include <utility>
#include <vector>

namespace galaxy {

/** SI speed of light (m/s). TPF bounce closure (Sec. IX, Eq. 31). */
constexpr double c = 299792458.0;

/** Solar mass [kg]; consistent with display_units / Earth–Moon SI literals. */
constexpr double kSolarMassKg = 1.98847e30;
/** Default SMBH mass for galaxy runs (~10^6 M_sun). */
constexpr double kDefaultBhMassKg = 1.0e6 * kSolarMassKg;
/** Default disk scale [m] ~10 kpc (galactic disk order). */
constexpr double kDefaultGalaxyRadiusM = 3.0e20;
/** Default inner disk radius [m] ~1 kpc. */
constexpr double kDefaultInnerRadiusM = 3.0e19;
/** Default Plummer-style softening [m] ~330 AU (numerical at galaxy scale). */
constexpr double kDefaultSofteningM = 1.0e16;
/** Generic engine default validation radius [m] (mode defaults may override with benchmark-specific values). */
constexpr double kDefaultValidationTwoBodyRadiusM = 5.0e18;
/** Earth–Moon benchmark defaults (SI). */
constexpr double kDefaultEarthMassKg = 5.972e24;
constexpr double kDefaultMoonMassKg = 7.348e22;
constexpr double kDefaultEarthMoonDistanceM = 3.844e8;
constexpr double kDefaultMoonTangentialSpeedMps = 1022.0;

// Simulation modes: galaxy + validation modes + TPFCore inspection/audit modes.
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

/** One key occurrence found in a config file (in file order). */
struct ConfigKeyOccurrence {
  int line_number = 0;
  std::string value;
};

// Load key=value pairs from a .cfg file into config. Returns true if file was read.
// Throws on malformed values (invalid number, etc.). Unknown keys are skipped.
bool load_config_file(const std::string& path, Config& config);

/** Apply one key=value (same rules as .cfg lines). Returns true if key was recognized. Throws on bad value. */
bool apply_config_kv(const std::string& key, const std::string& val, Config& config);

// Probe a config file for a single key. Returns value if found, else empty string.
// Returns empty if file does not exist or key not found.
std::string probe_config_key(const std::string& path, const std::string& key);

// Scan a config file and return all occurrences of `key` in file order, with line numbers.
// Returns empty if file does not exist or key is not present.
std::vector<ConfigKeyOccurrence> scan_config_key_occurrences(
    const std::string& path,
    const std::string& key);

// Find first existing run config path (root configs/ only). Returns empty if none found.
std::string find_run_config_path();

// If engine/configs/ exists: warn when using root configs, fail when run config was only there. Returns false to exit(1).
bool check_run_config_canonical(const std::string& run_config_path);

// Find package defaults path for given package name. Returns path if exists, else empty.
std::string find_package_defaults_path(const std::string& package_name);

/** Serialize all Config fields into stable key/value text pairs (enum values emitted as canonical strings). */
std::vector<std::pair<std::string, std::string>> serialize_config_kv(const Config& config);

struct Config {
  struct ExplicitOverrides {
    bool dt = false;
    bool softening = false;
    bool bh_mass = false;
    bool enable_star_star_gravity = false;
    bool tpfcore_source_softening = false;
    bool validation_two_body_radius = false;
    bool validation_two_body_speed_ratio = false;
    bool validation_earth_mass = false;
    bool validation_moon_mass = false;
    bool validation_earth_moon_distance = false;
    bool validation_moon_tangential_speed = false;
    bool validation_symmetric_include_bh = false;
    bool validation_symmetric_separation = false;
    bool validation_symmetric_speed = false;
    bool validation_small_n = false;
    bool validation_n_steps = false;
    bool validation_snapshot_every = false;
  } explicit_overrides;

  SimulationMode simulation_mode = SimulationMode::galaxy;

  int n_stars = 5000;
  double star_mass = kSolarMassKg;
  double bh_mass = kDefaultBhMassKg;

  double inner_radius = kDefaultInnerRadiusM;
  double outer_radius = kDefaultGalaxyRadiusM;
  /** Galaxy disk: stars sampled uniformly in area between 0.05 * galaxy_radius and galaxy_radius (m). */
  double galaxy_radius = kDefaultGalaxyRadiusM;

  double dt = 0.01;
  int n_steps = 50000;
  int snapshot_every = 50;

  double softening = 0.0;
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
   * - v11_weak_field_truncation: correspondence truncation implementation (Eq. 42-44 scalar superposition) using alpha_si.
   * - direct_tpf: principal-part implementation (Theta/I/kappa baseline, DeltaC omitted in current scope) with optional additive VDSG.
   */
  std::string tpf_dynamics_mode = "legacy_readout";
  /**
   * TPFCore correspondence-helper dynamics coupling alpha [SI] in Eq. (42)-(44): nabla^2 phi = alpha rho.
   * Used only by tpf_dynamics_mode=v11_weak_field_truncation (legacy weak_field_correspondence alias resolves there).
   * This is distinct from tpf_kappa (external key that maps to internal TPFCore kappa for direct_tpf and
   * also feeds derived-radial closure parameters where configured).
   */
  double tpf_weak_field_correspondence_alpha_si = -6.67430e-11;

  /**
   * TPFCore analysis/audit layer; does not replace dynamics.
   * none (default) | v11_weak_field_correspondence — static correspondence tensor audit outputs.
   */
  std::string tpf_analysis_mode = "none";

  /**
   * Static correspondence-audit geometry selector.
   * axis_monopole — +z axis tensor audit (default).
   * earth_moon_line_of_centers — Earth/Moon line-of-centers phi/acceleration CSV; no particle stepping.
   */
  std::string v11_weak_field_correspondence_benchmark = "axis_monopole";
  /** Paper Table II (v11): Earth mass (kg). */
  double v11_em_mass_earth_kg = 5.972e24;
  /** Paper Table II (v11): Moon mass (kg). */
  double v11_em_mass_moon_kg = 7.348e22;
  /** Paper Table II (v11): mean Earth–Moon center distance D (m). */
  double v11_em_mean_distance_m = 3.844e8;
  /**
   * Sidereal orbital period T (s) for Ω = 2π/T (paper Sec. XI D). If <= 0, use Kepler
   * T = 2π√(D³/(G(ME+MM))) with the same ME, MM, D and TPF_G_SI.
   */
  double v11_em_sidereal_period_s = 0.0;
  /** Calibration point x_cal (m from Earth center) for aTPF(x_cal) ≈ -calib_g (paper: Earth surface along line). */
  double v11_em_calib_surface_radius_m = 6.371e6;
  /** Target magnitude of surface acceleration (m/s²); paper uses ~9.81. */
  double v11_em_calib_surface_g_m_s2 = 9.81;

  /** TPFCore only: gate for legacy_readout dynamics (required for that path; ignored for direct_tpf). Nonzero VDSG can dominate final applied acceleration on legacy_readout. Default false. */
  bool tpfcore_enable_provisional_readout = false;
  /** TPFCore configured readout mode string (see provisional_readout.cpp). May appear in labels/metrics even when final acceleration includes other contributions. */
  std::string tpfcore_readout_mode = "tensor_radial_projection";
  /** TPFCore readout: magnitude scale factor. Default 1.0. */
  double tpfcore_readout_scale = 1.0;
  /** TPFCore derived-radial readout modes: theta_tt scale in diagnostics (derived closure; not added to ax, ay). Default 1.0. */
  double tpfcore_theta_tt_scale = 1.0;
  /** TPFCore derived-radial readout modes: theta_tr scale in diagnostics (not added to ax, ay). Default 1.0. */
  double tpfcore_theta_tr_scale = 1.0;
  /**
   * External config key for TPFCore coupling:
   * - external key name remains `tpf_kappa`
   * - maps to internal TPFCore `kappa`
   * - also feeds derived-radial closure parameters where configured
   * Default 1e32.
   */
  double tpf_kappa = 1.0e32;
  /**
   * TPFCore VDSG coupling λ: additive per-interaction velocity modifier on top of baseline acceleration.
   * Current implementation uses relative-speed scaling beta = |v_rel|/c with doppler_scale = 1 + λ_eff * beta.
   * λ_eff is mass-normalized (see tpf_vdsg_mass_baseline_kg).
   */
  double tpf_vdsg_coupling = 1.0e-20;
  /**
   * VDSG mass baseline M_ref (kg) for log normalization: λ_eff = λ · log10(M_ref) / log10(M_source).
   * If <= 0, uses star_mass at runtime (same units as simulation).
   */
  double tpf_vdsg_mass_baseline_kg = 0.0;
  /**
   * TPFCore: optional global |a| cap after acceleration assembly (fraction of |v|/dt per particle).
   * Not implied by tpf_vdsg_coupling. Default false.
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
   * TPFCore dynamical runs: fraction of n_steps using startup radial damping (1% per step).
   * Snapshots are not recorded during that interval. Ignored when physics_package != TPFCore.
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
  /** TPFCore inspection (single-source only): opt-in experimental planar Xi configuration-equation exterior solve dump. */
  bool tpf_xi_constraint_exterior_inspect = false;
  /** TPFCore inspection (single-source only): planar Xi exterior grid size N (NxN). */
  int tpf_xi_constraint_grid_n = 65;
  /** TPFCore inspection (single-source only): planar Xi exterior half-extent L for domain [-L, +L]. */
  double tpf_xi_constraint_half_extent = 10.0;
  /** TPFCore inspection (single-source only): inner excision/pinned radius for exterior solve. */
  double tpf_xi_constraint_inner_radius = 0.5;
  /** TPFCore inspection (single-source only): max iterations for experimental planar Xi solver. */
  int tpf_xi_constraint_max_iters = 250;
  /** TPFCore inspection (single-source only): convergence tolerance for experimental planar Xi solver. */
  double tpf_xi_constraint_tol = 1e-8;

  /** Galaxy mode template selector (see galaxy_init.hpp). May provide template defaults when related parameters remain neutral; explicit user-set values are not overwritten. Noise is controlled separately by galaxy_init_* noise keys and galaxy_init_master_chaos. */
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
   * Isotropic Cartesian velocity noise (fraction of v_circ per component) when all galaxy_init_* noise controls are zero.
   * If any of position/angle/magnitude noise (after master_chaos) is > 0, those controls determine perturbations.
   */
  double velocity_noise = 0.05;
  double initial_velocity_scale = 1.0;

  bool save_snapshots = true;
  bool save_run_info = true;
  /** Post-process hint for plot_cpp_run.py: smooth animation viewport vs per-frame target (velocity-gated). */
  bool plot_animation_dynamic_zoom = false;
  /** plot_cpp_compare.py: per-frame pooled-star coverage fraction for fixed shared compare viewport. */
  double plot_compare_smart_zoom_coverage = 0.95;
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
  /**
   * Display-unit controls for postprocess (plot_cpp_run.py, render.py, diagnostics.py).
   * Internal simulation, integration, snapshots, and numeric run_info values remain SI.
   * Supported distance: auto | m | km | AU | ly | pc | kpc
   * Supported time: auto | s | min | hr | day | yr | kyr | Myr
   * Supported velocity: auto | m/s | km/s
   */
  std::string display_distance_unit = "auto";
  std::string display_time_unit = "auto";
  std::string display_velocity_unit = "auto";
  /** Include active display units in overlays/captions when true. */
  bool display_units_in_overlay = true;
  /** Show a compact display-unit reference block when true. */
  bool display_show_unit_reference = true;

  // Validation-only
  double validation_two_body_radius = kDefaultValidationTwoBodyRadiusM;
  double validation_two_body_speed_ratio = 1.0;
  double validation_earth_mass = kDefaultEarthMassKg;
  double validation_moon_mass = kDefaultMoonMassKg;
  double validation_earth_moon_distance = kDefaultEarthMoonDistanceM;
  double validation_moon_tangential_speed = kDefaultMoonTangentialSpeedMps;
  bool validation_symmetric_include_bh = true;
  /** Half-separation along x (|x| of each star); full separation = 2 * this [m]. ~0.5 AU default. */
  double validation_symmetric_separation = 7.48e10;
  /** Tangential speed scale [m/s] (~30 km/s for AU-scale equal-mass binary order). */
  double validation_symmetric_speed = 3.0e4;
  int validation_small_n = 5;
  int validation_n_steps = 5000;
  int validation_snapshot_every = 5;

  std::string output_dir = "outputs";
  std::string run_id;
};

}  // namespace galaxy

#endif
