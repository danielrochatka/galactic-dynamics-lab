#ifndef GALAXY_GALAXY_INIT_HPP
#define GALAXY_GALAXY_INIT_HPP

#include "config.hpp"
#include "types.hpp"
#include <string>
#include <utility>
#include <vector>

namespace galaxy {

/**
 * Named galaxy IC recipes (galaxy mode only). Each enables a subset of structured
 * density seeds (m2/m3/bar/spiral/clumps). **Template-specific defaults** are applied
 * when structured/noise parameters are still at neutral Config defaults so the template
 * name matches visible structure (see apply_galaxy_init_template_defaults). User-set
 * nonzero values are never overwritten.
 */
enum class GalaxyInitTemplate {
  symmetric_disk,
  symmetric_disk_noisy,
  clumpy_disk,
  weak_m2,
  weak_m3,
  weak_bar,
  preformed_spiral
};

/** Parse template name (case-insensitive). Throws on unknown. */
GalaxyInitTemplate parse_galaxy_init_template(const std::string& s);

std::string galaxy_init_template_to_string(GalaxyInitTemplate t);

/** Log lines from apply_galaxy_init_template_defaults (for audit / stderr). */
struct GalaxyInitTemplateDefaultsLog {
  std::vector<std::string> applied;
  std::vector<std::string> warnings;
};

/** Filled during initialize_galaxy_disk for run_info and diagnostics. */
struct GalaxyInitAudit {
  bool valid = false;
  std::string template_name;
  unsigned seed = 0;
  double master_chaos = 1.0;
  bool master_scales_position_noise = false;
  bool master_scales_velocity_angle_noise = false;
  bool master_scales_velocity_magnitude_noise = false;

  double raw_position_noise = 0.0;
  double raw_velocity_angle_noise = 0.0;
  double raw_velocity_magnitude_noise = 0.0;
  double eff_position_noise = 0.0;
  double eff_velocity_angle_noise_rad = 0.0;
  double eff_velocity_magnitude_noise = 0.0;

  bool used_new_state_noise = false;
  bool used_legacy_velocity_noise = false;

  bool structured_m2 = false;
  bool structured_m3 = false;
  bool structured_bar = false;
  bool structured_spiral = false;
  bool structured_clumps = false;

  /** Raw config file values (before template defaults). */
  double galaxy_init_clumpiness_raw = 0.0;
  int galaxy_init_num_clumps_raw = 0;
  double galaxy_init_clump_radius_fraction_raw = 0.0;
  double galaxy_init_m2_amplitude_raw = 0.0;
  double galaxy_init_m3_amplitude_raw = 0.0;
  double galaxy_init_bar_amplitude_raw = 0.0;
  double galaxy_init_bar_axis_ratio_raw = 1.0;
  double galaxy_init_spiral_amplitude_raw = 0.0;
  double galaxy_init_spiral_winding_raw = 0.0;
  double galaxy_init_spiral_phase_raw = 0.0;
  double velocity_noise_raw = 0.05;

  /** Effective values used for placement (after template defaults). */
  double galaxy_init_clumpiness = 0.0;
  int galaxy_init_num_clumps = 0;
  double galaxy_init_clump_radius_fraction = 0.0;
  double galaxy_init_m2_amplitude = 0.0;
  double galaxy_init_m3_amplitude = 0.0;
  double galaxy_init_bar_amplitude = 0.0;
  double galaxy_init_bar_axis_ratio = 1.0;
  double galaxy_init_spiral_amplitude = 0.0;
  double galaxy_init_spiral_winding = 0.0;
  double galaxy_init_spiral_phase = 0.0;
  /** Legacy velocity_noise after resolution (often 0 when new-style noise is used). */
  double velocity_noise_effective = 0.05;
  /** galaxy_init_* noise after template defaults (for sync to Config / run_info). */
  double galaxy_init_position_noise_resolved = 0.0;
  double galaxy_init_velocity_angle_noise_resolved = 0.0;
  double galaxy_init_velocity_magnitude_noise_resolved = 0.0;

  bool template_defaults_used = false;
  GalaxyInitTemplateDefaultsLog template_defaults_log;

  double weight_w_max = 1.0;
  int rejection_fallbacks = 0;

  std::vector<std::pair<double, double>> clump_centers_xy;
};

/** Last galaxy initialization audit (galaxy mode). Invalid until initialize_galaxy_disk runs. */
const GalaxyInitAudit& last_galaxy_init_audit();

/**
 * Apply template-specific defaults to `effective` when structured/noise fields are still
 * at neutral built-in defaults. Does not modify `symmetric_disk`. Logs human-readable lines
 * to `log_out` and may print warnings to stderr.
 */
void apply_galaxy_init_template_defaults(GalaxyInitTemplate tmpl, Config& effective,
                                         GalaxyInitTemplateDefaultsLog* log_out = nullptr);

/**
 * After initialize_galaxy_disk, copy effective IC parameters from the last audit into
 * `config` so run_info / render_manifest match what was used for placement.
 */
void sync_config_galaxy_init_from_last_audit(Config& config);

/**
 * Full galaxy disk initialization: templates, structured seeds, reproducible RNG,
 * physically sensible tangential speeds + optional perturbations.
 */
void initialize_galaxy_disk(const Config& config, State& state, GalaxyInitAudit* audit_out = nullptr);

/** Text + optional CSV of initial state statistics (galaxy mode). */
void write_galaxy_init_diagnostics(const std::string& output_dir,
                                   const State& state,
                                   const Config& config,
                                   const GalaxyInitAudit& audit);

}  // namespace galaxy

#endif
