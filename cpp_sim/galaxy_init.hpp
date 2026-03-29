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
 * density seeds; noise is controlled separately via galaxy_init_* noise keys and
 * galaxy_init_master_chaos (noise only — not m2/m3/bar/spiral/clump strengths).
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

  double weight_w_max = 1.0;
  int rejection_fallbacks = 0;

  std::vector<std::pair<double, double>> clump_centers_xy;
};

/** Last galaxy initialization audit (galaxy mode). Invalid until initialize_galaxy_disk runs. */
const GalaxyInitAudit& last_galaxy_init_audit();

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
