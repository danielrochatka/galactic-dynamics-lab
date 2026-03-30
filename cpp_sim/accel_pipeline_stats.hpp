#ifndef GALAXY_ACCEL_PIPELINE_STATS_HPP
#define GALAXY_ACCEL_PIPELINE_STATS_HPP

namespace galaxy {

/** Last-step TPFCore pipeline instrumentation (readout baseline + VDSG + optional |a| shunt). */
struct AccelPipelineStats {
  bool valid = false;
  /** Mean |a| from provisional readout baseline only (before VDSG add-on). */
  double mean_baseline_mag = 0.0;
  /** Mean |a| of the VDSG additive modifier alone (|dax,day|). */
  double mean_vdsg_mag = 0.0;
  /** mean_vdsg_mag / max(mean_baseline_mag, tiny); audit-only ratio. */
  double vdsg_over_baseline_ratio = 0.0;
  /** Shunt applications in the last integrator call (per-particle caps). */
  unsigned shunt_events_last_step = 0;
  /** shunt_events_last_step / n_particles (0 if n=0). */
  double frac_capped_last_step = 0.0;
  bool shunt_enabled = false;
  double shunt_fraction = 0.0;
};

}  // namespace galaxy

#endif
