#ifndef GALAXY_RENDER_AUDIT_HPP
#define GALAXY_RENDER_AUDIT_HPP

#include "config.hpp"
#include <string>

namespace galaxy {

struct GalaxyInitAudit;

/** TPF dynamics branch label: Newtonian, TPF_direct, TPF_legacy_readout[:mode], TPF_legacy_readout_plus_VDSG[:mode], or disabled. */
std::string compute_active_dynamics_branch(const Config& config);

/** Provisional readout / diagnostic identity (Theta→readout path); not necessarily ax,ay driver when VDSG on. */
std::string compute_active_metrics_branch(const Config& config);

/** C++ symbol-level path for audit_full overlay / manifest. */
std::string compute_acceleration_code_path(const Config& config);

/** JSON + text manifest beside run_info (galaxy mode recommended). */
void write_render_manifest(const std::string& output_dir,
                           const Config& config,
                           int n_steps_done,
                           int n_snapshots,
                           int n_particles,
                           const GalaxyInitAudit* galaxy_init_audit);

}  // namespace galaxy

#endif
