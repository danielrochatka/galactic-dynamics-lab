#ifndef GALAXY_OUTPUT_HPP
#define GALAXY_OUTPUT_HPP

#include "config.hpp"
#include "types.hpp"
#include <vector>
#include <string>

namespace galaxy {

struct GalaxyInitAudit;

// Write run_info.txt with key config (dt, n_steps, softening, etc.) and run summary.
// n_particles: actual particle count (if < 0, use config.n_stars).
// run_config_path, package_defaults_path: if non-empty, written for reproducibility.
// galaxy_init_audit: if non-null and valid, writes resolved galaxy IC block (galaxy mode).
void write_run_info(const std::string& output_dir,
                    const Config& config,
                    int n_steps_done,
                    int n_snapshots,
                    int n_particles = -1,
                    const std::string& run_config_path = "",
                    const std::string& package_defaults_path = "",
                    const GalaxyInitAudit* galaxy_init_audit = nullptr);

// Write each snapshot to output_dir/snapshot_<step>.csv (columns: i,x,y,vx,vy,mass).
void write_snapshots(const std::string& output_dir,
                     const std::vector<Snapshot>& snapshots);

}  // namespace galaxy

#endif
