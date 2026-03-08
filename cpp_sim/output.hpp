#ifndef GALAXY_OUTPUT_HPP
#define GALAXY_OUTPUT_HPP

#include "config.hpp"
#include "types.hpp"
#include <vector>
#include <string>

namespace galaxy {

// Write run_info.txt with key config (dt, n_steps, softening, etc.) and run summary.
void write_run_info(const std::string& output_dir,
                    const Config& config,
                    int n_steps_done,
                    int n_snapshots);

// Write each snapshot to output_dir/snapshot_<step>.csv (columns: i,x,y,vx,vy,mass).
void write_snapshots(const std::string& output_dir,
                     const std::vector<Snapshot>& snapshots);

}  // namespace galaxy

#endif
