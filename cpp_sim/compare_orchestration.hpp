#ifndef GALAXY_COMPARE_ORCHESTRATION_HPP
#define GALAXY_COMPARE_ORCHESTRATION_HPP

#include <string>

namespace galaxy {

bool should_run_compare_parallel(bool compare_parallel_flag,
                                 bool compare_mode_requested,
                                 bool compare_same_package,
                                 bool platform_supports_process_parallelism);

#ifndef _WIN32
bool child_exit_ok(int status, const char* side, std::string* detail);
#endif

}  // namespace galaxy

#endif
