#include "compare_orchestration.hpp"

#ifndef _WIN32
#include <sys/wait.h>
#include <sstream>
#endif

namespace galaxy {

bool should_run_compare_parallel(bool compare_parallel_flag,
                                 bool compare_mode_requested,
                                 bool compare_same_package,
                                 bool platform_supports_process_parallelism) {
  return compare_parallel_flag && compare_mode_requested && !compare_same_package &&
         platform_supports_process_parallelism;
}

#ifndef _WIN32
bool child_exit_ok(int status, const char* side, std::string* detail) {
  std::ostringstream os;
  if (WIFEXITED(status)) {
    const int code = WEXITSTATUS(status);
    if (code == 0) return true;
    os << side << " child exited with code " << code;
  } else if (WIFSIGNALED(status)) {
    os << side << " child terminated by signal " << WTERMSIG(status);
  } else {
    os << side << " child ended with unexpected wait status " << status;
  }
  if (detail) *detail = os.str();
  return false;
}
#endif

}  // namespace galaxy
