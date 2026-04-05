#include "compare_orchestration.hpp"
#include "doctest.h"

TEST_CASE("compare orchestration chooses process-parallel only when all gates are true") {
  CHECK(galaxy::should_run_compare_parallel(true, true, false, true));
  CHECK_FALSE(galaxy::should_run_compare_parallel(false, true, false, true));
  CHECK_FALSE(galaxy::should_run_compare_parallel(true, false, false, true));
  CHECK_FALSE(galaxy::should_run_compare_parallel(true, true, true, true));
  CHECK_FALSE(galaxy::should_run_compare_parallel(true, true, false, false));
}

#ifndef _WIN32
TEST_CASE("child exit status decoding propagates failures") {
  std::string detail;
  CHECK(galaxy::child_exit_ok(0, "left", &detail));

  detail.clear();
  const int fail_status = (7 << 8);  // exited with status 7
  CHECK_FALSE(galaxy::child_exit_ok(fail_status, "right", &detail));
  CHECK(detail.find("right child exited with code 7") != std::string::npos);
}
#endif
