#include "doctest.h"
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

/* Regression: snapshot_*.csv must not use default float formatting (6 digits), which collapses
 * sub-50 m motion at ~1e8 m scales and corrupts long-run two-body diagnostics. */

TEST_CASE("snapshot CSV scientific format: low precision hides EM-scale motion") {
  const double x0 = 3.843663e8;
  const double x1 = x0 + 30.0;  // 30 m — typical between snapshots is smaller; still rounds away
  std::ostringstream a, b;
  a << std::scientific << std::setprecision(6) << x0;
  b << std::scientific << std::setprecision(6) << x1;
  CHECK(a.str() == b.str());
}

TEST_CASE("snapshot CSV scientific format: max_digits10 preserves sub-meter deltas at 1e8 m") {
  const double x0 = 3.843663e8;
  const double x1 = x0 + 0.25;
  std::ostringstream u, v;
  u << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << x0;
  v << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << x1;
  CHECK(u.str() != v.str());
}
