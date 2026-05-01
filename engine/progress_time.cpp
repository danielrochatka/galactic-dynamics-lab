#include "progress_time.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace {
constexpr double kDay = 86400.0;
constexpr double kYear = 365.25 * kDay;
constexpr double kKyr = 1e3 * kYear;
constexpr double kMyr = 1e6 * kYear;

std::pair<double, const char*> choose_unit(double sim_time_s, const std::string& preferred_unit) {
  if (preferred_unit == "min") return {60.0, "min"};
  if (preferred_unit == "hr") return {3600.0, "hr"};
  if (preferred_unit == "day") return {kDay, "day"};
  if (preferred_unit == "yr") return {kYear, "yr"};
  if (preferred_unit == "kyr") return {kKyr, "kyr"};
  if (preferred_unit == "Myr") return {kMyr, "Myr"};
  if (preferred_unit == "s") return {1.0, "s"};

  const double abs_t = std::fabs(sim_time_s);
  if (abs_t >= 0.5 * kMyr) return {kMyr, "Myr"};
  if (abs_t >= 0.5 * kKyr) return {kKyr, "kyr"};
  if (abs_t >= 2.0 * kYear) return {kYear, "yr"};
  if (abs_t >= 2.0 * kDay) return {kDay, "day"};
  if (abs_t >= 2.0 * 3600.0) return {3600.0, "hr"};
  if (abs_t >= 2.0 * 60.0) return {60.0, "min"};
  return {1.0, "s"};
}
}  // namespace

namespace galaxy {

std::string format_sim_time_for_progress(double sim_time_s, const std::string& preferred_unit) {
  const std::pair<double, const char*> unit = choose_unit(sim_time_s, preferred_unit);
  std::ostringstream os;
  os << std::fixed << std::setprecision(2) << (sim_time_s / unit.first) << " " << unit.second;
  return os.str();
}

}  // namespace galaxy
