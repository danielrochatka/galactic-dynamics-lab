#pragma once

#include <string>

namespace galaxy {
std::string format_sim_time_for_progress(double sim_time_s, const std::string& preferred_unit);
}
