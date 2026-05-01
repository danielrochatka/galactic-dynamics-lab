#include "doctest.h"

#include "progress_time.hpp"

TEST_CASE("format_sim_time_for_progress honors preferred units") {
  CHECK(galaxy::format_sim_time_for_progress(3.15576e13, "Myr") == "1.00 Myr");
  CHECK(galaxy::format_sim_time_for_progress(3.15576e10, "kyr") == "1.00 kyr");
  CHECK(galaxy::format_sim_time_for_progress(3.15576e7, "yr") == "1.00 yr");
}

TEST_CASE("format_sim_time_for_progress auto-selects large units") {
  const std::string rendered = galaxy::format_sim_time_for_progress(3.0e15, "auto");
  CHECK(rendered.find("Myr") != std::string::npos);
}
