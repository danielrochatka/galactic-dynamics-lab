#include "config.hpp"
#include "doctest.h"
#include "render_audit.hpp"

TEST_CASE("render audit labels tpf_xi_theta_v1 as active route") {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";

  const std::string dyn = galaxy::compute_active_dynamics_branch(c);
  const std::string met = galaxy::compute_active_metrics_branch(c);
  const std::string acc = galaxy::compute_acceleration_code_path(c);

  CHECK(dyn.find("tpf_xi_theta_v1") != std::string::npos);
  CHECK(met.find("tpf_xi_theta_v1") != std::string::npos);
  CHECK(acc.find("Xi_total") != std::string::npos);
}

TEST_CASE("render audit marks non-v1 modes as unsupported on this branch") {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";

  CHECK(galaxy::compute_active_dynamics_branch(c).find("unsupported_non_v1") != std::string::npos);
  CHECK(galaxy::compute_active_metrics_branch(c).find("unsupported_non_v1") != std::string::npos);
  CHECK(galaxy::compute_acceleration_code_path(c).find("unsupported_non_v1") != std::string::npos);
}
