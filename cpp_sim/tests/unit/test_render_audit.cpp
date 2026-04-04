#include "config.hpp"
#include "doctest.h"
#include "render_audit.hpp"

using galaxy::Config;

TEST_CASE("compute_active_dynamics_branch: Newtonian") {
  Config c;
  c.physics_package = "Newtonian";
  CHECK(galaxy::compute_active_dynamics_branch(c) == "Newtonian_pairwise_G_SI");
}

TEST_CASE("compute_active_dynamics_branch: unknown package is tagged non-TPFCore") {
  Config c;
  c.physics_package = "CustomLabModel";
  CHECK(galaxy::compute_active_dynamics_branch(c) == "CustomLabModel (non-TPFCore)");
  CHECK(galaxy::compute_active_metrics_branch(c) == "unknown");
  CHECK(galaxy::compute_acceleration_code_path(c) == "unknown_package");
}
