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

TEST_CASE("compute_active_dynamics_branch: direct_tpf is canonical low-order truncated entry with strict offs") {
  Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";
  CHECK(galaxy::compute_active_dynamics_branch(c) ==
        "TPF_direct_tpf_canonical_entry_using_v11_low_order_static_quasistatic_truncation_DeltaC_omitted_"
        "VDSG_off_provisional_readout_off_shunt_off_cooling_off");
  CHECK(galaxy::compute_active_metrics_branch(c) ==
        "direct_tpf_metrics_v11_low_order_static_quasistatic_truncation_DeltaC_omitted_"
        "VDSG_off_provisional_readout_off_shunt_off_cooling_off");
}
