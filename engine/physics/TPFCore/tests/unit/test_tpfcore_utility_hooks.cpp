#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"

TEST_CASE("TPFCore utility hooks advertise supported utility modes") {
  galaxy::TPFCorePackage pkg;
  CHECK(pkg.supports_utility_mode(galaxy::SimulationMode::tpf_single_source_inspect));
  CHECK(pkg.supports_utility_mode(galaxy::SimulationMode::tpf_symmetric_pair_inspect));
  CHECK(pkg.supports_utility_mode(galaxy::SimulationMode::tpf_source_field_benchmark));
  CHECK(pkg.supports_utility_mode(galaxy::SimulationMode::tpf_4d_static_residual_benchmark));
  CHECK(pkg.supports_utility_mode(galaxy::SimulationMode::tpf_4d_static_motion_readout_benchmark));
  CHECK(pkg.supports_utility_mode(galaxy::SimulationMode::tpf_4d_xi_motion_probe_benchmark));
  CHECK(pkg.supports_utility_mode(galaxy::SimulationMode::tpf_weak_field_calibration));
  CHECK(pkg.supports_utility_mode(galaxy::SimulationMode::tpf_newtonian_force_compare));
  CHECK(pkg.supports_utility_mode(galaxy::SimulationMode::tpf_diagnostic_consistency_audit));
  CHECK_FALSE(pkg.supports_utility_mode(galaxy::SimulationMode::galaxy));
}

TEST_CASE("TPFCore post-run diagnostics hook is safe for empty snapshots") {
  galaxy::TPFCorePackage pkg;
  galaxy::Config cfg;
  cfg.physics_package = "TPFCore";
  cfg.output_dir = ".";
  std::vector<galaxy::Snapshot> snapshots;
  CHECK(pkg.write_post_run_diagnostics(snapshots, cfg, cfg.output_dir));
}
