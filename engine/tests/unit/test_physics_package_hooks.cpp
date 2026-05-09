#include "config.hpp"
#include "doctest.h"
#include "physics/Newtonian/newtonian.hpp"
#include "physics/physics_package.hpp"

namespace {

class MinimalHookPackage final : public galaxy::PhysicsPackage {
 public:
  const char* name() const override { return "MinimalHookPackage"; }

  void compute_accelerations(const galaxy::State& state,
                             double,
                             double,
                             bool,
                             std::vector<double>& ax,
                             std::vector<double>& ay) const override {
    ax.assign(state.n(), 0.0);
    ay.assign(state.n(), 0.0);
  }
};

galaxy::State two_body_state() {
  galaxy::State s;
  s.resize(2);
  s.x[0] = 1.0;
  s.y[0] = 0.0;
  s.mass[0] = 1.0;
  s.x[1] = -1.0;
  s.y[1] = 0.0;
  s.mass[1] = 2.0;
  return s;
}

}  // namespace

TEST_CASE("physics package default hooks are no-op/empty") {
  MinimalHookPackage pkg;
  galaxy::Config cfg;

  CHECK(pkg.display_name() == "MinimalHookPackage");
  CHECK(pkg.runtime_metadata().empty());
  CHECK_FALSE(pkg.supports_utility_mode(galaxy::SimulationMode::galaxy));
  CHECK_FALSE(pkg.run_utility_mode(cfg, ""));
  CHECK(pkg.run_info_metadata(cfg).empty());
  CHECK(pkg.render_metadata(cfg).empty());
  CHECK(pkg.package_config_metadata().empty());
}

TEST_CASE("newtonian package name resolves through registry") {
  galaxy::PhysicsPackage* newtonian = galaxy::get_physics_package("Newtonian");
  REQUIRE(newtonian != nullptr);
  CHECK(std::string(newtonian->name()) == "Newtonian");
}

TEST_CASE("calling generic hooks does not change newtonian acceleration behavior") {
  galaxy::NewtonianPackage pkg;
  galaxy::State s = two_body_state();

  std::vector<double> ax_before, ay_before;
  pkg.compute_accelerations(s, 5.0, 0.01, true, ax_before, ay_before);

  galaxy::Config cfg;
  (void)pkg.display_name();
  (void)pkg.runtime_metadata();
  (void)pkg.supports_utility_mode(galaxy::SimulationMode::galaxy);
  (void)pkg.run_utility_mode(cfg, "");
  (void)pkg.run_info_metadata(cfg);
  (void)pkg.render_metadata(cfg);
  (void)pkg.package_config_metadata();

  std::vector<double> ax_after, ay_after;
  pkg.compute_accelerations(s, 5.0, 0.01, true, ax_after, ay_after);

  REQUIRE(ax_before.size() == ax_after.size());
  REQUIRE(ay_before.size() == ay_after.size());
  for (std::size_t i = 0; i < ax_before.size(); ++i) {
    CHECK(ax_before[i] == doctest::Approx(ax_after[i]));
    CHECK(ay_before[i] == doctest::Approx(ay_after[i]));
  }
}
