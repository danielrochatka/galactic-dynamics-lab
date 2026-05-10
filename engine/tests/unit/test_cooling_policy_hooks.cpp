#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "physics/physics_package.hpp"

TEST_CASE("default package cooling hook is inactive with zero steps") {
  galaxy::Config c;
  galaxy::PhysicsPackage* pkg = galaxy::get_physics_package("Newtonian");
  REQUIRE(pkg != nullptr);
  CHECK_FALSE(pkg->cooling_active(c));
  CHECK(pkg->cooling_steps(c, 100) == 0);
}

TEST_CASE("TPFCore cooling hook is inactive when fraction <= 0") {
  galaxy::Config c;
  c.tpf_cooling_fraction = 0.0;
  galaxy::TPFCorePackage pkg;
  CHECK_FALSE(pkg.cooling_active(c));
  CHECK(pkg.cooling_steps(c, 100) == 0);

  c.tpf_cooling_fraction = -0.5;
  CHECK_FALSE(pkg.cooling_active(c));
  CHECK(pkg.cooling_steps(c, 100) == 0);
}

TEST_CASE("TPFCore cooling hook returns clamped step count when fraction > 0") {
  galaxy::Config c;
  galaxy::TPFCorePackage pkg;

  c.tpf_cooling_fraction = 0.25;
  CHECK(pkg.cooling_active(c));
  CHECK(pkg.cooling_steps(c, 100) == 25);

  c.tpf_cooling_fraction = 2.0;
  CHECK(pkg.cooling_active(c));
  CHECK(pkg.cooling_steps(c, 100) == 100);
}
