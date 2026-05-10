#include "config.hpp"
#include "doctest.h"
#include "physics/physics_package.hpp"

TEST_CASE("default package cooling hook is inactive with zero steps") {
  galaxy::Config c;
  galaxy::PhysicsPackage* pkg = galaxy::get_physics_package("Newtonian");
  REQUIRE(pkg != nullptr);
  CHECK_FALSE(pkg->cooling_active(c));
  CHECK(pkg->cooling_steps(c, 100) == 0);
}
