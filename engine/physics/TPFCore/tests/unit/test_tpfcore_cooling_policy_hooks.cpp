#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"

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

TEST_CASE("TPFCore cooling hooks preserve damping and snapshot suppression semantics") {
  galaxy::Config c;
  c.tpf_cooling_fraction = 0.5;
  galaxy::TPFCorePackage pkg;

  galaxy::State s;
  s.resize(1);
  s.x[0] = 1.0;
  s.y[0] = 0.0;
  s.vx[0] = 1.0;
  s.vy[0] = 0.0;
  s.mass[0] = 1.0;

  const int cooling_steps = pkg.cooling_steps(c, 10);
  REQUIRE(cooling_steps == 5);
  pkg.apply_cooling_step(s, c, 1, cooling_steps);
  CHECK(s.vx[0] == doctest::Approx(0.99));
  CHECK(pkg.suppress_snapshot_for_cooling(c, 1, cooling_steps));
  CHECK_FALSE(pkg.suppress_snapshot_for_cooling(c, cooling_steps, cooling_steps));
}
