#include "config.hpp"
#include "doctest.h"
#include "output.hpp"
#include "resolved_scenario.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace {

std::string slurp(const std::string& path) {
  std::ifstream f(path);
  std::ostringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

}  // namespace

TEST_CASE("resolve earth_moon_benchmark centralizes effective conditions") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::earth_moon_benchmark;

  galaxy::ResolvedScenario r = galaxy::resolve_scenario(c);

  CHECK(r.initializer_used == "init_two_body");
  CHECK(r.config.bh_mass == doctest::Approx(0.0));
  CHECK(r.config.enable_star_star_gravity == true);
  CHECK(r.effective_n_steps == 5000);
  CHECK(r.effective_snapshot_every == 5);
  REQUIRE(r.initial_state.n() == 2);
  CHECK(r.initial_state.mass[0] == doctest::Approx(galaxy::kDefaultEarthMassKg));
  CHECK(r.initial_state.mass[1] == doctest::Approx(galaxy::kDefaultMoonMassKg));
}

TEST_CASE("resolve bh_orbit_validation and symmetric_pair mode defaults") {
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::bh_orbit_validation;
    c.enable_star_star_gravity = true;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.initializer_used == "init_two_body_star_around_bh");
    CHECK(r.config.enable_star_star_gravity == false);
    CHECK(r.effective_n_steps == 5000);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::symmetric_pair;
    c.validation_symmetric_include_bh = false;
    c.bh_mass = 42.0;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.initializer_used == "init_symmetric_pair");
    CHECK(r.config.bh_mass == doctest::Approx(0.0));
    CHECK(r.effective_snapshot_every == 5);
  }
}

TEST_CASE("user overrides survive scenario resolution") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::earth_moon_benchmark;
  c.validation_n_steps = 111;
  c.validation_snapshot_every = 7;
  c.bh_mass = 123.0;

  auto r = galaxy::resolve_scenario(c);
  CHECK(r.config.bh_mass == doctest::Approx(123.0));
  CHECK(r.effective_n_steps == 111);
  CHECK(r.effective_snapshot_every == 7);
}

TEST_CASE("resolved scenario artifact writer outputs expected keys") {
  galaxy::Config c;
  c.output_dir = "outputs/test_resolved_scenario_artifact";
  c.simulation_mode = galaxy::SimulationMode::earth_moon_benchmark;
  c.physics_package = "Newtonian";
  galaxy::ResolvedScenario r = galaxy::resolve_scenario(c);

#ifdef _WIN32
  _mkdir("outputs");
  _mkdir(c.output_dir.c_str());
#else
  const int mk_ok = std::system((std::string("mkdir -p ") + c.output_dir).c_str());
  (void)mk_ok;
#endif

  galaxy::write_resolved_scenario_artifacts(c.output_dir, r);
  const std::string txt = slurp(c.output_dir + "/resolved_scenario.txt");
  const std::string js = slurp(c.output_dir + "/resolved_scenario.json");

  CHECK(txt.find("simulation_mode\tearth_moon_benchmark") != std::string::npos);
  CHECK(txt.find("effective_bh_mass\t0") != std::string::npos);
  CHECK(js.find("\"initializer_used\": \"init_two_body\"") != std::string::npos);
  CHECK(js.find("\"n_steps\": 5000") != std::string::npos);
}
