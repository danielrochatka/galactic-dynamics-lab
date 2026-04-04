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
  CHECK(r.config.dt == doctest::Approx(3600.0));
  CHECK(r.config.softening == doctest::Approx(0.0));
  CHECK(r.effective_n_steps == 1440);
  CHECK(r.effective_snapshot_every == 6);
  CHECK(r.effective_total_sim_time == doctest::Approx(3600.0 * 1440.0));
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
    CHECK(r.config.dt == doctest::Approx(10000.0));
    CHECK(r.config.softening == doctest::Approx(0.0));
    CHECK(r.config.validation_two_body_radius == doctest::Approx(1.0e13));
    CHECK(r.effective_n_steps == 6000);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::symmetric_pair;
    c.validation_symmetric_include_bh = false;
    c.bh_mass = 42.0;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.initializer_used == "init_symmetric_pair");
    CHECK(r.config.bh_mass == doctest::Approx(0.0));
    CHECK(r.config.dt == doctest::Approx(3600.0));
    CHECK(r.config.softening == doctest::Approx(0.0));
    CHECK(r.effective_snapshot_every == 6);
  }
}

TEST_CASE("resolve small_n_conservation and timestep_convergence use explicit mode timing/softening") {
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::small_n_conservation;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.config.dt == doctest::Approx(1.0e-4));
    CHECK(r.config.softening == doctest::Approx(0.0));
    CHECK(r.config.bh_mass == doctest::Approx(1.0e18));
    CHECK(r.effective_n_steps == 20000);
    CHECK(r.effective_snapshot_every == 20);
  }
  {
    galaxy::Config c;
    c.simulation_mode = galaxy::SimulationMode::timestep_convergence;
    auto r = galaxy::resolve_scenario(c);
    CHECK(r.config.dt == doctest::Approx(10000.0));
    CHECK(r.config.softening == doctest::Approx(0.0));
    CHECK(r.config.validation_two_body_radius == doctest::Approx(1.0e13));
    CHECK(r.effective_n_steps == 6000);
    CHECK(r.effective_snapshot_every == 10);
  }
}

TEST_CASE("user overrides survive scenario resolution") {
  galaxy::Config c;
  c.simulation_mode = galaxy::SimulationMode::earth_moon_benchmark;
  c.validation_n_steps = 111;
  c.validation_snapshot_every = 7;
  c.bh_mass = 123.0;
  c.dt = 12.5;
  c.softening = 2.0;

  auto r = galaxy::resolve_scenario(c);
  CHECK(r.config.bh_mass == doctest::Approx(123.0));
  CHECK(r.config.dt == doctest::Approx(12.5));
  CHECK(r.config.softening == doctest::Approx(2.0));
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
  CHECK(txt.find("effective_dt\t3600") != std::string::npos);
  CHECK(txt.find("effective_softening\t0") != std::string::npos);
  CHECK(txt.find("effective_total_sim_time\t") != std::string::npos);
  CHECK(txt.find("timing_policy\tearth_moon_hourly_step_60d_horizon") != std::string::npos);
  CHECK(js.find("\"initializer_used\": \"init_two_body\"") != std::string::npos);
  CHECK(js.find("\"n_steps\": 1440") != std::string::npos);
  CHECK(js.find("\"total_sim_time\": 5184000") != std::string::npos);
  CHECK(js.find("\"softening\": 0") != std::string::npos);
  CHECK(js.find("\"timing\": \"earth_moon_hourly_step_60d_horizon\"") != std::string::npos);
}
