#include "config.hpp"
#include "doctest.h"

using galaxy::apply_config_kv;
using galaxy::Config;

TEST_CASE("config defaults") {
  Config c;
  CHECK(c.simulation_mode == galaxy::SimulationMode::galaxy);
  CHECK(c.physics_package == "Newtonian");
  CHECK(c.tpf_vdsg_coupling == doctest::Approx(1.0e-20));
  CHECK(c.render_overlay_mode == "none");
}

TEST_CASE("legacy alias tpf_gdd_coupling maps to tpf_vdsg_coupling") {
  Config c;
  CHECK(apply_config_kv("tpf_gdd_coupling", "3.125e-10", c));
  CHECK(c.tpf_vdsg_coupling == doctest::Approx(3.125e-10));
}

TEST_CASE("render_overlay_mode parsing") {
  Config c;
  CHECK(apply_config_kv("render_overlay_mode", "minimal", c));
  CHECK(c.render_overlay_mode == "minimal");
  CHECK(apply_config_kv("render_overlay_mode", "audit_full", c));
  CHECK(c.render_overlay_mode == "audit_full");
}

TEST_CASE("galaxy init keys") {
  Config c;
  CHECK(apply_config_kv("galaxy_init_template", "preformed_spiral", c));
  CHECK(c.galaxy_init_template == "preformed_spiral");
  CHECK(apply_config_kv("galaxy_init_seed", "424242", c));
  CHECK(c.galaxy_init_seed == 424242u);
  CHECK(apply_config_kv("galaxy_init_master_chaos", "2.5", c));
  CHECK(c.galaxy_init_master_chaos == doctest::Approx(2.5));
}

TEST_CASE("explicit tpf_vdsg_coupling still works") {
  Config c;
  CHECK(apply_config_kv("tpf_vdsg_coupling", "1e-99", c));
  CHECK(c.tpf_vdsg_coupling == doctest::Approx(1e-99));
}
