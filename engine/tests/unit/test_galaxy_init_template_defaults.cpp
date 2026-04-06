#include "doctest.h"
#include "galaxy_init.hpp"

using galaxy::apply_galaxy_init_template_defaults;
using galaxy::Config;
using galaxy::GalaxyInitTemplate;
using galaxy::GalaxyInitTemplateDefaultsLog;

TEST_CASE("weak_m2 default amplitude when zero") {
  Config c;
  c.galaxy_init_m2_amplitude = 0.0;
  GalaxyInitTemplateDefaultsLog log;
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::weak_m2, c, &log);
  CHECK(c.galaxy_init_m2_amplitude > 0.0);
  CHECK(!log.applied.empty());
}

TEST_CASE("weak_m3 default amplitude when zero") {
  Config c;
  c.galaxy_init_m3_amplitude = 0.0;
  GalaxyInitTemplateDefaultsLog log;
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::weak_m3, c, &log);
  CHECK(c.galaxy_init_m3_amplitude > 0.0);
}

TEST_CASE("weak_bar default bar amp and axis ratio") {
  Config c;
  c.galaxy_init_bar_amplitude = 0.0;
  c.galaxy_init_bar_axis_ratio = 1.0;
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::weak_bar, c, nullptr);
  CHECK(c.galaxy_init_bar_amplitude > 0.0);
  CHECK(c.galaxy_init_bar_axis_ratio > 1.01);
}

TEST_CASE("preformed_spiral default spiral parameters") {
  Config c;
  c.galaxy_init_spiral_amplitude = 0.0;
  c.galaxy_init_spiral_winding = 1.0;
  c.galaxy_init_spiral_phase = 0.0;
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::preformed_spiral, c, nullptr);
  CHECK(c.galaxy_init_spiral_amplitude > 0.0);
  CHECK(c.galaxy_init_spiral_winding > 1.01);
}

TEST_CASE("clumpy_disk default clumpiness") {
  Config c;
  c.galaxy_init_clumpiness = 0.0;
  c.galaxy_init_num_clumps = 8;
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::clumpy_disk, c, nullptr);
  CHECK(c.galaxy_init_clumpiness > 0.0);
  CHECK(c.galaxy_init_num_clumps >= 8);
}

TEST_CASE("symmetric_disk_noisy default noise when all zero") {
  Config c;
  c.galaxy_init_position_noise = 0.0;
  c.galaxy_init_velocity_angle_noise = 0.0;
  c.galaxy_init_velocity_magnitude_noise = 0.0;
  c.velocity_noise = 0.05;
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::symmetric_disk_noisy, c, nullptr);
  CHECK(c.galaxy_init_position_noise > 0.0);
  CHECK(c.galaxy_init_velocity_angle_noise > 0.0);
  CHECK(c.velocity_noise == doctest::Approx(0.0));
}

TEST_CASE("explicit m2 amplitude preserved") {
  Config c;
  c.galaxy_init_m2_amplitude = 0.12;
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::weak_m2, c, nullptr);
  CHECK(c.galaxy_init_m2_amplitude == doctest::Approx(0.12));
}

TEST_CASE("symmetric_disk leaves amplitudes zero") {
  Config c;
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::symmetric_disk, c, nullptr);
  CHECK(c.galaxy_init_m2_amplitude == doctest::Approx(0.0));
}

TEST_CASE("weak_m2 default differs from symmetric_disk neutral m2 amplitude") {
  Config sym;
  Config m2;
  m2.galaxy_init_m2_amplitude = 0.0;
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::symmetric_disk, sym, nullptr);
  apply_galaxy_init_template_defaults(GalaxyInitTemplate::weak_m2, m2, nullptr);
  CHECK(sym.galaxy_init_m2_amplitude == doctest::Approx(0.0));
  CHECK(m2.galaxy_init_m2_amplitude > sym.galaxy_init_m2_amplitude);
}
