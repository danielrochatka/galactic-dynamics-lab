#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "types.hpp"

#include <vector>

namespace {
galaxy::State one_body_state(double x, double y, double m) {
  galaxy::State s;
  s.resize(1);
  s.x[0] = x; s.y[0] = y; s.vx[0] = 0.0; s.vy[0] = 0.0; s.mass[0] = m;
  return s;
}
}

TEST_CASE("tpf_xi_theta_v1 is accepted and computes acceleration") {
  galaxy::Config c;
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";
  c.tpf_4d_xi_kernel_mode = "off";
  galaxy::TPFCorePackage p;
  p.init_from_config(c);
  galaxy::State s = one_body_state(10.0, 0.0, 1.0);
  std::vector<double> ax, ay;
  p.compute_accelerations(s, 100.0, 0.1, false, ax, ay);
  CHECK(ax.size() == 1);
  CHECK(ay.size() == 1);
}

TEST_CASE("non-v1 dynamics route names are rejected") {
  const char* modes[] = {
      "geodesic_correspondence",
      "direct_tpf",
      "v11_weak_field_truncation",
      "legacy_readout",
      "xi_kernel_deformed",
  };
  for (const char* m : modes) {
    galaxy::Config c;
    CHECK_THROWS(galaxy::apply_config_kv("tpf_dynamics_mode", m, c));
  }
}
