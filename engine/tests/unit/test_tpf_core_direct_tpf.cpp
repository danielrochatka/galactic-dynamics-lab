#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"

TEST_CASE("direct_tpf runtime route is rejected under strict tpf_xi_theta_v1 branch policy") {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "direct_tpf";

  galaxy::TPFCorePackage pkg;
  CHECK_THROWS_WITH(pkg.init_from_config(c),
                    "TPFCore on this branch supports only tpf_dynamics_mode=tpf_xi_theta_v1.");
}
