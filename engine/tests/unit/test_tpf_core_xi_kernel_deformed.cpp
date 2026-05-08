#include "config.hpp"
#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"

TEST_CASE("xi_kernel_deformed runtime route is rejected under strict tpf_xi_theta_v1 branch policy") {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "xi_kernel_deformed";

  galaxy::TPFCorePackage pkg;
  CHECK_THROWS_WITH(pkg.init_from_config(c),
                    "TPFCore on this branch supports only tpf_dynamics_mode=tpf_xi_theta_v1.");
}
