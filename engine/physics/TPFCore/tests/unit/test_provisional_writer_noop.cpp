#include "doctest.h"
#include "physics/TPFCore/tpf_core_package.hpp"
#include "physics/physics_package.hpp"
#include "config.hpp"
#include "types.hpp"

#include <cstdio>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {
bool exists(const std::string& path) {
  struct stat st;
  return ::stat(path.c_str(), &st) == 0;
}
}

TEST_CASE("closure diagnostics writer is no-op on this branch") {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";
  c.tpfcore_enable_provisional_readout = true;
  c.tpfcore_readout_mode = "derived_tpf_radial_readout";

  galaxy::TPFCorePackage p;
  p.init_from_config(c);

  galaxy::State s;
  s.resize(1);
  s.x[0] = 1.0; s.y[0] = 0.0; s.vx[0] = s.vy[0] = 0.0; s.mass[0] = 1.0;

  std::vector<double> ax, ay;
  p.compute_accelerations(s, 1.0e6, 1.0e-3, false, ax, ay);

  galaxy::Snapshot snap;
  snap.step = 0;
  snap.time = 0.0;
  snap.state = s;
  std::vector<galaxy::Snapshot> snaps(1, snap);

  char dir_template[] = "/tmp/tpf_closure_noop_XXXXXX";
  char* out_c = ::mkdtemp(dir_template);
  REQUIRE(out_c != nullptr);
  const std::string out_dir(out_c);

  p.write_closure_diagnostics(snaps, c, out_dir);

  CHECK_FALSE(exists(out_dir + "/tpf_closure_diagnostics.csv"));
  CHECK_FALSE(exists(out_dir + "/tpf_closure_diagnostics.txt"));

  ::rmdir(out_dir.c_str());
}

TEST_CASE("step0 orbit audit writer is inactive for v1 route on this branch") {
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";
  c.simulation_mode = galaxy::SimulationMode::bh_orbit_validation;
  c.tpfcore_enable_provisional_readout = true;

  galaxy::TPFCorePackage p;
  p.init_from_config(c);

  galaxy::State s;
  s.resize(1);
  s.x[0] = 1.0e8; s.y[0] = 0.0; s.vx[0] = 0.0; s.vy[0] = 0.0; s.mass[0] = 1.0;

  galaxy::Snapshot snap;
  snap.step = 0;
  snap.time = 0.0;
  snap.state = s;
  std::vector<galaxy::Snapshot> snaps(1, snap);

  char dir_template[] = "/tmp/tpf_step0_unavail_XXXXXX";
  char* out_c = ::mkdtemp(dir_template);
  REQUIRE(out_c != nullptr);
  const std::string out_dir(out_c);

  p.write_step0_orbit_audit(snaps, c, out_dir);
  CHECK_FALSE(exists(out_dir + "/tpf_step0_orbit_audit.txt"));
  ::rmdir(out_dir.c_str());
}
