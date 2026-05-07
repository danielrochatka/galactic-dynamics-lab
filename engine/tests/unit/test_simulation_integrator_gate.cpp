#include "config.hpp"
#include "doctest.h"
#include "integrator.hpp"
#include "physics/physics_package.hpp"
#include "simulation.hpp"

namespace {

class HarmonicPhysicsPackage final : public galaxy::PhysicsPackage {
 public:
  const char* name() const override { return "HarmonicTest"; }
  void init_from_config(const galaxy::Config&) override {}

  void compute_accelerations(const galaxy::State& state,
                             double,
                             double,
                             bool,
                             std::vector<double>& ax,
                             std::vector<double>& ay) const override {
    ax.assign(state.n(), 0.0);
    ay.assign(state.n(), 0.0);
    for (int i = 0; i < state.n(); ++i) {
      ax[i] = -state.x[i];
      ay[i] = -state.y[i];
    }
  }

  double compute_potential_energy(const galaxy::State&, double, double, bool) const override { return 0.0; }
};

galaxy::State one_body_state() {
  galaxy::State s;
  s.resize(1);
  s.x[0] = 1.0;
  s.y[0] = 0.0;
  s.vx[0] = 0.4;
  s.vy[0] = 0.0;
  s.mass[0] = 1.0;
  return s;
}

}  // namespace

TEST_CASE("run_simulation uses semi-implicit Euler for active tpf_xi_theta_v1 Xi-kernel deformation") {
  HarmonicPhysicsPackage physics;
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";
  c.tpf_4d_xi_kernel_mode = "metric_velocity";
  c.tpf_4d_xi_kernel_coupling = 1.0;
  c.dt = 0.1;

  galaxy::State initial = one_body_state();
  std::vector<galaxy::Snapshot> out = galaxy::run_simulation(c, initial, &physics, 1, 1, nullptr, 0);
  REQUIRE(out.size() == 2);

  galaxy::State expected = one_body_state();
  std::vector<double> ax, ay;
  galaxy::semi_implicit_euler_step(expected, &physics, c.bh_mass, c.softening, c.enable_star_star_gravity, c.dt, ax, ay);

  CHECK(out[1].state.x[0] == doctest::Approx(expected.x[0]));
  CHECK(out[1].state.vx[0] == doctest::Approx(expected.vx[0]));
}

TEST_CASE("run_simulation uses Velocity-Verlet for tpf_xi_theta_v1 when Xi-kernel mode is off") {
  HarmonicPhysicsPackage physics;
  galaxy::Config c;
  c.physics_package = "TPFCore";
  c.tpf_dynamics_mode = "tpf_xi_theta_v1";
  c.tpf_4d_xi_kernel_mode = "off";
  c.tpf_4d_xi_kernel_coupling = 1.0;
  c.dt = 0.1;

  galaxy::State initial = one_body_state();
  std::vector<galaxy::Snapshot> out = galaxy::run_simulation(c, initial, &physics, 1, 1, nullptr, 0);
  REQUIRE(out.size() == 2);

  galaxy::State expected = one_body_state();
  std::vector<double> ax, ay;
  galaxy::velocity_verlet_step(expected, &physics, c.bh_mass, c.softening, c.enable_star_star_gravity, c.dt, ax, ay);

  CHECK(out[1].state.x[0] == doctest::Approx(expected.x[0]));
  CHECK(out[1].state.vx[0] == doctest::Approx(expected.vx[0]));
}
