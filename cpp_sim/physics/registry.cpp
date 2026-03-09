#include "physics_package.hpp"
#include "Newtonian/newtonian.hpp"
#include "TPF/tpf_package.hpp"
#include <stdexcept>

namespace galaxy {

double compute_kinetic_energy(const State& state) {
  double ke = 0.0;
  for (int i = 0; i < state.n(); ++i)
    ke += 0.5 * state.mass[i] * (state.vx[i] * state.vx[i] + state.vy[i] * state.vy[i]);
  return ke;
}

namespace {

NewtonianPackage s_newtonian;
TPFPackage s_tpf;

PhysicsPackage* s_packages[] = {
  &s_newtonian,
  &s_tpf,
};
const int s_num_packages = sizeof(s_packages) / sizeof(s_packages[0]);

}  // namespace

PhysicsPackage* get_physics_package(const std::string& name) {
  for (int i = 0; i < s_num_packages; ++i) {
    if (name == s_packages[i]->name())
      return s_packages[i];
  }
  return nullptr;
}

bool has_physics_package(const std::string& name) {
  return get_physics_package(name) != nullptr;
}

}  // namespace galaxy
