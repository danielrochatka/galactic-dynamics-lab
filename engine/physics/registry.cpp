#include "physics_package.hpp"
#include <map>
#include <utility>

namespace galaxy {

double compute_kinetic_energy(const State& state) {
  double ke = 0.0;
  for (int i = 0; i < state.n(); ++i)
    ke += 0.5 * state.mass[i] * (state.vx[i] * state.vx[i] + state.vy[i] * state.vy[i]);
  return ke;
}

namespace {

std::map<std::string, PhysicsPackageFactory>& package_factories() {
  static std::map<std::string, PhysicsPackageFactory> factories;
  return factories;
}

std::map<std::string, std::unique_ptr<PhysicsPackage>>& package_instances() {
  static std::map<std::string, std::unique_ptr<PhysicsPackage>> instances;
  return instances;
}

}  // namespace

PhysicsPackage* get_physics_package(const std::string& name) {
  auto& instances = package_instances();
  auto it = instances.find(name);
  if (it != instances.end()) return it->second.get();

  auto fit = package_factories().find(name);
  if (fit == package_factories().end()) return nullptr;

  std::unique_ptr<PhysicsPackage> pkg = fit->second();
  if (!pkg) return nullptr;
  const std::string canonical = pkg->name() ? std::string(pkg->name()) : std::string();
  if (canonical.empty() || canonical != name) return nullptr;
  PhysicsPackage* out = pkg.get();
  instances.emplace(name, std::move(pkg));
  return out;
}

bool has_physics_package(const std::string& name) {
  return package_factories().find(name) != package_factories().end();
}

PackageModeInfo resolve_package_mode_token(const std::string& package_name,
                                           const std::string& mode_token) {
  PhysicsPackage* pkg = get_physics_package(package_name);
  if (!pkg) return {};
  return pkg->resolve_mode_token(mode_token);
}

bool register_physics_package_factory(const std::string& name, PhysicsPackageFactory factory) {
  if (name.empty() || !factory) return false;
  auto& factories = package_factories();
  if (factories.find(name) != factories.end()) {
    return false;
  }
  factories.emplace(name, std::move(factory));
  return true;
}

}  // namespace galaxy
