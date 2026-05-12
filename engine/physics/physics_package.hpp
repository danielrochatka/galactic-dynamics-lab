#ifndef GALAXY_PHYSICS_PACKAGE_HPP
#define GALAXY_PHYSICS_PACKAGE_HPP

/**
 * Physics package interface.
 * Every physics package (Newtonian, custom, etc.) must implement this interface.
 * The simulator and integrator depend only on this interface, not on any specific implementation.
 */

#include "../config.hpp"
#include "../types.hpp"
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace galaxy {

struct PackageMetadataEntry {
  enum class RenderPlacement {
    AfterAccelerationCodePath,
    AfterDynamicsMode
  };
  std::string key;
  std::string value;
  enum class ValueType { String, Number, Bool };
  ValueType value_type = ValueType::String;
  double number_value = 0.0;
  bool bool_value = false;
  RenderPlacement render_placement = RenderPlacement::AfterAccelerationCodePath;
};


struct PackageModeInfo {
  bool recognized = false;
  std::string canonical_token;
  bool requires_output_dir = true;
  bool is_utility_mode = false;
};

struct RunInfoContext {
  const void* package_runtime_stats = nullptr;
  const char* package_runtime_stats_kind = nullptr;
};

enum class PackageRunInfoSection {
  BranchDiagnosticsSupplement,
  PostCodeProvenance,
  LateRuntimeDetails
};

class PhysicsPackage {
 public:
  virtual ~PhysicsPackage() = default;

  /** Package identifier (e.g. "Newtonian"). Must match config physics_package name. */
  virtual const char* name() const = 0;

  /**
   * Compute accelerations for all particles. Required.
   * Fills ax[i], ay[i] for each particle i. Same signature used by the integrator.
   */
  virtual void compute_accelerations(const State& state,
                                     double bh_mass,
                                     double softening,
                                     bool star_star,
                                     std::vector<double>& ax,
                                     std::vector<double>& ay) const = 0;

  /**
   * Optional: total gravitational potential energy for diagnostics/validation.
   * star_star: include star-star contributions (must match compute_accelerations).
   * Default returns 0. Override in packages that define a potential.
   */
  virtual double compute_potential_energy(const State& state,
                                         double bh_mass,
                                         double softening,
                                         bool star_star = true) const {
    (void)state;
    (void)bh_mass;
    (void)softening;
    (void)star_star;
    return 0.0;
  }

  /** Optional: called once before simulation (e.g. load tables). Default no-op. */
  virtual void init() {}

  /** Optional: called with config before run. Packages may cache package-specific params. Default no-op. */
  virtual void init_from_config(const Config&) {}

  /** Optional: validation/reporting name or hook. Default returns name(). */
  virtual const char* validation_name() const { return name(); }

  /** Optional: human-facing package name. Default returns name(). */
  virtual std::string display_name() const { return name(); }

  /** Optional: runtime metadata entries to append to run-level reporting. Default empty. */
  virtual std::vector<PackageMetadataEntry> runtime_metadata() const { return {}; }

  virtual PackageModeInfo resolve_mode_token(const std::string& mode_token) const {
    (void)mode_token;
    return {};
  }

  /** Optional: synchronize package-owned options into any package-owned compatibility state. Default no-op. */
  virtual void sync_config_from_package_options(Config& config,
                                                const std::string& changed_key) const {
    (void)config;
    (void)changed_key;
  }

  /** Optional: package capability check for utility modes. Default unsupported. */
  virtual bool supports_utility_mode(SimulationMode) const { return false; }

  /** Optional: utility mode dispatch. Default no-op/false (not handled). */
  virtual bool run_utility_mode(const Config&, const std::string&) { return false; }
  virtual bool should_auto_plot_utility_mode(const Config& config, const PackageModeInfo& mode_info) const {
    (void)config;
    (void)mode_info;
    return false;
  }

  /** Optional: package run-info metadata. Default empty. */
  virtual std::vector<PackageMetadataEntry> run_info_metadata(const Config&, const std::string&) const { return {}; }

  /** Optional: package run-info supplement metadata emitted in branch/physics supplements section. */
  virtual std::vector<PackageMetadataEntry> run_info_supplement_metadata(const Config&) const { return {}; }

  /** Optional: package run-info metadata that depends on runtime context/stats. */
  virtual std::vector<PackageMetadataEntry> run_info_runtime_metadata(const Config&,
                                                                      const RunInfoContext&) const {
    return {};
  }
  virtual void write_run_info_section(std::ostream& out,
                                      const Config& config,
                                      const RunInfoContext& context,
                                      PackageRunInfoSection section) const {
    (void)out;
    (void)config;
    (void)context;
    (void)section;
  }

  /** Optional: package render/report metadata. Default empty. */
  virtual std::vector<PackageMetadataEntry> render_metadata(const Config&) const { return {}; }
  virtual std::string render_active_dynamics_branch(const Config&) const { return {}; }
  virtual std::string render_active_metrics_branch(const Config&) const { return {}; }
  virtual std::string render_acceleration_code_path(const Config&) const { return {}; }

  /** Optional: package config schema/default metadata. Default empty. */
  virtual std::vector<PackageMetadataEntry> package_config_metadata() const { return {}; }

  /** Optional: package-owned resolved/effective runtime metadata. Default empty. */
  virtual std::vector<PackageMetadataEntry> resolved_runtime_metadata(const Config& config,
                                                                      SimulationMode mode) const {
    (void)config;
    (void)mode;
    return {};
  }

  /** Optional: package-owned post-run diagnostics emission. Default no-op/success. */
  virtual bool write_post_run_diagnostics(const std::vector<Snapshot>& snapshots,
                                          const Config& config,
                                          const std::string& output_dir) {
    (void)snapshots;
    (void)config;
    (void)output_dir;
    return true;
  }

  /** Optional: package-owned cooling policy activation. Default inactive. */
  virtual bool cooling_active(const Config& config) const {
    (void)config;
    return false;
  }

  /** Optional: package-owned cooling policy step count. Default zero. */
  virtual int cooling_steps(const Config& config, int n_steps) const {
    (void)config;
    (void)n_steps;
    return 0;
  }

  /** Optional: package-owned per-step cooling application. Default no-op. */
  virtual void apply_cooling_step(State& state, const Config& config, int step, int cooling_steps) const {
    (void)state;
    (void)config;
    (void)step;
    (void)cooling_steps;
  }

  /** Optional: package-owned cooling snapshot suppression policy. Default false (keep snapshots). */
  virtual bool suppress_snapshot_for_cooling(const Config& config, int step, int cooling_steps) const {
    (void)config;
    (void)step;
    (void)cooling_steps;
    return false;
  }
};

/** Generic kinetic energy (same for all packages): 0.5 * sum(m_i * v_i^2). */
double compute_kinetic_energy(const State& state);

/**
 * Registry: get a physics package by name.
 * Returns nullptr if unknown. Caller must not delete the pointer (static/registry-owned).
 */
PhysicsPackage* get_physics_package(const std::string& name);

/** Returns true if the named package is available. */
bool has_physics_package(const std::string& name);

/** Factory function type used by the package registry. */
using PhysicsPackageFactory = std::function<std::unique_ptr<PhysicsPackage>()>;

/**
 * Register a package factory under a stable package name.
 * Returns false if name is empty, factory is null, or the name is already registered.
 */
bool register_physics_package_factory(const std::string& name, PhysicsPackageFactory factory);

PackageModeInfo resolve_package_mode_token(const std::string& package_name,
                                           const std::string& mode_token);

}  // namespace galaxy

#endif
