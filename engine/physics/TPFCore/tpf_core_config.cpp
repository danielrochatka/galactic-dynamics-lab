#include "tpf_core_config.hpp"

#include "../../config.hpp"
#include "../../package_discovery.hpp"

#include <stdexcept>

namespace galaxy {

namespace {

bool get_opt(const Config& c, const std::string& k, std::string& out) {
  auto it = c.package_options.find(k);
  if (it == c.package_options.end()) return false;
  out = it->second;
  return true;
}

}

void sync_tpfcore_legacy_fields_from_package_options(Config& config) {
  std::string v;
  if (get_opt(config, "tpf_dynamics_mode", v)) {
    if (v != "tpf_xi_theta_v1") {
      throw std::runtime_error("tpf_dynamics_mode must be tpf_xi_theta_v1 on this branch, got: " + v);
    }
    config.tpf_dynamics_mode = v;
  }
  if (get_opt(config, "tpf_analysis_mode", v)) {
    if (v != "none" && v != "v11_weak_field_correspondence") {
      throw std::runtime_error("tpf_analysis_mode must be none or v11_weak_field_correspondence, got: " + v);
    }
    config.tpf_analysis_mode = v;
  }
  if (get_opt(config, "tpfcore_readout_mode", v)) config.tpfcore_readout_mode = v;
  if (get_opt(config, "tpfcore_readout_scale", v)) config.tpfcore_readout_scale = std::stod(v);
  if (get_opt(config, "tpfcore_theta_tt_scale", v)) config.tpfcore_theta_tt_scale = std::stod(v);
  if (get_opt(config, "tpfcore_theta_tr_scale", v)) config.tpfcore_theta_tr_scale = std::stod(v);
  if (get_opt(config, "tpf_kappa", v)) config.tpf_kappa = std::stod(v);
  if (get_opt(config, "tpfcore_closure_kappa", v)) config.tpf_kappa = std::stod(v);
  if (get_opt(config, "tpf_vdsg_coupling", v)) config.tpf_vdsg_coupling = std::stod(v);
  if (get_opt(config, "tpf_global_accel_shunt_enable", v)) config.tpf_global_accel_shunt_enable = (v == "1" || v == "true" || v == "yes");
  if (get_opt(config, "tpf_global_accel_shunt_fraction", v)) config.tpf_global_accel_shunt_fraction = std::stod(v);
  if (get_opt(config, "tpf_cooling_fraction", v)) config.tpf_cooling_fraction = std::stod(v);
  if (get_opt(config, "tpf_accel_pipeline_diagnostics_csv", v)) config.tpf_accel_pipeline_diagnostics_csv = (v == "1" || v == "true" || v == "yes");
  if (get_opt(config, "tpf_poisson_bins", v)) config.tpf_poisson_bins = std::stoi(v);
  if (get_opt(config, "tpf_poisson_max_radius", v)) config.tpf_poisson_max_radius = std::stod(v);
  if (get_opt(config, "tpfcore_dump_readout_debug", v)) config.tpfcore_dump_readout_debug = (v == "1" || v == "true" || v == "yes");
  if (get_opt(config, "tpfcore_live_orbit_force_audit", v)) config.tpfcore_live_orbit_force_audit = (v == "1" || v == "true" || v == "yes");
}

bool tpfcore_owns_mode_token(const std::string& mode_token) {
  const PackageMetadata* md = find_package_metadata("TPFCore");
  return md && md->modes.count(mode_token) > 0;
}

} // namespace galaxy
