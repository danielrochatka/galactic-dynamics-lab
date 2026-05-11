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

bool parse_bool_opt(const std::string& v) {
  return v == "1" || v == "true" || v == "yes";
}

}

TPFCoreConfig build_tpfcore_config(const Config& config) {
  TPFCoreConfig out;
  std::string v;
  out.tpf_dynamics_mode = config.tpf_dynamics_mode;
  out.tpfcore_enable_provisional_readout = config.tpfcore_enable_provisional_readout;
  out.tpfcore_readout_mode = config.tpfcore_readout_mode;
  out.tpfcore_readout_scale = config.tpfcore_readout_scale;
  out.tpfcore_theta_tt_scale = config.tpfcore_theta_tt_scale;
  out.tpfcore_theta_tr_scale = config.tpfcore_theta_tr_scale;
  out.tpf_kappa = config.tpf_kappa;
  out.tpf_weak_field_correspondence_alpha_si = config.tpf_weak_field_correspondence_alpha_si;
  out.tpf_vdsg_coupling = config.tpf_vdsg_coupling;
  out.tpf_vdsg_mass_baseline_kg = config.tpf_vdsg_mass_baseline_kg;
  out.tpf_vdsg_mode = config.tpf_vdsg_mode;
  out.tpf_cooling_fraction = config.tpf_cooling_fraction;
  out.tpf_global_accel_shunt_enable = config.tpf_global_accel_shunt_enable;
  out.tpf_global_accel_shunt_fraction = config.tpf_global_accel_shunt_fraction;
  out.tpf_poisson_bins = config.tpf_poisson_bins;
  out.tpf_poisson_max_radius = config.tpf_poisson_max_radius;

  if (get_opt(config, "tpf_dynamics_mode", v)) out.tpf_dynamics_mode = v;
  if (get_opt(config, "tpfcore_enable_provisional_readout", v)) out.tpfcore_enable_provisional_readout = parse_bool_opt(v);
  if (get_opt(config, "tpfcore_readout_mode", v)) out.tpfcore_readout_mode = v;
  if (get_opt(config, "tpfcore_readout_scale", v)) out.tpfcore_readout_scale = std::stod(v);
  if (get_opt(config, "tpfcore_theta_tt_scale", v)) out.tpfcore_theta_tt_scale = std::stod(v);
  if (get_opt(config, "tpfcore_theta_tr_scale", v)) out.tpfcore_theta_tr_scale = std::stod(v);
  if (get_opt(config, "tpf_kappa", v)) out.tpf_kappa = std::stod(v);
  if (get_opt(config, "tpf_weak_field_correspondence_alpha_si", v)) out.tpf_weak_field_correspondence_alpha_si = std::stod(v);
  if (get_opt(config, "tpf_vdsg_coupling", v)) out.tpf_vdsg_coupling = std::stod(v);
  if (get_opt(config, "tpf_vdsg_mass_baseline_kg", v)) out.tpf_vdsg_mass_baseline_kg = std::stod(v);
  if (get_opt(config, "tpf_vdsg_mode", v)) out.tpf_vdsg_mode = v;
  if (get_opt(config, "tpf_cooling_fraction", v)) out.tpf_cooling_fraction = std::stod(v);
  if (get_opt(config, "tpf_global_accel_shunt_enable", v)) out.tpf_global_accel_shunt_enable = parse_bool_opt(v);
  if (get_opt(config, "tpf_global_accel_shunt_fraction", v)) out.tpf_global_accel_shunt_fraction = std::stod(v);
  if (get_opt(config, "tpf_poisson_bins", v)) out.tpf_poisson_bins = std::stoi(v);
  if (get_opt(config, "tpf_poisson_max_radius", v)) out.tpf_poisson_max_radius = std::stod(v);
  return out;
}

void sync_tpfcore_legacy_fields_from_package_options(Config& config, bool validate_dynamics_mode) {
  std::string v;
  if (get_opt(config, "tpf_dynamics_mode", v)) {
    if (validate_dynamics_mode && v != "tpf_xi_theta_v1") {
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
  if (get_opt(config, "tpfcore_enable_provisional_readout", v)) config.tpfcore_enable_provisional_readout = (v == "1" || v == "true" || v == "yes");
  if (get_opt(config, "tpf_weak_field_correspondence_alpha_si", v)) {
    config.tpf_weak_field_correspondence_alpha_si = std::stod(v);
    config.tpf_weak_field_correspondence_alpha_si_explicitly_set = true;
  }
  if (get_opt(config, "tpfcore_readout_mode", v)) config.tpfcore_readout_mode = v;
  if (get_opt(config, "tpfcore_readout_scale", v)) config.tpfcore_readout_scale = std::stod(v);
  if (get_opt(config, "tpfcore_theta_tt_scale", v)) config.tpfcore_theta_tt_scale = std::stod(v);
  if (get_opt(config, "tpfcore_theta_tr_scale", v)) config.tpfcore_theta_tr_scale = std::stod(v);
  if (get_opt(config, "tpf_kappa", v)) config.tpf_kappa = std::stod(v);
  if (get_opt(config, "tpfcore_closure_kappa", v)) config.tpf_kappa = std::stod(v);
  if (get_opt(config, "tpf_vdsg_coupling", v)) config.tpf_vdsg_coupling = std::stod(v);
  if (get_opt(config, "tpf_vdsg_mass_baseline_kg", v)) config.tpf_vdsg_mass_baseline_kg = std::stod(v);
  if (get_opt(config, "tpf_vdsg_mode", v)) config.tpf_vdsg_mode = v;
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
