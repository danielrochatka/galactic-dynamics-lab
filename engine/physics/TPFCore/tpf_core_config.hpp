#ifndef GALAXY_TPF_CORE_CONFIG_HPP
#define GALAXY_TPF_CORE_CONFIG_HPP

#include <string>

namespace galaxy {
struct Config;

struct TPFCoreConfig {
  std::string tpf_dynamics_mode = "tpf_xi_theta_v1";
  bool tpfcore_enable_provisional_readout = false;
  std::string tpfcore_readout_mode = "tensor_radial_projection";
  double tpfcore_readout_scale = 1.0;
  double tpfcore_theta_tt_scale = 1.0;
  double tpfcore_theta_tr_scale = 1.0;
  double tpf_kappa = 1.0e32;
  double tpf_weak_field_correspondence_alpha_si = -6.67430e-11;
  double tpf_vdsg_coupling = 0.0;
  double tpf_vdsg_mass_baseline_kg = 0.0;
  std::string tpf_vdsg_mode = "legacy_speed";
  double tpf_cooling_fraction = 0.0;
  bool tpf_global_accel_shunt_enable = false;
  double tpf_global_accel_shunt_fraction = 0.001;
  int tpf_poisson_bins = 100;
  double tpf_poisson_max_radius = 0.0;
};

TPFCoreConfig build_tpfcore_config(const Config& config);

void sync_tpfcore_legacy_fields_from_package_options(Config& config);

bool tpfcore_owns_mode_token(const std::string& mode_token);

} // namespace galaxy

#endif
