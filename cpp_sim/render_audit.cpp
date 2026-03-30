#include "render_audit.hpp"
#include "galaxy_init.hpp"
#include "physics/TPFCore/derived_tpf_radial.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace galaxy {

namespace {

void json_esc(std::ostream& o, const std::string& s) {
  o << '"';
  for (char c : s) {
    if (c == '\\' || c == '"')
      o << '\\' << c;
    else if (c == '\n')
      o << "\\n";
    else if (c == '\r')
      o << "\\r";
    else if (c == '\t')
      o << "\\t";
    else
      o << c;
  }
  o << '"';
}

void json_kv(std::ostream& j, bool& first, const char* key, const std::string& val) {
  if (!first) j << ',';
  first = false;
  j << '\n' << "  \"" << key << "\": ";
  json_esc(j, val);
}

void json_kv_num(std::ostream& j, bool& first, const char* key, double v) {
  if (!first) j << ',';
  first = false;
  j << '\n' << "  \"" << key << "\": " << std::scientific << std::setprecision(17) << v;
}

void json_kv_int(std::ostream& j, bool& first, const char* key, long long v) {
  if (!first) j << ',';
  first = false;
  j << '\n' << "  \"" << key << "\": " << v;
}

void json_kv_bool(std::ostream& j, bool& first, const char* key, bool v) {
  if (!first) j << ',';
  first = false;
  j << '\n' << "  \"" << key << "\": " << (v ? "true" : "false");
}

std::string run_id_from_paths(const std::string& output_dir, const Config& config) {
  if (!config.run_id.empty()) return config.run_id;
  /* basename of output_dir */
  std::string p = output_dir;
  while (!p.empty() && (p.back() == '/' || p.back() == '\\')) p.pop_back();
  auto slash = p.find_last_of("/\\");
  if (slash == std::string::npos) return p;
  return p.substr(slash + 1);
}

}  // namespace

std::string compute_active_dynamics_branch(const Config& config) {
  if (config.physics_package == "Newtonian") return "Newtonian_pairwise_G_SI";
  if (config.physics_package != "TPFCore") return config.physics_package + " (non-TPFCore)";
  if (!config.tpfcore_enable_provisional_readout)
    return "TPFCore_dynamics_DISABLED (provisional_readout off)";
  /* Baseline readout + optional additive VDSG modifier; branch label is readout identity (not VDSG-only). */
  return "TPF_readout_acceleration:" + config.tpfcore_readout_mode;
}

std::string compute_active_metrics_branch(const Config& config) {
  if (config.physics_package == "Newtonian") return "none";
  if (config.physics_package == "TPFCore") {
    if (config.tpfcore_enable_provisional_readout)
      return "tpfcore_readout:" + config.tpfcore_readout_mode;
    return "TPFCore_metrics_n/a (provisional_readout off)";
  }
  return "unknown";
}

std::string compute_acceleration_code_path(const Config& config) {
  if (config.physics_package == "Newtonian") return "NewtonianPackage::compute_accelerations";
  if (config.physics_package != "TPFCore") return "unknown_package";
  if (!config.tpfcore_enable_provisional_readout)
    return "TPFCorePackage::compute_accelerations (throws without provisional readout)";
  std::string base;
  if (tpfcore::is_derived_tpf_radial_readout_mode(config.tpfcore_readout_mode))
    base = "TPFCorePackage::compute_provisional_readout_acceleration + derived_tpf_radial_profile";
  else
    base = "TPFCorePackage::compute_provisional_readout_acceleration (" + config.tpfcore_readout_mode + ")";
  if (config.tpf_vdsg_coupling != 0.0) return base + " + accumulate_vdsg_velocity_modifier";
  return base;
}

void write_render_manifest(const std::string& output_dir,
                           const Config& config,
                           int n_steps_done,
                           int n_snapshots,
                           int n_particles,
                           const GalaxyInitAudit* galaxy_init_audit) {
  const std::string run_id = run_id_from_paths(output_dir, config);
  const int n_star = (n_particles >= 0) ? n_particles : config.n_stars;
  const std::string dyn = compute_active_dynamics_branch(config);
  const std::string met = compute_active_metrics_branch(config);
  const std::string acc = compute_acceleration_code_path(config);
  const bool cooling_on =
      (config.physics_package == "TPFCore" && config.tpf_cooling_fraction > 0.0 &&
       config.simulation_mode == SimulationMode::galaxy);

  std::ostringstream json_path;
  json_path << output_dir << "/render_manifest.json";
  std::ofstream jf(json_path.str());
  if (jf) {
    bool first = true;
    jf << '{';
    json_kv(jf, first, "schema", "galaxy_render_manifest_v1");
    json_kv(jf, first, "run_id", run_id);
    json_kv(jf, first, "output_dir", output_dir);
    json_kv(jf, first, "physics_package", config.physics_package);
    json_kv(jf, first, "simulation_mode", mode_to_string(config.simulation_mode));
    json_kv_int(jf, first, "n_stars", n_star);
    json_kv_int(jf, first, "n_steps", n_steps_done);
    json_kv_int(jf, first, "n_snapshots", n_snapshots);
    json_kv_num(jf, first, "dt", config.dt);
    json_kv_num(jf, first, "softening", config.softening);
    json_kv_num(jf, first, "bh_mass", config.bh_mass);
    json_kv_num(jf, first, "star_mass", config.star_mass);
    json_kv_num(jf, first, "galaxy_radius", config.galaxy_radius);
    json_kv(jf, first, "active_dynamics_branch", dyn);
    json_kv(jf, first, "active_metrics_branch", met);
    json_kv(jf, first, "acceleration_code_path", acc);
    json_kv_num(jf, first, "tpf_vdsg_coupling", config.tpf_vdsg_coupling);
    json_kv_num(jf, first, "tpf_kappa", config.tpf_kappa);
    json_kv_num(jf, first, "tpf_cooling_fraction", config.tpf_cooling_fraction);
    json_kv_bool(jf, first, "tpf_cooling_active_this_run", cooling_on);
    json_kv_bool(jf, first, "tpfcore_enable_provisional_readout", config.tpfcore_enable_provisional_readout);
    json_kv(jf, first, "tpfcore_readout_mode", config.tpfcore_readout_mode);
    json_kv(jf, first, "render_overlay_mode", config.render_overlay_mode);
    json_kv(jf, first, "galaxy_init_template", config.galaxy_init_template);
    json_kv_int(jf, first, "galaxy_init_seed", static_cast<long long>(config.galaxy_init_seed));
    json_kv_num(jf, first, "galaxy_init_master_chaos", config.galaxy_init_master_chaos);
    json_kv_num(jf, first, "galaxy_init_position_noise", config.galaxy_init_position_noise);
    json_kv_num(jf, first, "galaxy_init_velocity_angle_noise", config.galaxy_init_velocity_angle_noise);
    json_kv_num(jf, first, "galaxy_init_velocity_magnitude_noise", config.galaxy_init_velocity_magnitude_noise);
    json_kv_num(jf, first, "galaxy_init_clumpiness", config.galaxy_init_clumpiness);
    json_kv_int(jf, first, "galaxy_init_num_clumps", config.galaxy_init_num_clumps);
    json_kv_num(jf, first, "galaxy_init_clump_radius_fraction", config.galaxy_init_clump_radius_fraction);
    json_kv_num(jf, first, "galaxy_init_m2_amplitude", config.galaxy_init_m2_amplitude);
    json_kv_num(jf, first, "galaxy_init_m3_amplitude", config.galaxy_init_m3_amplitude);
    json_kv_num(jf, first, "galaxy_init_bar_amplitude", config.galaxy_init_bar_amplitude);
    json_kv_num(jf, first, "galaxy_init_bar_axis_ratio", config.galaxy_init_bar_axis_ratio);
    json_kv_num(jf, first, "galaxy_init_spiral_amplitude", config.galaxy_init_spiral_amplitude);
    json_kv_num(jf, first, "galaxy_init_spiral_winding", config.galaxy_init_spiral_winding);
    json_kv_num(jf, first, "galaxy_init_spiral_phase", config.galaxy_init_spiral_phase);
    json_kv_bool(jf, first, "enable_star_star_gravity", config.enable_star_star_gravity);
    json_kv_num(jf, first, "tpf_vdsg_mass_baseline_kg", config.tpf_vdsg_mass_baseline_kg);
    json_kv(jf, first, "config_key_aliases_note",
            "Legacy key tpf_gdd_coupling maps to tpf_vdsg_coupling; parser accepts both (canonical: tpf_vdsg_coupling).");
    jf << ",\n  \"config_key_aliases\": [\n";
    jf << "    {\"legacy\": \"tpf_gdd_coupling\", \"canonical\": \"tpf_vdsg_coupling\", "
          "\"accepted_in_config_parser\": true}\n  ]";
    jf << "\n}\n";
  }

  std::ostringstream txt_path;
  txt_path << output_dir << "/render_manifest.txt";
  std::ofstream tf(txt_path.str());
  if (tf) {
    tf << "=== render_manifest.txt (galaxy render / audit; mirrors render_manifest.json) ===\n";
    tf << "run_id\t" << run_id << "\n";
    tf << "output_dir\t" << output_dir << "\n";
    tf << "active_dynamics_branch\t" << dyn << "\n";
    tf << "active_metrics_branch\t" << met << "\n";
    tf << "acceleration_code_path\t" << acc << "\n";
    tf << "physics_package\t" << config.physics_package << "\n";
    tf << "simulation_mode\t" << mode_to_string(config.simulation_mode) << "\n";
    tf << "n_stars\t" << n_star << "\n";
    tf << "n_steps\t" << n_steps_done << "\n";
    tf << "n_snapshots\t" << n_snapshots << "\n";
    tf << "dt\t" << std::scientific << config.dt << "\n";
    tf << "softening\t" << config.softening << "\n";
    tf << "bh_mass\t" << config.bh_mass << "\n";
    tf << "star_mass\t" << config.star_mass << "\n";
    tf << "galaxy_radius\t" << config.galaxy_radius << "\n";
    tf << "enable_star_star_gravity\t" << (config.enable_star_star_gravity ? 1 : 0) << "\n";
    tf << "tpf_vdsg_coupling\t" << config.tpf_vdsg_coupling << "\n";
    tf << "tpf_kappa\t" << config.tpf_kappa << "\n";
    tf << "tpf_cooling_fraction\t" << config.tpf_cooling_fraction << "\n";
    tf << "tpf_cooling_active_this_run\t" << (cooling_on ? 1 : 0) << "\n";
    tf << "tpfcore_enable_provisional_readout\t" << (config.tpfcore_enable_provisional_readout ? 1 : 0) << "\n";
    tf << "tpfcore_readout_mode\t" << config.tpfcore_readout_mode << "\n";
    tf << "render_overlay_mode\t" << config.render_overlay_mode << "\n";
    tf << "galaxy_init_template\t" << config.galaxy_init_template << "\n";
    tf << "galaxy_init_seed\t" << config.galaxy_init_seed << "\n";
    tf << "galaxy_init_master_chaos\t" << config.galaxy_init_master_chaos << "\n";
    tf << "galaxy_init_position_noise\t" << config.galaxy_init_position_noise << "\n";
    tf << "galaxy_init_velocity_angle_noise\t" << config.galaxy_init_velocity_angle_noise << "\n";
    tf << "galaxy_init_velocity_magnitude_noise\t" << config.galaxy_init_velocity_magnitude_noise << "\n";
    tf << "galaxy_init_clumpiness\t" << config.galaxy_init_clumpiness << "\n";
    tf << "galaxy_init_num_clumps\t" << config.galaxy_init_num_clumps << "\n";
    tf << "galaxy_init_clump_radius_fraction\t" << config.galaxy_init_clump_radius_fraction << "\n";
    tf << "galaxy_init_m2_amplitude\t" << config.galaxy_init_m2_amplitude << "\n";
    tf << "galaxy_init_m3_amplitude\t" << config.galaxy_init_m3_amplitude << "\n";
    tf << "galaxy_init_bar_amplitude\t" << config.galaxy_init_bar_amplitude << "\n";
    tf << "galaxy_init_bar_axis_ratio\t" << config.galaxy_init_bar_axis_ratio << "\n";
    tf << "galaxy_init_spiral_amplitude\t" << config.galaxy_init_spiral_amplitude << "\n";
    tf << "galaxy_init_spiral_winding\t" << config.galaxy_init_spiral_winding << "\n";
    tf << "galaxy_init_spiral_phase\t" << config.galaxy_init_spiral_phase << "\n";
    if (galaxy_init_audit && galaxy_init_audit->valid) {
      tf << "galaxy_init_audit_template\t" << galaxy_init_audit->template_name << "\n";
      tf << "galaxy_init_audit_used_legacy_velocity_noise\t"
         << (galaxy_init_audit->used_legacy_velocity_noise ? 1 : 0) << "\n";
    }
    tf << "config_key_aliases\ttpf_gdd_coupling -> tpf_vdsg_coupling (legacy alias; canonical tpf_vdsg_coupling)\n";
    tf << "=== end render_manifest.txt ===\n";
  }
}

}  // namespace galaxy
