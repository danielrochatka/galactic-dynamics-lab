#include "render_audit.hpp"
#include "galaxy_init.hpp"
#include "git_provenance.hpp"
#include "physics/TPFCore/derived_tpf_radial.hpp"
#include <cmath>
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

namespace {

bool tpf_vdsg_active_for_audit(const Config& config) {
  return std::isfinite(config.tpf_vdsg_coupling) && (config.tpf_vdsg_coupling != 0.0);
}

bool tpf_v11_weak_field_truncation_active(const Config& config) {
  return config.tpf_dynamics_mode == "v11_weak_field_truncation" ||
         config.tpf_dynamics_mode == "weak_field_correspondence";
}

}  // namespace

std::string compute_active_dynamics_branch(const Config& config) {
  if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
    if (config.v11_weak_field_correspondence_benchmark == "earth_moon_line_of_centers") {
      return "TPF_v11_weak_field_correspondence_audit_earth_moon_line_benchmark (correspondence-only; not "
             "legacy_readout; not direct_tpf dynamics; not orbit integration)";
    }
    return "TPF_v11_weak_field_correspondence_audit (correspondence-only; not legacy_readout; not direct_tpf "
           "dynamics)";
  }
  if (config.physics_package == "Newtonian") return "Newtonian_pairwise_G_SI";
  if (config.physics_package != "TPFCore") return config.physics_package + " (non-TPFCore)";
  if (tpf_v11_weak_field_truncation_active(config)) {
    return "TPF_v11_weak_field_truncation_weak_field_correspondence_helper_alpha_si_path";
  }
  if (config.tpf_dynamics_mode == "direct_tpf") {
    return std::string("TPF_direct_tpf_tensor_principal_part_DeltaC_omitted_") +
           (tpf_vdsg_active_for_audit(config) ? "VDSG_on" : "VDSG_off") +
           "_provisional_readout_off_shunt_off_cooling_off";
  }
  /* legacy_readout (default when key omitted) */
  if (!config.tpfcore_enable_provisional_readout)
    return "TPFCore_PROVISIONAL_legacy_readout_DISABLED (provisional_readout off)";
  const std::string& mode = config.tpfcore_readout_mode;
  if (tpf_vdsg_active_for_audit(config)) return "TPF_PROVISIONAL_legacy_readout_plus_EXPLORATORY_VDSG:" + mode;
  return "TPF_PROVISIONAL_legacy_readout:" + mode;
}

std::string compute_active_metrics_branch(const Config& config) {
  if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
    if (config.v11_weak_field_correspondence_benchmark == "earth_moon_line_of_centers") {
      return "v11_earth_moon_line_benchmark (phi Eq.44-45 vs Newtonian Eq.46 CSV; correspondence; DeltaC omitted)";
    }
    return "v11_paper_tensors (Xi,Theta,I,C_principal per Eq.(10) minus DeltaC; DeltaC omitted per v11 scope)";
  }
  if (config.physics_package == "Newtonian") return "none";
  if (config.physics_package == "TPFCore") {
    if (tpf_v11_weak_field_truncation_active(config))
      return "v11_weak_field_truncation_metrics_weak_field_correspondence_helper_alpha_si_path";
    if (config.tpf_dynamics_mode == "direct_tpf") {
      return std::string("direct_tpf_metrics_tensor_principal_part_DeltaC_omitted_") +
             (tpf_vdsg_active_for_audit(config) ? "VDSG_on" : "VDSG_off") +
             "_provisional_readout_off_shunt_off_cooling_off";
    }
    if (config.tpfcore_enable_provisional_readout)
      return "tpfcore_readout:" + config.tpfcore_readout_mode;
    return "TPFCore_metrics_n/a (provisional_readout off)";
  }
  return "unknown";
}

std::string compute_acceleration_code_path(const Config& config) {
  if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
    return "none (v11_weak_field_correspondence audit-only; no particle acceleration from this path)";
  }
  if (config.physics_package == "Newtonian") return "NewtonianPackage::compute_accelerations";
  if (config.physics_package != "TPFCore") return "unknown_package";
  if (tpf_v11_weak_field_truncation_active(config)) {
    return "TPFCorePackage::compute_v11_weak_field_truncation_accelerations (v11 Eq.42-44 weak-field correspondence helper; "
           "alpha_si correspondence path; no VDSG/readout/shunt/cooling)";
  }
  if (config.tpf_dynamics_mode == "direct_tpf") {
    return std::string("TPFCorePackage::compute_direct_tpf_accelerations "
                       "(tensor principal-part route: field_evaluation -> Theta3D -> principal_Cij -> tensor_projection; "
                       "Theta/I/kappa baseline; DeltaC omitted in current implementation scope; readout/shunt/cooling rejected)")
           + (tpf_vdsg_active_for_audit(config)
                  ? std::string(" + accumulate_vdsg_velocity_modifier (optional additive VDSG extension)")
                  : std::string(" + accumulate_vdsg_velocity_modifier (continuous zero contribution at tpf_vdsg_coupling == 0)"));
  }
  if (!config.tpfcore_enable_provisional_readout)
    return "TPFCorePackage::compute_accelerations (legacy_readout; throws without provisional readout)";
  std::string base;
  if (tpfcore::is_derived_tpf_radial_readout_mode(config.tpfcore_readout_mode))
    base = "TPFCorePackage::compute_provisional_readout_acceleration + derived_tpf_radial_profile";
  else
    base = "TPFCorePackage::compute_provisional_readout_acceleration (" + config.tpfcore_readout_mode + ")";
  std::string tail = " + accumulate_vdsg_velocity_modifier";
  if (config.tpf_global_accel_shunt_enable)
    tail += " + apply_global_accel_magnitude_shunt (when tpf_global_accel_shunt_enable)";
  else
    tail += " (global |a| shunt OFF — clean readout+VDSG path without velocity cap)";
  return base + tail;
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
  const GitProvenance gp = resolve_git_provenance();
  const bool cooling_on =
      (config.physics_package == "TPFCore" && config.tpf_cooling_fraction > 0.0 &&
       config.simulation_mode == SimulationMode::galaxy &&
       !tpf_v11_weak_field_truncation_active(config));

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
    json_kv(jf, first, "physics_package_compare", config.physics_package_compare);
    json_kv(jf, first, "simulation_mode", mode_to_string(config.simulation_mode));
    json_kv(jf, first, "artifact_naming_convention", "<mode>__<physics>__<scope>__<quantity>__<stage>.<ext>");
    json_kv(jf, first, "artifact_scope_note_primary", "primary = main interpretation outputs for selected mode");
    json_kv(jf, first, "artifact_scope_note_secondary",
            "secondary = lab/origin-radial diagnostics when pair-frame diagnostics are primary");
    json_kv(jf, first, "tpf_analysis_mode", config.tpf_analysis_mode);
    if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
      json_kv_bool(jf, first, "v11_delta_c_computed", false);
      json_kv(jf, first, "v11_weak_field_correspondence_benchmark", config.v11_weak_field_correspondence_benchmark);
      if (config.v11_weak_field_correspondence_benchmark == "earth_moon_line_of_centers") {
        json_kv(jf, first, "v11_weak_field_audit_scope",
                "earth_moon_line_of_centers_phi_Eq44_acceleration_Eq45_vs_Newtonian_Eq46; not_full_many_body_tpf; "
                "DeltaC_omitted");
        json_kv(jf, first, "v11_earth_moon_benchmark_png_note",
                "auto-generated when plot_v11_earth_moon_line_benchmark.py succeeds; requires python3 numpy matplotlib");
        json_kv(jf, first, "v11_earth_moon_benchmark_png_files",
                "tpf_v11_earth_moon_line_correspondence_benchmark_compare.png;"
                "tpf_v11_earth_moon_line_correspondence_benchmark_difference.png;"
                "tpf_v11_earth_moon_line_correspondence_benchmark_earth_zoom.png;"
                "tpf_v11_earth_moon_line_correspondence_benchmark_normalized_shape.png");
      } else {
        json_kv(jf, first, "v11_weak_field_audit_scope",
                "static_axis_benchmark_correspondence_only; DeltaC_omitted_per_manuscript_v11");
        json_kv(jf, first, "v11_eq10_C_principal_scaling",
                "C_mu_nu = kappa * (principal_bracket_including_minus_half_gI); C00_SI = kappa * I/2 when Theta0mu=0");
      }
    }
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
    if (config.physics_package == "TPFCore" &&
        config.simulation_mode != SimulationMode::tpf_v11_weak_field_correspondence) {
      const bool is_direct = (config.tpf_dynamics_mode == "direct_tpf");
      const bool is_v11_alias = tpf_v11_weak_field_truncation_active(config);
      json_kv(jf, first, "tpf_core_law_mode",
              is_direct ? "direct_tpf_canonical_entry"
                        : (is_v11_alias ? "v11_weak_field_truncation_compat_alias" : "legacy_readout_provisional"));
      json_kv(jf, first, "tpf_truncation_mode",
              is_direct
                  ? "direct_tpf_tensor_principal_part_Theta_I_kappa_baseline_DeltaC_omitted"
                  : (is_v11_alias
                         ? "v11_weak_field_correspondence_helper_alpha_si_path"
                         : "none"));
      json_kv(jf, first, "tpf_extension_mode",
              is_v11_alias ? "none_vdsg_rejected"
                           : (tpf_vdsg_active_for_audit(config) ? "vdsg" : "none"));
      json_kv(jf, first, "tpf_stabilizer_mode",
              (is_direct || is_v11_alias)
                  ? "none_shunt_and_cooling_rejected"
                  : ((config.tpf_global_accel_shunt_enable || config.tpf_cooling_fraction > 0.0)
                         ? "shunt_or_cooling_configured"
                         : "none"));
    }
    json_kv_num(jf, first, "tpf_vdsg_coupling", config.tpf_vdsg_coupling);
    json_kv_num(jf, first, "tpfcore_closure_kappa", config.tpf_kappa);
    json_kv_num(jf, first, "tpf_kappa", config.tpf_kappa);
    json_kv_num(jf, first, "tpf_cooling_fraction", config.tpf_cooling_fraction);
    json_kv_bool(jf, first, "tpf_cooling_active_this_run", cooling_on);
    json_kv_bool(jf, first, "tpfcore_enable_provisional_readout", config.tpfcore_enable_provisional_readout);
    if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
      json_kv_bool(jf, first, "v11_weak_field_correspondence_audit_only", true);
      json_kv_bool(jf, first, "tpfcore_enable_provisional_readout_operative_for_this_run", false);
      json_kv(jf, first, "tpf_dynamics_mode_operative_for_this_run", "none_audit_only");
    }
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
    json_kv_bool(jf, first, "tpf_global_accel_shunt_enable", config.tpf_global_accel_shunt_enable);
    json_kv_num(jf, first, "tpf_global_accel_shunt_fraction", config.tpf_global_accel_shunt_fraction);
    json_kv_bool(jf, first, "tpf_accel_pipeline_diagnostics_csv", config.tpf_accel_pipeline_diagnostics_csv);
    json_kv(jf, first, "git_commit_full", gp.git_commit_full);
    json_kv(jf, first, "git_commit_short", gp.git_commit_short);
    json_kv(jf, first, "git_branch", gp.git_branch);
    json_kv(jf, first, "git_tag", gp.git_tag);
    json_kv_bool(jf, first, "git_dirty", gp.git_dirty);
    json_kv(jf, first, "code_version_label", gp.code_version_label);
    json_kv(jf, first, "config_key_aliases_note",
            "Legacy key tpf_gdd_coupling maps to tpf_vdsg_coupling; parser accepts both (canonical: tpf_vdsg_coupling). "
            "Legacy key tpf_kappa maps to tpfcore_closure_kappa; parser accepts both (canonical: tpfcore_closure_kappa).");
    jf << ",\n  \"config_key_aliases\": [\n";
    jf << "    {\"legacy\": \"tpf_gdd_coupling\", \"canonical\": \"tpf_vdsg_coupling\", "
          "\"accepted_in_config_parser\": true},\n";
    jf << "    {\"legacy\": \"tpf_kappa\", \"canonical\": \"tpfcore_closure_kappa\", "
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
    if (config.physics_package == "TPFCore" &&
        config.simulation_mode != SimulationMode::tpf_v11_weak_field_correspondence) {
      const bool is_direct = (config.tpf_dynamics_mode == "direct_tpf");
      const bool is_v11_alias = tpf_v11_weak_field_truncation_active(config);
      tf << "tpf_core_law_mode\t"
         << (is_direct ? "direct_tpf_canonical_entry"
                       : (is_v11_alias ? "v11_weak_field_truncation_compat_alias" : "legacy_readout_provisional"))
         << "\n";
      tf << "tpf_truncation_mode\t"
         << (is_direct
                 ? "direct_tpf_tensor_principal_part_Theta_I_kappa_baseline_DeltaC_omitted"
                 : (is_v11_alias
                        ? "v11_weak_field_correspondence_helper_alpha_si_path"
                        : "none"))
         << "\n";
      tf << "tpf_extension_mode\t"
         << (is_v11_alias ? "none_vdsg_rejected"
                          : (tpf_vdsg_active_for_audit(config) ? "vdsg" : "none"))
         << "\n";
      tf << "tpf_stabilizer_mode\t"
         << ((is_direct || is_v11_alias)
                 ? "none_shunt_and_cooling_rejected"
                 : ((config.tpf_global_accel_shunt_enable || config.tpf_cooling_fraction > 0.0)
                        ? "shunt_or_cooling_configured"
                        : "none"))
         << "\n";
    }
    tf << "physics_package\t" << config.physics_package << "\n";
    tf << "physics_package_compare\t" << config.physics_package_compare << "\n";
    tf << "simulation_mode\t" << mode_to_string(config.simulation_mode) << "\n";
    tf << "artifact_naming_convention\t<mode>__<physics>__<scope>__<quantity>__<stage>.<ext>\n";
    tf << "artifact_scope_note_primary\tprimary=main interpretation outputs for selected mode\n";
    tf << "artifact_scope_note_secondary\tsecondary=lab/origin-radial diagnostics when pair frame is primary\n";
    tf << "tpf_analysis_mode\t" << config.tpf_analysis_mode << "\n";
    if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
      tf << "v11_delta_c_computed\t0\n";
      tf << "v11_weak_field_correspondence_benchmark\t" << config.v11_weak_field_correspondence_benchmark << "\n";
      if (config.v11_weak_field_correspondence_benchmark == "earth_moon_line_of_centers") {
        tf << "v11_weak_field_audit_scope\tearth_moon_line_phi_Eq44_Eq45_vs_Newtonian_Eq46_not_full_many_body_DeltaC_omitted\n";
        tf << "v11_earth_moon_benchmark_png_note\tauto_generated_when_plot_v11_earth_moon_line_benchmark_py_succeeds\n";
        tf << "v11_earth_moon_benchmark_png_files\ttpf_v11_earth_moon_line_correspondence_benchmark_compare.png;"
              "tpf_v11_earth_moon_line_correspondence_benchmark_difference.png;"
              "tpf_v11_earth_moon_line_correspondence_benchmark_earth_zoom.png;"
              "tpf_v11_earth_moon_line_correspondence_benchmark_normalized_shape.png\n";
      } else {
        tf << "v11_weak_field_audit_scope\tstatic_axis_benchmark_correspondence_only_DeltaC_omitted_per_v11\n";
        tf << "v11_eq10_C_principal_scaling\tC_mu_nu=kappa*(principal_bracket_including_minus_half_gI);C00_SI=kappa*I/2_when_Theta0mu=0\n";
      }
    }
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
    tf << "tpfcore_closure_kappa\t" << config.tpf_kappa << "\n";
    tf << "tpf_kappa\t" << config.tpf_kappa << "\n";
    tf << "tpf_cooling_fraction\t" << config.tpf_cooling_fraction << "\n";
    tf << "tpf_cooling_active_this_run\t" << (cooling_on ? 1 : 0) << "\n";
    if (config.simulation_mode == SimulationMode::tpf_v11_weak_field_correspondence) {
      tf << "v11_weak_field_correspondence_audit_only\t1\n";
      tf << "v11_audit_tpfcore_dynamics_note\tno particle integration; configured legacy_readout/direct_tpf fields not operative\n";
      tf << "tpf_dynamics_mode_configured\t" << config.tpf_dynamics_mode << "\n";
      tf << "tpf_dynamics_mode_operative_for_this_run\tnone_audit_only\n";
      tf << "tpfcore_enable_provisional_readout_configured\t" << (config.tpfcore_enable_provisional_readout ? 1 : 0)
         << "\n";
      tf << "tpfcore_enable_provisional_readout_operative_for_this_run\t0\n";
    } else {
      tf << "tpfcore_enable_provisional_readout\t" << (config.tpfcore_enable_provisional_readout ? 1 : 0) << "\n";
      tf << "tpf_dynamics_mode\t" << config.tpf_dynamics_mode << "\n";
    }
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
    tf << "tpf_global_accel_shunt_enable\t" << (config.tpf_global_accel_shunt_enable ? 1 : 0) << "\n";
    tf << "tpf_global_accel_shunt_fraction\t" << config.tpf_global_accel_shunt_fraction << "\n";
    tf << "tpf_accel_pipeline_diagnostics_csv\t" << (config.tpf_accel_pipeline_diagnostics_csv ? 1 : 0) << "\n";
    tf << "git_commit_full\t" << gp.git_commit_full << "\n";
    tf << "git_commit_short\t" << gp.git_commit_short << "\n";
    tf << "git_branch\t" << gp.git_branch << "\n";
    tf << "git_tag\t" << gp.git_tag << "\n";
    tf << "git_dirty\t" << (gp.git_dirty ? 1 : 0) << "\n";
    tf << "code_version_label\t" << gp.code_version_label << "\n";
    if (galaxy_init_audit && galaxy_init_audit->valid) {
      tf << "galaxy_init_audit_template\t" << galaxy_init_audit->template_name << "\n";
      tf << "galaxy_init_audit_used_legacy_velocity_noise\t"
         << (galaxy_init_audit->used_legacy_velocity_noise ? 1 : 0) << "\n";
    }
    tf << "config_key_aliases\ttpf_gdd_coupling -> tpf_vdsg_coupling (legacy alias; canonical tpf_vdsg_coupling)\n";
    tf << "config_key_aliases\ttpf_kappa -> tpfcore_closure_kappa (legacy alias; canonical tpfcore_closure_kappa)\n";
    tf << "=== end render_manifest.txt ===\n";
  }
}

}  // namespace galaxy
