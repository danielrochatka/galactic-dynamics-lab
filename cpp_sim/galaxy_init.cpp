#include "galaxy_init.hpp"
#include "physics/TPFCore/derived_tpf_radial.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <numeric>
#include <random>
#include <stdexcept>

namespace galaxy {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;
/** Inner edge of disk annulus as fraction of galaxy_radius (matches legacy init_galaxy_disk). */
constexpr double kDiskInnerFraction = 0.05;
constexpr int kMaxRejectionPerStar = 8000;
constexpr double kWeightFloor = 0.12;
constexpr double kBarInnerScale = 0.28;  // exp(-(r/(kBarInnerScale*R))^2) concentration
constexpr double kSpiralLogOffset = 1e-9;

GalaxyInitAudit g_last_audit;

struct TemplateFlags {
  bool m2 = false;
  bool m3 = false;
  bool bar = false;
  bool spiral = false;
  bool clumps = false;
};

/** Template-specific defaults when structured/noise keys are still at Config neutrals (see config.hpp). */
constexpr double kTemplateDefaultM2Amp = 0.08;
constexpr double kTemplateDefaultM3Amp = 0.08;
constexpr double kTemplateDefaultBarAmp = 0.08;
constexpr double kTemplateDefaultBarAxisRatio = 1.35;
constexpr double kTemplateDefaultSpiralAmp = 0.06;
constexpr double kTemplateDefaultSpiralWinding = 1.2;
constexpr double kTemplateDefaultClumpiness = 0.28;
constexpr int kTemplateDefaultNumClumps = 12;
/** "Warm" noise preset (configs/galaxy_init_presets.cfg). */
constexpr double kTemplateNoisyPos = 0.01;
constexpr double kTemplateNoisyAng = 0.08;
constexpr double kTemplateNoisyMag = 0.05;

TemplateFlags flags_for_template(GalaxyInitTemplate t) {
  TemplateFlags f;
  switch (t) {
    case GalaxyInitTemplate::symmetric_disk:
    case GalaxyInitTemplate::symmetric_disk_noisy:
      break;
    case GalaxyInitTemplate::clumpy_disk:
      f.clumps = true;
      break;
    case GalaxyInitTemplate::weak_m2:
      f.m2 = true;
      break;
    case GalaxyInitTemplate::weak_m3:
      f.m3 = true;
      break;
    case GalaxyInitTemplate::weak_bar:
      f.bar = true;
      break;
    case GalaxyInitTemplate::preformed_spiral:
      f.spiral = true;
      break;
  }
  return f;
}

double clamp_axis_ratio(double a) {
  if (!(a > 0.0) || !std::isfinite(a)) return 1.0;
  if (a < 0.05) return 0.05;
  if (a > 25.0) return 25.0;
  return a;
}

double structured_weight(const TemplateFlags& tf, double r, double theta, double R, double r_min,
                         const Config& cfg) {
  double w = 1.0;
  if (tf.m2)
    w *= 1.0 + cfg.galaxy_init_m2_amplitude * std::cos(2.0 * theta);
  if (tf.m3)
    w *= 1.0 + cfg.galaxy_init_m3_amplitude * std::cos(3.0 * theta);
  if (tf.bar) {
    w *= 1.0 + cfg.galaxy_init_bar_amplitude * std::cos(2.0 * theta);
    double rr = r / std::max(R * kBarInnerScale, 1e-30);
    w *= 1.0 + 0.5 * cfg.galaxy_init_bar_amplitude * std::exp(-rr * rr);
  }
  if (tf.spiral) {
    double lr = std::log(r / r_min + kSpiralLogOffset);
    double phi = 2.0 * theta + cfg.galaxy_init_spiral_winding * lr + cfg.galaxy_init_spiral_phase;
    w *= 1.0 + cfg.galaxy_init_spiral_amplitude * std::cos(phi);
  }
  return std::max(w, kWeightFloor);
}

double structured_w_max(const TemplateFlags& tf, const Config& cfg) {
  double w = 1.0;
  if (tf.m2) w *= 1.0 + std::abs(cfg.galaxy_init_m2_amplitude);
  if (tf.m3) w *= 1.0 + std::abs(cfg.galaxy_init_m3_amplitude);
  if (tf.bar) {
    w *= 1.0 + std::abs(cfg.galaxy_init_bar_amplitude);
    w *= 1.0 + 0.5 * std::abs(cfg.galaxy_init_bar_amplitude);
  }
  if (tf.spiral) w *= 1.0 + std::abs(cfg.galaxy_init_spiral_amplitude);
  return std::max(w, kWeightFloor);
}

const char* template_debug_name(GalaxyInitTemplate t) {
  switch (t) {
    case GalaxyInitTemplate::symmetric_disk:
      return "symmetric_disk";
    case GalaxyInitTemplate::symmetric_disk_noisy:
      return "symmetric_disk_noisy";
    case GalaxyInitTemplate::clumpy_disk:
      return "clumpy_disk";
    case GalaxyInitTemplate::weak_m2:
      return "weak_m2";
    case GalaxyInitTemplate::weak_m3:
      return "weak_m3";
    case GalaxyInitTemplate::weak_bar:
      return "weak_bar";
    case GalaxyInitTemplate::preformed_spiral:
      return "preformed_spiral";
  }
  return "symmetric_disk";
}

void log_warn_ic(GalaxyInitTemplateDefaultsLog* log_out, const std::string& msg) {
  if (log_out) log_out->warnings.push_back(msg);
  std::cerr << "galaxy_init: " << msg << "\n";
}

void log_applied_ic(GalaxyInitTemplateDefaultsLog* log_out, const std::string& msg) {
  if (log_out) log_out->applied.push_back(msg);
}

void validate_structured_amplitudes(GalaxyInitTemplate tmpl, const Config& g, GalaxyInitTemplateDefaultsLog* log_out) {
  TemplateFlags tf = flags_for_template(tmpl);
  bool need = tf.m2 || tf.m3 || tf.bar || tf.spiral;
  if (!need) return;
  double w = structured_w_max(tf, g);
  if (w <= 1.0 + 1e-12) {
    log_warn_ic(log_out, std::string("template '") + template_debug_name(tmpl) +
                             "' has structured flags but rejection weight cap ~1 (check amplitudes / template defaults).");
  }
}

void apply_galaxy_init_template_defaults_impl(GalaxyInitTemplate tmpl, Config& g,
                                               GalaxyInitTemplateDefaultsLog* log_out) {
  switch (tmpl) {
    case GalaxyInitTemplate::symmetric_disk:
      break;
    case GalaxyInitTemplate::symmetric_disk_noisy: {
      bool all_new_zero = (g.galaxy_init_position_noise == 0.0 && g.galaxy_init_velocity_angle_noise == 0.0 &&
                           g.galaxy_init_velocity_magnitude_noise == 0.0);
      if (all_new_zero) {
        if (g.velocity_noise > 0.0) {
          log_warn_ic(log_out,
                      "symmetric_disk_noisy: applying default galaxy_init_* noise; setting legacy velocity_noise "
                      "to 0 (was " +
                          std::to_string(g.velocity_noise) + ") so the new noise pipeline is used.");
        }
        g.galaxy_init_position_noise = kTemplateNoisyPos;
        g.galaxy_init_velocity_angle_noise = kTemplateNoisyAng;
        g.galaxy_init_velocity_magnitude_noise = kTemplateNoisyMag;
        g.velocity_noise = 0.0;
        log_applied_ic(log_out, "symmetric_disk_noisy: position_noise=" + std::to_string(kTemplateNoisyPos) +
                                    ", velocity_angle_noise=" + std::to_string(kTemplateNoisyAng) +
                                    ", velocity_magnitude_noise=" + std::to_string(kTemplateNoisyMag) +
                                    ", velocity_noise=0");
      }
      break;
    }
    case GalaxyInitTemplate::clumpy_disk: {
      bool had_zero_clump = (g.galaxy_init_clumpiness == 0.0);
      if (had_zero_clump) {
        g.galaxy_init_clumpiness = kTemplateDefaultClumpiness;
        log_applied_ic(log_out, "clumpy_disk: galaxy_init_clumpiness=" + std::to_string(kTemplateDefaultClumpiness));
        if (g.galaxy_init_num_clumps == 8) {
          g.galaxy_init_num_clumps = kTemplateDefaultNumClumps;
          log_applied_ic(log_out,
                         "clumpy_disk: galaxy_init_num_clumps=" + std::to_string(kTemplateDefaultNumClumps));
        }
      }
      break;
    }
    case GalaxyInitTemplate::weak_m2:
      if (g.galaxy_init_m2_amplitude == 0.0) {
        g.galaxy_init_m2_amplitude = kTemplateDefaultM2Amp;
        log_applied_ic(log_out, "weak_m2: galaxy_init_m2_amplitude=" + std::to_string(kTemplateDefaultM2Amp));
      }
      break;
    case GalaxyInitTemplate::weak_m3:
      if (g.galaxy_init_m3_amplitude == 0.0) {
        g.galaxy_init_m3_amplitude = kTemplateDefaultM3Amp;
        log_applied_ic(log_out, "weak_m3: galaxy_init_m3_amplitude=" + std::to_string(kTemplateDefaultM3Amp));
      }
      break;
    case GalaxyInitTemplate::weak_bar:
      if (g.galaxy_init_bar_amplitude == 0.0) {
        g.galaxy_init_bar_amplitude = kTemplateDefaultBarAmp;
        log_applied_ic(log_out, "weak_bar: galaxy_init_bar_amplitude=" + std::to_string(kTemplateDefaultBarAmp));
      }
      if (g.galaxy_init_bar_axis_ratio == 1.0) {
        g.galaxy_init_bar_axis_ratio = kTemplateDefaultBarAxisRatio;
        log_applied_ic(log_out,
                       "weak_bar: galaxy_init_bar_axis_ratio=" + std::to_string(kTemplateDefaultBarAxisRatio));
      }
      break;
    case GalaxyInitTemplate::preformed_spiral:
      if (g.galaxy_init_spiral_amplitude == 0.0) {
        g.galaxy_init_spiral_amplitude = kTemplateDefaultSpiralAmp;
        log_applied_ic(log_out,
                       "preformed_spiral: galaxy_init_spiral_amplitude=" + std::to_string(kTemplateDefaultSpiralAmp));
      }
      if (g.galaxy_init_spiral_winding == 1.0 && g.galaxy_init_spiral_phase == 0.0) {
        g.galaxy_init_spiral_winding = kTemplateDefaultSpiralWinding;
        log_applied_ic(log_out,
                       "preformed_spiral: galaxy_init_spiral_winding=" + std::to_string(kTemplateDefaultSpiralWinding) +
                           " (neutral 1.0 + phase 0)");
      }
      break;
    default:
      break;
  }

  if (tmpl == GalaxyInitTemplate::clumpy_disk && g.galaxy_init_clumpiness == 0.0) {
    log_warn_ic(log_out, "clumpy_disk: effective clumpiness is still zero after resolution.");
  }
  validate_structured_amplitudes(tmpl, g, log_out);
}

void apply_bar_axis_stretch(double& x, double& y, double axis_ratio, double r_min, double R) {
  if (std::abs(axis_ratio - 1.0) < 1e-9) return;
  double s = std::sqrt(axis_ratio);
  double invs = 1.0 / s;
  x *= s;
  y *= invs;
  double r = std::hypot(x, y);
  if (r < 1e-300) return;
  double scale = 1.0;
  if (r < r_min) scale = r_min / r;
  else if (r > R) scale = R / r;
  x *= scale;
  y *= scale;
}

void apply_position_jitter_polar(double& x, double& y, double R, double r_min,
                                 double eff_pos_noise, std::mt19937& rng,
                                 std::normal_distribution<double>& normal) {
  if (eff_pos_noise <= 0.0 || !std::isfinite(eff_pos_noise)) return;
  double r = std::hypot(x, y);
  double th = std::atan2(y, x);
  double span_r = std::max(R - r_min, 1e-30);
  r += eff_pos_noise * span_r * normal(rng);
  th += eff_pos_noise * kPi * normal(rng);
  if (r < r_min) r = r_min;
  if (r > R) r = R;
  x = r * std::cos(th);
  y = r * std::sin(th);
}

void rotate_velocity_tangent_to_radial_mix(double& vx, double& vy, double cos_t, double sin_t,
                                           double delta_theta) {
  /* tangential CCW: t = (-sin_t, cos_t); radial: (cos_t, sin_t). Small delta mixes them. */
  double tx = -sin_t;
  double ty = cos_t;
  double rx = cos_t;
  double ry = sin_t;
  double cd = std::cos(delta_theta);
  double sd = std::sin(delta_theta);
  double ex = cd * tx + sd * rx;
  double ey = cd * ty + sd * ry;
  double v = std::hypot(vx, vy);
  vx = v * ex;
  vy = v * ey;
}

}  // namespace

GalaxyInitTemplate parse_galaxy_init_template(const std::string& s) {
  size_t a = 0;
  while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  size_t b = s.size();
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  std::string t;
  for (size_t i = a; i < b; ++i) {
    char c = s[i];
    if (c >= 'A' && c <= 'Z') t += static_cast<char>(c - 'A' + 'a');
    else t += c;
  }
  if (t == "symmetric_disk") return GalaxyInitTemplate::symmetric_disk;
  if (t == "symmetric_disk_noisy") return GalaxyInitTemplate::symmetric_disk_noisy;
  if (t == "clumpy_disk") return GalaxyInitTemplate::clumpy_disk;
  if (t == "weak_m2") return GalaxyInitTemplate::weak_m2;
  if (t == "weak_m3") return GalaxyInitTemplate::weak_m3;
  if (t == "weak_bar") return GalaxyInitTemplate::weak_bar;
  if (t == "preformed_spiral") return GalaxyInitTemplate::preformed_spiral;
  throw std::runtime_error("Unknown galaxy_init_template: " + s);
}

std::string galaxy_init_template_to_string(GalaxyInitTemplate t) {
  switch (t) {
    case GalaxyInitTemplate::symmetric_disk: return "symmetric_disk";
    case GalaxyInitTemplate::symmetric_disk_noisy: return "symmetric_disk_noisy";
    case GalaxyInitTemplate::clumpy_disk: return "clumpy_disk";
    case GalaxyInitTemplate::weak_m2: return "weak_m2";
    case GalaxyInitTemplate::weak_m3: return "weak_m3";
    case GalaxyInitTemplate::weak_bar: return "weak_bar";
    case GalaxyInitTemplate::preformed_spiral: return "preformed_spiral";
  }
  return "symmetric_disk";
}

const GalaxyInitAudit& last_galaxy_init_audit() { return g_last_audit; }

void apply_galaxy_init_template_defaults(GalaxyInitTemplate tmpl, Config& effective,
                                         GalaxyInitTemplateDefaultsLog* log_out) {
  apply_galaxy_init_template_defaults_impl(tmpl, effective, log_out);
}

void sync_config_galaxy_init_from_last_audit(Config& config) {
  const GalaxyInitAudit& a = g_last_audit;
  if (!a.valid) return;
  config.galaxy_init_clumpiness = a.galaxy_init_clumpiness;
  config.galaxy_init_num_clumps = a.galaxy_init_num_clumps;
  config.galaxy_init_clump_radius_fraction = a.galaxy_init_clump_radius_fraction;
  config.galaxy_init_m2_amplitude = a.galaxy_init_m2_amplitude;
  config.galaxy_init_m3_amplitude = a.galaxy_init_m3_amplitude;
  config.galaxy_init_bar_amplitude = a.galaxy_init_bar_amplitude;
  config.galaxy_init_bar_axis_ratio = a.galaxy_init_bar_axis_ratio;
  config.galaxy_init_spiral_amplitude = a.galaxy_init_spiral_amplitude;
  config.galaxy_init_spiral_winding = a.galaxy_init_spiral_winding;
  config.galaxy_init_spiral_phase = a.galaxy_init_spiral_phase;
  config.galaxy_init_position_noise = a.galaxy_init_position_noise_resolved;
  config.galaxy_init_velocity_angle_noise = a.galaxy_init_velocity_angle_noise_resolved;
  config.galaxy_init_velocity_magnitude_noise = a.galaxy_init_velocity_magnitude_noise_resolved;
  config.velocity_noise = a.velocity_noise_effective;
}

void initialize_galaxy_disk(const Config& config, State& state, GalaxyInitAudit* audit_out) {
  GalaxyInitAudit audit;
  audit.valid = true;
  GalaxyInitTemplate tmpl = parse_galaxy_init_template(config.galaxy_init_template);
  audit.template_name = galaxy_init_template_to_string(tmpl);

  audit.galaxy_init_m2_amplitude_raw = config.galaxy_init_m2_amplitude;
  audit.galaxy_init_m3_amplitude_raw = config.galaxy_init_m3_amplitude;
  audit.galaxy_init_bar_amplitude_raw = config.galaxy_init_bar_amplitude;
  audit.galaxy_init_bar_axis_ratio_raw = config.galaxy_init_bar_axis_ratio;
  audit.galaxy_init_spiral_amplitude_raw = config.galaxy_init_spiral_amplitude;
  audit.galaxy_init_spiral_winding_raw = config.galaxy_init_spiral_winding;
  audit.galaxy_init_spiral_phase_raw = config.galaxy_init_spiral_phase;
  audit.galaxy_init_clumpiness_raw = config.galaxy_init_clumpiness;
  audit.galaxy_init_num_clumps_raw = config.galaxy_init_num_clumps;
  audit.galaxy_init_clump_radius_fraction_raw = config.galaxy_init_clump_radius_fraction;
  audit.velocity_noise_raw = config.velocity_noise;

  Config effective = config;
  GalaxyInitTemplateDefaultsLog tlog;
  apply_galaxy_init_template_defaults_impl(tmpl, effective, &tlog);
  audit.template_defaults_used = !tlog.applied.empty();
  audit.template_defaults_log = std::move(tlog);

  audit.seed = effective.galaxy_init_seed;
  audit.master_chaos = effective.galaxy_init_master_chaos;
  if (!(audit.master_chaos > 0.0) || !std::isfinite(audit.master_chaos)) audit.master_chaos = 1.0;

  audit.raw_position_noise = config.galaxy_init_position_noise;
  audit.raw_velocity_angle_noise = config.galaxy_init_velocity_angle_noise;
  audit.raw_velocity_magnitude_noise = config.galaxy_init_velocity_magnitude_noise;

  audit.galaxy_init_position_noise_resolved = effective.galaxy_init_position_noise;
  audit.galaxy_init_velocity_angle_noise_resolved = effective.galaxy_init_velocity_angle_noise;
  audit.galaxy_init_velocity_magnitude_noise_resolved = effective.galaxy_init_velocity_magnitude_noise;
  audit.velocity_noise_effective = effective.velocity_noise;

  audit.eff_position_noise = effective.galaxy_init_position_noise * audit.master_chaos;
  audit.eff_velocity_angle_noise_rad = effective.galaxy_init_velocity_angle_noise * audit.master_chaos;
  audit.eff_velocity_magnitude_noise = effective.galaxy_init_velocity_magnitude_noise * audit.master_chaos;

  audit.master_scales_position_noise = (effective.galaxy_init_position_noise > 0.0);
  audit.master_scales_velocity_angle_noise = (effective.galaxy_init_velocity_angle_noise > 0.0);
  audit.master_scales_velocity_magnitude_noise = (effective.galaxy_init_velocity_magnitude_noise > 0.0);

  audit.used_new_state_noise =
      (audit.eff_position_noise > 0.0 || audit.eff_velocity_angle_noise_rad > 0.0 ||
       audit.eff_velocity_magnitude_noise > 0.0);
  audit.used_legacy_velocity_noise =
      !audit.used_new_state_noise && (effective.velocity_noise > 0.0);

  TemplateFlags tf = flags_for_template(tmpl);
  audit.structured_m2 = tf.m2;
  audit.structured_m3 = tf.m3;
  audit.structured_bar = tf.bar;
  audit.structured_spiral = tf.spiral;
  audit.structured_clumps = tf.clumps;

  audit.galaxy_init_clumpiness = effective.galaxy_init_clumpiness;
  audit.galaxy_init_num_clumps = std::max(1, effective.galaxy_init_num_clumps);
  audit.galaxy_init_clump_radius_fraction = effective.galaxy_init_clump_radius_fraction;
  audit.galaxy_init_m2_amplitude = effective.galaxy_init_m2_amplitude;
  audit.galaxy_init_m3_amplitude = effective.galaxy_init_m3_amplitude;
  audit.galaxy_init_bar_amplitude = effective.galaxy_init_bar_amplitude;
  audit.galaxy_init_bar_axis_ratio = effective.galaxy_init_bar_axis_ratio;
  audit.galaxy_init_spiral_amplitude = effective.galaxy_init_spiral_amplitude;
  audit.galaxy_init_spiral_winding = effective.galaxy_init_spiral_winding;
  audit.galaxy_init_spiral_phase = effective.galaxy_init_spiral_phase;

  const int n = config.n_stars;
  state.resize(n);

  std::mt19937 rng(effective.galaxy_init_seed);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);

  const double R = effective.galaxy_radius;
  const double r_min = kDiskInnerFraction * R;
  const double bh_mass = effective.bh_mass;
  const double star_mass = effective.star_mass;

  double w_max = structured_w_max(tf, effective);
  audit.weight_w_max = w_max;

  /* Clump centers (disk annulus, uniform area per clump center draw) */
  audit.clump_centers_xy.clear();
  int n_clumps = tf.clumps ? audit.galaxy_init_num_clumps : 0;
  double clump_p = tf.clumps ? std::max(0.0, std::min(1.0, effective.galaxy_init_clumpiness)) : 0.0;
  double sigma_clump = std::max(effective.galaxy_init_clump_radius_fraction, 1e-6) * R;

  for (int c = 0; c < n_clumps; ++c) {
    double u = u01(rng);
    double r_sq = r_min * r_min + u * (R * R - r_min * r_min);
    double rc = std::sqrt(r_sq);
    double thc = kTwoPi * u01(rng);
    audit.clump_centers_xy.push_back({rc * std::cos(thc), rc * std::sin(thc)});
  }

  std::vector<double> radii(n);
  std::vector<double> theta(n);

  audit.rejection_fallbacks = 0;

  for (int i = 0; i < n; ++i) {
    double x = 0.0, y = 0.0;
    bool placed_clump = false;

    if (tf.clumps && clump_p > 0.0 && u01(rng) < clump_p && n_clumps > 0) {
      int k = static_cast<int>(u01(rng) * n_clumps);
      if (k >= n_clumps) k = n_clumps - 1;
      double cx = audit.clump_centers_xy[k].first;
      double cy = audit.clump_centers_xy[k].second;
      for (int attempt = 0; attempt < kMaxRejectionPerStar; ++attempt) {
        x = cx + sigma_clump * normal(rng);
        y = cy + sigma_clump * normal(rng);
        double r = std::hypot(x, y);
        if (r >= r_min && r <= R) {
          placed_clump = true;
          break;
        }
      }
      if (!placed_clump) {
        ++audit.rejection_fallbacks;
        double u = u01(rng);
        double r_sq = r_min * r_min + u * (R * R - r_min * r_min);
        double r = std::sqrt(r_sq);
        double th = kTwoPi * u01(rng);
        x = r * std::cos(th);
        y = r * std::sin(th);
      }
    } else {
      bool need_structure = tf.m2 || tf.m3 || tf.bar || tf.spiral;
      bool placed_struct = false;
      if (!need_structure) {
        double u = u01(rng);
        double r_sq = r_min * r_min + u * (R * R - r_min * r_min);
        double r = std::sqrt(r_sq);
        double th = kTwoPi * u01(rng);
        x = r * std::cos(th);
        y = r * std::sin(th);
      } else {
        for (int attempt = 0; attempt < kMaxRejectionPerStar; ++attempt) {
          double u = u01(rng);
          double r_sq = r_min * r_min + u * (R * R - r_min * r_min);
          double r = std::sqrt(r_sq);
          double th = kTwoPi * u01(rng);
          double w = structured_weight(tf, r, th, R, r_min, effective);
          if (u01(rng) < w / w_max) {
            x = r * std::cos(th);
            y = r * std::sin(th);
            placed_struct = true;
            break;
          }
        }
        if (!placed_struct) {
          ++audit.rejection_fallbacks;
          double u = u01(rng);
          double r_sq = r_min * r_min + u * (R * R - r_min * r_min);
          double r = std::sqrt(r_sq);
          double th = kTwoPi * u01(rng);
          x = r * std::cos(th);
          y = r * std::sin(th);
        }
      }
    }

    if (tf.bar) {
      apply_bar_axis_stretch(x, y, clamp_axis_ratio(effective.galaxy_init_bar_axis_ratio), r_min, R);
    }

    apply_position_jitter_polar(x, y, R, r_min, audit.eff_position_noise, rng, normal);

    state.x[i] = x;
    state.y[i] = y;
    state.mass[i] = star_mass;
    radii[i] = std::hypot(x, y);
    theta[i] = std::atan2(y, x);
  }

  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&radii](int a, int b) { return radii[a] < radii[b]; });
  std::vector<int> n_inside(n);
  for (int k = 0; k < n; ++k) n_inside[order[k]] = k;

  const bool use_derived_tpf_init =
      (effective.physics_package == "TPFCore" && effective.tpfcore_enable_provisional_readout &&
       tpfcore::is_derived_tpf_radial_readout_mode(effective.tpfcore_readout_mode));

  if (use_derived_tpf_init) {
    const double eps =
        (effective.tpfcore_source_softening > 0.0) ? effective.tpfcore_source_softening : effective.softening;
    tpfcore::DerivedTpfPoissonConfig dcfg;
    dcfg.kappa = effective.tpf_kappa;
    dcfg.bins = effective.tpf_poisson_bins;
    dcfg.max_radius = effective.tpf_poisson_max_radius;
    dcfg.galaxy_radius = effective.galaxy_radius;
    tpfcore::TpfRadialGravityProfile profile =
        tpfcore::build_tpf_gravity_profile(state, bh_mass, dcfg, eps);
    for (int i = 0; i < n; ++i) {
      double th = theta[i];
      double cos_t = std::cos(th);
      double sin_t = std::sin(th);
      double r_cyl = std::hypot(state.x[i], state.y[i]);
      double a_s =
          tpfcore::radial_acceleration_scalar_derived(state, bh_mass, profile, r_cyl, eps);
      double v_circ = std::sqrt(std::max(0.0, std::abs(a_s) * r_cyl));
      double vx = -sin_t * v_circ;
      double vy = cos_t * v_circ;

      if (audit.used_new_state_noise) {
        double dth = audit.eff_velocity_angle_noise_rad * normal(rng);
        rotate_velocity_tangent_to_radial_mix(vx, vy, cos_t, sin_t, dth);
        double mag_scale = 1.0 + audit.eff_velocity_magnitude_noise * normal(rng);
        if (mag_scale < 0.05) mag_scale = 0.05;
        vx *= mag_scale;
        vy *= mag_scale;
      } else if (audit.used_legacy_velocity_noise) {
        double scale = effective.velocity_noise * v_circ;
        vx += scale * normal(rng);
        vy += scale * normal(rng);
      }

      state.vx[i] = effective.initial_velocity_scale * vx;
      state.vy[i] = effective.initial_velocity_scale * vy;
    }
  } else {
    constexpr double G_SI = 6.6743e-11;
    for (int i = 0; i < n; ++i) {
      double r = radii[i];
      double th = theta[i];
      double cos_t = std::cos(th);
      double sin_t = std::sin(th);
      double enclosed_mass = bh_mass + n_inside[i] * star_mass;
      double v_circ = std::sqrt((G_SI * enclosed_mass) / std::max(r, 1e-30));
      double vx = -sin_t * v_circ;
      double vy = cos_t * v_circ;

      if (audit.used_new_state_noise) {
        double dth = audit.eff_velocity_angle_noise_rad * normal(rng);
        rotate_velocity_tangent_to_radial_mix(vx, vy, cos_t, sin_t, dth);
        double mag_scale = 1.0 + audit.eff_velocity_magnitude_noise * normal(rng);
        if (mag_scale < 0.05) mag_scale = 0.05;
        vx *= mag_scale;
        vy *= mag_scale;
      } else if (audit.used_legacy_velocity_noise) {
        double scale = effective.velocity_noise * v_circ;
        vx += scale * normal(rng);
        vy += scale * normal(rng);
      }

      state.vx[i] = effective.initial_velocity_scale * vx;
      state.vy[i] = effective.initial_velocity_scale * vy;
    }
  }

  g_last_audit = audit;
  if (audit_out) *audit_out = audit;
}

void write_galaxy_init_diagnostics(const std::string& output_dir,
                                   const State& state,
                                   const Config& config,
                                   const GalaxyInitAudit& audit) {
  const int n = state.n();
  if (n <= 0) return;

  std::vector<double> r(n), vmag(n), lz(n);
  double lz_sum = 0.0;
  for (int i = 0; i < n; ++i) {
    r[i] = std::hypot(state.x[i], state.y[i]);
    vmag[i] = std::hypot(state.vx[i], state.vy[i]);
    double lzi = state.mass[i] * (state.x[i] * state.vy[i] - state.y[i] * state.vx[i]);
    lz[i] = lzi;
    lz_sum += lzi;
  }

  auto mean = [](const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return s / std::max(1, static_cast<int>(v.size()));
  };

  double r_mean = mean(r);
  double v_mean = mean(vmag);

  std::vector<double> rs = r;
  std::sort(rs.begin(), rs.end());
  std::vector<double> vs = vmag;
  std::sort(vs.begin(), vs.end());

  std::ostringstream diag_path;
  diag_path << output_dir << "/galaxy_init_diagnostics.txt";
  std::ofstream txt(diag_path.str());
  if (!txt) return;

  txt << std::scientific << std::setprecision(16);
  txt << "=== Galaxy initial-condition diagnostics ===\n";
  txt << "template\t" << audit.template_name << "\n";
  txt << "galaxy_init_seed\t" << audit.seed << "\n";
  txt << "n_stars\t" << n << "\n";
  txt << "galaxy_radius\t" << config.galaxy_radius << "\n";
  txt << "disk_annulus_inner_fraction\t" << kDiskInnerFraction << "\n";
  txt << "--- Radial r = hypot(x,y) (m) ---\n";
  txt << "r_min\t" << rs.front() << "\n";
  txt << "r_max\t" << rs.back() << "\n";
  txt << "r_mean\t" << r_mean << "\n";
  txt << "r_median\t" << rs[rs.size() / 2] << "\n";
  if (!rs.empty()) {
    size_t i10 = (rs.size() - 1) * 10 / 100;
    size_t i90 = (rs.size() - 1) * 90 / 100;
    txt << "r_p10\t" << rs[i10] << "\n";
    txt << "r_p90\t" << rs[i90] << "\n";
  }
  txt << "--- Speed |v| (m/s) ---\n";
  txt << "v_min\t" << vs.front() << "\n";
  txt << "v_max\t" << vs.back() << "\n";
  txt << "v_mean\t" << v_mean << "\n";
  txt << "v_median\t" << vs[vs.size() / 2] << "\n";
  txt << "--- Angular momentum (z) ---\n";
  txt << "L_z_total\t" << lz_sum << "  (sum_i m_i (x_i v_y_i - y_i v_x_i))\n";
  txt << "L_z_mean_per_star\t" << (lz_sum / static_cast<double>(n)) << "\n";
  txt << "--- IC template resolution (raw config vs effective) ---\n";
  txt << "template_defaults_used\t" << (audit.template_defaults_used ? 1 : 0) << "\n";
  txt << "galaxy_init_m2_amplitude_raw\t" << audit.galaxy_init_m2_amplitude_raw << "\n";
  txt << "galaxy_init_m3_amplitude_raw\t" << audit.galaxy_init_m3_amplitude_raw << "\n";
  txt << "galaxy_init_bar_amplitude_raw\t" << audit.galaxy_init_bar_amplitude_raw << "\n";
  txt << "galaxy_init_bar_axis_ratio_raw\t" << audit.galaxy_init_bar_axis_ratio_raw << "\n";
  txt << "galaxy_init_spiral_amplitude_raw\t" << audit.galaxy_init_spiral_amplitude_raw << "\n";
  txt << "galaxy_init_spiral_winding_raw\t" << audit.galaxy_init_spiral_winding_raw << "\n";
  txt << "galaxy_init_spiral_phase_raw\t" << audit.galaxy_init_spiral_phase_raw << "\n";
  txt << "galaxy_init_clumpiness_raw\t" << audit.galaxy_init_clumpiness_raw << "\n";
  txt << "galaxy_init_num_clumps_raw\t" << audit.galaxy_init_num_clumps_raw << "\n";
  txt << "galaxy_init_clump_radius_fraction_raw\t" << audit.galaxy_init_clump_radius_fraction_raw << "\n";
  txt << "velocity_noise_raw\t" << audit.velocity_noise_raw << "\n";
  txt << "galaxy_init_m2_amplitude_effective\t" << audit.galaxy_init_m2_amplitude << "\n";
  txt << "galaxy_init_m3_amplitude_effective\t" << audit.galaxy_init_m3_amplitude << "\n";
  txt << "galaxy_init_bar_amplitude_effective\t" << audit.galaxy_init_bar_amplitude << "\n";
  txt << "galaxy_init_bar_axis_ratio_effective\t" << audit.galaxy_init_bar_axis_ratio << "\n";
  txt << "galaxy_init_spiral_amplitude_effective\t" << audit.galaxy_init_spiral_amplitude << "\n";
  txt << "galaxy_init_spiral_winding_effective\t" << audit.galaxy_init_spiral_winding << "\n";
  txt << "galaxy_init_spiral_phase_effective\t" << audit.galaxy_init_spiral_phase << "\n";
  txt << "galaxy_init_clumpiness_effective\t" << audit.galaxy_init_clumpiness << "\n";
  txt << "galaxy_init_num_clumps_effective\t" << audit.galaxy_init_num_clumps << "\n";
  txt << "galaxy_init_clump_radius_fraction_effective\t" << audit.galaxy_init_clump_radius_fraction << "\n";
  txt << "velocity_noise_effective\t" << audit.velocity_noise_effective << "\n";
  txt << "galaxy_init_position_noise_resolved\t" << audit.galaxy_init_position_noise_resolved << "\n";
  txt << "galaxy_init_velocity_angle_noise_resolved\t" << audit.galaxy_init_velocity_angle_noise_resolved
      << "\n";
  txt << "galaxy_init_velocity_magnitude_noise_resolved\t" << audit.galaxy_init_velocity_magnitude_noise_resolved
      << "\n";
  for (size_t i = 0; i < audit.template_defaults_log.applied.size(); ++i) {
    txt << "template_default_applied_" << i << "\t" << audit.template_defaults_log.applied[i] << "\n";
  }
  for (size_t i = 0; i < audit.template_defaults_log.warnings.size(); ++i) {
    txt << "template_warning_" << i << "\t" << audit.template_defaults_log.warnings[i] << "\n";
  }
  txt << "--- Audit ---\n";
  txt << "used_new_state_noise\t" << (audit.used_new_state_noise ? 1 : 0) << "\n";
  txt << "used_legacy_velocity_noise\t" << (audit.used_legacy_velocity_noise ? 1 : 0) << "\n";
  txt << "master_chaos\t" << audit.master_chaos << "\n";
  txt << "eff_position_noise\t" << audit.eff_position_noise << "\n";
  txt << "eff_velocity_angle_noise_rad\t" << audit.eff_velocity_angle_noise_rad << "\n";
  txt << "eff_velocity_magnitude_noise\t" << audit.eff_velocity_magnitude_noise << "\n";
  txt << "rejection_fallbacks\t" << audit.rejection_fallbacks << "\n";
  txt << "structured_weight_cap\t" << audit.weight_w_max << "\n";
  if (!audit.clump_centers_xy.empty()) {
    txt << "clump_centers_xy_count\t" << audit.clump_centers_xy.size() << "\n";
    for (size_t i = 0; i < audit.clump_centers_xy.size(); ++i) {
      txt << "clump_center_" << i << "_x\t" << audit.clump_centers_xy[i].first << "\n";
      txt << "clump_center_" << i << "_y\t" << audit.clump_centers_xy[i].second << "\n";
    }
  }
  txt << "=== End galaxy_init_diagnostics ===\n";

  if (!config.save_snapshots) {
    std::ostringstream csv_path;
    csv_path << output_dir << "/galaxy_init_snapshot.csv";
    std::ofstream csv(csv_path.str());
    if (csv) {
      csv << "# step,0,time,0.000000e+00\n";
      csv << "i,x,y,vx,vy,mass\n";
      csv << std::scientific << std::setprecision(16);
      for (int i = 0; i < n; ++i) {
        csv << i << "," << state.x[i] << "," << state.y[i] << "," << state.vx[i] << "," << state.vy[i] << ","
            << state.mass[i] << "\n";
      }
    }
  }
}

}  // namespace galaxy
