#include "softening_policy.hpp"
#include <cmath>
#include <stdexcept>
namespace galaxy {
namespace { constexpr double kSolarRadiusM = 6.957e8; }

double plummer_softening_scale(double r_sq, double softening) {
  if (softening <= 0.0) return 1.0;
  const double eps2 = softening * softening;
  const double softened = r_sq + eps2;
  if (r_sq <= 0.0 || softened <= 0.0) return 0.0;
  const double ratio = r_sq / softened;
  return ratio * std::sqrt(ratio);
}

void apply_plummer_softening(double dx,
                             double dy,
                             double softening,
                             double& ax,
                             double& ay) {
  const double scale = plummer_softening_scale(dx * dx + dy * dy, softening);
  ax *= scale;
  ay *= scale;
}

ResolvedSoftening resolve_softening(const Config& cfg, const State& state) {
  ResolvedSoftening r{}; r.mode = cfg.softening_mode; r.profile = cfg.softening_auto_profile;
  if (r.mode == "off") { r.source = "softening_mode=off"; return r; }
  if (r.mode == "manual") { r.effective_softening = std::max(0.0, cfg.softening); r.source = "softening_mode=manual(config.softening)"; return r; }
  if (r.mode != "auto") throw std::runtime_error("invalid softening_mode");
  if (r.profile.empty()) r.profile = "collisionless";
  if (r.profile == "collisionless") {
    const int dim = (cfg.auto_softening_dimension > 0) ? cfg.auto_softening_dimension : 2;
    if (dim <= 0) throw std::runtime_error("auto_softening_dimension must be > 0");
    const double fac = (cfg.auto_softening_factor > 0.0) ? cfg.auto_softening_factor : 1.8;
    if (!(fac > 0.0)) throw std::runtime_error("auto_softening_factor must be > 0 for collisionless");
    const double rin = (cfg.inner_radius > 0.0) ? cfg.inner_radius : 0.05 * cfg.galaxy_radius;
    const double rout = (cfg.outer_radius > rin) ? cfg.outer_radius : cfg.galaxy_radius;
    const int n = (state.n() > 0) ? state.n() : cfg.n_stars;
    if (n <= 0) throw std::runtime_error("auto softening requires positive particle count");
    if (!(rout > 0.0) || !(rout > rin)) throw std::runtime_error("collisionless auto requires outer_radius > inner_radius and outer_radius > 0");
    const double nball_outer = (dim == 2) ? (M_PI * rout * rout) : ((M_PI * M_PI * rout * rout * rout * rout) / 2.0);
    const double nball_inner = (dim == 2) ? (M_PI * rin * rin) : ((M_PI * M_PI * rin * rin * rin * rin) / 2.0);
    const double volume = nball_outer - nball_inner;
    if (!(volume > 0.0)) throw std::runtime_error("collisionless auto annulus volume must be > 0");
    const double mean_sep = std::pow(volume / static_cast<double>(n), 1.0 / static_cast<double>(dim));
    r.mean_separation = mean_sep; r.radius_inner_used = rin; r.radius_outer_used = rout; r.dimension = dim; r.factor = fac;
    r.effective_softening = mean_sep / fac; r.source = "auto:collisionless(mean_separation/factor)";
  } else if (r.profile == "stellar_physical" || r.profile == "nuclear_cluster") {
    const double def = (r.profile == "stellar_physical") ? 10.0 : 3.0;
    const double fac = (cfg.auto_softening_factor > 0.0) ? cfg.auto_softening_factor : def;
    r.factor = fac; r.dimension = 0; r.mean_separation = 0.0;
    r.effective_softening = fac * kSolarRadiusM;
    r.source = "auto:" + r.profile + "(solar_radius_multiplier)";
  } else throw std::runtime_error("invalid softening_auto_profile");
  double max_cap = 0.0;
  if (cfg.auto_softening_max > 0.0) {
    max_cap = cfg.auto_softening_max;
    r.max_cap_source = "explicit";
  } else if (r.mode == "auto" && r.profile == "collisionless" && r.radius_outer_used > 0.0) {
    max_cap = 0.05 * r.radius_outer_used;
    r.max_cap_source = "contextual";
  }
  if (max_cap > 0.0) {
    const double prev = r.effective_softening;
    r.effective_softening = std::min(r.effective_softening, max_cap);
    r.max_capped = (r.effective_softening < prev);
    if (r.max_capped) r.source += "|max_cap=" + r.max_cap_source;
  }
  if (cfg.auto_softening_min > 0.0) {
    const double prev = r.effective_softening;
    r.effective_softening = std::max(r.effective_softening, cfg.auto_softening_min);
    r.min_floored = (r.effective_softening > prev);
    if (r.min_floored) r.source += "|min_floor=explicit";
  }
  if (!std::isfinite(r.effective_softening) || r.effective_softening < 0.0) throw std::runtime_error("resolved softening invalid");
  return r;
}
}
