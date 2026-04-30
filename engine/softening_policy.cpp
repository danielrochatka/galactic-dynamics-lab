#include "softening_policy.hpp"
#include <cmath>
#include <stdexcept>
namespace galaxy {
namespace { constexpr double kSolarRadiusM = 6.957e8; }
ResolvedSoftening resolve_softening(const Config& cfg, const State& state) {
  ResolvedSoftening r{}; r.mode = cfg.softening_mode; r.profile = cfg.softening_auto_profile;
  if (r.mode == "off") { r.source = "softening_mode=off"; return r; }
  if (r.mode == "manual") { r.effective_softening = std::max(0.0, cfg.softening); r.source = "softening_mode=manual(config.softening)"; return r; }
  if (r.mode != "auto") throw std::runtime_error("invalid softening_mode");
  const bool tpf_like = (cfg.physics_package == "TPFCore" || cfg.tpf_dynamics_mode == "xi_kernel_deformed" || cfg.simulation_mode == SimulationMode::tpf_4d_xi_motion_probe_benchmark);
  if (r.profile.empty()) r.profile = tpf_like ? "stellar_physical" : "collisionless";
  if (r.profile == "collisionless") {
    const int dim = (cfg.auto_softening_dimension > 0) ? cfg.auto_softening_dimension : 2;
    if (dim != 2) throw std::runtime_error("collisionless auto currently supports auto_softening_dimension=2 only");
    const double fac = (cfg.auto_softening_factor > 0.0) ? cfg.auto_softening_factor : 1.8;
    if (!(fac > 0.0)) throw std::runtime_error("auto_softening_factor must be > 0 for collisionless");
    const double rin = (cfg.inner_radius > 0.0) ? cfg.inner_radius : 0.05 * cfg.galaxy_radius;
    const double rout = (cfg.outer_radius > rin) ? cfg.outer_radius : cfg.galaxy_radius;
    const int n = (state.n() > 0) ? state.n() : cfg.n_stars;
    if (n <= 0) throw std::runtime_error("auto softening requires positive particle count");
    if (!(rout > 0.0) || !(rout > rin)) throw std::runtime_error("collisionless auto requires outer_radius > inner_radius and outer_radius > 0");
    const double area = M_PI * (rout * rout - rin * rin);
    if (!(area > 0.0)) throw std::runtime_error("collisionless auto annulus area must be > 0");
    const double mean_sep = std::sqrt(area / static_cast<double>(n));
    r.mean_separation = mean_sep; r.radius_inner_used = rin; r.radius_outer_used = rout; r.dimension = dim; r.factor = fac;
    r.effective_softening = mean_sep / fac; r.source = "auto:collisionless(mean_separation/factor)";
  } else if (r.profile == "stellar_physical" || r.profile == "nuclear_cluster") {
    const double def = (r.profile == "stellar_physical") ? 10.0 : 3.0;
    const double fac = (cfg.auto_softening_factor > 0.0) ? cfg.auto_softening_factor : def;
    r.factor = fac; r.dimension = 0; r.mean_separation = 0.0;
    r.effective_softening = fac * kSolarRadiusM;
    r.source = "auto:" + r.profile + "(solar_radius_multiplier)";
  } else throw std::runtime_error("invalid softening_auto_profile");
  if (cfg.auto_softening_min > 0.0) r.effective_softening = std::max(r.effective_softening, cfg.auto_softening_min);
  if (cfg.auto_softening_max > 0.0) r.effective_softening = std::min(r.effective_softening, cfg.auto_softening_max);
  if (!std::isfinite(r.effective_softening) || r.effective_softening < 0.0) throw std::runtime_error("resolved softening invalid");
  return r;
}
}
