#include "preflight.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

namespace galaxy {

GalaxyPreflightSummary build_galaxy_preflight_summary(const Config& config,
                                                       const ResolvedScenario& resolved) {
  GalaxyPreflightSummary out;
  out.disk_mass = static_cast<double>(config.n_stars) * config.star_mass;
  out.disk_to_bh_mass_ratio = out.disk_mass / config.bh_mass;
  out.total_sim_time = config.dt * resolved.effective_n_steps;
  const double G_SI = 6.6743e-11;
  const double pi = std::acos(-1.0);
  out.outer_orbital_period = 2.0 * pi * std::sqrt(std::pow(config.outer_radius, 3.0) /
                                                   (G_SI * (out.disk_mass + config.bh_mass)));
  out.simulated_outer_orbits = out.total_sim_time / out.outer_orbital_period;
  out.softening_to_outer_radius = resolved.config.softening / config.outer_radius;

  if (out.simulated_outer_orbits < 0.25)
    out.warnings.push_back("simulated_outer_orbits < 0.25");
  if (out.softening_to_outer_radius > 0.05)
    out.warnings.push_back("softening_to_outer_radius > 0.05");
  if (out.disk_to_bh_mass_ratio < 0.1)
    out.warnings.push_back("disk_to_bh_mass_ratio < 0.1");
  if (out.disk_to_bh_mass_ratio > 1000.0)
    out.warnings.push_back("disk_to_bh_mass_ratio > 1000");
  if (config.star_mass > config.bh_mass && config.n_stars > 10)
    out.warnings.push_back("star_mass > bh_mass and n_stars > 10");

  return out;
}

void append_galaxy_preflight_to_run_info(const std::string& output_dir,
                                         const GalaxyPreflightSummary& summary) {
  std::ofstream f(output_dir + "/run_info.txt", std::ios::app);
  if (!f) return;
  f << "\n=== Galaxy preflight sanity check ===\n";
  f << std::scientific << std::setprecision(17);
  f << "preflight_disk_mass\t" << summary.disk_mass << "\n";
  f << "preflight_disk_to_bh_mass_ratio\t" << summary.disk_to_bh_mass_ratio << "\n";
  f << "preflight_total_sim_time\t" << summary.total_sim_time << "\n";
  f << "preflight_outer_orbital_period\t" << summary.outer_orbital_period << "\n";
  f << "preflight_simulated_outer_orbits\t" << summary.simulated_outer_orbits << "\n";
  f << "preflight_softening_to_outer_radius\t" << summary.softening_to_outer_radius << "\n";
  f << "preflight_warning_count\t" << summary.warnings.size() << "\n";
  for (std::size_t i = 0; i < summary.warnings.size(); ++i)
    f << "preflight_warning_" << i << "\t" << summary.warnings[i] << "\n";
}

}
