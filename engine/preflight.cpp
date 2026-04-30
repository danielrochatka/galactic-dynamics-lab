#include "preflight.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace galaxy {
namespace {
constexpr double kG = 6.6743e-11;
constexpr double kMetersPerKpc = 3.0856775814913673e19;
constexpr double kSecondsPerMyr = 31557600.0 * 1.0e6;
}

GalaxyPreflightSummary build_galaxy_preflight_summary(const Config& config,
                                                       const ResolvedScenario& resolved) {
  GalaxyPreflightSummary out;
  out.n_stars = config.n_stars;
  out.star_mass = config.star_mass;
  out.star_mass_solar = config.star_mass / kSolarMassKg;
  out.disk_mass = static_cast<double>(config.n_stars) * config.star_mass;
  out.disk_mass_solar = out.disk_mass / kSolarMassKg;
  out.bh_mass = config.bh_mass;
  out.bh_mass_solar = config.bh_mass / kSolarMassKg;
  out.disk_to_bh_mass_ratio = out.disk_mass / config.bh_mass;
  out.inner_radius_m = config.inner_radius;
  out.inner_radius_kpc = config.inner_radius / kMetersPerKpc;
  out.outer_radius_m = config.outer_radius;
  out.outer_radius_kpc = config.outer_radius / kMetersPerKpc;
  out.galaxy_radius_m = config.galaxy_radius;
  out.galaxy_radius_kpc = config.galaxy_radius / kMetersPerKpc;
  out.effective_softening_m = resolved.config.softening;
  out.effective_softening_kpc = resolved.config.softening / kMetersPerKpc;
  out.softening_to_outer_radius = resolved.config.softening / config.outer_radius;
  out.dt = config.dt;
  out.n_steps = resolved.effective_n_steps;
  out.total_sim_time_s = config.dt * resolved.effective_n_steps;
  out.total_sim_time_myr = out.total_sim_time_s / kSecondsPerMyr;
  const double pi = std::acos(-1.0);
  const double total_mass = out.disk_mass + out.bh_mass;
  out.outer_orbital_period_s = 2.0 * pi * std::sqrt(std::pow(config.outer_radius, 3.0) / (kG * total_mass));
  out.outer_orbital_period_myr = out.outer_orbital_period_s / kSecondsPerMyr;
  out.simulated_outer_orbits = out.total_sim_time_s / out.outer_orbital_period_s;

  if (out.simulated_outer_orbits < 0.25) out.warnings.push_back("simulated_outer_orbits < 0.25");
  if (out.softening_to_outer_radius > 0.05) out.warnings.push_back("softening_to_outer_radius > 0.05");
  if (out.disk_to_bh_mass_ratio < 0.1) out.warnings.push_back("disk_to_bh_mass_ratio < 0.1");
  if (out.disk_to_bh_mass_ratio > 1000.0) out.warnings.push_back("disk_to_bh_mass_ratio > 1000");
  if (config.star_mass > config.bh_mass && config.n_stars > 10) out.warnings.push_back("star_mass > bh_mass and n_stars > 10");
  return out;
}

void print_galaxy_preflight_summary(const GalaxyPreflightSummary& s) {
  std::cout << "Galaxy setup summary:\n" << std::scientific << std::setprecision(6)
            << "  n_stars: " << s.n_stars << "\n"
            << "  star_mass: " << s.star_mass << " kg (" << s.star_mass_solar << " M_sun)\n"
            << "  disk_mass: " << s.disk_mass << " kg (" << s.disk_mass_solar << " M_sun)\n"
            << "  bh_mass: " << s.bh_mass << " kg (" << s.bh_mass_solar << " M_sun)\n"
            << "  disk_to_bh_mass_ratio: " << s.disk_to_bh_mass_ratio << "\n"
            << "  inner_radius: " << s.inner_radius_m << " m (" << s.inner_radius_kpc << " kpc)\n"
            << "  outer_radius: " << s.outer_radius_m << " m (" << s.outer_radius_kpc << " kpc)\n"
            << "  galaxy_radius: " << s.galaxy_radius_m << " m (" << s.galaxy_radius_kpc << " kpc)\n"
            << "  effective_softening: " << s.effective_softening_m << " m (" << s.effective_softening_kpc << " kpc)\n"
            << "  softening_to_outer_radius: " << s.softening_to_outer_radius << "\n"
            << "  dt: " << s.dt << "\n"
            << "  n_steps: " << s.n_steps << "\n"
            << "  total_sim_time: " << s.total_sim_time_s << " s (" << s.total_sim_time_myr << " Myr)\n"
            << "  estimated_outer_orbital_period: " << s.outer_orbital_period_s << " s (" << s.outer_orbital_period_myr << " Myr)\n"
            << "  simulated_outer_orbits: " << s.simulated_outer_orbits << "\n";
}

void append_galaxy_preflight_to_run_info(const std::string& output_dir,
                                         const GalaxyPreflightSummary& s) {
  std::ofstream f(output_dir + "/run_info.txt", std::ios::app);
  if (!f) return;
  f << "\n=== Galaxy preflight sanity check ===\n" << std::scientific << std::setprecision(17)
    << "preflight_n_stars\t" << s.n_stars << "\n"
    << "preflight_star_mass_kg\t" << s.star_mass << "\n"
    << "preflight_star_mass_msun\t" << s.star_mass_solar << "\n"
    << "preflight_disk_mass_kg\t" << s.disk_mass << "\n"
    << "preflight_disk_mass_msun\t" << s.disk_mass_solar << "\n"
    << "preflight_bh_mass_kg\t" << s.bh_mass << "\n"
    << "preflight_bh_mass_msun\t" << s.bh_mass_solar << "\n"
    << "preflight_disk_to_bh_mass_ratio\t" << s.disk_to_bh_mass_ratio << "\n"
    << "preflight_inner_radius_m\t" << s.inner_radius_m << "\n"
    << "preflight_inner_radius_kpc\t" << s.inner_radius_kpc << "\n"
    << "preflight_outer_radius_m\t" << s.outer_radius_m << "\n"
    << "preflight_outer_radius_kpc\t" << s.outer_radius_kpc << "\n"
    << "preflight_galaxy_radius_m\t" << s.galaxy_radius_m << "\n"
    << "preflight_galaxy_radius_kpc\t" << s.galaxy_radius_kpc << "\n"
    << "preflight_effective_softening_m\t" << s.effective_softening_m << "\n"
    << "preflight_effective_softening_kpc\t" << s.effective_softening_kpc << "\n"
    << "preflight_softening_to_outer_radius\t" << s.softening_to_outer_radius << "\n"
    << "preflight_dt\t" << s.dt << "\n"
    << "preflight_n_steps\t" << s.n_steps << "\n"
    << "preflight_total_sim_time_s\t" << s.total_sim_time_s << "\n"
    << "preflight_total_sim_time_myr\t" << s.total_sim_time_myr << "\n"
    << "preflight_outer_orbital_period_s\t" << s.outer_orbital_period_s << "\n"
    << "preflight_outer_orbital_period_myr\t" << s.outer_orbital_period_myr << "\n"
    << "preflight_simulated_outer_orbits\t" << s.simulated_outer_orbits << "\n"
    << "galaxy_preflight_disk_to_bh_mass_ratio\t" << s.disk_to_bh_mass_ratio << "\n"
    << "galaxy_preflight_softening_to_outer_radius\t" << s.softening_to_outer_radius << "\n"
    << "galaxy_preflight_estimated_outer_orbital_period\t" << s.outer_orbital_period_s << "\n"
    << "galaxy_preflight_simulated_outer_orbits\t" << s.simulated_outer_orbits << "\n"
    << "preflight_warning_count\t" << s.warnings.size() << "\n"
    << "galaxy_preflight_warning_count\t" << s.warnings.size() << "\n";
  std::string joined;
  for (std::size_t i = 0; i < s.warnings.size(); ++i) {
    if (i) joined += " | ";
    joined += s.warnings[i];
    f << "preflight_warning_" << i << "\t" << s.warnings[i] << "\n";
  }
  f << "galaxy_preflight_warnings\t" << joined << "\n";
}

}
