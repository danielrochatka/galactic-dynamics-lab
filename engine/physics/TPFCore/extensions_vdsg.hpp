#ifndef GALAXY_PHYSICS_TPFCORE_EXTENSIONS_VDSG_HPP
#define GALAXY_PHYSICS_TPFCORE_EXTENSIONS_VDSG_HPP

#include "../../types.hpp"
#include <string>
#include <vector>

namespace galaxy {
namespace tpfcore {

double vdsg_effective_coupling(double lambda0, double source_mass_kg, double baseline_mass_kg);

struct VdsgDiagnosticsSummary {
  bool mode_is_radial = false;
  unsigned long long pairs_evaluated = 0;
  unsigned long long beta_positive = 0;
  unsigned long long beta_negative = 0;
  unsigned long long beta_near_zero = 0;
  unsigned long long x_clamped = 0;
  double min_beta_rad = 0.0, max_beta_rad = 0.0, sum_abs_beta_rad = 0.0;
  double min_x_raw = 0.0, max_x_raw = 0.0;
  double min_x_clamped = 0.0, max_x_clamped = 0.0;
  double min_factor = 1.0, max_factor = 1.0, sum_factor = 0.0;
  double max_abs_delta_accel = 0.0, sum_abs_delta_accel = 0.0;
};

bool is_valid_vdsg_mode(const std::string& mode);

void accumulate_vdsg_velocity_modifier(const State& state,
                                       double bh_mass,
                                       double softening,
                                       bool star_star,
                                       double vdsg_coupling,
                                       double vdsg_mass_baseline_kg,
                                       const std::string& vdsg_mode,
                                       double mass_gate_m0_kg,
                                       double mass_gate_alpha,
                                       double x_clamp,
                                       bool weak_field_gate_enable,
                                       double weak_field_a0,
                                       double weak_field_power,
                                       double bounded_amplitude,
                                       std::vector<double>& ax,
                                       std::vector<double>& ay,
                                       VdsgDiagnosticsSummary* diagnostics_out = nullptr);

}  // namespace tpfcore
}  // namespace galaxy

#endif
