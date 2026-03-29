#ifndef GALAXY_PHYSICS_TPFCORE_DERIVED_TPF_RADIAL_HPP
#define GALAXY_PHYSICS_TPFCORE_DERIVED_TPF_RADIAL_HPP

/**
 * Hybrid TPF radial: Sec. IX bounce metric on baryons + macroscopic Poisson ledger (κ * I) with bounce
 * regularization on ρ_eff, spherical shells → M_eff_enc(r).
 */

#include "source_ansatz.hpp"
#include "../../types.hpp"
#include <string>
#include <vector>

namespace galaxy {
namespace tpfcore {

/** SI Newton constant (matches manuscript numerical value). */
constexpr double TPF_G_SI = 6.6743e-11;

/** Config for derived radial Poisson grid (κ, bins, radial extent). */
struct DerivedTpfPoissonConfig {
  /** Coupling κ: ρ_raw = κ * I (SI). */
  double kappa = 1.0e32;
  int bins = 100;
  /** If <= 0, use galaxy_radius. */
  double max_radius = 0.0;
  double galaxy_radius = 50.0;

  double max_radius_resolved() const {
    return max_radius > 0.0 ? max_radius : galaxy_radius;
  }
};

Theta3D evaluate_derived_theta(double mass_kg, double dx, double dy, double dz, double eps);

/** Frobenius sum Θ_ij Θ_ij only; not manuscript I = Θ_μν Θ^μν − λ Θ². κ-ledger uses compute_invariant_I. */
double derived_invariant_I_contracted(const Theta3D& theta);

Theta3D sum_derived_theta_at_point(const State& state, double bh_mass, double px, double py, double pz,
                                   double eps);

/** R6 bounce enclosed mass (SI): m(r) = M_total * r^6 / (r^6 + r_s^6), r_s = 2 G M_total / c^2. */
double get_tpf_mass_at_r(double M_total, double r);

struct TpfRadialGravityProfile {
  int bins = 0;
  double max_radius = 0.0;
  double delta_r = 0.0;
  std::vector<double> r_outer;
  /** Cumulative effective mass from Poisson shells (κ I, bounce-regularized ρ_eff). */
  std::vector<double> M_eff_enc;

  /** Cumulative M_eff enclosed within cylindrical radius r. */
  double M_eff_at_cylindrical_r(double r_cyl) const;
  /** Alias for the same interpolation (hybrid ledger). */
  double get_effective_mass_at(double r_cyl) const { return M_eff_at_cylindrical_r(r_cyl); }
};

TpfRadialGravityProfile build_tpf_gravity_profile(const State& state, double bh_mass, double max_radius,
                                                  int bins, double tpf_kappa, double eps);

inline TpfRadialGravityProfile build_tpf_gravity_profile(const State& state, double bh_mass,
                                                         const DerivedTpfPoissonConfig& cfg, double eps) {
  return build_tpf_gravity_profile(state, bh_mass, cfg.max_radius_resolved(), cfg.bins, cfg.kappa, eps);
}

double enclosed_stellar_mass_cyl(const State& state, double r_cyl);

/**
 * a_r = -G * (M_baryon_bounced + M_eff) / r_soft^2,
 * M_baryon_bounced = get_tpf_mass_at_r(M_BH + M_stars_enc, r_cyl), M_eff from profile.
 */
double radial_acceleration_scalar_derived(const State& state, double bh_mass,
                                          const TpfRadialGravityProfile& profile, double r_cyl,
                                          double eps);

inline bool is_derived_tpf_radial_readout_mode(const std::string& mode) {
  return mode == "tr_coherence_readout" || mode == "derived_tpf_radial_readout";
}

}  // namespace tpfcore
}  // namespace galaxy

#endif
