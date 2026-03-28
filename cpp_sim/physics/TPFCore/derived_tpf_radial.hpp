#ifndef GALAXY_PHYSICS_TPFCORE_DERIVED_TPF_RADIAL_HPP
#define GALAXY_PHYSICS_TPFCORE_DERIVED_TPF_RADIAL_HPP

/**
 * Derived TPF radial gravity: manuscript Hessian superposition + TPF bounce closure (Sec. IX, Eq. 31).
 * Enclosed mass m(r) = M * r^3 / (r^3 + r_s^3) with r_s = 2 G M / c^2 (SI); radial acceleration uses m(r).
 * Optional 1D radial grid (legacy diagnostics): cumulative M_eff is zero; dynamics do not use Poisson shells.
 */

#include "source_ansatz.hpp"
#include "../../types.hpp"
#include <string>
#include <vector>

namespace galaxy {
namespace tpfcore {

/** SI Newton constant (matches manuscript numerical value). */
constexpr double TPF_G_SI = 6.6743e-11;

/** Config for derived radial grid (diagnostics only; bounce law uses no tunable couplings). */
struct DerivedTpfPoissonConfig {
  int bins = 100;
  /** If <= 0, use galaxy_radius. */
  double max_radius = 0.0;
  double galaxy_radius = 50.0;

  double max_radius_resolved() const {
    return max_radius > 0.0 ? max_radius : galaxy_radius;
  }
};

/**
 * Hessian Theta_ab = G M (3 d_a d_b - delta_ab d^2) / d^5 with d = displacement from source to field point.
 * Softening: d^2 = dx^2+dy^2+dz^2 + eps^2 (isotropic on distance magnitude only).
 */
Theta3D evaluate_derived_theta(double mass_kg, double dx, double dy, double dz, double eps);

/** Full contraction I = Theta_ab Theta^ab (Euclidean Frobenius squared). */
double derived_invariant_I_contracted(const Theta3D& theta);

/** Linear superposition of derived Theta at field point (px, py, pz) from BH + all disk stars (z=0). */
Theta3D sum_derived_theta_at_point(const State& state, double bh_mass, double px, double py, double pz,
                                   double eps);

/** TPF bounce enclosed mass (Eq. 31, SI): m(r) = M_total * r^3 / (r^3 + r_s^3), r_s = 2 G M_total / c^2. */
double get_tpf_mass_at_r(double M_total, double r);

/** 1D radial grid for optional diagnostics; M_eff_enc is identically zero. */
struct TpfRadialGravityProfile {
  int bins = 0;
  double max_radius = 0.0;
  double delta_r = 0.0;
  /** Outer radius of shell k is (k+1) * delta_r. */
  std::vector<double> r_outer;
  /** Legacy cumulative effective mass (always zero with bounce closure). */
  std::vector<double> M_eff_enc;

  double M_eff_at_cylindrical_r(double r_cyl) const;
};

TpfRadialGravityProfile build_tpf_gravity_profile(const State& state, double bh_mass, double max_radius,
                                                  int bins, double eps);

inline TpfRadialGravityProfile build_tpf_gravity_profile(const State& state, double bh_mass,
                                                         const DerivedTpfPoissonConfig& cfg, double eps) {
  return build_tpf_gravity_profile(state, bh_mass, cfg.max_radius_resolved(), cfg.bins, eps);
}

/** Stellar mass with sqrt(x^2+y^2) <= r_cyl. */
double enclosed_stellar_mass_cyl(const State& state, double r_cyl);

/**
 * a_r (outward radial scalar) = -G * m(r_cyl) / r_soft^2 with m from get_tpf_mass_at_r(M_enc, r_cyl),
 * M_enc = M_BH + enclosed stars, r_soft^2 = r_cyl^2 + eps^2.
 */
double radial_acceleration_scalar_derived(const State& state, double bh_mass, double r_cyl, double eps);

inline bool is_derived_tpf_radial_readout_mode(const std::string& mode) {
  return mode == "tr_coherence_readout" || mode == "derived_tpf_radial_readout";
}

}  // namespace tpfcore
}  // namespace galaxy

#endif
