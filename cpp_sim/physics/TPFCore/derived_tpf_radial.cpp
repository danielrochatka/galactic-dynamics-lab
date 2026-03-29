/**
 * Derived TPF radial: hybrid bounce baryons + κ Poisson ledger.
 */

#include "derived_tpf_radial.hpp"
#include "../../config.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace galaxy {
namespace tpfcore {

namespace {

constexpr double kPi = 3.14159265358979323846;

Theta3D add_theta(const Theta3D& a, const Theta3D& b) {
  Theta3D o;
  o.xx = a.xx + b.xx;
  o.xy = a.xy + b.xy;
  o.xz = a.xz + b.xz;
  o.yy = a.yy + b.yy;
  o.yz = a.yz + b.yz;
  o.zz = a.zz + b.zz;
  return o;
}

}  // namespace

Theta3D evaluate_derived_theta(double mass_kg, double dx, double dy, double dz, double eps) {
  Theta3D t;
  double d_sq = dx * dx + dy * dy + dz * dz + eps * eps;
  double d = std::sqrt(d_sq);
  if (d < 1e-30) d = 1e-30;
  double d5 = d * d * d * d * d;
  double gm = TPF_G_SI * mass_kg;
  t.xx = gm * (3.0 * dx * dx - d_sq) / d5;
  t.yy = gm * (3.0 * dy * dy - d_sq) / d5;
  t.zz = gm * (3.0 * dz * dz - d_sq) / d5;
  t.xy = gm * (3.0 * dx * dy) / d5;
  t.xz = gm * (3.0 * dx * dz) / d5;
  t.yz = gm * (3.0 * dy * dz) / d5;
  return t;
}

double derived_invariant_I_contracted(const Theta3D& theta) {
  return theta.xx * theta.xx + theta.yy * theta.yy + theta.zz * theta.zz +
         2.0 * (theta.xy * theta.xy + theta.xz * theta.xz + theta.yz * theta.yz);
}

Theta3D sum_derived_theta_at_point(const State& state, double bh_mass, double px, double py, double pz,
                                   double eps) {
  Theta3D sum;
  sum.xx = sum.xy = sum.xz = sum.yy = sum.yz = sum.zz = 0.0;
  if (bh_mass > 0.0) {
    Theta3D bh = evaluate_derived_theta(bh_mass, px - 0.0, py - 0.0, pz - 0.0, eps);
    sum = add_theta(sum, bh);
  }
  const int n = state.n();
  for (int j = 0; j < n; ++j) {
    double mj = state.mass[j];
    if (mj <= 0.0) continue;
    Theta3D tj =
        evaluate_derived_theta(mj, px - state.x[j], py - state.y[j], pz - 0.0, eps);
    sum = add_theta(sum, tj);
  }
  return sum;
}

double enclosed_stellar_mass_cyl(const State& state, double r_cyl) {
  double M = 0.0;
  const int n = state.n();
  for (int j = 0; j < n; ++j) {
    double rj = std::hypot(state.x[j], state.y[j]);
    if (rj <= r_cyl) M += state.mass[j];
  }
  return M;
}

double get_tpf_mass_at_r(double M_total, double r) {
  if (M_total <= 0.0) return 0.0;
  double rr = std::abs(r);
  if (rr <= 0.0) return 0.0;

  double r_s = (2.0 * TPF_G_SI * M_total) / (c * c);
  double r2 = rr * rr;
  double r6 = r2 * r2 * r2;
  double rs2 = r_s * r_s;
  double rs6 = rs2 * rs2 * rs2;

  return M_total * (r6 / (r6 + rs6));
}

/**
 * Data pipeline (enforced):
 * 1. Here: R6 bounce on ρ ledger; cumulative_M_eff += shell_mass only if finite; profile.M_eff_enc[i] set each bin.
 * 2. radial_acceleration_scalar_derived: M_eff = profile.get_effective_mass_at(r_cyl) (interpolates M_eff_enc).
 */
TpfRadialGravityProfile build_tpf_gravity_profile(const State& state, double bh_mass, double max_radius,
                                                  int bins, double tpf_kappa, double eps) {
  static int warning_count = 0;
  static int stability_alert_prints = 0;
  const int MAX_WARNINGS = 5;

  TpfRadialGravityProfile profile;
  profile.bins = std::max(1, bins);
  profile.max_radius = std::max(max_radius, 1e-30);
  profile.delta_r = profile.max_radius / static_cast<double>(profile.bins);
  const size_t num_bins = static_cast<size_t>(profile.bins);
  profile.r_outer.resize(num_bins);
  profile.M_eff_enc.resize(num_bins);

  double r_s_bh = 0.0;
  if (bh_mass > 0.0) r_s_bh = (2.0 * TPF_G_SI * bh_mass) / (c * c);

  const double dr = profile.delta_r;
  double cumulative_M_eff = 0.0;

  for (size_t i = 0; i < num_bins; ++i) {
    double R_b = (static_cast<double>(i) + 0.5) * dr;
    double r = R_b;
    double px = R_b;
    double py = 0.0;
    double pz = 0.0;
    Theta3D theta_tot = sum_derived_theta_at_point(state, bh_mass, px, py, pz, eps);
    double I_total = derived_invariant_I_contracted(theta_tot);

    // --- HIGHER-ORDER SPATIAL MEMORY TERM (CORRECTED) ---
    // The delta-Gamma variation dictates that the central structural constraint
    // spreads macroscopically. In 3D space, this spatial memory dilutes
    // with the surface area (1/r^2), transitioning the invariant from
    // a local 1/r^6 drop-off to a galactic 1/r^2 halo.

    // Extract the base structural intensity at the core scale
    double theta_mag = std::sqrt(std::abs(I_total));

    // Project the constraint outward using a geometric 1/r^2 spatial memory
    // We use a macroscopic scale (L_MACRO) to balance the dimensions
    const double L_MACRO = 1.0e20;
    double spatial_memory_invariant = theta_mag * (L_MACRO / r) * (L_MACRO / r);

    // Inject the macroscopic spatial inertia into the total invariant
    I_total += spatial_memory_invariant;
    // ----------------------------------------------------

    double rho_raw = std::abs(I_total) * std::abs(tpf_kappa);
    if (!std::isfinite(rho_raw) || rho_raw > 1e300) {
      ++warning_count;
      if (stability_alert_prints < MAX_WARNINGS) {
        std::cerr << "![STABILITY ALERT]: Numerical singularity at r=" << r << std::endl;
        ++stability_alert_prints;
      }
      rho_raw = 0.0;
    }

    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double rs2 = r_s_bh * r_s_bh;
    double rs6 = rs2 * rs2 * rs2;
    double bounce_ratio = r6 / (r6 + rs6);

    double rho_eff = rho_raw * bounce_ratio;
    if (rho_eff < 0.0) rho_eff = 0.0;

    double shell_mass = 4.0 * kPi * r * r * rho_eff * dr;
    if (std::isfinite(shell_mass)) cumulative_M_eff += shell_mass;

    profile.r_outer[i] = (static_cast<double>(i) + 1.0) * dr;
    profile.M_eff_enc[i] = cumulative_M_eff;
  }

  /* Print once per process: κ and M_BH are arguments tpf_kappa, bh_mass (not Config). */
  static bool diagnostics_printed = false;
  if (!diagnostics_printed) {
    std::cout << "\n--- TPF MACROSCOPIC LEDGER (INITIAL CALIBRATION) ---" << std::endl;
    std::cout << "Coupling Constant (kappa): " << tpf_kappa << std::endl;
    std::cout << "Baryonic Black Hole Mass:  " << bh_mass << " kg" << std::endl;
    std::cout << "Total TPF Effective Mass:  " << cumulative_M_eff << " kg" << std::endl;
    if (bh_mass > 0.0)
      std::cout << "Ratio (TPF / BH):          " << (cumulative_M_eff / bh_mass) << std::endl;
    else
      std::cout << "Ratio (TPF / BH):          (undefined, bh_mass=0)" << std::endl;
    std::cout << "Numerical Stability Shunts triggered: " << warning_count << " bins." << std::endl;
    if (warning_count > MAX_WARNINGS) {
      std::cout << "![WARNING]: Core gravity may be under-resolved (singularity shunts exceeded "
                << MAX_WARNINGS << ")." << std::endl;
    }
    std::cout << "----------------------------------------------------\n" << std::endl;
    std::cout << std::flush;
    diagnostics_printed = true;
  }

  return profile;
}

double TpfRadialGravityProfile::M_eff_at_cylindrical_r(double r_cyl) const {
  if (r_cyl <= 0.0 || bins <= 0 || r_outer.empty()) return 0.0;
  if (r_cyl >= r_outer.back()) return M_eff_enc.back();
  for (size_t k = 0; k < r_outer.size(); ++k) {
    if (r_cyl <= r_outer[k]) {
      double r_lo = (k == 0) ? 0.0 : r_outer[k - 1];
      double M_lo = (k == 0) ? 0.0 : M_eff_enc[k - 1];
      double r_hi = r_outer[k];
      double M_hi = M_eff_enc[k];
      double t = (r_cyl - r_lo) / std::max(r_hi - r_lo, 1e-30);
      return M_lo + t * (M_hi - M_lo);
    }
  }
  return M_eff_enc.back();
}

/** Uses profile built by build_tpf_gravity_profile; M_eff is never ignored when this overload is used. */
double radial_acceleration_scalar_derived(const State& state, double bh_mass,
                                          const TpfRadialGravityProfile& profile, double r_cyl,
                                          double eps) {
  double r_soft_sq = r_cyl * r_cyl + eps * eps;
  if (r_soft_sq < 1e-60) r_soft_sq = 1e-60;
  double M_stars_enc = enclosed_stellar_mass_cyl(state, r_cyl);
  double M_baryon_bounced = get_tpf_mass_at_r(bh_mass + M_stars_enc, r_cyl);
  double M_eff = profile.get_effective_mass_at(r_cyl);
  return -TPF_G_SI * (M_baryon_bounced + M_eff) / r_soft_sq;
}

}  // namespace tpfcore
}  // namespace galaxy
