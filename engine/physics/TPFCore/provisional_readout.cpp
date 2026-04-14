/**
 * PROVISIONAL motion/readout layer for TPFCore.
 *
 * EXPLORATORY closures downstream of the ansatz (see readout_closure.hpp).
 * VDSG dynamics path bypasses readout for integrator ax, ay — see TPFCorePackage::compute_accelerations.
 */

#include "provisional_readout.hpp"

#include "readout_model_families.hpp"
#include "regime_diagnostics.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

namespace galaxy {
namespace tpfcore {

static double effective_eps(double source_softening, double global_softening) {
  return (source_softening > 0.0) ? source_softening : global_softening;
}

static bool is_negated_mode(const std::string& mode) {
  return mode == "tensor_radial_projection_negated";
}

void compute_provisional_readout_acceleration(const State& state,
                                               int i,
                                               double bh_mass,
                                               bool star_star,
                                               double softening,
                                               double source_softening,
                                               const std::string& readout_mode,
                                               double readout_scale,
                                               double theta_tt_scale,
                                               double theta_tr_scale,
                                               double& ax,
                                               double& ay,
                                               const DerivedTpfPoissonConfig* derived_poisson,
                                               const TpfRadialGravityProfile* derived_profile) {
  const double eps = effective_eps(source_softening, softening);

  if (is_derived_tpf_radial_readout_mode(readout_mode)) {
    static const DerivedTpfPoissonConfig kDefaultDerivedPoisson;
    const DerivedTpfPoissonConfig& dcfg = derived_poisson ? *derived_poisson : kDefaultDerivedPoisson;
    readout_models::apply_derived_radial_readout(state, i, bh_mass, eps, dcfg, derived_profile, readout_scale,
                                                 theta_tt_scale, theta_tr_scale, ax, ay, nullptr);
    return;
  }

  if (readout_mode == "experimental_radial_r_scaling") {
    readout_models::apply_experimental_radial_scaling_readout(state, i, bh_mass, star_star, eps, readout_scale, ax,
                                                              ay, nullptr);
    return;
  }

  if (readout_mode != "tensor_radial_projection" && readout_mode != "tensor_radial_projection_negated") {
    ax = ay = 0.0;
    return;
  }

  readout_models::apply_tensor_radial_projection_readout(state, i, bh_mass, star_star, eps, readout_scale,
                                                         is_negated_mode(readout_mode), ax, ay, nullptr);
}

void compute_provisional_readout_with_diagnostics(const State& state,
                                                   int i,
                                                   double bh_mass,
                                                   bool star_star,
                                                   double softening,
                                                   double source_softening,
                                                   const std::string& readout_mode,
                                                   double readout_scale,
                                                   double theta_tt_scale,
                                                   double theta_tr_scale,
                                                   double& ax,
                                                   double& ay,
                                                   ReadoutDiagnostics& diag,
                                                   const DerivedTpfPoissonConfig* derived_poisson,
                                                   const TpfRadialGravityProfile* derived_profile) {
  const double eps = effective_eps(source_softening, softening);

  if (is_derived_tpf_radial_readout_mode(readout_mode)) {
    static const DerivedTpfPoissonConfig kDefaultDerivedPoisson;
    const DerivedTpfPoissonConfig& dcfg = derived_poisson ? *derived_poisson : kDefaultDerivedPoisson;
    readout_models::apply_derived_radial_readout(state, i, bh_mass, eps, dcfg, derived_profile, readout_scale,
                                                 theta_tt_scale, theta_tr_scale, ax, ay, &diag);
    return;
  }

  if (readout_mode == "experimental_radial_r_scaling") {
    readout_models::apply_experimental_radial_scaling_readout(state, i, bh_mass, star_star, eps, readout_scale, ax,
                                                              ay, &diag);
    return;
  }

  if (readout_mode != "tensor_radial_projection" && readout_mode != "tensor_radial_projection_negated") {
    ax = ay = 0.0;
    diag.theta_xx = diag.theta_xy = diag.theta_yy = diag.theta_trace = diag.invariant_I = 0.0;
    return;
  }

  readout_models::apply_tensor_radial_projection_readout(state, i, bh_mass, star_star, eps, readout_scale,
                                                         is_negated_mode(readout_mode), ax, ay, &diag);
}

// --- Debug CSV: owned by readout module (column layout by mode) ---
void write_readout_debug_csv(const std::vector<Snapshot>& snapshots,
                             const std::string& output_dir,
                             double softening,
                             double bh_mass,
                             bool star_star,
                             double source_softening,
                             const std::string& readout_mode,
                             double readout_scale,
                             double theta_tt_scale,
                             double theta_tr_scale,
                             const DerivedTpfPoissonConfig& derived_poisson) {
  if (snapshots.empty()) return;

  const double eps = effective_eps(source_softening, softening);
  std::ofstream f(output_dir + "/tpf_readout_debug.csv");
  if (!f) return;

  const bool tr_style = is_derived_tpf_radial_readout_mode(readout_mode);
  const bool experimental_r_scaling = (readout_mode == "experimental_radial_r_scaling");
  const bool use_tr_style_columns = tr_style || experimental_r_scaling;
  /* residual_available=0 for multi-source (no analytic residual); residual_norm=0 when not available */
  if (use_tr_style_columns) {
    f << "time,particle,x,y,vx,vy,radius,theta_rr,theta_tt,theta_tr,theta_rr_plus_theta_tt,"
      << "provisional_radial_readout,provisional_tangential_readout,ax,ay,a_radial,a_inward,a_tangential,"
      << "theta_norm,invariant_I,regime,residual_available,residual_norm\n";
  } else {
    f << "time,particle,x,y,vx,vy,ax,ay,radius,radial_unit_x,radial_unit_y,"
      << "a_radial,a_inward,a_tangential,theta_xx,theta_xy,theta_yy,theta_trace,invariant_I,"
      << "theta_norm,regime,residual_available,residual_norm\n";
  }

  for (const auto& snap : snapshots) {
    const State& s = snap.state;
    const double t = snap.time;
    TpfRadialGravityProfile derived_prof;
    const TpfRadialGravityProfile* derived_prof_ptr = nullptr;
    if (is_derived_tpf_radial_readout_mode(readout_mode)) {
      derived_prof = build_tpf_gravity_profile(s, bh_mass, derived_poisson, eps);
      derived_prof_ptr = &derived_prof;
    }
    for (int i = 0; i < s.n(); ++i) {
      double x = s.x[i], y = s.y[i];
      double vx = s.vx[i], vy = s.vy[i];

      double ax = 0, ay = 0;
      ReadoutDiagnostics diag;
      compute_provisional_readout_with_diagnostics(
          s, i, bh_mass, star_star, softening, source_softening,
          readout_mode, readout_scale,
          theta_tt_scale, theta_tr_scale, ax, ay, diag,
          is_derived_tpf_radial_readout_mode(readout_mode) ? &derived_poisson : nullptr,
          derived_prof_ptr);

      double r2 = x * x + y * y + eps * eps;
      double r = std::sqrt(r2);
      double radial_unit_x = (r > 1e-30) ? (x / r) : 1.0;
      double radial_unit_y = (r > 1e-30) ? (y / r) : 0.0;

      double a_radial;
      double a_inward;
      double tangential_unit_x = -radial_unit_y;
      double tangential_unit_y = radial_unit_x;
      double a_tangential;

      /* Hybrid ledger: CSV must match radial_acceleration_scalar_derived (M_baryon_bounced + M_eff). */
      if (tr_style && derived_prof_ptr != nullptr) {
        const double r_cyl = std::hypot(x, y);
        double r_soft_sq = r_cyl * r_cyl + eps * eps;
        if (r_soft_sq < 1e-60) r_soft_sq = 1e-60;
        const double M_stars_enc = enclosed_stellar_mass_cyl(s, r_cyl);
        const double M_baryon_bounced = get_tpf_mass_at_r(bh_mass + M_stars_enc, r_cyl);
        const double M_eff = std::abs(derived_prof.get_effective_mass_at(r_cyl));
        const double a_s = -TPF_G_SI * (M_baryon_bounced + M_eff) / r_soft_sq;
        if (r < 1e-30) {
          ax = ay = 0.0;
        } else {
          ax = a_s * (x / r);
          ay = a_s * (y / r);
        }
        a_radial = ax * radial_unit_x + ay * radial_unit_y;
        a_inward = -a_radial;
        a_tangential = ax * tangential_unit_x + ay * tangential_unit_y;
        diag.provisional_radial_readout = a_s;
      } else {
        a_radial = ax * radial_unit_x + ay * radial_unit_y;
        a_inward = -a_radial;
        a_tangential = ax * tangential_unit_x + ay * tangential_unit_y;
      }

      const std::string regime_out =
          (tr_style && !diag.regime.empty()) ? diag.regime
                                           : std::string(regime_label_from_theta_norm(diag.theta_norm));
      if (use_tr_style_columns) {
        f << std::scientific << t << "," << i << "," << x << "," << y << "," << vx << "," << vy << ","
          << r << "," << diag.theta_rr << "," << diag.theta_tt << "," << diag.theta_tr << ","
          << diag.theta_rr_plus_theta_tt << "," << diag.provisional_radial_readout << ","
          << diag.provisional_tangential_readout << "," << ax << "," << ay << ","
          << a_radial << "," << a_inward << "," << a_tangential << ","
          << diag.theta_norm << "," << diag.invariant_I << "," << regime_out << ",0,0\n";
      } else {
        f << std::scientific << t << "," << i << "," << x << "," << y << "," << vx << "," << vy << ","
          << ax << "," << ay << "," << r << "," << radial_unit_x << "," << radial_unit_y << ","
          << a_radial << "," << a_inward << "," << a_tangential << ","
          << diag.theta_xx << "," << diag.theta_xy << "," << diag.theta_yy << ","
          << diag.theta_trace << "," << diag.invariant_I << ","
          << diag.theta_norm << "," << regime_out << ",0,0\n";
      }
    }
  }
}

}  // namespace tpfcore
}  // namespace galaxy
