/**
 * PROVISIONAL motion/readout layer for TPFCore.
 *
 * EXPLORATORY: This is NOT the full derived TPF dynamics.
 * Closures are downstream of the ansatz (see readout_closure.hpp).
 */

#include "provisional_readout.hpp"
#include "derived_tpf_radial.hpp"
#include "field_evaluation.hpp"
#include "readout_closure.hpp"
#include "regime_diagnostics.hpp"
#include "source_ansatz.hpp"
#include "../../types.hpp"
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

// --- Closure: tensor_radial (per-source Theta·r_hat, superposed; optional negated) ---
static void apply_tensor_radial_closure(const State& state,
                                 int i,
                                 double bh_mass,
                                 bool star_star,
                                 double eps,
                                 double readout_scale,
                                 double& ax,
                                 double& ay,
                                 Theta3D* theta_sum,
                                 bool* has_theta_sum) {
  ax = 0.0;
  ay = 0.0;
  if (theta_sum) {
    theta_sum->xx = theta_sum->xy = theta_sum->xz = theta_sum->yy = theta_sum->yz = theta_sum->zz = 0.0;
    *has_theta_sum = false;
  }

  const double x = state.x[i];
  const double y = state.y[i];
  const int n = state.n();

  auto add_contribution = [&](double xs, double ys, double m) {
    if (m <= 0.0) return;
    double dx = x - xs;
    double dy = y - ys;
    double r2 = dx * dx + dy * dy + eps * eps;
    double r = std::sqrt(r2);
    if (r < 1e-30) return;
    double rx = dx / r;
    double ry = dy / r;

    FieldAtPoint field = evaluate_provisional_field_single_source(xs, ys, m, x, y, eps);
    const Theta3D& theta = field.theta;
    double ax_contrib = theta.xx * rx + theta.xy * ry;
    double ay_contrib = theta.xy * rx + theta.yy * ry;

    ax += readout_scale * ax_contrib;
    ay += readout_scale * ay_contrib;

    if (theta_sum) {
      theta_sum->xx += theta.xx;
      theta_sum->xy += theta.xy;
      theta_sum->xz += theta.xz;
      theta_sum->yy += theta.yy;
      theta_sum->yz += theta.yz;
      theta_sum->zz += theta.zz;
      *has_theta_sum = true;
    }
  };

  if (bh_mass > 0.0) add_contribution(0.0, 0.0, bh_mass);
  if (star_star) {
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      add_contribution(state.x[j], state.y[j], state.mass[j]);
    }
  }
}

// --- Superposition: field at particle from all sources (used by tr_coherence closure) ---
static void compute_theta_sum(const State& state,
                              int i,
                              double bh_mass,
                              bool star_star,
                              double eps,
                              double /*x*/, double /*y*/,
                              Theta3D& theta_sum) {
  FieldAtPoint field = evaluate_provisional_field_multi_source(state, i, bh_mass, star_star, eps);
  theta_sum = field.theta;
}

// --- Closure: hybrid derived TPF radial (bounce baryons + κ ledger + Hessian superposition). ---
static void apply_derived_tpf_radial_readout_closure(const State& state,
                                                     int i,
                                                     double bh_mass,
                                                     double eps,
                                                     const DerivedTpfPoissonConfig& dcfg,
                                                     const TpfRadialGravityProfile* profile_in,
                                                     double readout_scale,
                                                     double theta_tt_scale,
                                                     double theta_tr_scale,
                                                     double& ax,
                                                     double& ay,
                                                     ReadoutDiagnostics* diag) {
  const double x = state.x[i];
  const double y = state.y[i];
  double r2 = x * x + y * y + eps * eps;
  double r = std::sqrt(r2);
  if (r < 1e-30) {
    ax = ay = 0.0;
    if (diag) {
      diag->theta_rr = diag->theta_tt = diag->theta_tr = diag->theta_rr_plus_theta_tt =
        diag->provisional_radial_readout = diag->provisional_tangential_readout = 0.0;
      diag->regime.clear();
    }
    return;
  }

  TpfRadialGravityProfile profile_storage;
  const TpfRadialGravityProfile* profile = profile_in;
  if (!profile) {
    profile_storage = build_tpf_gravity_profile(state, bh_mass, dcfg, eps);
    profile = &profile_storage;
  }

  double r_cyl = std::hypot(x, y);
  double a_s = radial_acceleration_scalar_derived(state, bh_mass, *profile, r_cyl, eps);
  ax = a_s * (x / r);
  ay = a_s * (y / r);

  Theta3D theta_sum = sum_derived_theta_at_point(state, bh_mass, x, y, 0.0, eps);
  double rx = x / r;
  double ry = y / r;
  double tx = -ry;
  double ty = rx;
  double theta_rr = rx * rx * theta_sum.xx + 2.0 * rx * ry * theta_sum.xy + ry * ry * theta_sum.yy;
  double theta_tt = theta_tt_scale * (-theta_rr);
  double theta_rr_plus_theta_tt = theta_rr + theta_tt;
  double theta_tr =
      tx * (theta_sum.xx * rx + theta_sum.xy * ry) + ty * (theta_sum.xy * rx + theta_sum.yy * ry);
  double provisional_tangential = readout_scale * theta_tr_scale * theta_tr;

  if (diag) {
    diag->theta_rr = theta_rr;
    diag->theta_tt = theta_tt;
    diag->theta_tr = theta_tr;
    diag->theta_rr_plus_theta_tt = theta_rr_plus_theta_tt;
    diag->provisional_radial_readout = a_s;
    diag->provisional_tangential_readout = provisional_tangential;
    diag->theta_xx = theta_sum.xx;
    diag->theta_xy = theta_sum.xy;
    diag->theta_yy = theta_sum.yy;
    diag->theta_trace = theta_sum.trace();
    diag->invariant_I = derived_invariant_I_contracted(theta_sum);
    diag->theta_norm = theta_frobenius_norm(theta_sum);
    diag->regime = "derived-tpf-radial";
  }
}

// --- EXPERIMENTAL closure: radial only, with r-scaling (Theta ~ 1/r^3 -> a_rad ~ 1/r^2). ---
static void apply_experimental_radial_r_scaling_closure(const State& state,
                                                        int i,
                                                        double bh_mass,
                                                        bool star_star,
                                                        double eps,
                                                        double readout_scale,
                                                        double& ax,
                                                        double& ay,
                                                        ReadoutDiagnostics* diag) {
  const double x = state.x[i];
  const double y = state.y[i];
  double r2 = x * x + y * y + eps * eps;
  double r = std::sqrt(r2);
  if (r < 1e-30) {
    ax = ay = 0.0;
    if (diag) diag->theta_rr = diag->theta_tt = diag->theta_tr = diag->theta_rr_plus_theta_tt =
      diag->provisional_radial_readout = diag->provisional_tangential_readout = 0.0;
    return;
  }
  double rx = x / r;
  double ry = y / r;
  double tx = -ry;
  double ty = rx;

  Theta3D theta_sum;
  compute_theta_sum(state, i, bh_mass, star_star, eps, x, y, theta_sum);

  double theta_rr = rx * rx * theta_sum.xx + 2.0 * rx * ry * theta_sum.xy + ry * ry * theta_sum.yy;
  double theta_tr = tx * (theta_sum.xx * rx + theta_sum.xy * ry) + ty * (theta_sum.xy * rx + theta_sum.yy * ry);

  // Inward radial: magnitude = readout_scale * (-theta_rr) * r (so effective ~ 1/r^2 when theta_rr ~ 1/r^3).
  // Apply as acceleration = -magnitude * r_hat (inward). r_hat = (rx, ry) points outward; Newtonian a is inward.
  double provisional_radial = readout_scale * (-theta_rr) * r;
  double provisional_tangential = 0.0;

  ax = -provisional_radial * rx;
  ay = -provisional_radial * ry;

  if (diag) {
    diag->theta_rr = theta_rr;
    diag->theta_tt = 0.0;
    diag->theta_tr = theta_tr;
    diag->theta_rr_plus_theta_tt = theta_rr;
    diag->provisional_radial_readout = provisional_radial;
    diag->provisional_tangential_readout = provisional_tangential;
    diag->theta_xx = theta_sum.xx;
    diag->theta_xy = theta_sum.xy;
    diag->theta_yy = theta_sum.yy;
    diag->theta_trace = theta_sum.trace();
    diag->invariant_I = compute_invariant_I(theta_sum);
    diag->theta_norm = theta_frobenius_norm(theta_sum);
  }
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
    apply_derived_tpf_radial_readout_closure(state, i, bh_mass, eps, dcfg, derived_profile,
                                             readout_scale, theta_tt_scale, theta_tr_scale, ax, ay,
                                             nullptr);
    return;
  }

  if (readout_mode == "experimental_radial_r_scaling") {
    apply_experimental_radial_r_scaling_closure(state, i, bh_mass, star_star, eps,
                                                readout_scale, ax, ay, nullptr);
    return;
  }

  if (readout_mode != "tensor_radial_projection" && readout_mode != "tensor_radial_projection_negated")
    return;

  apply_tensor_radial_closure(state, i, bh_mass, star_star, eps, readout_scale, ax, ay, nullptr,
                              nullptr);

  if (is_negated_mode(readout_mode)) {
    ax = -ax;
    ay = -ay;
  }
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
    apply_derived_tpf_radial_readout_closure(state, i, bh_mass, eps, dcfg, derived_profile,
                                             readout_scale, theta_tt_scale, theta_tr_scale, ax, ay,
                                             &diag);
    diag.ax = ax;
    diag.ay = ay;
    return;
  }

  if (readout_mode == "experimental_radial_r_scaling") {
    apply_experimental_radial_r_scaling_closure(state, i, bh_mass, star_star, eps,
                                                readout_scale, ax, ay, &diag);
    diag.ax = ax;
    diag.ay = ay;
    return;
  }

  if (readout_mode != "tensor_radial_projection" && readout_mode != "tensor_radial_projection_negated") {
    ax = ay = 0.0;
    diag.theta_xx = diag.theta_xy = diag.theta_yy = diag.theta_trace = diag.invariant_I = 0.0;
    return;
  }

  Theta3D theta_sum;
  bool has_theta;
  apply_tensor_radial_closure(state, i, bh_mass, star_star, eps, readout_scale, ax, ay, &theta_sum, &has_theta);

  if (is_negated_mode(readout_mode)) {
    ax = -ax;
    ay = -ay;
  }

  diag.ax = ax;
  diag.ay = ay;
  if (has_theta) {
    diag.theta_xx = theta_sum.xx;
    diag.theta_xy = theta_sum.xy;
    diag.theta_yy = theta_sum.yy;
    diag.theta_trace = theta_sum.trace();
    diag.invariant_I = compute_invariant_I(theta_sum);
    diag.theta_norm = theta_frobenius_norm(theta_sum);
  } else {
    diag.theta_xx = diag.theta_xy = diag.theta_yy = diag.theta_trace = diag.invariant_I = diag.theta_norm = 0.0;
  }
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

      double a_radial = ax * radial_unit_x + ay * radial_unit_y;
      double a_inward = -a_radial;
      double tangential_unit_x = -radial_unit_y;
      double tangential_unit_y = radial_unit_x;
      double a_tangential = ax * tangential_unit_x + ay * tangential_unit_y;

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
