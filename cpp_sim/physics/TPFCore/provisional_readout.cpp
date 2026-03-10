/**
 * PROVISIONAL motion/readout layer for TPFCore.
 *
 * EXPLORATORY: This is NOT the full derived TPF dynamics.
 * Closures are downstream of the ansatz (see readout_closure.hpp).
 */

#include "provisional_readout.hpp"
#include "field_evaluation.hpp"
#include "readout_closure.hpp"
#include "regime_diagnostics.hpp"
#include "source_ansatz.hpp"
#include "../../types.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>

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
                                 double c,
                                 double readout_scale,
                                 double& ax,
                                 double& ay,
                                 Theta2D* theta_sum,
                                 bool* has_theta_sum) {
  ax = 0.0;
  ay = 0.0;
  if (theta_sum) {
    theta_sum->xx = theta_sum->xy = theta_sum->yy = 0.0;
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

    FieldAtPoint field = evaluate_provisional_field_single_source(xs, ys, m, x, y, eps, c);
    const Theta2D& theta = field.theta;
    double ax_contrib = theta.xx * rx + theta.xy * ry;
    double ay_contrib = theta.xy * rx + theta.yy * ry;

    ax += readout_scale * ax_contrib;
    ay += readout_scale * ay_contrib;

    if (theta_sum) {
      theta_sum->xx += theta.xx;
      theta_sum->xy += theta.xy;
      theta_sum->yy += theta.yy;
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
                              double c,
                              double /*x*/, double /*y*/,
                              Theta2D& theta_sum) {
  FieldAtPoint field = evaluate_provisional_field_multi_source(state, i, bh_mass, star_star, eps, c);
  theta_sum = field.theta;
}

// --- Closure: tr_coherence (superposed Theta -> Theta_rr/Theta_tt/Theta_tr formula). EXPLORATORY. ---
static void apply_tr_coherence_closure(const State& state,
                                         int i,
                                         double bh_mass,
                                         bool star_star,
                                         double eps,
                                         double c,
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
    if (diag) diag->theta_rr = diag->theta_tt = diag->theta_tr = diag->theta_rr_plus_theta_tt =
      diag->provisional_radial_readout = diag->provisional_tangential_readout = 0.0;
    return;
  }
  double rx = x / r;
  double ry = y / r;
  double tx = -ry;
  double ty = rx;

  Theta2D theta_sum;
  compute_theta_sum(state, i, bh_mass, star_star, eps, c, x, y, theta_sum);

  double theta_rr = rx * rx * theta_sum.xx + 2.0 * rx * ry * theta_sum.xy + ry * ry * theta_sum.yy;
  double theta_tt = theta_tt_scale * (-theta_rr);
  double theta_rr_plus_theta_tt = theta_rr + theta_tt;
  double theta_tr = tx * (theta_sum.xx * rx + theta_sum.xy * ry) + ty * (theta_sum.xy * rx + theta_sum.yy * ry);

  double provisional_radial = readout_scale * theta_rr_plus_theta_tt;
  double provisional_tangential = readout_scale * theta_tr_scale * theta_tr;

  ax = provisional_radial * rx + provisional_tangential * tx;
  ay = provisional_radial * ry + provisional_tangential * ty;

  if (diag) {
    diag->theta_rr = theta_rr;
    diag->theta_tt = theta_tt;
    diag->theta_tr = theta_tr;
    diag->theta_rr_plus_theta_tt = theta_rr_plus_theta_tt;
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
                                               double c,
                                               const std::string& readout_mode,
                                               double readout_scale,
                                               double theta_tt_scale,
                                               double theta_tr_scale,
                                               double& ax,
                                               double& ay) {
  const double eps = effective_eps(source_softening, softening);

  if (readout_mode == "tr_coherence_readout") {
    apply_tr_coherence_closure(state, i, bh_mass, star_star, eps, c,
                               readout_scale, theta_tt_scale, theta_tr_scale, ax, ay, nullptr);
    return;
  }

  if (readout_mode != "tensor_radial_projection" && readout_mode != "tensor_radial_projection_negated")
    return;

  apply_tensor_radial_closure(state, i, bh_mass, star_star, eps, c, readout_scale, ax, ay, nullptr, nullptr);

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
                                                   double c,
                                                   const std::string& readout_mode,
                                                   double readout_scale,
                                                   double theta_tt_scale,
                                                   double theta_tr_scale,
                                                   double& ax,
                                                   double& ay,
                                                   ReadoutDiagnostics& diag) {
  const double eps = effective_eps(source_softening, softening);

  if (readout_mode == "tr_coherence_readout") {
    apply_tr_coherence_closure(state, i, bh_mass, star_star, eps, c,
                               readout_scale, theta_tt_scale, theta_tr_scale, ax, ay, &diag);
    diag.ax = ax;
    diag.ay = ay;
    return;
  }

  if (readout_mode != "tensor_radial_projection" && readout_mode != "tensor_radial_projection_negated") {
    ax = ay = 0.0;
    diag.theta_xx = diag.theta_xy = diag.theta_yy = diag.theta_trace = diag.invariant_I = 0.0;
    return;
  }

  Theta2D theta_sum;
  bool has_theta;
  apply_tensor_radial_closure(state, i, bh_mass, star_star, eps, c, readout_scale, ax, ay, &theta_sum, &has_theta);

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
                             double isotropic_c,
                             const std::string& readout_mode,
                             double readout_scale,
                             double theta_tt_scale,
                             double theta_tr_scale) {
  if (snapshots.empty()) return;

  const double eps = effective_eps(source_softening, softening);
  std::ofstream f(output_dir + "/tpf_readout_debug.csv");
  if (!f) return;

  const bool tr_coherence = (readout_mode == "tr_coherence_readout");
  /* residual_available=0 for multi-source (no analytic residual); residual_norm=0 when not available */
  if (tr_coherence) {
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
    for (int i = 0; i < s.n(); ++i) {
      double x = s.x[i], y = s.y[i];
      double vx = s.vx[i], vy = s.vy[i];

      double ax = 0, ay = 0;
      ReadoutDiagnostics diag;
      compute_provisional_readout_with_diagnostics(
          s, i, bh_mass, star_star, softening, source_softening,
          isotropic_c, readout_mode, readout_scale,
          theta_tt_scale, theta_tr_scale, ax, ay, diag);

      double r2 = x * x + y * y + eps * eps;
      double r = std::sqrt(r2);
      double radial_unit_x = (r > 1e-30) ? (x / r) : 1.0;
      double radial_unit_y = (r > 1e-30) ? (y / r) : 0.0;

      double a_radial = ax * radial_unit_x + ay * radial_unit_y;
      double a_inward = -a_radial;
      double tangential_unit_x = -radial_unit_y;
      double tangential_unit_y = radial_unit_x;
      double a_tangential = ax * tangential_unit_x + ay * tangential_unit_y;

      const char* regime = regime_label_from_theta_norm(diag.theta_norm);
      if (tr_coherence) {
        f << std::scientific << t << "," << i << "," << x << "," << y << "," << vx << "," << vy << ","
          << r << "," << diag.theta_rr << "," << diag.theta_tt << "," << diag.theta_tr << ","
          << diag.theta_rr_plus_theta_tt << "," << diag.provisional_radial_readout << ","
          << diag.provisional_tangential_readout << "," << ax << "," << ay << ","
          << a_radial << "," << a_inward << "," << a_tangential << ","
          << diag.theta_norm << "," << diag.invariant_I << "," << regime << ",0,0\n";
      } else {
        f << std::scientific << t << "," << i << "," << x << "," << y << "," << vx << "," << vy << ","
          << ax << "," << ay << "," << r << "," << radial_unit_x << "," << radial_unit_y << ","
          << a_radial << "," << a_inward << "," << a_tangential << ","
          << diag.theta_xx << "," << diag.theta_xy << "," << diag.theta_yy << ","
          << diag.theta_trace << "," << diag.invariant_I << ","
          << diag.theta_norm << "," << regime << ",0,0\n";
      }
    }
  }
}

}  // namespace tpfcore
}  // namespace galaxy
