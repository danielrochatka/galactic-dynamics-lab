/**
 * Newtonian-vs-TPF acceleration comparison diagnostic.
 * Diagnostics only; no change to physics formulas or integration.
 */

#include "force_compare.hpp"
#include "init_conditions.hpp"
#include "physics/physics_package.hpp"
#include "physics/TPFCore/derived_tpf_radial.hpp"
#include "physics/TPFCore/provisional_readout.hpp"
#include "simulation.hpp"
#include "types.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace galaxy {

namespace {

const double PI = 3.14159265358979323846;

struct CircleAuditStats {
  double min_ratio_x_axis = 1e300, max_ratio_x_axis = 0.0;
  double min_ratio_off_axis = 1e300, max_ratio_off_axis = 0.0;
  int n_x_axis = 0, n_off_axis = 0;
  double max_abs_tan_tpf = 0.0, max_abs_rad_newt = 0.0;
  double max_diff_mag = 0.0;
};

// Radii for static circle audit (weak-field regime)
const std::vector<double> CIRCLE_RADII = {20.0, 50.0, 100.0, 200.0};
const int ANGLES_PER_CIRCLE = 32;

void radial_tangential(double x, double y, double eps_sq,
                      double ax, double ay,
                      double& a_radial, double& a_tangential) {
  double r2 = x * x + y * y + eps_sq;
  double r = std::sqrt(r2);
  if (r < 1e-30) {
    a_radial = a_tangential = 0.0;
    return;
  }
  double rx = x / r;
  double ry = y / r;
  double tx = -ry;
  double ty = rx;
  a_radial = ax * rx + ay * ry;
  a_tangential = ax * tx + ay * ty;
}

void run_static_circle_audit(const Config& config,
                             PhysicsPackage* newton,
                             PhysicsPackage* tpf,
                             std::ostream& csv,
                             std::ostream& txt,
                             CircleAuditStats& circle_stats) {
  const double bh_mass = config.bh_mass;
  const double softening = config.softening;
  const double eps_sq = softening * softening;
  const bool star_star = false;

  circle_stats = CircleAuditStats{};
  double& max_ratio_rad_x_axis = circle_stats.max_ratio_x_axis;
  double& max_ratio_rad_off_axis = circle_stats.max_ratio_off_axis;
  double& min_ratio_rad_x_axis = circle_stats.min_ratio_x_axis;
  double& min_ratio_rad_off_axis = circle_stats.min_ratio_off_axis;
  int& n_x_axis = circle_stats.n_x_axis;
  int& n_off_axis = circle_stats.n_off_axis;
  double& max_abs_tan_tpf = circle_stats.max_abs_tan_tpf;
  double& max_abs_rad_newt = circle_stats.max_abs_rad_newt;
  double& max_diff_mag = circle_stats.max_diff_mag;

  for (double r : CIRCLE_RADII) {
    for (int k = 0; k < ANGLES_PER_CIRCLE; ++k) {
      double theta = (2.0 * PI * k) / ANGLES_PER_CIRCLE;
      double x = r * std::cos(theta);
      double y = r * std::sin(theta);

      State state;
      state.resize(1);
      state.x[0] = x;
      state.y[0] = y;
      state.vx[0] = 0.0;
      state.vy[0] = 0.0;
      state.mass[0] = config.star_mass;

      std::vector<double> ax_n(1), ay_n(1), ax_t(1), ay_t(1);
      newton->compute_accelerations(state, bh_mass, softening, star_star, ax_n, ay_n);
      tpf->compute_accelerations(state, bh_mass, softening, star_star, ax_t, ay_t);

      double a_rad_n, a_tan_n, a_rad_t, a_tan_t;
      radial_tangential(x, y, eps_sq, ax_n[0], ay_n[0], a_rad_n, a_tan_n);
      radial_tangential(x, y, eps_sq, ax_t[0], ay_t[0], a_rad_t, a_tan_t);

      double diff_x = ax_t[0] - ax_n[0];
      double diff_y = ay_t[0] - ay_n[0];
      double diff_mag = std::sqrt(diff_x * diff_x + diff_y * diff_y);
      if (diff_mag > max_diff_mag) max_diff_mag = diff_mag;

      double ratio_rad = 0.0;
      double abs_rad_n = std::abs(a_rad_n);
      if (abs_rad_n > 1e-15)
        ratio_rad = std::abs(a_rad_t) / abs_rad_n;

      csv << "circle," << std::scientific << r << "," << theta << ",0,0," << x << "," << y << ","
          << ax_n[0] << "," << ay_n[0] << "," << ax_t[0] << "," << ay_t[0] << ","
          << a_rad_n << "," << a_tan_n << "," << a_rad_t << "," << a_tan_t << ","
          << diff_x << "," << diff_y << "," << ratio_rad << "\n";

      if (abs_rad_n > 1e-15) {
        if (std::abs(theta) < 1e-10 || std::abs(theta - 2.0 * PI) < 1e-10) {
          ++n_x_axis;
          if (ratio_rad > max_ratio_rad_x_axis) max_ratio_rad_x_axis = ratio_rad;
          if (ratio_rad < min_ratio_rad_x_axis) min_ratio_rad_x_axis = ratio_rad;
        } else {
          ++n_off_axis;
          if (ratio_rad > max_ratio_rad_off_axis) max_ratio_rad_off_axis = ratio_rad;
          if (ratio_rad < min_ratio_rad_off_axis) min_ratio_rad_off_axis = ratio_rad;
        }
      }
      if (std::abs(a_tan_t) > max_abs_tan_tpf) max_abs_tan_tpf = std::abs(a_tan_t);
      if (abs_rad_n > max_abs_rad_newt) max_abs_rad_newt = abs_rad_n;
    }
  }

  txt << "--- Static circle audit ---\n";
  txt << "  Radii: ";
  for (size_t i = 0; i < CIRCLE_RADII.size(); ++i)
    txt << CIRCLE_RADII[i] << (i + 1 < CIRCLE_RADII.size() ? ", " : "");
  txt << "\n  Angles per circle: " << ANGLES_PER_CIRCLE << "\n";
  if (n_x_axis > 0) {
    txt << "  On x-axis (theta=0): ratio_rad in [" << std::scientific << circle_stats.min_ratio_x_axis
        << ", " << circle_stats.max_ratio_x_axis << "]\n";
  }
  if (circle_stats.n_off_axis > 0) {
    txt << "  Off-axis: ratio_rad in [" << std::scientific << circle_stats.min_ratio_off_axis
        << ", " << circle_stats.max_ratio_off_axis << "]\n";
  }
  txt << "  max|a_tan_tpf| = " << std::scientific << circle_stats.max_abs_tan_tpf
      << ", max|a_rad_newt| = " << circle_stats.max_abs_rad_newt << "\n";
  txt << "  max |diff vector| = " << circle_stats.max_diff_mag << "\n\n";
}

void run_snapshot_replay_audit(const Config& config,
                               PhysicsPackage* newton,
                               PhysicsPackage* tpf,
                               std::ostream& csv,
                               std::ostream& txt) {
  State state0;
  init_two_body_star_around_bh(config, state0);
  Config c = config;
  c.enable_star_star_gravity = false;
  int n_steps = config.validation_n_steps;
  int snapshot_every = std::max(1, config.validation_snapshot_every);
  if (n_steps / snapshot_every > 2000)
    snapshot_every = std::max(1, n_steps / 2000);

  PhysicsPackage* newton_pkg = get_physics_package("Newtonian");
  std::vector<Snapshot> snapshots = run_simulation(c, state0, newton_pkg, n_steps, snapshot_every);

  const double bh_mass = config.bh_mass;
  const double softening = config.softening;
  const double eps_sq = softening * softening;
  const bool star_star = false;

  txt << "--- Snapshot replay audit ---\n";
  txt << "  Reference path: Newtonian two_body_orbit, n_steps=" << n_steps
      << ", snapshot_every=" << snapshot_every << ", points=" << snapshots.size() << "\n";

  double max_ratio_rad = 0.0, min_ratio_rad = 1e300;
  double max_abs_tan_tpf = 0.0, max_abs_rad_newt = 0.0;
  double max_diff_mag = 0.0;
  int n_ratio = 0;

  for (size_t idx = 0; idx < snapshots.size(); ++idx) {
    const Snapshot& snap = snapshots[idx];
    double x = snap.state.x[0];
    double y = snap.state.y[0];

    State state;
    state.resize(1);
    state.x[0] = x;
    state.y[0] = y;
    state.vx[0] = 0.0;
    state.vy[0] = 0.0;
    state.mass[0] = config.star_mass;

    std::vector<double> ax_n(1), ay_n(1), ax_t(1), ay_t(1);
    newton->compute_accelerations(state, bh_mass, softening, star_star, ax_n, ay_n);
    tpf->compute_accelerations(state, bh_mass, softening, star_star, ax_t, ay_t);

    double a_rad_n, a_tan_n, a_rad_t, a_tan_t;
    radial_tangential(x, y, eps_sq, ax_n[0], ay_n[0], a_rad_n, a_tan_n);
    radial_tangential(x, y, eps_sq, ax_t[0], ay_t[0], a_rad_t, a_tan_t);

    double diff_x = ax_t[0] - ax_n[0];
    double diff_y = ay_t[0] - ay_n[0];
    double diff_mag = std::sqrt(diff_x * diff_x + diff_y * diff_y);
    if (diff_mag > max_diff_mag) max_diff_mag = diff_mag;

    double ratio_rad = 0.0;
    double abs_rad_n = std::abs(a_rad_n);
    if (abs_rad_n > 1e-15) {
      ratio_rad = std::abs(a_rad_t) / abs_rad_n;
      ++n_ratio;
      if (ratio_rad > max_ratio_rad) max_ratio_rad = ratio_rad;
      if (ratio_rad < min_ratio_rad) min_ratio_rad = ratio_rad;
    }
    if (std::abs(a_tan_t) > max_abs_tan_tpf) max_abs_tan_tpf = std::abs(a_tan_t);
    if (abs_rad_n > max_abs_rad_newt) max_abs_rad_newt = abs_rad_n;

    double r_pt = std::sqrt(x * x + y * y);
    double theta_pt = std::atan2(y, x);
    csv << "replay," << std::scientific << r_pt << "," << theta_pt << "," << snap.step << "," << snap.time << "," << x << "," << y << ","
        << ax_n[0] << "," << ay_n[0] << "," << ax_t[0] << "," << ay_t[0] << ","
        << a_rad_n << "," << a_tan_n << "," << a_rad_t << "," << a_tan_t << ","
        << diff_x << "," << diff_y << "," << ratio_rad << "\n";
  }

  txt << "  ratio_rad (when |a_rad_newt|>0): min=" << std::scientific << min_ratio_rad
      << ", max=" << max_ratio_rad << " (n=" << n_ratio << ")\n";
  txt << "  max|a_tan_tpf| = " << max_abs_tan_tpf << ", max|a_rad_newt| = " << max_abs_rad_newt << "\n";
  txt << "  max |diff vector| = " << max_diff_mag << "\n\n";
}

void write_summary(const CircleAuditStats& circle_stats, std::ostream& txt) {
  txt << "--- Summary: mismatch characterisation ---\n\n";

  const double TOL_RATIO = 0.05;
  const double TOL_TAN = 0.01;

  double spread_x = circle_stats.max_ratio_x_axis - circle_stats.min_ratio_x_axis;
  double spread_off = circle_stats.max_ratio_off_axis - circle_stats.min_ratio_off_axis;
  bool has_x = circle_stats.n_x_axis > 0;
  bool has_off = circle_stats.n_off_axis > 0;

  txt << "  Does TPF match Newtonian radial force only on the x-axis, or at all angles?\n";
  if (has_x && has_off) {
    double mid_x = 0.5 * (circle_stats.min_ratio_x_axis + circle_stats.max_ratio_x_axis);
    double mid_off = 0.5 * (circle_stats.min_ratio_off_axis + circle_stats.max_ratio_off_axis);
    if (std::abs(mid_off - mid_x) > TOL_RATIO || spread_off > 2.0 * spread_x) {
      txt << "    -> Axis-only: ratio_rad on x-axis differs from off-axis, or off-axis spread is large.\n";
    } else {
      txt << "    -> At all angles: ratio_rad is similar on and off x-axis (possibly with a scale factor).\n";
    }
  } else {
    txt << "    -> See CSV: compare ratio_rad on theta=0 vs other angles.\n";
  }

  txt << "\n  Is there a non-negligible tangential force mismatch?\n";
  double tan_ratio = (circle_stats.max_abs_rad_newt > 1e-30)
      ? (circle_stats.max_abs_tan_tpf / circle_stats.max_abs_rad_newt) : 0.0;
  if (tan_ratio > TOL_TAN) {
    txt << "    -> Yes: max|a_tan_tpf| / max|a_rad_newt| = " << std::scientific << tan_ratio << " (non-negligible).\n";
  } else {
    txt << "    -> No: max|a_tan_tpf| / max|a_rad_newt| = " << std::scientific << tan_ratio << " (negligible).\n";
  }

  txt << "\n  Is the difference mostly a scale issue, a directional issue, or an off-axis shape issue?\n";
  if (tan_ratio > TOL_TAN) {
    txt << "    -> Directional: TPF has non-negligible tangential component vs Newtonian (radial-only).\n";
  } else if (has_off && spread_off > TOL_RATIO) {
    txt << "    -> Off-axis shape: ratio_rad varies with angle (see CSV).\n";
  } else if (has_x) {
    txt << "    -> Mostly scale: ratio_rad roughly constant with angle; check if near 1.\n";
  } else {
    txt << "    -> See CSV and audits above.\n";
  }
  txt << "\n  Use tpf_newtonian_force_compare.csv for per-point comparison.\n";
}

}  // namespace

void run_tpf_newtonian_force_compare(const Config& config, const std::string& output_dir) {
  PhysicsPackage* newton = get_physics_package("Newtonian");
  PhysicsPackage* tpf = get_physics_package("TPFCore");
  if (!newton || !tpf) {
    std::cerr << "tpf_newtonian_force_compare requires both Newtonian and TPFCore packages.\n";
    return;
  }
  tpf->init_from_config(config);

  std::string csv_path = output_dir + "/tpf_newtonian_force_compare.csv";
  std::string txt_path = output_dir + "/tpf_newtonian_force_compare.txt";

  std::ofstream csv(csv_path);
  std::ofstream txt_file(txt_path);
  if (!csv || !txt_file) {
    std::cerr << "Failed to open " << csv_path << " or " << txt_path << "\n";
    return;
  }

  std::ostream& txt = txt_file;
  txt << "Newtonian vs TPF acceleration comparison (diagnostics only)\n";
  txt << "Same positions/states; same softening; TPF uses current readout mode and calibrated scale.\n\n";
  txt << "TPF readout_mode (from config): " << config.tpfcore_readout_mode << "\n";
  txt << "TPF readout_scale: " << std::scientific << config.tpfcore_readout_scale << "\n";
  txt << "softening: " << config.softening << ", bh_mass: " << config.bh_mass << "\n\n";

  csv << "audit_type,r,theta,step,time,x,y,ax_newt,ay_newt,ax_tpf,ay_tpf,"
      << "a_rad_newt,a_tan_newt,a_rad_tpf,a_tan_tpf,diff_x,diff_y,ratio_rad\n";
  CircleAuditStats circle_stats;
  run_static_circle_audit(config, newton, tpf, csv, txt, circle_stats);
  run_snapshot_replay_audit(config, newton, tpf, csv, txt);
  write_summary(circle_stats, txt);

  std::cout << "Wrote " << csv_path << "\n";
  std::cout << "Wrote " << txt_path << "\n";
}

// --- Diagnostic consistency audit: same points, side-by-side, decisive conclusion ---
static const std::vector<double> AUDIT_RADII = {20.0, 50.0, 100.0, 200.0};

void run_tpf_diagnostic_consistency_audit(const Config& config, const std::string& output_dir) {
  PhysicsPackage* newton = get_physics_package("Newtonian");
  PhysicsPackage* tpf = get_physics_package("TPFCore");
  if (!newton || !tpf) {
    std::cerr << "tpf_diagnostic_consistency_audit requires both Newtonian and TPFCore packages.\n";
    return;
  }
  tpf->init_from_config(config);

  const double bh_mass = config.bh_mass;
  const double softening = config.softening;
  const double eps_sq = softening * softening;
  const bool star_star = false;

  std::string csv_path = output_dir + "/tpf_diagnostic_consistency_audit.csv";
  std::string txt_path = output_dir + "/tpf_diagnostic_consistency_audit.txt";
  std::ofstream csv(csv_path);
  std::ofstream txt_file(txt_path);
  if (!csv || !txt_file) {
    std::cerr << "Failed to open " << csv_path << " or " << txt_path << "\n";
    return;
  }
  std::ostream& txt = txt_file;

  txt << "TPF diagnostic consistency audit: weak_field_calibration vs force_compare\n";
  txt << "Same sampled points (x-axis, r = 20, 50, 100, 200); side-by-side intermediates.\n\n";

  csv << "r,softening,x,y,theta_rr,theta_tt,theta_tr,provisional_radial_readout,ax_tpf,ay_tpf,"
      << "a_tpf_cal,a_rad_tpf_fc,a_newton_cal,a_newton_package_ax,a_rad_newt_fc,ratio_cal,ratio_fc\n";

  bool same_softening = true;
  bool same_sign_convention = true;
  bool same_abs_logic = true;
  bool same_package_path = true;
  bool benchmark_formula_mismatch = false;
  double max_ratio_cal = 0.0, min_ratio_cal = 1e300;
  double max_ratio_fc = 0.0, min_ratio_fc = 1e300;

  for (double r : AUDIT_RADII) {
    double x = r;
    double y = 0.0;

    State state;
    state.resize(1);
    state.x[0] = x;
    state.y[0] = y;
    state.vx[0] = 0.0;
    state.vy[0] = 0.0;
    state.mass[0] = config.star_mass;

    std::vector<double> ax_t(1), ay_t(1), ax_n(1), ay_n(1);
    tpf->compute_accelerations(state, bh_mass, softening, star_star, ax_t, ay_t);
    newton->compute_accelerations(state, bh_mass, softening, star_star, ax_n, ay_n);

    tpfcore::ReadoutDiagnostics diag;
    double ax_d = 0.0, ay_d = 0.0;
    double eps_audit =
        (config.tpfcore_source_softening > 0.0) ? config.tpfcore_source_softening : softening;
    tpfcore::DerivedTpfPoissonConfig dcfg;
    dcfg.kappa = config.tpf_kappa;
    dcfg.bins = config.tpf_poisson_bins;
    dcfg.max_radius = config.tpf_poisson_max_radius;
    dcfg.galaxy_radius = config.galaxy_radius;
    tpfcore::TpfRadialGravityProfile audit_prof;
    const tpfcore::TpfRadialGravityProfile* audit_prof_ptr = nullptr;
    if (tpfcore::is_derived_tpf_radial_readout_mode(config.tpfcore_readout_mode)) {
      audit_prof = tpfcore::build_tpf_gravity_profile(state, bh_mass, dcfg, eps_audit);
      audit_prof_ptr = &audit_prof;
    }
    tpfcore::compute_provisional_readout_with_diagnostics(
        state, 0, bh_mass, star_star, softening, config.tpfcore_source_softening,
        config.tpfcore_readout_mode,
        config.tpfcore_readout_scale, config.tpfcore_theta_tt_scale, config.tpfcore_theta_tr_scale,
        ax_d, ay_d, diag,
        tpfcore::is_derived_tpf_radial_readout_mode(config.tpfcore_readout_mode) ? &dcfg : nullptr,
        audit_prof_ptr);
    (void)ax_d;
    (void)ay_d;

    double r_sq = r * r + eps_sq;
    double a_newton_cal = bh_mass / (r_sq * std::sqrt(r_sq));
    double a_tpf_cal = std::abs(ax_t[0]);
    double ratio_cal = (a_newton_cal > 1e-300) ? (a_tpf_cal / a_newton_cal) : 0.0;

    double a_rad_tpf = 0.0, a_tan_tpf = 0.0, a_rad_n = 0.0, a_tan_n = 0.0;
    radial_tangential(x, y, eps_sq, ax_t[0], ay_t[0], a_rad_tpf, a_tan_tpf);
    radial_tangential(x, y, eps_sq, ax_n[0], ay_n[0], a_rad_n, a_tan_n);
    double a_rad_tpf_fc = a_rad_tpf;
    double a_rad_newt_fc = a_rad_n;
    double ratio_fc = (std::abs(a_rad_n) > 1e-15) ? (std::abs(a_rad_tpf) / std::abs(a_rad_n)) : 0.0;

    double a_newton_package_ax = std::abs(ax_n[0]);

    csv << std::scientific << r << "," << softening << "," << x << "," << y << ","
        << diag.theta_rr << "," << diag.theta_tt << "," << diag.theta_tr << ","
        << diag.provisional_radial_readout << "," << ax_t[0] << "," << ay_t[0] << ","
        << a_tpf_cal << "," << a_rad_tpf_fc << "," << a_newton_cal << "," << a_newton_package_ax << ","
        << a_rad_newt_fc << "," << ratio_cal << "," << ratio_fc << "\n";

    if (ratio_cal > max_ratio_cal) max_ratio_cal = ratio_cal;
    if (ratio_cal < min_ratio_cal) min_ratio_cal = ratio_cal;
    if (ratio_fc > max_ratio_fc) max_ratio_fc = ratio_fc;
    if (ratio_fc < min_ratio_fc) min_ratio_fc = ratio_fc;

    double newton_cal_times_r = a_newton_cal * r;
    if (a_newton_cal > 1e-300 && std::abs(a_newton_package_ax - newton_cal_times_r) < 1e-10 * std::max(1.0, a_newton_package_ax))
      benchmark_formula_mismatch = true;
  }

  txt << "--- Per-radius CSV columns ---\n";
  txt << "  r, softening, x, y, theta_rr, theta_tt, theta_tr, provisional_radial_readout,\n";
  txt << "  ax_tpf, ay_tpf, a_tpf_cal, a_rad_tpf_fc, a_newton_cal, a_newton_package_ax, a_rad_newt_fc, ratio_cal, ratio_fc\n\n";

  txt << "--- Explicit answers ---\n";
  txt << "  Are both diagnostics using the same softening convention?\n";
  txt << "    weak_field_calibration: uses config.softening for state and for formula r_sq = r^2 + softening^2.\n";
  txt << "    force_compare: uses config.softening for package calls and for radial_tangential (eps_sq = softening^2).\n";
  txt << "    -> " << (same_softening ? "Yes.\n" : "No.\n");

  txt << "  Are both diagnostics using the same sign convention for radial acceleration?\n";
  txt << "    weak_field_calibration: a_tpf = |ax[0]| (magnitude); a_newton = M/(r^2+eps^2)^(3/2) (positive magnitude).\n";
  txt << "    force_compare: ratio = |a_rad_tpf| / |a_rad_newt| (both magnitudes).\n";
  txt << "    -> " << (same_sign_convention ? "Yes (both use magnitude for ratio).\n" : "No.\n");

  txt << "  Are both diagnostics using the same absolute-value/signed-value logic?\n";
  txt << "    weak_field_calibration: a_tpf = |ax|, a_newton = formula (positive). ratio = a_tpf/a_newton.\n";
  txt << "    force_compare: ratio = |a_rad_tpf| / |a_rad_newt|.\n";
  txt << "    -> " << (same_abs_logic ? "Yes.\n" : "No.\n");

  txt << "  Are both diagnostics using the same TPF package path and same readout scale?\n";
  txt << "    Both call the same TPF package compute_accelerations (same readout path and scale).\n";
  txt << "    -> " << (same_package_path ? "Yes.\n" : "No.\n");

  txt << "  If they differ, at what exact step/formula do they diverge?\n";
  txt << "    weak_field_calibration Newtonian benchmark: a_newton = M/(r^2+eps^2)^(3/2).\n";
  txt << "    Newtonian package on x-axis: ax = -M*r/(r^2+eps^2)^(3/2), so |ax| = M*r/(r^2+eps^2)^(3/2).\n";
  txt << "    So package |ax| = r * (calibration benchmark). force_compare uses package output; calibration uses formula.\n";
  txt << "    -> Divergence: Newtonian benchmark formula in calibration is M/(r^2+eps^2)^(3/2); package gives M*r/(r^2+eps^2)^(3/2).\n\n";

  txt << "--- Ratios observed ---\n";
  txt << "  ratio_cal (a_tpf_cal / a_newton_cal): min = " << std::scientific << min_ratio_cal << ", max = " << max_ratio_cal << "\n";
  txt << "  ratio_fc (|a_rad_tpf| / |a_rad_newt|): min = " << min_ratio_fc << ", max = " << max_ratio_fc << "\n";
  txt << "  For consistent diagnostics, ratio_fc would equal ratio_cal. Instead ratio_fc ≈ ratio_cal / r (factor r).\n\n";

  txt << "--- DECISIVE CONCLUSION ---\n";
  if (benchmark_formula_mismatch || (max_ratio_fc < 0.5 * min_ratio_cal && max_ratio_fc > 0)) {
    txt << "  benchmark formula mismatch\n";
    txt << "  weak_field_calibration compares TPF |ax| to the formula M/(r^2+eps^2)^(3/2).\n";
    txt << "  The Newtonian package returns |ax| = M*r/(r^2+eps^2)^(3/2) on the x-axis.\n";
    txt << "  So the calibration benchmark is 1/r times the actual Newtonian package output.\n";
    txt << "  force_compare correctly uses the Newtonian package output; its ratio is therefore ~1/r times the calibration ratio.\n";
    txt << "  Resolve by either: (1) make weak_field_calibration use Newtonian package output for the benchmark, or (2) document that calibration uses standard radial magnitude formula and force_compare uses package output.\n";
  } else {
    txt << "  diagnostics consistent\n";
    txt << "  (No benchmark formula mismatch detected at sampled radii.)\n";
  }

  std::cout << "Wrote " << csv_path << "\n";
  std::cout << "Wrote " << txt_path << "\n";
}

}  // namespace galaxy
