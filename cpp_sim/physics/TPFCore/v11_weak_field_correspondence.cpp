/**
 * v11 weak-field correspondence audit (manuscript v11, static sector).
 * Correspondence-only: benchmark potentials / 1D fields as stated in the paper.
 * No ΔC_{μν}, no particle dynamics, no legacy readout, no VDSG.
 */

#include "v11_weak_field_correspondence.hpp"
#include "derived_tpf_radial.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace galaxy {

namespace {

/** Paper v11: λ = 1/n in n = 4 dimensions (fixed; not tunable). */
constexpr double kLambda4D = 0.25;

/** Paper v11 weak-field calibration: κ = 16πG (SI); appears in Eq. (10) quadratic Θ·Θ term only. */
double v11_kappa_si() { return 16.0 * M_PI * tpfcore::TPF_G_SI; }

/**
 * Softened monopole Φ = -G M / R with R = sqrt(z^2 + eps^2) on the +z axis (x=y=0).
 * Matches correspondence use of a static potential; eps is numerical regularization only (config softening).
 */
struct AxisPointMassPhi {
  double G;
  double M_kg;
  double eps_sq;

  double R(double z) const { return std::sqrt(z * z + eps_sq); }
  double Phi(double z) const {
    double r = R(z);
    if (r < 1e-300) return 0.0;
    return -G * M_kg / r;
  }
  /** ∂_z Φ */
  double dPhi_dz(double z) const {
    double r = R(z);
    if (r < 1e-300) return 0.0;
    return G * M_kg * z / (r * r * r);
  }
  /** ∂_zz Φ on axis */
  double d2Phi_dzz(double z) const {
    double r2 = z * z + eps_sq;
    if (r2 < 1e-300) return 0.0;
    double r = std::sqrt(r2);
    double r5 = r2 * r2 * r;
    return G * M_kg * (eps_sq - 2.0 * z * z) / r5;
  }
  /** ∂_xx Φ = ∂_yy Φ on axis (x=y=0) */
  double d2Phi_dxx_on_axis(double z) const {
    double r2 = z * z + eps_sq;
    if (r2 < 1e-300) return 0.0;
    double r = std::sqrt(r2);
    double r5 = r2 * r2 * r;
    return G * M_kg / r5 * r2;  /* GM / R^3 */
  }
  /** Θ = ∇²Φ (4D static: spatial Laplacian) on axis = 3 G M eps² / R^5, R² = z²+eps² */
  double theta_trace_on_axis(double z) const {
    double r2 = z * z + eps_sq;
    if (r2 < 1e-300) return 0.0;
    double r = std::sqrt(r2);
    double r5 = r2 * r2 * r;
    return 3.0 * G * M_kg * eps_sq / r5;
  }
  /** dΘ/dz on axis (analytic); used for Eq. (9) residual audit */
  double d_theta_trace_dz_on_axis(double z) const {
    double u = z * z + eps_sq;
    if (u < 1e-300) return 0.0;
    double r7 = u * u * u * std::sqrt(u);
    return -15.0 * G * M_kg * eps_sq * z / r7;
  }
};

/** Ω = 2π/T; if T <= 0, Kepler circular period from D and G(ME+MM). */
double em_orbital_omega_rad_s(const Config& c) {
  const double G = tpfcore::TPF_G_SI;
  const double ME = c.v11_em_mass_earth_kg;
  const double MM = c.v11_em_mass_moon_kg;
  const double D = c.v11_em_mean_distance_m;
  double T = c.v11_em_sidereal_period_s;
  if (T <= 0.0 || !std::isfinite(T)) {
    if (!(D > 0.0) || !std::isfinite(D)) throw std::runtime_error("em_orbital_omega: invalid D");
    const double mtot = ME + MM;
    if (!(mtot > 0.0)) throw std::runtime_error("em_orbital_omega: ME+MM must be positive");
    T = 2.0 * M_PI * std::sqrt((D * D * D) / (G * mtot));
  }
  if (!(T > 0.0) || !std::isfinite(T)) throw std::runtime_error("em_orbital_omega: invalid period T");
  return 2.0 * M_PI / T;
}

/**
 * Paper v11 Sec. XI: x-axis from Earth center (0) to Moon (D). Point-mass exterior φ from Eq. (44):
 *   φ(x) = α ( ME/|x| + MM/|D−x| ),  0 < x < D.
 * Eq. (45): aTPF(x) = −dφ/dx + Ω²(x − xb),  xb = D MM/(ME+MM).
 * Single calibration: aTPF(x_cal) = −g0 with x_cal = calib radius (paper: Earth surface along +x).
 */
void run_axis_monopole_audit(const Config& config, const std::string& output_dir) {
  AxisPointMassPhi P;
  P.G = tpfcore::TPF_G_SI;
  P.M_kg = config.bh_mass;
  P.eps_sq = config.softening * config.softening;

  const double zmin = config.tpfcore_probe_radius_min;
  const double zmax = config.tpfcore_probe_radius_max;
  const int ns = std::max(2, config.tpfcore_probe_samples);

  std::ostringstream csv_path;
  csv_path << output_dir << "/tpf_v11_weak_field_correspondence_profile.csv";
  std::ofstream csv(csv_path.str());
  if (!csv) throw std::runtime_error("failed to open " + csv_path.str());

  csv << "z_m,softening_m,Phi_SI,Xi_x_SI,Xi_y_SI,Xi_z_SI,Theta_xx,Theta_yy,Theta_zz,Theta_trace,"
         "Laplacian_trace_check,I_v11,"
         "eq10_inner_bracket_C00_before_kappa,C00_eq10_principal_SI_kappa_times_full_bracket,"
         "Cxx_eq10_principal_SI_kappa_times_full_bracket,Cxx_eq10_normalized_by_kappa,"
         "Cyy_eq10_principal_SI_kappa_times_full_bracket,Cyy_eq10_normalized_by_kappa,"
         "Czz_eq10_principal_SI_kappa_times_full_bracket,Czz_eq10_normalized_by_kappa,"
         "eq9_residual_Rx_SI,eq9_residual_Ry_SI,eq9_residual_Rz_SI\n";

  const double one_minus_lambda = 1.0 - kLambda4D;

  double z_last = 0.0;
  double I_last = 0.0;
  double max_abs_Rz = 0.0;
  double z_at_max_abs_Rz = 0.0;
  /* Laplacian in Cartesian = Theta_trace for this Phi */
  for (int i = 0; i < ns; ++i) {
    double t = (ns == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(ns - 1);
    double z = zmin + t * (zmax - zmin);
    if (z <= 0.0) continue; /* stay on exterior axis away from singular division */

    double theta_xx, theta_yy, theta_zz, theta_tr, I_s, c00, cxx, cyy, czz;
    theta_xx = P.d2Phi_dxx_on_axis(z);
    theta_yy = theta_xx;
    theta_zz = P.d2Phi_dzz(z);
    theta_tr = theta_xx + theta_yy + theta_zz;
    double sum_sq = theta_xx * theta_xx + theta_yy * theta_yy + theta_zz * theta_zz;
    I_s = sum_sq - kLambda4D * theta_tr * theta_tr;

    const double kap = v11_kappa_si();
    /* Manuscript v11 Eq. (10) as used in weak-field reduction: C_{mu nu} = kappa * ( ... - (1/2) g_{mu nu} I + ... ),
       i.e. kappa multiplies the entire principal bracket including the -(1/2) g I term. DeltaC omitted here.
       Static weak field, Theta_{0 mu}=0, g_{00}=-1: inner bracket C00 = -(1/2) g_{00} I = +I/2 => C00 = kappa * I/2. */
    const double inner_00 = 0.5 * I_s;
    const double inner_xx =
        theta_xx * theta_xx - kLambda4D * theta_tr * theta_xx - 0.5 * I_s;
    const double inner_yy =
        theta_yy * theta_yy - kLambda4D * theta_tr * theta_yy - 0.5 * I_s;
    const double inner_zz =
        theta_zz * theta_zz - kLambda4D * theta_tr * theta_zz - 0.5 * I_s;
    c00 = kap * inner_00;
    cxx = kap * inner_xx;
    cyy = kap * inner_yy;
    czz = kap * inner_zz;

    /* Eq. (9) flat static: partial_i Theta^i{}_j - lambda partial_j Theta = (1-lambda) partial_j Theta
       for this Phi ansatz (see summary). On +z axis: Rx=Ry=0, Rz=(1-lambda)*dTheta/dz. */
    const double dth_dz = P.d_theta_trace_dz_on_axis(z);
    const double rx = 0.0;
    const double ry = 0.0;
    const double rz = one_minus_lambda * dth_dz;
    const double arz = std::abs(rz);
    if (arz > max_abs_Rz) {
      max_abs_Rz = arz;
      z_at_max_abs_Rz = z;
    }

    csv << std::scientific << std::setprecision(17) << z << "," << config.softening << "," << P.Phi(z) << ","
        << 0.0 << "," << 0.0 << "," << P.dPhi_dz(z) << "," << theta_xx << "," << theta_yy << "," << theta_zz << ","
        << theta_tr << "," << theta_tr << "," << I_s << "," << inner_00 << "," << c00 << ","
        << cxx << "," << inner_xx << "," << cyy << "," << inner_yy << "," << czz << "," << inner_zz << ","
        << rx << "," << ry << "," << rz << "\n";
    z_last = z;
    I_last = I_s;
  }

  std::ostringstream txt_path;
  txt_path << output_dir << "/tpf_v11_weak_field_correspondence_summary.txt";
  std::ofstream txt(txt_path.str());
  if (!txt) throw std::runtime_error("failed to open " + txt_path.str());

  txt << "=== TPF manuscript v11 — weak-field correspondence audit (NOT dynamics) ===\n\n";
  txt << "Scope: static weak-field sector only. This run does not implement direct TPF particle integration,\n";
  txt << "legacy provisional readout accelerations, geodesic motion, or Einstein equations.\n\n";
  txt << "Benchmark: axis_monopole (this summary). Potential Phi (correspondence benchmark only; numerical softening eps = config softening):\n";
  txt << "  Phi(z) = -G*M / sqrt(z^2 + eps^2)  on the +z axis (point mass M = bh_mass at origin).\n";
  txt << "Paper construction (static, flat spatial background, Theta_{0mu} = 0):\n";
  txt << "  Xi_i = d_i Phi, Xi_0 = 0  =>  on axis Xi_z = dPhi/dz.\n";
  txt << "  Theta_ij = d_i d_j Phi  (Cartesian; covariant = partial in this limit).\n";
  txt << "  Theta = trace(Theta_ij) = Laplacian(Phi) in Cartesian coordinates.\n";
  txt << "  I = Theta_{mu nu} Theta^{mu nu} - lambda Theta^2  with lambda = 1/4 (4D).\n";
  txt << "  Here I reduces to sum_ij Theta_ij^2 - lambda * Theta^2 (spatial; Theta_{0mu}=0).\n\n";
  txt << "C_{mu nu} — principal part from manuscript Eq. (10) ONLY (Delta C_{mu nu} deliberately omitted):\n";
  txt << "  C_principal = kappa * ( Theta_mu^alpha Theta_{nu alpha} - lambda * Theta * Theta_{mu nu}"
         " - (1/2) * g_{mu nu} * I\n";
  txt << "  eta = diag(-1,1,1,1) for weak-field g; kappa = 16*pi*G (SI). The factor kappa multiplies the ENTIRE\n";
  txt << "  parenthesis (including the -(1/2) g I term), matching the paper’s weak-field display C00 ~ (kappa/2) I\n";
  txt << "  when Theta_{0 mu} ~ 0. CSV: C*_eq10_principal_SI_kappa_times_full_bracket = kappa * (inner bracket);\n";
  txt << "  C*_eq10_normalized_by_kappa = inner bracket = (principal / kappa) (same units as I; not a dynamics proof).\n";
  txt << "  Delta C: NOT computed (omitted per v11 paper scope; connection-variation terms deferred).\n\n";
  txt << "VDSG: must be off for this audit (tpf_vdsg_coupling == 0). This path does not use VDSG.\n\n";
  txt << "Correspondence reminder (paper): in the calibrated static weak-field sector, Theta_ij Theta^ij\n";
  txt << "is tied to coarse matter density rho via Poisson correspondence; exterior vacuum behavior is\n";
  txt << "a limiting case. This CSV is a tensor profile audit, not a claim of exact global equality.\n\n";
  txt << "Eq. (9) residual audit (correspondence construction only; NOT a proof of full TPF dynamics):\n";
  txt << "  Under flat static assumptions, the paper's Eq. (9) is nabla_mu (Theta^mu{}_nu - lambda delta^mu{}_nu Theta)=0.\n";
  txt << "  For smooth Phi with Theta_ij = d_i d_j Phi and Theta = nabla^2 Phi (static 4D trace), one has\n";
  txt << "    sum_i d_i Theta_ij - lambda d_j Theta = (1-lambda) d_j Theta  (algebraic identity for this ansatz).\n";
  txt << "  On the +z axis, symmetry gives Rx=Ry=0 and Rz=(1-lambda)*dTheta/dz with Theta=trace(Hessian Phi).\n";
  txt << "  Interpretation: nonzero Rz does NOT mean 'TPF failure' — this benchmark Phi is a softened monopole\n";
  txt << "  (numerical eps), not asserted to be an on-shell solution of the full TPF configuration problem.\n";
  txt << "  Exact harmonic exterior (nabla^2 Phi=0 everywhere): Theta=0 => all residuals vanish.\n";
  txt << "  With softening, Theta and dTheta/dz are generally nonzero near the regularization scale; residuals\n";
  txt << "  are expected unless eps<<z in the vacuum-like regime.\n\n";
  txt << "Eq. (9) residual summary (this run):\n";
  txt << "  sup_z |Rz| = " << std::scientific << max_abs_Rz << " SI at z = " << z_at_max_abs_Rz
         << " m (Rx=Ry=0 on axis by symmetry).\n";
  txt << "  Classification: approximate audit only — flat Cartesian, static, axis benchmark; softening breaks\n";
  txt << "  exact harmonic vacuum; nonzero Rz reflects the chosen regularized Phi, not a standalone 'TPF error'.\n\n";
  if (z_last > 0.0 && config.softening > 0.0 && z_last > 10.0 * config.softening) {
    const double gm = P.G * P.M_kg;
    const double I_asymptotic = 6.0 * (gm / (z_last * z_last * z_last)) * (gm / (z_last * z_last * z_last));
    const double ratio = (I_asymptotic > 0.0) ? (I_last / I_asymptotic) : 0.0;
    txt << "Scaling sanity (correspondence-level only; axis + exterior vacuum limit, eps << z):\n";
    txt << "  At largest tabulated z = " << std::scientific << z_last << " m, I_computed = " << I_last << "\n";
    txt << "  Leading vacuum axis form I ~ 6 (G M / z^3)^2 gives I_asymptotic = " << I_asymptotic
           << "  ratio I/I_asymptotic = " << ratio << "\n";
    txt << "  (Deviations expected when softening is not negligible vs z; not a proof of global field equations.)\n\n";
  }
  txt << "Outputs:\n";
  txt << "  " << csv_path.str() << "\n";
  txt << "  " << txt_path.str() << "\n";
}

void run_earth_moon_line_correspondence_benchmark(const Config& config, const std::string& output_dir) {
  const double G = tpfcore::TPF_G_SI;
  const double ME = config.v11_em_mass_earth_kg;
  const double MM = config.v11_em_mass_moon_kg;
  const double D = config.v11_em_mean_distance_m;
  const double RE = config.v11_em_calib_surface_radius_m;
  const double g0 = config.v11_em_calib_surface_g_m_s2;

  if (!(ME > 0.0) || !(MM > 0.0) || !(D > 0.0) || !(RE > 0.0) || !(RE < D) || !(g0 > 0.0)) {
    throw std::runtime_error(
        "earth_moon_line_of_centers: require positive ME, MM, D, g0 and 0 < RE < D (paper line-of-centers).");
  }

  const double Omega = em_orbital_omega_rad_s(config);
  const double xb = D * MM / (ME + MM);

  /* Denominator for calibration: d/dx (ME/x + MM/(D-x)) at x=RE is -ME/RE^2 + MM/(D-RE)^2 */
  const double dphi_prefactor_denom = ME / (RE * RE)-MM / ((D - RE) * (D - RE));
  if (!(std::isfinite(dphi_prefactor_denom)) || std::abs(dphi_prefactor_denom) < 1e-200) {
    throw std::runtime_error("earth_moon_line_of_centers: invalid calibration denominator at RE (check RE vs D).");
  }
  /* aTPF(RE) = alpha * (ME/RE^2 - MM/(D-RE)^2) + Omega^2(RE-xb) = -g0  (paper: ~9.81 m/s^2 along +x toward Earth) */
  const double alpha = (-g0 - Omega * Omega * (RE - xb)) / dphi_prefactor_denom;

  const double x_min = config.tpfcore_probe_radius_min;
  const double x_max = config.tpfcore_probe_radius_max;
  const int ns = std::max(2, config.tpfcore_probe_samples);

  if (!(x_min > 0.0) || !(x_max > x_min) || !(x_max < D)) {
    throw std::runtime_error(
        "earth_moon_line_of_centers: require 0 < tpfcore_probe_radius_min < tpfcore_probe_radius_max < "
        "v11_em_mean_distance_m (interior point-mass line; paper Sec. XI D).");
  }

  std::ostringstream csv_path;
  csv_path << output_dir << "/tpf_v11_earth_moon_line_correspondence_benchmark.csv";
  std::ofstream csv(csv_path.str());
  if (!csv) throw std::runtime_error("failed to open " + csv_path.str());

  csv << "x_m,D_m,M_Earth_kg,M_Moon_kg,Omega_rad_s,x_b_m,G_SI,alpha_calibrated_SI,"
         "phi_tpf_SI,dphi_dx_SI,neg_dphi_dx_SI,centrifugal_tpf_Omega2_x_minus_xb_m_s2,"
         "a_tpf_correspondence_Eq45_m_s2,a_newtonian_line_Eq46_m_s2,"
         "benchmark_difference_TPF_Eq45_minus_Newtonian_Eq46_m_s2\n";

  for (int i = 0; i < ns; ++i) {
    double t = (ns == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(ns - 1);
    double x = x_min + t * (x_max - x_min);
    if (!(x > 0.0) || !(x < D)) continue;

    /* Eq. (44): φ(x) = α ( ME/x + MM/(D-x) ) */
    const double phi = alpha * (ME / x + MM / (D - x));
    /* dφ/dx = α ( -ME/x² + MM/(D-x)² ) */
    const double dphi_dx = alpha * (-ME / (x * x) + MM / ((D - x) * (D - x)));
    const double neg_dphi_dx = -dphi_dx;

    /* Eq. (45): co-rotating frame about barycenter */
    const double cf_tpf = Omega * Omega * (x - xb);
    const double a_tpf = neg_dphi_dx + cf_tpf;

    /* Eq. (46): benchmark Newtonian line-of-centers in rotating frame (paper display; centrifugal −Ω² x) */
    const double a_n46 = -G * ME / (x * x) + G * MM / ((D - x) * (D - x)) - Omega * Omega * x;

    csv << std::scientific << std::setprecision(17) << x << "," << D << "," << ME << "," << MM << "," << Omega << ","
        << xb << "," << G << "," << alpha << "," << phi << "," << dphi_dx << "," << neg_dphi_dx << "," << cf_tpf << ","
        << a_tpf << "," << a_n46 << "," << (a_tpf - a_n46) << "\n";
  }

  const char* em_csv_name = "tpf_v11_earth_moon_line_correspondence_benchmark.csv";
  std::ostringstream gnu_script_path;
  gnu_script_path << output_dir << "/tpf_v11_earth_moon_line_correspondence_benchmark.gnu";
  std::ofstream gnu(gnu_script_path.str());
  bool wrote_gnu = false;
  if (gnu) {
    wrote_gnu = true;
    gnu << "# Correspondence benchmark only — not direct TPF dynamics.\n";
    gnu << "# Plots TPF weak-field correspondence (Eq. 45) vs benchmark Newtonian line-of-centers (Eq. 46).\n";
    gnu << "# Run from the output directory: gnuplot tpf_v11_earth_moon_line_correspondence_benchmark.gnu\n";
    gnu << "set datafile separator comma\n";
    gnu << "set key top right\n";
    gnu << "set xlabel \"x from Earth center (m) — line of centers toward Moon\"\n";
    gnu << "set ylabel \"a_x (m/s^2) — co-rotating frame (paper Sec. XI D)\"\n";
    gnu << "plot '" << em_csv_name << "' using 1:13 every ::1 with lines lw 2 title "
           "'TPF correspondence Eq.45 (calibrated phi Eq.44)', \\\n";
    gnu << "     '' using 1:14 every ::1 with lines lw 2 title 'Benchmark Newtonian Eq.46', \\\n";
    gnu << "     '' using 1:15 every ::1 with lines dt 2 title "
           "'Benchmark difference (TPF Eq.45 minus Newtonian Eq.46)'\n";
  }

  std::ostringstream txt_path;
  txt_path << output_dir << "/tpf_v11_earth_moon_line_correspondence_benchmark_summary.txt";
  std::ofstream txt(txt_path.str());
  if (!txt) throw std::runtime_error("failed to open " + txt_path.str());

  txt << "=== TPF manuscript v11 — Earth–Moon weak-field LINE-OF-CENTERS correspondence benchmark ===\n\n";
  txt << "NOT full TPF many-body dynamics. NOT direct_tpf. NOT a particle integrator. NOT GR completion.\n";
  txt << "NOT VDSG. Delta C_mu_nu omitted. Paper-faithful static / weak-field / correspondence mapping only.\n\n";
  txt << "Paper references (v11 Sec. XI C–D):\n";
  txt << "  (42) nabla^2 phi = alpha rho  (weak-field Poisson correspondence).\n";
  txt << "  (43)-(44) phi from 3D Green's function for compact sources: point-mass exterior form used here.\n";
  txt << "  (45) aTPF(x) = -dphi/dx + Omega^2 (x - xb)  (co-rotating frame; xb = D*MM/(ME+MM)).\n";
  txt << "  (46) benchmark Newtonian line-of-centers field in the same rotating frame (paper display):\n";
  txt << "       a(x) = -G*ME/x^2*sgn(x) + G*MM/(D-x)^2*sgn(D-x) - Omega^2 x  (interior: 0 < x < D).\n";
  txt << "  (47) interpretive identification aTPF ~ C_xx along the line in the calibrated weak-field limit (audit labels only).\n\n";
  txt << "Inputs (defaults from paper Table II where applicable):\n";
  txt << "  ME = " << std::scientific << ME << " kg, MM = " << MM << " kg, D = " << D << " m\n";
  txt << "  Omega = " << Omega << " rad/s";
  if (config.v11_em_sidereal_period_s <= 0.0)
    txt << " (Kepler T = 2*pi*sqrt(D^3/(G*(ME+MM))) because v11_em_sidereal_period_s <= 0)\n";
  else
    txt << " (T = " << config.v11_em_sidereal_period_s << " s from config)\n";
  txt << "  xb = " << xb << " m\n";
  txt << "  Calibration (paper): aTPF(x_cal) = -g0 with x_cal = RE = " << RE << " m, g0 = " << g0 << " m/s^2\n";
  txt << "  => alpha = " << std::scientific << alpha << " SI (coupling from Eq. (44); not asserted equal to -G)\n";
  txt << "  Compare: -G = " << -G << " SI (Newtonian gradient uses G explicitly in Eq. 46 column).\n\n";
  txt << "Sampling: x in [" << x_min << ", " << x_max << "] m, n = " << ns << " (tpfcore_probe_* keys).\n\n";
  txt << "Column semantics:\n";
  txt << "  a_tpf_correspondence_Eq45_m_s2 — correspondence benchmark from (44)-(45) only.\n";
  txt << "  a_newtonian_line_Eq46_m_s2 — standard rotating-frame line-of-centers from paper Eq. (46).\n";
  txt << "  benchmark_difference_* — algebraic gap between those two benchmark constructions; NOT 'new physics'.\n\n";
  txt << "Artifacts:\n";
  txt << "  " << csv_path.str() << "\n";
  txt << "  " << txt_path.str() << "\n";
  if (wrote_gnu) txt << "  " << gnu_script_path.str() << " (gnuplot; correspondence plot)\n";
}

}  // namespace

void run_v11_weak_field_correspondence_audit(const Config& config, const std::string& output_dir) {
  if (config.v11_weak_field_correspondence_benchmark == "earth_moon_line_of_centers") {
    run_earth_moon_line_correspondence_benchmark(config, output_dir);
    return;
  }
  run_axis_monopole_audit(config, output_dir);
}

}  // namespace galaxy
