/**
 * v11 weak-field correspondence audit (manuscript v11, static sector).
 * Correspondence-only: Φ is the paper's benchmark static potential for a point mass (softened for numerics).
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

/** Paper v11 weak-field calibration: κ = 16πG (SI) relating C_{00} sector to Poisson correspondence. */
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
};

}  // namespace

void run_v11_weak_field_correspondence_audit(const Config& config, const std::string& output_dir) {
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
         "Laplacian_trace_check,I_v11,C00_principal,Cxx_principal,Cyy_principal,Czz_principal\n";

  double z_last = 0.0;
  double I_last = 0.0;
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
    c00 = 0.5 * I_s;
    cxx = kap * (theta_xx * theta_xx) - kLambda4D * theta_tr * theta_xx - 0.5 * I_s;
    cyy = kap * (theta_yy * theta_yy) - kLambda4D * theta_tr * theta_yy - 0.5 * I_s;
    czz = kap * (theta_zz * theta_zz) - kLambda4D * theta_tr * theta_zz - 0.5 * I_s;

    csv << std::scientific << std::setprecision(17) << z << "," << config.softening << "," << P.Phi(z) << ","
        << 0.0 << "," << 0.0 << "," << P.dPhi_dz(z) << "," << theta_xx << "," << theta_yy << "," << theta_zz << ","
        << theta_tr << "," << theta_tr << "," << I_s << "," << c00 << "," << cxx << "," << cyy << "," << czz << "\n";
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
  txt << "Potential Phi (correspondence benchmark only; numerical softening eps = config softening):\n";
  txt << "  Phi(z) = -G*M / sqrt(z^2 + eps^2)  on the +z axis (point mass M = bh_mass at origin).\n";
  txt << "Paper construction (static, flat spatial background, Theta_{0mu} = 0):\n";
  txt << "  Xi_i = d_i Phi, Xi_0 = 0  =>  on axis Xi_z = dPhi/dz.\n";
  txt << "  Theta_ij = d_i d_j Phi  (Cartesian; covariant = partial in this limit).\n";
  txt << "  Theta = trace(Theta_ij) = Laplacian(Phi) in Cartesian coordinates.\n";
  txt << "  I = Theta_{mu nu} Theta^{mu nu} - lambda Theta^2  with lambda = 1/4 (4D).\n";
  txt << "  Here I reduces to sum_ij Theta_ij^2 - lambda * Theta^2 (spatial; Theta_{0mu}=0).\n\n";
  txt << "C_{mu nu} — principal part from manuscript Eq. (10) ONLY (Delta C_{mu nu} deliberately omitted):\n";
  txt << "  C_principal = kappa * Theta_mu^alpha Theta_{nu alpha} - lambda * Theta * Theta_{mu nu}"
         " - (1/2) * eta_{mu nu} * I\n";
  txt << "  eta = diag(-1,1,1,1). kappa = 16*pi*G (SI) per v11 weak-field calibration note.\n";
  txt << "  Delta C: NOT computed (omitted per v11 paper scope; connection-variation terms deferred).\n\n";
  txt << "VDSG: must be off for this audit (tpf_vdsg_coupling == 0). This path does not use VDSG.\n\n";
  txt << "Correspondence reminder (paper): in the calibrated static weak-field sector, Theta_ij Theta^ij\n";
  txt << "is tied to coarse matter density rho via Poisson correspondence; exterior vacuum behavior is\n";
  txt << "a limiting case. This CSV is a tensor profile audit, not a claim of exact global equality.\n\n";
  txt << "Eq. (9) configuration equation: NOT audited here (no covariant divergence check in this artifact).\n\n";
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

}  // namespace galaxy
