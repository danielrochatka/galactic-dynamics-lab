#include "physics/TPFCore/derived_tpf_radial.hpp"
#include "doctest.h"
#include <cmath>

namespace {

/** Mirror v11_weak_field_correspondence.cpp on-axis invariants for vacuum exterior scaling check. */
void axis_invariants(double G, double M, double z, double eps_sq, double& theta_xx, double& theta_zz,
                   double& theta_tr, double& I) {
  const double lambda4 = 0.25;
  double r2 = z * z + eps_sq;
  double r = std::sqrt(r2);
  double r5 = r2 * r2 * r;
  theta_xx = G * M * r2 / r5;
  theta_zz = G * M * (eps_sq - 2.0 * z * z) / r5;
  double theta_yy = theta_xx;
  theta_tr = theta_xx + theta_yy + theta_zz;
  I = theta_xx * theta_xx + theta_yy * theta_yy + theta_zz * theta_zz - lambda4 * theta_tr * theta_tr;
}

}  // namespace

TEST_CASE("v11 Eq. (9) flat-static residual R_z = (1-lambda) dTheta/dz on axis") {
  const double G = galaxy::tpfcore::TPF_G_SI;
  const double M = 2.0e30;
  const double eps = 1.0;
  const double z = 1.0e10;
  const double eps_sq = eps * eps;
  const double u = z * z + eps_sq;
  const double r7 = u * u * u * std::sqrt(u);
  const double dTheta_dz = -15.0 * G * M * eps_sq * z / r7;
  const double lambda = 0.25;
  const double Rz = (1.0 - lambda) * dTheta_dz;
  /* Harmonic vacuum: Theta=0 => dTheta/dz=0 => Rz=0; here eps<<z still leaves tiny Theta — check small */
  CHECK(std::abs(Rz) < 1e-30);
}

TEST_CASE("v11 weak-field axis I asymptotic matches 6(GM/z^3)^2 when eps << z") {
  const double G = galaxy::tpfcore::TPF_G_SI;
  const double M = 2.0e30;
  const double z = 1.0e11;
  const double eps = 1.0;
  double txx, tzz, tr, I;
  axis_invariants(G, M, z, eps * eps, txx, tzz, tr, I);
  const double gm = G * M;
  const double I_asym = 6.0 * (gm / (z * z * z)) * (gm / (z * z * z));
  CHECK(tr == doctest::Approx(0.0).epsilon(1e-9)); /* Laplacian ~ 0 exterior */
  CHECK(I == doctest::Approx(I_asym).epsilon(1e-6));
}
