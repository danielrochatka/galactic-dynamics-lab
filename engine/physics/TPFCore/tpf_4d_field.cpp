#include "physics/TPFCore/tpf_4d_field.hpp"

#include <stdexcept>

namespace galaxy {
namespace tpfcore {

namespace {

double& xi_component_ref(Xi4& xi, std::size_t mu) {
  switch (mu) {
    case 0: return xi.t;
    case 1: return xi.x;
    case 2: return xi.y;
    case 3: return xi.z;
    default: throw std::out_of_range("Xi4 index out of range");
  }
}

double theta_component_value(const Theta4& theta, std::size_t mu, std::size_t nu) {
  if (mu == 0 && nu == 0) return theta.tt;
  if (mu == 0 && nu == 1) return theta.tx;
  if (mu == 0 && nu == 2) return theta.ty;
  if (mu == 0 && nu == 3) return theta.tz;
  if (mu == 1 && nu == 0) return theta.xt;
  if (mu == 1 && nu == 1) return theta.xx;
  if (mu == 1 && nu == 2) return theta.xy;
  if (mu == 1 && nu == 3) return theta.xz;
  if (mu == 2 && nu == 0) return theta.yt;
  if (mu == 2 && nu == 1) return theta.yx;
  if (mu == 2 && nu == 2) return theta.yy;
  if (mu == 2 && nu == 3) return theta.yz;
  if (mu == 3 && nu == 0) return theta.zt;
  if (mu == 3 && nu == 1) return theta.zx;
  if (mu == 3 && nu == 2) return theta.zy;
  if (mu == 3 && nu == 3) return theta.zz;
  throw std::out_of_range("Theta4 index out of range");
}

double& theta_component_ref(Theta4& theta, std::size_t mu, std::size_t nu) {
  if (mu == 0 && nu == 0) return theta.tt;
  if (mu == 0 && nu == 1) return theta.tx;
  if (mu == 0 && nu == 2) return theta.ty;
  if (mu == 0 && nu == 3) return theta.tz;
  if (mu == 1 && nu == 0) return theta.xt;
  if (mu == 1 && nu == 1) return theta.xx;
  if (mu == 1 && nu == 2) return theta.xy;
  if (mu == 1 && nu == 3) return theta.xz;
  if (mu == 2 && nu == 0) return theta.yt;
  if (mu == 2 && nu == 1) return theta.yx;
  if (mu == 2 && nu == 2) return theta.yy;
  if (mu == 2 && nu == 3) return theta.yz;
  if (mu == 3 && nu == 0) return theta.zt;
  if (mu == 3 && nu == 1) return theta.zx;
  if (mu == 3 && nu == 2) return theta.zy;
  if (mu == 3 && nu == 3) return theta.zz;
  throw std::out_of_range("Theta4 index out of range");
}

}  // namespace

double& Xi4::component(std::size_t mu) {
  return xi_component_ref(*this, mu);
}

double Xi4::component(std::size_t mu) const {
  switch (mu) {
    case 0: return t;
    case 1: return x;
    case 2: return y;
    case 3: return z;
    default: throw std::out_of_range("Xi4 index out of range");
  }
}

double& Theta4::component(std::size_t mu, std::size_t nu) {
  return theta_component_ref(*this, mu, nu);
}

double Theta4::component(std::size_t mu, std::size_t nu) const {
  return theta_component_value(*this, mu, nu);
}

double metric_signature_sign(std::size_t mu) {
  if (mu == 0) return -1.0;
  if (mu < 4) return 1.0;
  throw std::out_of_range("metric index out of range");
}

double trace_contraction_4d(const Theta4& theta) {
  return metric_signature_sign(0) * theta.tt + metric_signature_sign(1) * theta.xx +
         metric_signature_sign(2) * theta.yy + metric_signature_sign(3) * theta.zz;
}

double theta_double_contraction_4d(const Theta4& theta) {
  double out = 0.0;
  for (std::size_t mu = 0; mu < 4; ++mu) {
    for (std::size_t nu = 0; nu < 4; ++nu) {
      const double sign = metric_signature_sign(mu) * metric_signature_sign(nu);
      const double value = theta.component(mu, nu);
      out += sign * value * value;
    }
  }
  return out;
}

double invariant_I_4d(const Theta4& theta) {
  const double trace = trace_contraction_4d(theta);
  return theta_double_contraction_4d(theta) - LAMBDA_4D * trace * trace;
}

}  // namespace tpfcore
}  // namespace galaxy
