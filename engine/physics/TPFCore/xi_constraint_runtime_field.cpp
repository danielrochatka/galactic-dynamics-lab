#include "xi_constraint_runtime_field.hpp"

#include <cmath>

namespace galaxy {
namespace tpfcore {
namespace {

struct PlanarNodeSample {
  double xi_x = 0.0;
  double xi_y = 0.0;
  double theta_xx = 0.0;
  double theta_xy = 0.0;
  double theta_yx = 0.0;
  double theta_yy = 0.0;
};

double node_value(const std::vector<double>& a, const PlanarXiGrid& g, int i, int j) {
  return a[g.index(i, j)];
}

double derivative_x_at_node(const std::vector<double>& a, const PlanarXiGrid& g, int i, int j) {
  if (g.dx <= 0.0) return 0.0;
  if (i > 0 && i + 1 < g.nx) {
    return (node_value(a, g, i + 1, j) - node_value(a, g, i - 1, j)) / (2.0 * g.dx);
  }
  if (i + 1 < g.nx) {
    return (node_value(a, g, i + 1, j) - node_value(a, g, i, j)) / g.dx;
  }
  if (i > 0) {
    return (node_value(a, g, i, j) - node_value(a, g, i - 1, j)) / g.dx;
  }
  return 0.0;
}

double derivative_y_at_node(const std::vector<double>& a, const PlanarXiGrid& g, int i, int j) {
  if (g.dx <= 0.0) return 0.0;
  if (j > 0 && j + 1 < g.ny) {
    return (node_value(a, g, i, j + 1) - node_value(a, g, i, j - 1)) / (2.0 * g.dx);
  }
  if (j + 1 < g.ny) {
    return (node_value(a, g, i, j + 1) - node_value(a, g, i, j)) / g.dx;
  }
  if (j > 0) {
    return (node_value(a, g, i, j) - node_value(a, g, i, j - 1)) / g.dx;
  }
  return 0.0;
}

PlanarNodeSample build_node_sample(const PlanarXiGrid& grid, int i, int j) {
  PlanarNodeSample s;
  s.xi_x = grid.xi_x[grid.index(i, j)];
  s.xi_y = grid.xi_y[grid.index(i, j)];
  s.theta_xx = derivative_x_at_node(grid.xi_x, grid, i, j);
  s.theta_xy = derivative_x_at_node(grid.xi_y, grid, i, j);
  s.theta_yx = derivative_y_at_node(grid.xi_x, grid, i, j);
  s.theta_yy = derivative_y_at_node(grid.xi_y, grid, i, j);
  return s;
}

double bilerp(double q00, double q10, double q01, double q11, double tx, double ty) {
  const double a = (1.0 - tx) * q00 + tx * q10;
  const double b = (1.0 - tx) * q01 + tx * q11;
  return (1.0 - ty) * a + ty * b;
}

}  // namespace

PlanarXiSample2D sample_planar_xi_runtime_field_bilinear(const PlanarXiGrid& grid, double x, double y) {
  PlanarXiSample2D out;
  if (grid.nx < 2 || grid.ny < 2 || grid.dx <= 0.0) return out;

  const double u = (x + grid.L) / grid.dx;
  const double v = (y + grid.L) / grid.dx;
  int i0 = static_cast<int>(std::floor(u));
  int j0 = static_cast<int>(std::floor(v));
  if (i0 < 0) i0 = 0;
  if (j0 < 0) j0 = 0;
  if (i0 > grid.nx - 2) i0 = grid.nx - 2;
  if (j0 > grid.ny - 2) j0 = grid.ny - 2;
  const int i1 = i0 + 1;
  const int j1 = j0 + 1;

  double tx = u - static_cast<double>(i0);
  double ty = v - static_cast<double>(j0);
  if (tx < 0.0) tx = 0.0;
  if (tx > 1.0) tx = 1.0;
  if (ty < 0.0) ty = 0.0;
  if (ty > 1.0) ty = 1.0;

  const PlanarNodeSample s00 = build_node_sample(grid, i0, j0);
  const PlanarNodeSample s10 = build_node_sample(grid, i1, j0);
  const PlanarNodeSample s01 = build_node_sample(grid, i0, j1);
  const PlanarNodeSample s11 = build_node_sample(grid, i1, j1);

  out.xi_x = bilerp(s00.xi_x, s10.xi_x, s01.xi_x, s11.xi_x, tx, ty);
  out.xi_y = bilerp(s00.xi_y, s10.xi_y, s01.xi_y, s11.xi_y, tx, ty);
  out.theta_xx = bilerp(s00.theta_xx, s10.theta_xx, s01.theta_xx, s11.theta_xx, tx, ty);
  out.theta_xy = bilerp(s00.theta_xy, s10.theta_xy, s01.theta_xy, s11.theta_xy, tx, ty);
  out.theta_yx = bilerp(s00.theta_yx, s10.theta_yx, s01.theta_yx, s11.theta_yx, tx, ty);
  out.theta_yy = bilerp(s00.theta_yy, s10.theta_yy, s01.theta_yy, s11.theta_yy, tx, ty);
  return out;
}

PlanarEq10UnsymBaseline2D compute_planar_eq10_unsym_baseline_2d(const PlanarXiSample2D& s,
                                                                 double kappa,
                                                                 double lambda) {
  PlanarEq10UnsymBaseline2D out;
  out.theta_trace_2d = s.theta_xx + s.theta_yy;
  out.invariant_I_2d = s.theta_xx * s.theta_xx + s.theta_xy * s.theta_xy + s.theta_yx * s.theta_yx +
                       s.theta_yy * s.theta_yy - lambda * out.theta_trace_2d * out.theta_trace_2d;

  const double b_xx = s.theta_xx * s.theta_xx + s.theta_yx * s.theta_xy;
  const double b_xy = s.theta_xx * s.theta_yx + s.theta_yx * s.theta_yy;
  const double b_yx = s.theta_xy * s.theta_xx + s.theta_yy * s.theta_xy;
  const double b_yy = s.theta_xy * s.theta_yx + s.theta_yy * s.theta_yy;

  out.c_xx = kappa * (b_xx - lambda * out.theta_trace_2d * s.theta_xx - 0.5 * out.invariant_I_2d);
  out.c_xy = kappa * (b_xy - lambda * out.theta_trace_2d * s.theta_xy);
  out.c_yx = kappa * (b_yx - lambda * out.theta_trace_2d * s.theta_yx);
  out.c_yy = kappa * (b_yy - lambda * out.theta_trace_2d * s.theta_yy - 0.5 * out.invariant_I_2d);
  return out;
}

PlanarEq10UnsymReadout2D compute_planar_eq10_unsym_readout_2d(const PlanarXiSample2D& sample,
                                                              const PlanarEq10UnsymBaseline2D& baseline,
                                                              double small_epsilon) {
  PlanarEq10UnsymReadout2D out;
  const double xi_norm = std::sqrt(sample.xi_x * sample.xi_x + sample.xi_y * sample.xi_y);
  if (xi_norm <= small_epsilon) return out;

  const double ux = sample.xi_x / xi_norm;
  const double uy = sample.xi_y / xi_norm;

  out.ax = -(baseline.c_xx * ux + baseline.c_xy * uy);
  out.ay = -(baseline.c_yx * ux + baseline.c_yy * uy);
  return out;
}

}  // namespace tpfcore
}  // namespace galaxy
