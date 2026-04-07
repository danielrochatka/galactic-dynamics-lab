#include "xi_constraint_exterior_solver.hpp"

#include <algorithm>
#include <cmath>

namespace galaxy {
namespace tpfcore {
namespace {

inline bool has_full_cross_stencil(const PlanarXiGrid& grid, int i, int j) {
  static const int kDi[8] = {1, -1, 0, 0, 1, 1, -1, -1};
  static const int kDj[8] = {0, 0, 1, -1, 1, -1, 1, -1};
  for (int k = 0; k < 8; ++k) {
    const int ni = i + kDi[k];
    const int nj = j + kDj[k];
    if (!grid.in_bounds(ni, nj)) return false;
    const int nidx = grid.index(ni, nj);
    if (!grid.is_exterior[nidx]) return false;
  }
  return true;
}

inline Xi2D ansatz_xi_at(const XiConstraintExteriorParams& params, double x, double y) {
  return provisional_point_source_field(0.0, 0.0, params.source_mass, x, y, params.softening).xi;
}

}  // namespace

PlanarXiGrid initialize_planar_xi_grid_from_ansatz(const XiConstraintExteriorParams& params) {
  PlanarXiGrid g;
  g.nx = params.n;
  g.ny = params.n;
  g.L = params.L;
  g.dx = (g.nx > 1) ? (2.0 * g.L / static_cast<double>(g.nx - 1)) : 0.0;

  const int ncell = g.nx * g.ny;
  g.xi_x.assign(ncell, 0.0);
  g.xi_y.assign(ncell, 0.0);
  g.is_exterior.assign(ncell, 0);
  g.is_pinned.assign(ncell, 0);

  for (int j = 0; j < g.ny; ++j) {
    for (int i = 0; i < g.nx; ++i) {
      const int idx = g.index(i, j);
      const double x = g.x_at(i);
      const double y = g.y_at(j);
      const double r = std::sqrt(x * x + y * y);
      g.is_exterior[idx] = (r >= params.r_inner) ? 1 : 0;
      Xi2D xi = ansatz_xi_at(params, x, y);
      g.xi_x[idx] = xi.x;
      g.xi_y[idx] = xi.y;

      const bool outer = (i == 0 || j == 0 || i == g.nx - 1 || j == g.ny - 1);
      if (outer) g.is_pinned[idx] = 1;
    }
  }

  for (int j = 1; j + 1 < g.ny; ++j) {
    for (int i = 1; i + 1 < g.nx; ++i) {
      const int idx = g.index(i, j);
      if (!g.is_exterior[idx]) continue;
      const bool adjacent_to_excision =
          (!g.is_exterior[g.index(i + 1, j)] || !g.is_exterior[g.index(i - 1, j)] ||
           !g.is_exterior[g.index(i, j + 1)] || !g.is_exterior[g.index(i, j - 1)]);
      if (adjacent_to_excision) g.is_pinned[idx] = 1;
    }
  }

  return g;
}

PlanarXiGradientUnsym2D compute_planar_xi_gradient_unsym2d(const PlanarXiGrid& grid, int i, int j) {
  PlanarXiGradientUnsym2D grad;
  if (!has_full_cross_stencil(grid, i, j) || grid.dx <= 0.0) return grad;

  const double inv2dx = 1.0 / (2.0 * grid.dx);
  grad.xx = (grid.xi_x[grid.index(i + 1, j)] - grid.xi_x[grid.index(i - 1, j)]) * inv2dx;
  grad.xy = (grid.xi_y[grid.index(i + 1, j)] - grid.xi_y[grid.index(i - 1, j)]) * inv2dx;
  grad.yx = (grid.xi_x[grid.index(i, j + 1)] - grid.xi_x[grid.index(i, j - 1)]) * inv2dx;
  grad.yy = (grid.xi_y[grid.index(i, j + 1)] - grid.xi_y[grid.index(i, j - 1)]) * inv2dx;
  return grad;
}

std::vector<PlanarConstraintResidual2D> compute_planar_constraint_residual_field(const PlanarXiGrid& grid) {
  const int ncell = grid.nx * grid.ny;
  std::vector<PlanarConstraintResidual2D> out(ncell);
  if (grid.dx <= 0.0) return out;

  const double h2 = grid.dx * grid.dx;
  const double inv_h2 = 1.0 / h2;

  for (int j = 1; j + 1 < grid.ny; ++j) {
    for (int i = 1; i + 1 < grid.nx; ++i) {
      const int idx = grid.index(i, j);
      if (!grid.is_exterior[idx]) continue;
      if (!has_full_cross_stencil(grid, i, j)) continue;

      const double xip = grid.xi_x[grid.index(i + 1, j)];
      const double xim = grid.xi_x[grid.index(i - 1, j)];
      const double xjp = grid.xi_x[grid.index(i, j + 1)];
      const double xjm = grid.xi_x[grid.index(i, j - 1)];
      const double x0 = grid.xi_x[idx];

      const double yip = grid.xi_y[grid.index(i + 1, j)];
      const double yim = grid.xi_y[grid.index(i - 1, j)];
      const double yjp = grid.xi_y[grid.index(i, j + 1)];
      const double yjm = grid.xi_y[grid.index(i, j - 1)];
      const double y0 = grid.xi_y[idx];

      const double dxx_x = (xip - 2.0 * x0 + xim) * inv_h2;
      const double dyy_x = (xjp - 2.0 * x0 + xjm) * inv_h2;
      const double dxx_y = (yip - 2.0 * y0 + yim) * inv_h2;
      const double dyy_y = (yjp - 2.0 * y0 + yjm) * inv_h2;

      const double dxy_y = (grid.xi_y[grid.index(i + 1, j + 1)] - grid.xi_y[grid.index(i + 1, j - 1)] -
                            grid.xi_y[grid.index(i - 1, j + 1)] + grid.xi_y[grid.index(i - 1, j - 1)]) /
                           (4.0 * h2);
      const double dxy_x = (grid.xi_x[grid.index(i + 1, j + 1)] - grid.xi_x[grid.index(i + 1, j - 1)] -
                            grid.xi_x[grid.index(i - 1, j + 1)] + grid.xi_x[grid.index(i - 1, j - 1)]) /
                           (4.0 * h2);

      PlanarConstraintResidual2D r;
      r.Rx = (1.0 - LAMBDA_4D) * dxx_x + dyy_x - LAMBDA_4D * dxy_y;
      r.Ry = dxx_y + (1.0 - LAMBDA_4D) * dyy_y - LAMBDA_4D * dxy_x;
      r.norm = std::sqrt(r.Rx * r.Rx + r.Ry * r.Ry);
      out[idx] = r;
    }
  }

  return out;
}

double compute_max_residual_norm_exterior(const PlanarXiGrid& grid,
                                          const std::vector<PlanarConstraintResidual2D>& residuals) {
  double max_norm = 0.0;
  const int ncell = grid.nx * grid.ny;
  for (int idx = 0; idx < ncell; ++idx) {
    if (!grid.is_exterior[idx]) continue;
    max_norm = std::max(max_norm, residuals[idx].norm);
  }
  return max_norm;
}

double compute_max_abs_xi_difference_vs_ansatz(const PlanarXiGrid& grid,
                                                const XiConstraintExteriorParams& params) {
  double max_abs = 0.0;
  for (int j = 0; j < grid.ny; ++j) {
    for (int i = 0; i < grid.nx; ++i) {
      const int idx = grid.index(i, j);
      if (!grid.is_exterior[idx]) continue;
      const Xi2D a = ansatz_xi_at(params, grid.x_at(i), grid.y_at(j));
      max_abs = std::max(max_abs, std::fabs(grid.xi_x[idx] - a.x));
      max_abs = std::max(max_abs, std::fabs(grid.xi_y[idx] - a.y));
    }
  }
  return max_abs;
}

XiConstraintExteriorSolveResult solve_xi_constraint_exterior(const XiConstraintExteriorParams& params) {
  // Experimental benchmark solve: planar exterior configuration equation only.
  // This is intentionally isolated and is not the full Eq. (10) direct dynamics solve.
  XiConstraintExteriorSolveResult result;
  result.grid = initialize_planar_xi_grid_from_ansatz(params);

  std::vector<PlanarConstraintResidual2D> residual0 = compute_planar_constraint_residual_field(result.grid);
  result.initial_max_residual_norm = compute_max_residual_norm_exterior(result.grid, residual0);

  const double h2 = result.grid.dx * result.grid.dx;
  const double a = 1.0 - LAMBDA_4D;
  const double denom = 2.0 * (a + 1.0);

  std::vector<double> next_x = result.grid.xi_x;
  std::vector<double> next_y = result.grid.xi_y;

  for (int iter = 0; iter < params.max_iterations; ++iter) {
    for (int j = 1; j + 1 < result.grid.ny; ++j) {
      for (int i = 1; i + 1 < result.grid.nx; ++i) {
        const int idx = result.grid.index(i, j);
        if (!result.grid.is_exterior[idx] || result.grid.is_pinned[idx]) continue;
        if (!has_full_cross_stencil(result.grid, i, j)) continue;

        const double xip = result.grid.xi_x[result.grid.index(i + 1, j)];
        const double xim = result.grid.xi_x[result.grid.index(i - 1, j)];
        const double xjp = result.grid.xi_x[result.grid.index(i, j + 1)];
        const double xjm = result.grid.xi_x[result.grid.index(i, j - 1)];

        const double yip = result.grid.xi_y[result.grid.index(i + 1, j)];
        const double yim = result.grid.xi_y[result.grid.index(i - 1, j)];
        const double yjp = result.grid.xi_y[result.grid.index(i, j + 1)];
        const double yjm = result.grid.xi_y[result.grid.index(i, j - 1)];

        const double cross_y = (result.grid.xi_y[result.grid.index(i + 1, j + 1)] -
                                result.grid.xi_y[result.grid.index(i + 1, j - 1)] -
                                result.grid.xi_y[result.grid.index(i - 1, j + 1)] +
                                result.grid.xi_y[result.grid.index(i - 1, j - 1)]) /
                               4.0;
        const double cross_x = (result.grid.xi_x[result.grid.index(i + 1, j + 1)] -
                                result.grid.xi_x[result.grid.index(i + 1, j - 1)] -
                                result.grid.xi_x[result.grid.index(i - 1, j + 1)] +
                                result.grid.xi_x[result.grid.index(i - 1, j - 1)]) /
                               4.0;

        next_x[idx] = (a * (xip + xim) + (xjp + xjm) - LAMBDA_4D * cross_y) / denom;
        next_y[idx] = ((yip + yim) + a * (yjp + yjm) - LAMBDA_4D * cross_x) / denom;
      }
    }

    result.grid.xi_x.swap(next_x);
    result.grid.xi_y.swap(next_y);

    const auto residuals = compute_planar_constraint_residual_field(result.grid);
    const double max_res = compute_max_residual_norm_exterior(result.grid, residuals);

    result.iterations = iter + 1;
    if (max_res <= params.tolerance) {
      result.converged = true;
      break;
    }

    if (h2 <= 0.0) break;
  }

  result.residuals = compute_planar_constraint_residual_field(result.grid);
  result.final_max_residual_norm = compute_max_residual_norm_exterior(result.grid, result.residuals);
  return result;
}

}  // namespace tpfcore
}  // namespace galaxy
