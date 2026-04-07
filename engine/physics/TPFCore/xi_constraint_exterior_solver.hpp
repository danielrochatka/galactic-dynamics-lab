#ifndef GALAXY_PHYSICS_TPFCORE_XI_CONSTRAINT_EXTERIOR_SOLVER_HPP
#define GALAXY_PHYSICS_TPFCORE_XI_CONSTRAINT_EXTERIOR_SOLVER_HPP

#include "source_ansatz.hpp"

#include <vector>

namespace galaxy {
namespace tpfcore {

/** Experimental 2D planar Xi grid for exterior configuration-equation benchmark solves. */
struct PlanarXiGrid {
  int nx = 0;
  int ny = 0;
  double L = 0.0;
  double dx = 0.0;
  std::vector<double> xi_x;
  std::vector<double> xi_y;
  std::vector<unsigned char> is_exterior;
  std::vector<unsigned char> is_pinned;

  int index(int i, int j) const { return j * nx + i; }
  bool in_bounds(int i, int j) const { return i >= 0 && i < nx && j >= 0 && j < ny; }
  double x_at(int i) const { return -L + dx * static_cast<double>(i); }
  double y_at(int j) const { return -L + dx * static_cast<double>(j); }
};

/** Unsymmetrized planar gradient Theta_{mu nu} = nabla_mu Xi_nu in 2D. */
struct PlanarXiGradientUnsym2D {
  double xx = 0.0;
  double xy = 0.0;
  double yx = 0.0;
  double yy = 0.0;
};

/** Planar configuration-equation residual R = (Rx, Ry). */
struct PlanarConstraintResidual2D {
  double Rx = 0.0;
  double Ry = 0.0;
  double norm = 0.0;
};

/** Parameters for experimental exterior Xi constraint solve on z=0 plane. */
struct XiConstraintExteriorParams {
  int n = 65;
  double L = 10.0;
  double r_inner = 0.5;
  double source_mass = 1.0;
  double softening = 0.1;
  int max_iterations = 250;
  double tolerance = 1e-8;
};

/** Solver result bundle for diagnostics and tests. */
struct XiConstraintExteriorSolveResult {
  PlanarXiGrid grid;
  std::vector<PlanarConstraintResidual2D> residuals;
  double initial_max_residual_norm = 0.0;
  double final_max_residual_norm = 0.0;
  int iterations = 0;
  bool converged = false;
};

PlanarXiGrid initialize_planar_xi_grid_from_ansatz(const XiConstraintExteriorParams& params);
PlanarXiGradientUnsym2D compute_planar_xi_gradient_unsym2d(const PlanarXiGrid& grid, int i, int j);
std::vector<PlanarConstraintResidual2D> compute_planar_constraint_residual_field(const PlanarXiGrid& grid);
double compute_max_residual_norm_exterior(const PlanarXiGrid& grid,
                                          const std::vector<PlanarConstraintResidual2D>& residuals);
double compute_max_abs_xi_difference_vs_ansatz(const PlanarXiGrid& grid,
                                                const XiConstraintExteriorParams& params);
XiConstraintExteriorSolveResult solve_xi_constraint_exterior(const XiConstraintExteriorParams& params);

}  // namespace tpfcore
}  // namespace galaxy

#endif
