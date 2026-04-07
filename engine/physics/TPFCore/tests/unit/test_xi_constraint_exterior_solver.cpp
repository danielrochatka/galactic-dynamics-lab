#include "doctest.h"
#include "physics/TPFCore/source_ansatz.hpp"
#include "physics/TPFCore/xi_constraint_exterior_solver.hpp"

#include <cmath>

TEST_CASE("xi exterior solver pins outer and excision-adjacent cells to ansatz") {
  galaxy::tpfcore::XiConstraintExteriorParams p;
  p.n = 51;
  p.L = 8.0;
  p.r_inner = 0.8;
  p.source_mass = 5.0;
  p.softening = 0.25;

  auto g = galaxy::tpfcore::initialize_planar_xi_grid_from_ansatz(p);

  for (int j = 0; j < g.ny; ++j) {
    for (int i = 0; i < g.nx; ++i) {
      const int idx = g.index(i, j);
      if (!g.is_exterior[idx] || !g.is_pinned[idx]) continue;
      auto pf = galaxy::tpfcore::provisional_point_source_field(0.0, 0.0, p.source_mass, g.x_at(i), g.y_at(j), p.softening);
      CHECK(g.xi_x[idx] == doctest::Approx(pf.xi.x));
      CHECK(g.xi_y[idx] == doctest::Approx(pf.xi.y));
    }
  }
}

TEST_CASE("xi exterior solver residual max norm does not increase") {
  galaxy::tpfcore::XiConstraintExteriorParams p;
  p.n = 65;
  p.L = 10.0;
  p.r_inner = 0.9;
  p.source_mass = 10.0;
  p.softening = 0.3;
  p.max_iterations = 100;
  p.tolerance = 1e-9;

  auto out = galaxy::tpfcore::solve_xi_constraint_exterior(p);
  CHECK(out.final_max_residual_norm <= out.initial_max_residual_norm + 1e-14);
}

TEST_CASE("xi exterior solver stays close to provisional ansatz under ansatz BC and initialization") {
  galaxy::tpfcore::XiConstraintExteriorParams p;
  p.n = 61;
  p.L = 12.0;
  p.r_inner = 1.0;
  p.source_mass = 7.0;
  p.softening = 0.4;
  p.max_iterations = 60;

  auto out = galaxy::tpfcore::solve_xi_constraint_exterior(p);
  double max_abs_diff = galaxy::tpfcore::compute_max_abs_xi_difference_vs_ansatz(out.grid, p);
  CHECK(max_abs_diff < 1.0);
}
