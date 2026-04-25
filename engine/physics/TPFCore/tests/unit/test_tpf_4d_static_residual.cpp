#include "doctest.h"
#include "physics/TPFCore/tpf_4d_static_residual.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

using galaxy::tpfcore::evaluate_static_configuration_residual_4d;
using galaxy::tpfcore::StaticResidualGridConfig;
using galaxy::tpfcore::StaticSourcePoint4D;

namespace {

StaticResidualGridConfig centered_grid(std::size_t n, double spacing, double exclusion_radius) {
  const double min_coord = -0.5 * spacing * static_cast<double>(n - 1);
  StaticResidualGridConfig cfg{};
  cfg.nx = n;
  cfg.ny = n;
  cfg.nz = n;
  cfg.spacing = spacing;
  cfg.origin_x = min_coord;
  cfg.origin_y = min_coord;
  cfg.origin_z = min_coord;
  cfg.field_softening_eps = 1e-3;
  cfg.source_exclusion_radius = exclusion_radius;
  return cfg;
}

}  // namespace

TEST_CASE("4D static residual monopole summary is finite and near-core residuals are larger") {
  const std::vector<StaticSourcePoint4D> sources{StaticSourcePoint4D{1.0, 0.0, 0.0, 0.0}};

  const StaticResidualGridConfig no_excl_cfg = centered_grid(11, 0.2, 0.0);
  const StaticResidualGridConfig excl_cfg = centered_grid(11, 0.2, 0.35);

  const galaxy::tpfcore::StaticResidualGridResult no_excl =
      evaluate_static_configuration_residual_4d(sources, no_excl_cfg);
  const galaxy::tpfcore::StaticResidualGridResult excl =
      evaluate_static_configuration_residual_4d(sources, excl_cfg);

  CHECK(std::isfinite(excl.summary.max_residual_spatial_norm));
  CHECK(std::isfinite(excl.summary.mean_residual_spatial_norm));
  CHECK(std::isfinite(excl.summary.median_residual_spatial_norm));
  CHECK(std::isfinite(excl.summary.max_normalized_residual));
  CHECK(std::isfinite(excl.summary.mean_normalized_residual));

  CHECK(excl.summary.free_cells_used > 0);
  CHECK(no_excl.summary.max_residual_spatial_norm > excl.summary.max_residual_spatial_norm);
}

TEST_CASE("4D static residual max decreases as source exclusion radius increases") {
  const std::vector<StaticSourcePoint4D> sources{StaticSourcePoint4D{1.0, 0.0, 0.0, 0.0}};

  StaticResidualGridConfig small_excl = centered_grid(11, 0.2, 0.15);
  StaticResidualGridConfig large_excl = centered_grid(11, 0.2, 0.45);

  const galaxy::tpfcore::StaticResidualGridResult small =
      evaluate_static_configuration_residual_4d(sources, small_excl);
  const galaxy::tpfcore::StaticResidualGridResult large =
      evaluate_static_configuration_residual_4d(sources, large_excl);

  CHECK(small.summary.free_cells_used > 0);
  CHECK(large.summary.free_cells_used > 0);
  CHECK(large.summary.max_residual_spatial_norm < small.summary.max_residual_spatial_norm);
}

TEST_CASE("4D static residual includes z-derivative sampling") {
  const std::vector<StaticSourcePoint4D> sources{StaticSourcePoint4D{1.0, 0.2, -0.1, 0.35}};

  StaticResidualGridConfig cfg = centered_grid(9, 0.25, 0.1);
  const galaxy::tpfcore::StaticResidualGridResult out =
      evaluate_static_configuration_residual_4d(sources, cfg);

  CHECK(out.summary.z_derivative_samples == out.summary.free_cells_used);
  CHECK(out.summary.z_derivative_samples > 0);

  bool found_nonzero_residual_z = false;
  for (std::size_t i = 0; i < out.points.size(); ++i) {
    if (out.points[i].used_in_summary && std::fabs(out.points[i].residual_z) > 1e-12) {
      found_nonzero_residual_z = true;
      break;
    }
  }
  CHECK(found_nonzero_residual_z);
}

TEST_CASE("4D static residual boundary cells are excluded from summary") {
  const std::vector<StaticSourcePoint4D> sources{StaticSourcePoint4D{1.0, 0.0, 0.0, 0.0}};

  StaticResidualGridConfig cfg = centered_grid(7, 0.3, 0.2);
  const galaxy::tpfcore::StaticResidualGridResult out =
      evaluate_static_configuration_residual_4d(sources, cfg);

  const std::size_t total = cfg.nx * cfg.ny * cfg.nz;
  const std::size_t interior = (cfg.nx - 2) * (cfg.ny - 2) * (cfg.nz - 2);

  CHECK(out.summary.total_grid_cells == total);
  CHECK(out.summary.interior_cells == interior);
  CHECK(out.summary.excluded_boundary_count == total - interior);
}

TEST_CASE("4D static residual multi-source smoke test") {
  const std::vector<StaticSourcePoint4D> sources{
      StaticSourcePoint4D{1.0, -0.4, 0.0, 0.2},
      StaticSourcePoint4D{0.8, 0.5, -0.3, -0.25},
  };

  StaticResidualGridConfig cfg = centered_grid(11, 0.2, 0.2);
  const galaxy::tpfcore::StaticResidualGridResult out =
      evaluate_static_configuration_residual_4d(sources, cfg);

  CHECK(out.summary.free_cells_used > 10);
  CHECK(out.summary.free_cells_used <= out.summary.interior_cells);
  CHECK(std::isfinite(out.summary.max_residual_spatial_norm));
  CHECK(std::isfinite(out.summary.mean_residual_spatial_norm));
  CHECK(std::isfinite(out.summary.max_normalized_residual));
}
