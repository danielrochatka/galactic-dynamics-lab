#ifndef GALAXY_PHYSICS_TPFCORE_TPF_4D_STATIC_RESIDUAL_HPP
#define GALAXY_PHYSICS_TPFCORE_TPF_4D_STATIC_RESIDUAL_HPP

#include <cstddef>
#include <vector>

#include "physics/TPFCore/tpf_4d_static_field.hpp"

namespace galaxy {
namespace tpfcore {

struct StaticResidualGridConfig {
  std::size_t nx;
  std::size_t ny;
  std::size_t nz;
  double spacing;
  double origin_x;
  double origin_y;
  double origin_z;
  double field_softening_eps;
  double source_exclusion_radius;
};

struct StaticResidualAtPoint {
  std::size_t i;
  std::size_t j;
  std::size_t k;
  double x;
  double y;
  double z;

  bool is_boundary;
  bool is_near_source;
  bool used_in_summary;

  double residual_t;
  double residual_x;
  double residual_y;
  double residual_z;
  double residual_spatial_norm;
  double residual_4_norm_like;
  double theta_spatial_frobenius_norm;
  double normalized_residual;
};

struct StaticResidualSummary {
  std::size_t total_grid_cells;
  std::size_t interior_cells;
  std::size_t excluded_boundary_count;
  std::size_t excluded_near_source_count;
  std::size_t free_cells_used;

  std::size_t z_derivative_samples;

  double max_residual_spatial_norm;
  double mean_residual_spatial_norm;
  double median_residual_spatial_norm;
  double max_normalized_residual;
  double mean_normalized_residual;
};

struct StaticResidualGridResult {
  StaticResidualGridConfig config;
  std::vector<StaticResidualAtPoint> points;
  StaticResidualSummary summary;
};

StaticResidualGridResult evaluate_static_configuration_residual_4d(
    const std::vector<StaticSourcePoint4D>& sources,
    const StaticResidualGridConfig& config);

}  // namespace tpfcore
}  // namespace galaxy

#endif
