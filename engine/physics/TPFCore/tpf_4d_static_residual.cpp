#include "physics/TPFCore/tpf_4d_static_residual.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "physics/TPFCore/source_ansatz.hpp"
#include "physics/TPFCore/tpf_4d_field.hpp"

namespace galaxy {
namespace tpfcore {

namespace {

struct ResidualTensorA {
  double component[4][4];
};

std::size_t flatten_index(const StaticResidualGridConfig& cfg, std::size_t i, std::size_t j, std::size_t k) {
  return (k * cfg.ny + j) * cfg.nx + i;
}

bool is_boundary_cell(const StaticResidualGridConfig& cfg, std::size_t i, std::size_t j, std::size_t k) {
  return i == 0 || j == 0 || k == 0 || i + 1 == cfg.nx || j + 1 == cfg.ny || k + 1 == cfg.nz;
}

double distance_to_source(const StaticSourcePoint4D& source, double x, double y, double z) {
  const double dx = x - source.x;
  const double dy = y - source.y;
  const double dz = z - source.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

bool is_near_any_source(const std::vector<StaticSourcePoint4D>& sources,
                        double x,
                        double y,
                        double z,
                        double exclusion_radius) {
  if (exclusion_radius <= 0.0) {
    return false;
  }
  for (std::vector<StaticSourcePoint4D>::const_iterator it = sources.begin(); it != sources.end(); ++it) {
    if (distance_to_source(*it, x, y, z) <= exclusion_radius) {
      return true;
    }
  }
  return false;
}

ResidualTensorA compute_A_mixed(const Field4DAtPoint& field) {
  ResidualTensorA out{};
  const double theta_trace = trace_contraction_4d(field.theta);
  for (std::size_t mu = 0; mu < 4; ++mu) {
    for (std::size_t nu = 0; nu < 4; ++nu) {
      const double theta_mixed = metric_signature_sign(mu) * field.theta.component(mu, nu);
      const double delta = (mu == nu) ? 1.0 : 0.0;
      out.component[mu][nu] = theta_mixed - LAMBDA_4D * delta * theta_trace;
    }
  }
  return out;
}

double theta_spatial_frobenius_norm(const Field4DAtPoint& field) {
  double s = 0.0;
  for (std::size_t i = 1; i < 4; ++i) {
    for (std::size_t j = 1; j < 4; ++j) {
      const double v = field.theta.component(i, j);
      s += v * v;
    }
  }
  return std::sqrt(s);
}

StaticResidualSummary summarize_points(const std::vector<StaticResidualAtPoint>& points,
                                       const StaticResidualGridConfig& config,
                                       std::size_t z_derivative_samples) {
  StaticResidualSummary summary{};
  summary.total_grid_cells = config.nx * config.ny * config.nz;
  summary.interior_cells = 0;
  if (config.nx >= 3 && config.ny >= 3 && config.nz >= 3) {
    summary.interior_cells = (config.nx - 2) * (config.ny - 2) * (config.nz - 2);
  }

  summary.excluded_boundary_count = summary.total_grid_cells - summary.interior_cells;
  summary.excluded_near_source_count = 0;
  summary.free_cells_used = 0;
  summary.z_derivative_samples = z_derivative_samples;

  std::vector<double> spatial_norms;
  double sum_spatial = 0.0;
  double sum_normed = 0.0;
  summary.max_residual_spatial_norm = 0.0;
  summary.max_normalized_residual = 0.0;

  for (std::vector<StaticResidualAtPoint>::const_iterator it = points.begin(); it != points.end(); ++it) {
    if (it->is_near_source && !it->is_boundary) {
      ++summary.excluded_near_source_count;
    }
    if (!it->used_in_summary) {
      continue;
    }
    ++summary.free_cells_used;
    spatial_norms.push_back(it->residual_spatial_norm);
    sum_spatial += it->residual_spatial_norm;
    sum_normed += it->normalized_residual;
    if (it->residual_spatial_norm > summary.max_residual_spatial_norm) {
      summary.max_residual_spatial_norm = it->residual_spatial_norm;
    }
    if (it->normalized_residual > summary.max_normalized_residual) {
      summary.max_normalized_residual = it->normalized_residual;
    }
  }

  if (summary.free_cells_used == 0) {
    summary.mean_residual_spatial_norm = 0.0;
    summary.median_residual_spatial_norm = 0.0;
    summary.mean_normalized_residual = 0.0;
    return summary;
  }

  summary.mean_residual_spatial_norm = sum_spatial / static_cast<double>(summary.free_cells_used);
  summary.mean_normalized_residual = sum_normed / static_cast<double>(summary.free_cells_used);

  std::sort(spatial_norms.begin(), spatial_norms.end());
  const std::size_t n = spatial_norms.size();
  if (n % 2 == 1) {
    summary.median_residual_spatial_norm = spatial_norms[n / 2];
  } else {
    summary.median_residual_spatial_norm = 0.5 * (spatial_norms[n / 2 - 1] + spatial_norms[n / 2]);
  }
  return summary;
}

}  // namespace

StaticResidualGridResult evaluate_static_configuration_residual_4d(
    const std::vector<StaticSourcePoint4D>& sources,
    const StaticResidualGridConfig& config) {
  if (config.nx < 3 || config.ny < 3 || config.nz < 3) {
    throw std::invalid_argument("static residual grid requires nx, ny, nz >= 3");
  }
  if (!(config.spacing > 0.0)) {
    throw std::invalid_argument("static residual grid spacing must be > 0");
  }

  const std::size_t total = config.nx * config.ny * config.nz;
  std::vector<Field4DAtPoint> fields(total);
  std::vector<ResidualTensorA> tensors_a(total);
  std::vector<StaticResidualAtPoint> points(total);

  for (std::size_t k = 0; k < config.nz; ++k) {
    for (std::size_t j = 0; j < config.ny; ++j) {
      for (std::size_t i = 0; i < config.nx; ++i) {
        const double x = config.origin_x + config.spacing * static_cast<double>(i);
        const double y = config.origin_y + config.spacing * static_cast<double>(j);
        const double z = config.origin_z + config.spacing * static_cast<double>(k);
        const std::size_t idx = flatten_index(config, i, j, k);

        const StaticProbePoint4D probe{x, y, z};
        fields[idx] = evaluate_static_sources_field_4d(sources, probe, config.field_softening_eps);
        tensors_a[idx] = compute_A_mixed(fields[idx]);

        StaticResidualAtPoint point{};
        point.i = i;
        point.j = j;
        point.k = k;
        point.x = x;
        point.y = y;
        point.z = z;
        point.is_boundary = is_boundary_cell(config, i, j, k);
        point.is_near_source = is_near_any_source(
            sources, x, y, z, config.source_exclusion_radius);
        point.used_in_summary = false;
        point.residual_t = 0.0;
        point.residual_x = 0.0;
        point.residual_y = 0.0;
        point.residual_z = 0.0;
        point.residual_spatial_norm = 0.0;
        point.residual_4_norm_like = 0.0;
        point.theta_spatial_frobenius_norm = theta_spatial_frobenius_norm(fields[idx]);
        point.normalized_residual = 0.0;
        points[idx] = point;
      }
    }
  }

  const double inv_2h = 1.0 / (2.0 * config.spacing);
  std::size_t z_derivative_samples = 0;

  for (std::size_t k = 1; k + 1 < config.nz; ++k) {
    for (std::size_t j = 1; j + 1 < config.ny; ++j) {
      for (std::size_t i = 1; i + 1 < config.nx; ++i) {
        const std::size_t idx = flatten_index(config, i, j, k);
        StaticResidualAtPoint& point = points[idx];

        if (point.is_near_source) {
          continue;
        }

        const std::size_t ip = flatten_index(config, i + 1, j, k);
        const std::size_t im = flatten_index(config, i - 1, j, k);
        const std::size_t jp = flatten_index(config, i, j + 1, k);
        const std::size_t jm = flatten_index(config, i, j - 1, k);
        const std::size_t kp = flatten_index(config, i, j, k + 1);
        const std::size_t km = flatten_index(config, i, j, k - 1);

        for (std::size_t nu = 0; nu < 4; ++nu) {
          const double d_dx = (tensors_a[ip].component[1][nu] - tensors_a[im].component[1][nu]) * inv_2h;
          const double d_dy = (tensors_a[jp].component[2][nu] - tensors_a[jm].component[2][nu]) * inv_2h;
          const double d_dz = (tensors_a[kp].component[3][nu] - tensors_a[km].component[3][nu]) * inv_2h;
          const double residual_nu = d_dx + d_dy + d_dz;

          if (nu == 0) point.residual_t = residual_nu;
          if (nu == 1) point.residual_x = residual_nu;
          if (nu == 2) point.residual_y = residual_nu;
          if (nu == 3) point.residual_z = residual_nu;
        }

        point.residual_spatial_norm = std::sqrt(point.residual_x * point.residual_x +
                                                point.residual_y * point.residual_y +
                                                point.residual_z * point.residual_z);
        point.residual_4_norm_like = std::sqrt(
            std::fabs(-point.residual_t * point.residual_t +
                      point.residual_x * point.residual_x +
                      point.residual_y * point.residual_y +
                      point.residual_z * point.residual_z));

        const double scale = std::max(point.theta_spatial_frobenius_norm, 1e-30);
        point.normalized_residual = point.residual_spatial_norm / scale;
        point.used_in_summary = true;
        ++z_derivative_samples;
      }
    }
  }

  StaticResidualGridResult result{};
  result.config = config;
  result.points.swap(points);
  result.summary = summarize_points(result.points, config, z_derivative_samples);
  return result;
}

}  // namespace tpfcore
}  // namespace galaxy
