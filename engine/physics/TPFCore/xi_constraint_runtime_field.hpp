#ifndef GALAXY_PHYSICS_TPFCORE_XI_CONSTRAINT_RUNTIME_FIELD_HPP
#define GALAXY_PHYSICS_TPFCORE_XI_CONSTRAINT_RUNTIME_FIELD_HPP

#include "xi_constraint_exterior_solver.hpp"

namespace galaxy {
namespace tpfcore {

/** Experimental planar runtime sample from solved Xi grid. */
struct PlanarXiSample2D {
  double xi_x = 0.0;
  double xi_y = 0.0;
  double theta_xx = 0.0;
  double theta_xy = 0.0;
  double theta_yx = 0.0;
  double theta_yy = 0.0;
};

/** Experimental planar Eq.(10)-style unsym principal-part baseline (DeltaC omitted). */
struct PlanarEq10UnsymBaseline2D {
  double theta_trace_2d = 0.0;
  double invariant_I_2d = 0.0;
  double c_xx = 0.0;
  double c_xy = 0.0;
  double c_yx = 0.0;
  double c_yy = 0.0;
};

/** Experimental Xi-directed readout from unsym planar C_ij. */
struct PlanarEq10UnsymReadout2D {
  double ax = 0.0;
  double ay = 0.0;
};

PlanarXiSample2D sample_planar_xi_runtime_field_bilinear(const PlanarXiGrid& grid, double x, double y);

PlanarEq10UnsymBaseline2D compute_planar_eq10_unsym_baseline_2d(const PlanarXiSample2D& sample,
                                                                 double kappa,
                                                                 double lambda);

PlanarEq10UnsymReadout2D compute_planar_eq10_unsym_readout_2d(const PlanarXiSample2D& sample,
                                                              const PlanarEq10UnsymBaseline2D& baseline,
                                                              double small_epsilon);

}  // namespace tpfcore
}  // namespace galaxy

#endif
