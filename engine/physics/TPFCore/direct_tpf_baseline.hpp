#ifndef GALAXY_PHYSICS_TPFCORE_DIRECT_TPF_BASELINE_HPP
#define GALAXY_PHYSICS_TPFCORE_DIRECT_TPF_BASELINE_HPP

#include "field_evaluation.hpp"

namespace galaxy {
namespace tpfcore {

/** Eq. (10) symbol mapping for direct_tpf baseline (DeltaC omitted placeholder). */
struct Eq10XiDisplacement {
  Xi2D xi;
};

struct Eq10ThetaTensor {
  Theta3D theta;
};

struct Eq10ThetaTrace {
  double value;
};

struct Eq10InvariantI {
  double value;
};

/** Direct baseline keeps DeltaC explicit but zero in current implementation scope. */
struct Eq10DeltaCPlaceholder {
  double xx;
  double xy;
  double yy;
  double zz;
};

/** Principal-part C_ij used by Xi-directed operational readout. */
struct Eq10PrincipalCij {
  double c_xx;
  double c_xy;
  double c_yy;
  double c_zz;
};

struct XiDirectedReadoutResult {
  double ax;
  double ay;
};

struct DirectTpfBaselineAccelerationResult {
  double ax;
  double ay;
};

struct DirectTpfBaselineArtifacts {
  Eq10XiDisplacement xi;
  Eq10ThetaTensor theta;
  Eq10ThetaTrace theta_trace;
  Eq10InvariantI invariant_I;
  Eq10DeltaCPlaceholder delta_c;
  Eq10PrincipalCij principal_cij;
};

Eq10ThetaTrace compute_theta_trace(const Eq10ThetaTensor& theta_tensor);
Eq10InvariantI compute_invariant_I(const Eq10ThetaTensor& theta_tensor);
Eq10DeltaCPlaceholder compute_deltaC_placeholder_zero();
Eq10PrincipalCij compute_principal_Cij_from_eq10_baseline(const Eq10ThetaTensor& theta_tensor,
                                                          const Eq10ThetaTrace& theta_trace,
                                                          const Eq10InvariantI& invariant_I,
                                                          const Eq10DeltaCPlaceholder& delta_c,
                                                          double kappa,
                                                          double lambda);
XiDirectedReadoutResult compute_xi_directed_tensor_readout(const Eq10XiDisplacement& xi,
                                                           const Eq10PrincipalCij& principal_cij);
DirectTpfBaselineArtifacts compute_direct_tpf_baseline_artifacts(const FieldAtPoint& field,
                                                                 double kappa,
                                                                 double lambda);
DirectTpfBaselineArtifacts compute_direct_tpf_baseline_artifacts(const CanonicalFieldObjects& field,
                                                                 double kappa,
                                                                 double lambda);
DirectTpfBaselineAccelerationResult compute_direct_tpf_baseline_acceleration(const FieldAtPoint& field,
                                                                             double kappa,
                                                                             double lambda);
DirectTpfBaselineAccelerationResult compute_direct_tpf_baseline_acceleration(const CanonicalFieldObjects& field,
                                                                             double kappa,
                                                                             double lambda);

}  // namespace tpfcore
}  // namespace galaxy

#endif
