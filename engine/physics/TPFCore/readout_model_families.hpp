#ifndef GALAXY_PHYSICS_TPFCORE_READOUT_MODEL_FAMILIES_HPP
#define GALAXY_PHYSICS_TPFCORE_READOUT_MODEL_FAMILIES_HPP

#include "derived_tpf_radial.hpp"
#include "../../types.hpp"

namespace galaxy {
namespace tpfcore {

struct ReadoutDiagnostics;

namespace readout_models {

void apply_tensor_radial_projection_readout(const State& state,
                                            int i,
                                            double bh_mass,
                                            bool star_star,
                                            double eps,
                                            double readout_scale,
                                            bool negated,
                                            double& ax,
                                            double& ay,
                                            ReadoutDiagnostics* diag);

void apply_derived_radial_readout(const State& state,
                                  int i,
                                  double bh_mass,
                                  double eps,
                                  const DerivedTpfPoissonConfig& dcfg,
                                  const TpfRadialGravityProfile* profile_in,
                                  double readout_scale,
                                  double theta_tt_scale,
                                  double theta_tr_scale,
                                  double& ax,
                                  double& ay,
                                  ReadoutDiagnostics* diag);

void apply_experimental_radial_scaling_readout(const State& state,
                                               int i,
                                               double bh_mass,
                                               bool star_star,
                                               double eps,
                                               double readout_scale,
                                               double& ax,
                                               double& ay,
                                               ReadoutDiagnostics* diag);

}  // namespace readout_models
}  // namespace tpfcore
}  // namespace galaxy

#endif
