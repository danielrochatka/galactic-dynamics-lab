#ifndef GALAXY_PHYSICS_TPFCORE_READOUT_CLOSURE_HPP
#define GALAXY_PHYSICS_TPFCORE_READOUT_CLOSURE_HPP

/**
 * Provisional readout closure boundary.
 *
 * Readout closures are DOWNSTREAM of the source ansatz: they take Theta (or
 * per-source Theta) and produce a provisional acceleration. They are NOT the
 * source theory; they are exploratory motion/readout formulas.
 *
 * Implementations (in provisional_readout.cpp):
 * - tensor_radial: per-source Theta·r_hat, superposed (with optional negated).
 * - derived_tpf_radial_readout and tr_coherence_readout: same path — hybrid radial a_s r̂ from
 *   bounced baryons + κ–I ledger (see derived_tpf_radial.*). Theta_rr/Theta_tt/Theta_tr are diagnostics;
 *   tangential readout is not added to ax, ay.
 * - experimental_radial_r_scaling: superposed Theta; radial-only a = readout_scale*(-theta_rr)*r*r_hat.
 *
 * Full scope / tiers: TPF_PAPER_V11_SCOPE.md in this directory.
 *
 * The public API is in provisional_readout.hpp (compute_provisional_readout_*).
 * This header documents the closure boundary only; no new symbols required by callers.
 */

#endif
