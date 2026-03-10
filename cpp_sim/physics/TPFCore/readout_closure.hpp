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
 * - tr_coherence: superposed Theta at particle, then Theta_rr/Theta_tt/Theta_tr formula.
 *
 * The public API is in provisional_readout.hpp (compute_provisional_readout_*).
 * This header documents the closure boundary only; no new symbols required by callers.
 */

#endif
