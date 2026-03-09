# TPFCore: Honest primitive TPF structure package

**This is NOT the removed weak-field Newtonian-like TPF package.** TPFCore implements the bare primitive structure from the paper. The previous isotropic placeholder ansatz has been **replaced** with a Hessian-based weak-field point-source ansatz.

## Source ansatz (provisional weak-field)

The implementation uses the paper's weak-field point-source construction:

- **Phi** = -M / sqrt(r² + eps²) — softened point-source scalar (numerical regularization)
- **Xi_i** = partial_i Phi — displacement field
- **Theta_ij** = partial_i partial_j Phi — configuration gradient (Hessian of Phi)

r² = dx² + dy², dx = x - xs, dy = y - ys. eps is the source softening (`tpfcore_source_softening` or global `softening`).

This is a **provisional weak-field point-source ansatz**, NOT the full nonlinear TPF source law.

## What is directly paper-derived

- **Xi^mu** — configuration displacement field
- **Theta_{mu nu} = nabla_mu Xi_nu** — configuration gradient tensor
- **I = Theta_{mu nu} Theta^{mu nu} - lambda Theta^2** — scalar invariant
- **lambda = 1/4** in 4D (fixed)

## Where the ansatz lives

- **`source_ansatz.hpp`** — Declaration, `PointSourceField`, `provisional_point_source_field`
- **`source_ansatz.cpp`** — Closed-form Phi, Xi, Theta derivatives

## CSV output columns

**theta_profile.csv** (single source: no axis column; symmetric pair: `axis` = x or y)

| Column | Meaning |
|--------|---------|
| axis | Probe axis (x or y); omitted for single-source |
| radius | Distance along probe axis |
| x, y | Field point coordinates |
| xi_x, xi_y | Displacement field components |
| theta_xx, theta_xy, theta_yy | Theta tensor components |
| theta_trace | Theta_xx + Theta_yy |
| invariant_I | I = Theta_mn Theta^mn - λ Theta² |

**invariant_profile.csv**

| Column | Meaning |
|--------|---------|
| radius | (or axis, radius for symmetric pair) |
| invariant_I | Scalar invariant |
| theta_trace | Trace of Theta |

## Symmetry behavior on +x probe line

For **single-source** at origin, probing along +x (y=0):

- **theta_xy ≈ 0** — by symmetry (y=0 is a principal axis)
- **theta_xx ≠ theta_yy** — radial vs transverse structure (anisotropy)
- **invariant_I** — decays smoothly with radius

`field_summary.txt` reports whether these expectations hold numerically.

## Config (defaults.cfg)

- `tpfcore_enable_provisional_readout` — false (dynamics not available)
- `tpfcore_probe_radius_min`, `tpfcore_probe_radius_max`, `tpfcore_probe_samples`
- `tpfcore_dump_theta_profile`, `tpfcore_dump_invariant_profile`
- `tpfcore_source_softening` — softening for Phi. If ≤ 0, use global `softening`.

## How to run

```bash
# In config: physics_package = TPFCore

# Single source at origin, probe along +x
./galaxy_sim tpf_single_source_inspect

# Symmetric pair at (±d,0), probe along +x and +y
./galaxy_sim tpf_symmetric_pair_inspect
```

## What is NOT implemented

- Full nonlinear/dynamic TPF field equations
- Galaxy-scale dynamics
- Acceleration readout (no -grad(phi) dynamics)
- Newtonian-style substitution

## Warning

TPFCore is inspection-first, not a completed galaxy-dynamics solver. The ansatz is provisional where the paper remains underspecified.
