# TPFCore: Honest primitive TPF structure package

**This is NOT the removed weak-field Newtonian-like TPF package.** TPFCore implements the bare primitive structure from the paper. The previous isotropic placeholder ansatz has been **replaced** with a Hessian-based weak-field point-source ansatz.

## Source ansatz (provisional weak-field)

The implementation uses the paper's weak-field point-source construction:

- **Phi** = -M / sqrt(r² + eps²) — softened point-source scalar (numerical regularization)
- **Xi_i** = partial_i Phi — displacement field (unchanged)
- **Theta_ij** = Hess_ij(Phi) + B(r) δ_ij — configuration gradient with optional isotropic correction
- **B(r)** = c·M / (r² + eps²)^(3/2) — exploratory isotropic correction (c configurable, default 0)

r² = dx² + dy², dx = x - xs, dy = y - ys. eps is the source softening (`tpfcore_source_softening` or global `softening`).

With c = 0, Theta is the pure Hessian; c ≠ 0 adds a minimal isotropic tensor term to test whether near-source residuals decrease. This is an **exploratory correction to improve structural consistency**, NOT a fully derived final TPF source law. The ansatz remains provisional.

## What is directly paper-derived

- **Xi^mu** — configuration displacement field
- **Theta_{mu nu} = nabla_mu Xi_nu** — configuration gradient tensor
- **I = Theta_{mu nu} Theta^{mu nu} - lambda Theta^2** — scalar invariant
- **lambda = 1/4** in 4D (fixed)

## Where the ansatz lives

- **`source_ansatz.hpp`** — Declaration, `PointSourceField`, `provisional_point_source_field`
- **`source_ansatz.cpp`** — Closed-form Phi, Xi, Theta derivatives

## Field-equation residual

TPFCore inspection modes compute the **configuration-equation residual** from the paper’s weak-field form:

> R_ν = ∂_i (Θ^i_ν − λ δ^i_ν Θ)

with λ = 1/4 and Θ = Θ_xx + Θ_yy. This is the spatial part of nabla_μ(Θ^μ_ν − λ δ^μ_ν Θ).

- **Computation**: Analytic. Third derivatives of Phi are closed-form, so the residual is evaluated directly from the ansatz.
- **Meaning of near-zero**: If the ansatz satisfied the paper’s configuration equation, R would be zero. A small residual indicates the ansatz is structurally closer; a large residual indicates a mismatch.
- **Scope**: This is a **structural field-equation sanity check**, not proof of full TPF dynamics. It helps judge whether the Hessian-based point-source ansatz approximates the configuration equation in the weak-field limit.

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
| residual_x, residual_y, residual_norm | Configuration-equation residual components and norm |

**invariant_profile.csv**

| Column | Meaning |
|--------|---------|
| radius | (or axis, radius for symmetric pair) |
| invariant_I | Scalar invariant |
| theta_trace | Trace of Theta |
| residual_norm | ||R|| = sqrt(R_x² + R_y²) |

## Symmetry behavior on probe lines

For **single-source** at origin, probing along +x (y=0):

- **theta_xy ≈ 0** — by symmetry (y=0 is a principal axis)
- **residual_y ≈ 0** — by symmetry (configuration equation)
- **theta_xx ≠ theta_yy** — radial vs transverse structure (anisotropy)
- **invariant_I** — decays smoothly with radius

For **symmetric pair** at (±d, 0):

- **residual_y ≈ 0** on +x axis (y=0 symmetry)
- **residual_x ≈ 0** on +y axis (x=0 symmetry)

`field_summary.txt` reports max residual magnitudes and whether these symmetry expectations hold.

## Config (defaults.cfg)

- `tpfcore_enable_provisional_readout` — false (dynamics not available)
- `tpfcore_probe_radius_min`, `tpfcore_probe_radius_max`, `tpfcore_probe_samples`
- `tpfcore_dump_theta_profile`, `tpfcore_dump_invariant_profile`
- `tpfcore_source_softening` — softening for Phi. If ≤ 0, use global `softening`.
- `tpfcore_residual_step` — step size for numerical residual (default 1e-6); not used when analytic.
- `tpfcore_isotropic_correction_c` — dimensionless coefficient for B(r) = c·M/(r²+eps²)^(3/2). Default 0.0 (pure Hessian). Try e.g. 0.1, 0.25, 0.5, 1.0 to test residual reduction.

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
