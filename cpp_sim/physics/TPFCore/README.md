# TPFCore: Honest primitive TPF structure package

**This is NOT the removed weak-field Newtonian-like TPF package.** TPFCore implements the bare primitive structure from the paper. The source uses the **3D** Hessian of a softened Coulomb potential evaluated on the **z = 0** simulation plane (source and field point at z = 0).

## Source ansatz (provisional weak-field)

- **Phi** = -M / R with **R²** = dx² + dy² + eps² — same as 3D distance with z = z_s = 0 and isotropic softening eps² in the radical
- **Xi** = (∂_x Phi, ∂_y Phi) — in-plane displacement; ∂_z Phi = 0 on the plane
- **Theta** = full symmetric **3×3** Hessian of Phi; on z = 0: Theta_xz = Theta_yz = 0, Theta_zz = m/R³

dx = x − x_s, dy = y − y_s. eps is the source softening (`tpfcore_source_softening` or global `softening`). The field-equation residual uses the **full 3D** spatial divergence (including ∂_z Theta_{xz}, ∂_z Theta_{yz}).

## What is directly paper-derived

- **Xi^mu** — configuration displacement field
- **Theta_{mu nu} = nabla_mu Xi_nu** — configuration gradient tensor
- **I = Theta_{mu nu} Theta^{mu nu} - lambda Theta^2** — scalar invariant
- **lambda = 1/4** in 4D (fixed)

## Where the code lives

- **`source_ansatz.hpp`**, **`source_ansatz.cpp`** — Primitive field: Phi, Xi, Theta (closed-form derivatives)
- **`provisional_readout.hpp`**, **`provisional_readout.cpp`** — Provisional motion layer: maps Theta → acceleration (exploratory; isolated from ansatz)

## Field-equation residual

TPFCore inspection modes compute the **configuration-equation residual** from the paper’s weak-field form:

> R_ν = ∂_i (Θ^i_ν − λ δ^i_ν Θ)

with λ = 1/4 and **Θ = tr(Theta) = Θ_xx + Θ_yy + Θ_zz**, and **i** summed over **x, y, z** (including contributions from ∂_z even on z = 0). This matches the spatial part of nabla_μ(Θ^μ_ν − λ δ^μ_ν Θ) for ν ∈ {x, y} on the plane.

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
| theta_trace | Theta_xx + Theta_yy + Theta_zz (full 3D trace on slice) |
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

- `tpfcore_enable_provisional_readout` — false. Set true to enable tensor-driven motion (exploratory).
- `tpfcore_readout_mode` — `tensor_radial_projection`, `tensor_radial_projection_negated`, or `tr_coherence_readout` (see below).
- `tpfcore_readout_scale` — scale factor for readout magnitude (default 1.0)
- `tpfcore_theta_tt_scale` — (tr_coherence_readout only) Theta_tt balancing companion scale (default 1.0)
- `tpfcore_theta_tr_scale` — (tr_coherence_readout only) Theta_tr mixed coupling scale (default 1.0)
- `tpfcore_dump_readout_debug` — write `tpf_readout_debug.csv` for dynamical runs (default true)
- `tpfcore_probe_radius_min`, `tpfcore_probe_radius_max`, `tpfcore_probe_samples`
- `tpfcore_dump_theta_profile`, `tpfcore_dump_invariant_profile`
- `tpfcore_source_softening` — softening for Phi. If ≤ 0, use global `softening`.
- `tpfcore_residual_step` — step size for numerical residual (default 1e-6); not used when analytic.

## How to run

```bash
# In config: physics_package = TPFCore

# Single source at origin, probe along +x
./galaxy_sim tpf_single_source_inspect

# Symmetric pair at (±d,0), probe along +x and +y
./galaxy_sim tpf_symmetric_pair_inspect
```

## Provisional motion/readout layer

A **provisional**, **exploratory** motion/readout layer maps the local tensor field (Theta) into an acceleration vector for the simulator. This is **NOT** the full derived TPF dynamics. It exists solely to test particle motion driven by the tensor field without reverting to Newtonian scalar-potential dynamics. All readout modes are explicitly exploratory; no silent fallback to Newtonian.

**Design:**
- Motion derived from local TPF tensor structure (Theta), **not** from Phi or −grad(Φ)
- Explicitly isolated in `provisional_readout.cpp`; easy to swap out later

**Enabling:** Set `tpfcore_enable_provisional_readout = true` in config.

### Readout mode: `tensor_radial_projection` (spatial)

- Per source: compute Theta_s at particle location
- r_hat = unit vector from source to particle
- Contribution: `scale × (Theta_s · r_hat)` — tensor-vector contraction
- Superpose over BH (if any) and other particles (when star_star)

**Note:** This spatial projection is exploratory only. It **did not produce bound two-body motion** (particle moved outward continuously in two_body_orbit). It is conceptually misaligned with the paper’s orbit/coherence discussion, which uses time–radial structure (Theta_tt, Theta_rr, Theta_tr).

**tensor_radial_projection_negated**: Same projection, opposite sign. For debugging sign errors only. Not a final theory result.

### Readout mode: `tr_coherence_readout` (paper-aligned t–r structure)

**Exploratory** readout inspired by the paper’s weak-field t–r orbit/coherence discussion (Theta_tt, Theta_rr, Theta_tr). It is **not** the final derived TPF motion law; it moves the provisional motion layer closer to the paper than the spatial projection.

- **Theta_rr**: spatial tensor projected onto local radial: r_hat^T Theta r_hat (from superposed Theta at particle).
- **Theta_tt**: configurable balancing companion term: `tpfcore_theta_tt_scale × (−Theta_rr)` so that Theta_rr + Theta_tt can be tuned (e.g. toward balance).
- **Theta_tr**: provisional mixed time–radial coupling: t_hat^T Theta r_hat (tangential unit t_hat ⊥ r_hat).
- **Provisional radial readout**: `readout_scale × (Theta_rr + Theta_tt)` along r_hat.
- **Provisional tangential readout**: `readout_scale × tpfcore_theta_tr_scale × Theta_tr` along t_hat.
- **Acceleration**: radial component + tangential component (no Phi gradient).

**Config:** `tpfcore_readout_mode = tr_coherence_readout`, `tpfcore_theta_tt_scale` (default 1.0), `tpfcore_theta_tr_scale` (default 1.0).

**Supported simulation modes** (with any readout enabled):
- `two_body_orbit`, `symmetric_pair`, `small_n_conservation`, `galaxy`

### Readout debug CSV (`tpf_readout_debug.csv`, when `tpfcore_dump_readout_debug=true`)

**For `tensor_radial_projection` (and negated):**  
`time`, `particle`, `x`, `y`, `vx`, `vy`, `ax`, `ay`, `radius`, `radial_unit_x`, `radial_unit_y`, `a_radial`, `a_inward`, `a_tangential`, `theta_xx`, `theta_xy`, `theta_yy`, `theta_trace`, `invariant_I`.

**For `tr_coherence_readout`:**  
`time`, `particle`, `x`, `y`, `vx`, `vy`, `radius`, `theta_rr`, `theta_tt`, `theta_tr`, `theta_rr_plus_theta_tt`, `provisional_radial_readout`, `provisional_tangential_readout`, `ax`, `ay`, `a_radial`, `a_inward`, `a_tangential`.

- **a_radial** = a · r_hat; **a_inward** = −a_radial (positive = toward origin); **a_tangential** = component perpendicular to r_hat.

Used to diagnose runaway vs bound behavior. For bound orbits, a_inward should be positive (attraction).

**run_info.txt** includes `tpfcore_readout_mode`, `tpfcore_readout_scale`, `tpfcore_theta_tt_scale`, `tpfcore_theta_tr_scale` (when applicable), `tpfcore_dump_readout_debug`.

## What is NOT implemented

- Full nonlinear/dynamic TPF field equations (the readout is exploratory, not derived)
- Newtonian-style substitution (−grad(Phi) acceleration)

## Warning

TPFCore is inspection-first. The provisional readout is exploratory and not the final TPF dynamics. The ansatz remains provisional where the paper remains underspecified.
