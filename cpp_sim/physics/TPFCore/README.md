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

## Where the code lives

- **`source_ansatz.hpp`**, **`source_ansatz.cpp`** — Primitive field: Phi, Xi, Theta (closed-form derivatives)
- **`provisional_readout.hpp`**, **`provisional_readout.cpp`** — Provisional motion layer: maps Theta → acceleration (exploratory; isolated from ansatz)

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

- `tpfcore_enable_provisional_readout` — false. Set true to enable tensor-driven motion (exploratory).
- `tpfcore_readout_mode` — `tensor_radial_projection` (Theta·r_hat per source)
- `tpfcore_readout_scale` — scale factor for readout magnitude (default 1.0)
- `tpfcore_probe_radius_min`, `tpfcore_probe_radius_max`, `tpfcore_probe_samples`
- `tpfcore_dump_theta_profile`, `tpfcore_dump_invariant_profile`
- `tpfcore_source_softening` — softening for Phi. If ≤ 0, use global `softening`.
- `tpfcore_residual_step` — step size for numerical residual (default 1e-6); not used when analytic.
- `tpfcore_isotropic_correction_c` — dimensionless coefficient for B(r) = c·M/(r²+eps²)^(3/2). Default 0.0 (pure Hessian). Try e.g. 0.1, 0.25, 0.5, 1.0 to test residual reduction.
- `tpfcore_c_sweep_min`, `tpfcore_c_sweep_max`, `tpfcore_c_sweep_steps` — c-sweep range and resolution (for `tpf_single_source_optimize_c`).
- `tpfcore_c_objective` — objective to minimize: `max_residual_norm`, `mean_residual_norm`, or `l2_residual_norm`.

## c-sweep utility (exploratory ansatz-tuning)

**`tpf_single_source_optimize_c`** is an exploratory ansatz-tuning tool. It numerically fits c against the field-equation residual of the current provisional ansatz. The optimal c may be irrational or otherwise non-obvious, so manual guessing is not enough—this utility sweeps c over a configurable range and selects the value that minimizes the chosen objective.

**Important:** The fitted c is **NOT** a final paper-derived constant. It is a numerically tuned value for this ansatz and geometry. Do not present it as a derived TPF parameter.

### How to run the c sweep

```bash
# In config: physics_package = TPFCore

./galaxy_sim tpf_single_source_optimize_c
```

Config keys: `tpfcore_c_sweep_min`, `tpfcore_c_sweep_max`, `tpfcore_c_sweep_steps`, `tpfcore_c_objective`.

### Output files

- **`c_sweep.csv`** — per-c row: `c`, `max_residual_norm`, `mean_residual_norm`, `l2_residual_norm`
- **`c_sweep_summary.txt`** — sweep range, number of steps, chosen objective, best c, best objective value

### Objective metrics (minimized)

| Objective | Meaning |
|-----------|---------|
| `max_residual_norm` | max over probe points of ‖R‖. Penalizes worst-case residual (e.g. near source). |
| `mean_residual_norm` | mean of ‖R‖ over probe points. Balances overall fit. |
| `l2_residual_norm` | sqrt(sum ‖R‖²). Stronger penalty on large local residuals. |

Lower values indicate the ansatz better satisfies the configuration equation for that c.

## How to run

```bash
# In config: physics_package = TPFCore

# Single source at origin, probe along +x
./galaxy_sim tpf_single_source_inspect

# Symmetric pair at (±d,0), probe along +x and +y
./galaxy_sim tpf_symmetric_pair_inspect

# Numerically fit c against residual (exploratory; fitted c is not a final constant)
./galaxy_sim tpf_single_source_optimize_c
```

## Provisional motion/readout layer

A **provisional**, **exploratory** motion/readout layer maps the local tensor field (Theta) into an acceleration vector for the simulator. This is **NOT** the full derived TPF dynamics. It exists solely to test particle motion driven by the tensor field without reverting to Newtonian scalar-potential dynamics.

**Design:**
- Motion derived from local TPF tensor structure (Theta), **not** from Phi or −grad(Phi)
- Tensor-driven: Theta_ij contracted with radial unit vector r_hat
- Explicitly isolated in `provisional_readout.cpp`; easy to swap out later

**Enabling:** Set `tpfcore_enable_provisional_readout = true` in config.

**Readout mode** `tensor_radial_projection`:
- Per source: compute Theta_s at particle location
- r_hat = unit vector from source to particle
- Contribution: `scale × (Theta_s · r_hat)` — tensor-vector contraction
- Superpose over BH (if any) and other particles (when star_star)

**Config:**
- `tpfcore_readout_mode` — `tensor_radial_projection` (default)
- `tpfcore_readout_scale` — magnitude scale (default 1.0)

**Supported simulation modes** (with readout enabled):
- `two_body_orbit`, `symmetric_pair`, `small_n_conservation`, `galaxy`

**run_info.txt** includes `tpfcore_readout_mode` and `tpfcore_readout_scale`.

## What is NOT implemented

- Full nonlinear/dynamic TPF field equations (the readout is exploratory, not derived)
- Newtonian-style substitution (−grad(Phi) acceleration)

## Warning

TPFCore is inspection-first. The provisional readout is exploratory and not the final TPF dynamics. The ansatz remains provisional where the paper remains underspecified.
