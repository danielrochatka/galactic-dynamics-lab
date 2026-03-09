# TPF weak-field correspondence package

**This package implements the paper's weak-field correspondence sector only.** It is not a numerical implementation of the full variational TPF field equations.

## What this package implements

- **Linearized weak-field / quasi-static sector** of the TPF framework:
  - Curl-free displacement field: **xi = grad(phi)**
  - Configuration gradient: **Theta_ij = partial_i partial_j phi**
  - Scalar invariant and trace: **Theta = nabla² phi**
  - Linearized trace equation: **nabla² phi = alpha * rho(x)**

- **Source model**: Particle masses as source masses. The scalar field phi(x) is built from superposition of Green-function contributions (no grid), equivalent to solving the Poisson equation for pointlike or softened sources.

- **Operational acceleration readout**: The simulator requires per-particle accelerations. This package derives them from the weak-field layer: **a = -grad(phi)**. This is an *operational readout* used by the engine, not the primitive TPF ontology. The primitive objects are the displacement field and configuration tensor; accelerations are a derived layer for integration.

- **Lambda**: Fixed at 1/4 in 4D (not exposed as a parameter).

## What this package does NOT implement

- Full variational TPF field equations
- Nonlinear / dynamic TPF evolution
- Arbitrary vector-force summation as the primitive law
- Earth–Moon heuristic finite-difference collapse analogue
- Disk-galaxy predictions from the full TPF equations
- Any physics beyond the weak-field Poisson-type reduction

## Mapping to the simulator

1. **Phi from sources**: phi(x) = -(alpha / 4π) × Σ m_j / r_ij (softened)
2. **Acceleration readout**: a_i = -grad(phi)|_i computed via superposition
3. Same BH + optional star–star structure as Newtonian; softening for regularity

## Configuration

Package-specific options live in **`physics/TPF/defaults.cfg`**. Users can override them from the run config (`configs/my.local.cfg`).

| Key | Meaning | Default |
|-----|---------|---------|
| `physics_package` | `TPF` to select this package | `Newtonian` |
| `tpf_alpha` | Weak-field coupling in nabla² phi = alpha×rho | (see below) |
| `tpf_match_newtonian_scale` | If true, use alpha=4π to match Newtonian magnitudes | `true` |
| `tpf_softening` | Softening for TPF Green-function. If ≤0, use global `softening` | `0` |

When `tpf_match_newtonian_scale = true`, alpha is set to 4π so the weak-field readout reproduces Newtonian-scale accelerations for comparison runs. When false, `tpf_alpha` is used directly (must be > 0).

See the main README for the layered config system (built-in → package defaults → run config).

## Output

- `run_info.txt` includes: `physics_package`, `tpf_alpha`, `tpf_match_newtonian_scale`, `tpf_softening`, and the loaded config paths
- Console indicates when the TPF weak-field package is in use

## Warning

This is **not** yet the full nonlinear, dynamic, or galactic TPF solver. It is a conservative, paper-faithful implementation of the weak-field correspondence sector for testing and comparison.
