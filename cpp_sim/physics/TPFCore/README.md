# TPFCore: Honest primitive TPF structure package

**This is NOT the removed weak-field Newtonian-like TPF package.** That package was misleading and has been deleted. TPFCore implements the bare primitive structure from the paper as far as honestly possible.

## What is directly paper-derived

- **Xi^mu** — configuration displacement field (primitive ontology)
- **Theta_{mu nu} = nabla_mu Xi_nu** — configuration gradient tensor
- **I = Theta_{mu nu} Theta^{mu nu} - lambda Theta^2** — scalar invariant
- **lambda = 1/4** in 4D (fixed, not tunable)

## What is provisional

- **Point-source ansatz** — The paper does not fully specify the point-source primitive field. The implementation in **`source_ansatz.hpp`** / **`source_ansatz.cpp`** is a clearly labeled PROVISIONAL ansatz. All assumptions are marked in code comments.
- **Acceleration readout** — Not implemented. TPFCore is inspection-first. A provisional readout layer could be added later (and must be explicitly enabled); it would NOT use scalar-phi or -grad(phi).

## Where the provisional ansatz lives

- **`source_ansatz.hpp`** — Declaration and documentation
- **`source_ansatz.cpp`** — Implementation of `provisional_point_source_theta` and `compute_invariant_I`

## What inspection outputs exist

For `tpf_single_source_inspect` and `tpf_symmetric_pair_inspect`:

| File | Contents |
|------|----------|
| `theta_profile.csv` | radius, theta_xx, theta_xy, theta_yy, theta_trace, invariant_I |
| `invariant_profile.csv` | radius, invariant_I |
| `field_summary.txt` | Geometry, source positions, provisional-ansatz flag |

## What is NOT yet implemented

- Full nonlinear/dynamic TPF field equations
- Variational evolution
- Galaxy-scale dynamics
- Acceleration readout for integration (provisional or otherwise)
- Any Newtonian-style scalar potential or -grad(phi) substitution

## How to run

```bash
# Single source at origin, probe along +x
physics_package = TPFCore
./galaxy_sim tpf_single_source_inspect

# Symmetric pair at (±d,0), probe along +x
physics_package = TPFCore
./galaxy_sim tpf_symmetric_pair_inspect
```

## Config (defaults.cfg)

- `tpfcore_enable_provisional_readout` — false (dynamics not available)
- `tpfcore_probe_radius_min`, `tpfcore_probe_radius_max`, `tpfcore_probe_samples`
- `tpfcore_dump_theta_profile`, `tpfcore_dump_invariant_profile`

## Warning

TPFCore is incomplete and provisional where the paper remains underspecified. It is an honest skeleton for inspection, not a completed galaxy-dynamics solver.
