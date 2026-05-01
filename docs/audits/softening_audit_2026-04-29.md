# Softening audit — 2026-04-29

## Scope
- Newtonian implementation: `engine/physics/Newtonian/newtonian.cpp`.
- TPF Xi-off runtime implementation: `engine/physics/TPFCore/tpf_core_package.cpp`.
- Unit tests expanded in:
  - `engine/physics/Newtonian/tests/unit/test_newtonian_forces.cpp`
  - `engine/physics/TPFCore/tests/unit/test_xi_kernel_runtime_audit.cpp`

## Formulas checked
1. Newtonian BH acceleration:
   - Verified implementation is `a = -G M r_vec / (r^2 + eps^2)^(3/2)`.
2. Newtonian star-star acceleration:
   - Verified implementation is `a_i += G m_j (x_j - x_i)/(|x_j-x_i|^2 + eps^2)^(3/2)`.
3. Newtonian potential:
   - Verified BH and pair terms use `U = -G M m / sqrt(r^2 + eps^2)` and softened pair distance consistently.
4. Xi-off readout correspondence:
   - For `tpf_dynamics_mode=xi_kernel_deformed`, `tpf_4d_xi_kernel_mode=off`, acceleration is `a = -K_xi Xi` with softened `1/r^3` kernel.

## Definite bug found and fixed
- **Bug:** Xi-kernel runtime route used global call softening unconditionally and ignored `tpfcore_source_softening` override.
- **Fix:** Xi-kernel route now uses effective source softening `eps = (tpfcore_source_softening > 0 ? tpfcore_source_softening : softening)`.
- **Why this matters:** The config truth contract for `tpfcore_source_softening` was not honored on Xi-kernel runtime path.

## Tests added

### Newtonian softening correctness tests
- Analytic BH softened acceleration equality test.
- Finite-difference potential-gradient vs acceleration consistency (BH-only).
- Finite-difference potential-gradient vs acceleration consistency (star-star pair).
- Pairwise momentum symmetry (`m_i a_i + m_j a_j = 0`) for unequal masses and `eps>0`.
- Softening-limit behaviors:
  - `eps=0` recovers inverse-square magnitude.
  - increasing `eps` monotonically decreases |a| at fixed r.
  - very large `eps` yields finite small acceleration.
  - `r=0, eps>0` is finite and non-NaN/Inf.

### TPF Xi-off correspondence + softening-source truth tests
- Xi-off with `K_xi=G_SI` matches Newtonian acceleration for:
  - BH-only (`star_star=false`),
  - star-star enabled (`star_star=true`),
  - multiple `eps` values.
- `tpfcore_source_softening` truth test:
  - positive override value is used instead of global softening;
  - zero value falls back to global softening.

## Audit conclusion
- Newtonian softening implementation is internally consistent between acceleration and potential.
- Newtonian pairwise softened forces conserve momentum within numerical tolerance in tested configurations.
- Xi-off path matches Newtonian when sharing softening and setting `K_xi = G_SI`.
- Similar morphology/clumping between Newtonian and Xi-off is therefore less likely to be from a shared softening formula bug in these audited paths; remaining clumping concern is more likely physical/IC/integration/regime dependent than a direct softening inconsistency.
