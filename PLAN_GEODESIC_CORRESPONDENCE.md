# Plan: geodesic_correspondence Acceleration Route
## Date: 2026-05-02
## Status: Investigation — read-only, no code modified

## Executive Summary

The proposed `geodesic_correspondence` route implements the chain matter → Ξ → ρ_Ξ → Poisson →
Φ → a = −∇Φ per TPF_FOUNDATIONS.md §4. The density step uses the linear flux correspondence
`ρ_Ξ = (1/4π) ∂_i Ξ^i` (Paper A v14.4 Appendix D), which recovers the original point mass via
Gauss's theorem with no softening calibration required. All necessary infrastructure either already
exists in the codebase or can be borrowed unchanged from the Newtonian package. Behavior in
non-trivial test cases is to be determined empirically once implemented.

---

## A. Code Paths Already in Place

### A.1 Θ Computation

**Location:** `engine/physics/TPFCore/source_ansatz.cpp`, function
`provisional_point_source_field()` (lines 12–37), and accumulated over all sources in
`field_evaluation.cpp::evaluate_provisional_field_multi_source()` (lines 40–69).

**Verified by reading function body.** The Theta components are:

```
Theta.xx = m * (r² − 3 dx²) / R⁵
Theta.xy = −3 m dx dy / R⁵
Theta.zz = m r² / R⁵      (= m / R³ at z = 0)
```

This is the 3D Hessian of Φ = −m/R (with R² = dx² + dy² + eps²), i.e., Θ_ij = ∂_i ∂_j Φ.
**This matches Paper A Eqs. (20)–(25).** The note in the header says it explicitly:
"Theta_ij = Hess_ij(Phi)" and "3D Hessian of Phi = -m/R at z = 0".

Multi-source superposition: `evaluate_provisional_field_multi_source` linearly sums per-source
Xi and Theta components, then recomputes the trace and invariant I from the summed Theta. This is
the correct procedure (linearity of the Hessian).

### A.2 C_μν Principal-Part Computation

**Location:** `engine/physics/TPFCore/direct_tpf_baseline.cpp`,
`compute_principal_Cij_from_eq10_baseline()` (lines 29–50).

**Verified by reading function body.** Computes:

```
C_xx = κ (Θ_xα Θ^α_x − λ tr(Θ) Θ_xx − ½ I)
C_xy = κ (Θ_xα Θ^α_y − λ tr(Θ) Θ_xy)
C_yy = κ (Θ_yα Θ^α_y − λ tr(Θ) Θ_yy − ½ I)
C_zz = κ (Θ_zα Θ^α_z − λ tr(Θ) Θ_zz − ½ I)
```

with λ = LAMBDA_4D = 0.25 (fixed, not tunable). This is correct as a field-equation tensor (the
analog of G_μν). The misuse is in `compute_xi_directed_tensor_readout()` (lines 52–66), which then
projects C_ij along the Xi direction to produce acceleration — that step is the wrong move
(C_ij is not an acceleration; it is a field equation tensor). The C_μν computation itself is
correct; only its interpretation as direct acceleration is wrong.

### A.3 Invariant I

**Location:** `source_ansatz.cpp::compute_invariant_I()` (lines 74–80); also
`derived_tpf_radial.cpp::derived_invariant_I_contracted()` (lines 48–52).

`compute_invariant_I` computes `I = Θ_μν Θ^μν − λ Θ²` (Frobenius norm squared minus λ times
trace-squared) — the field-equation invariant used in C_μν. `derived_invariant_I_contracted`
computes only `Θ_ij Θ^ij` (raw Frobenius squared, no λ correction) — its header comment says
"Frobenius Θ_ij Θ_ij only — not compute_invariant_I". **Neither function is the density proxy
for `geodesic_correspondence`.** The density step uses the linear flux divergence
`ρ_Ξ = (1/4π) ∂_i Ξ^i` (Paper A v14.4 Appendix D), not a Θ contraction. `Θ_ij Θ^ij` belongs
to the nonlinear energy-ledger sector (v14.4 §IV.A.d) and is documented here only for completeness.

### A.4 Existing Poisson Infrastructure

**Location:** `engine/physics/TPFCore/derived_tpf_radial.cpp`,
`build_tpf_gravity_profile()` (lines 102–182).

**Verified by reading function body.** This function:

1. Loops over radial bins.
2. At each bin center evaluates `sum_derived_theta_at_point()` (the summed Theta from all
   sources evaluated at a 1D probe point on the +x axis).
3. Computes `I_total = |compute_invariant_I(theta_tot)|`.
4. Forms `rho_raw = I_total * |kappa|`.
5. Applies an R6 bounce suppression near the BH Schwarzschild radius.
6. Integrates `shell_mass = 4π r² ρ_eff dr` cumulatively.

The result `M_eff_enc` is used by `radial_acceleration_scalar_derived()` which produces
`a_r = −G (M_baryon_bounced + M_eff) / r_soft²`. This **is** a Poisson solve — but over a
radially-binned 1D profile, not the full 2D plane, and uses `compute_invariant_I` (which includes
the λ correction) rather than the raw `Θ_ij Θ^ij` contraction. It is a pre-existing approximation
serving the `derived_tpf_radial_readout` / `tr_coherence_readout` legacy modes.

`derived_poisson_cfg_` (declared in `tpf_core_package.hpp` line 174) configures this: kappa,
bins, max_radius, galaxy_radius. It is populated from config keys `tpf_kappa`,
`tpf_poisson_bins`, `tpf_poisson_max_radius`, `galaxy_radius` (lines 124–127 of
`tpf_core_package.cpp`).

There is no full 2D Poisson solver in TPFCore. The existing infrastructure is strictly 1D radial.

### A.5 Acceleration Dispatch

**Location:** `tpf_core_package.cpp::TPFCorePackage::compute_accelerations()` (lines 802–987).

**Verified by reading function body.** The dispatch is an `if` / `if` / `if` chain (not if/else):

1. `tpf_dynamics_mode_ == "xi_kernel_deformed"` → `compute_xi_kernel_deformed_accelerations()`
   then `return`.
2. `tpf_dynamics_mode_ == "direct_tpf"` → `compute_direct_tpf_accelerations()` + optional VDSG
   additive then `return`.
3. `tpf_dynamics_mode_ == "v11_weak_field_truncation"` →
   `compute_v11_weak_field_correspondence_accelerations()` then `return`.
4. Fall-through: requires `provisional_readout_` to be true, else throws. Enters `eval_accel_pipeline()`.

The config key is `tpf_dynamics_mode` (`engine/config.hpp` line 161, default
`"xi_kernel_deformed"`).

The `v11_weak_field_truncation` path calls
`compute_v11_weak_field_correspondence_accelerations()` (verified in
`v11_weak_field_correspondence.cpp` lines 415–456). That function computes
`a = alpha_si * M_j * (x_i − x_j) / |r|³` — i.e., it is a Newtonian-gradient formula with
calibrated coefficient alpha_si. At default `alpha_si = −G_SI` (config line 168) this is exactly
Newtonian gravity with sign flipped. When properly configured, it is
`a = −G Σ M_j (x − x_j) / |x − x_j|³`, which is exactly the gradient of Φ_Newton.

---

## B. Code Paths Needed for geodesic_correspondence

### B.1 Dispatch Registration

**Where:** Rename the existing `"v11_weak_field_truncation"` branch in `compute_accelerations()`
(`tpf_core_package.cpp` line ~845) to match `"geodesic_correspondence"`. Add a second branch
catching the old string `"v11_weak_field_truncation"` as a deprecation alias — emit a warning
and forward to the renamed path when `alpha_si == −G_SI`, or run free-parameter as-is with a
notice otherwise. No new private method is required; `compute_v11_weak_field_truncation_accelerations()`
is renamed to `compute_geodesic_correspondence_accelerations()` in place.

Also update `"geodesic_correspondence"` in any validation or manifest enumerations. The
`init_from_config()` block at line 113 already handles a similar alias rename for
`"weak_field_correspondence"` and serves as the pattern to follow.

### B.2 Step 3 — ρ from linear flux divergence of Ξ (Paper A v14.4 Appendix D)

**Supersedes the earlier Θ_ij Θ^ij approach.** Paper A v14.4 Appendix D establishes the
correct correspondence as the linear flux divergence:

```
ρ_Ξ = (1/4π) ∂_i Ξ^i
```

For the point-source ansatz `Ξ^i = M x^i / r³` (already computed in
`provisional_point_source_field()`, verified: `xi.x = m * dx / R³`, `xi.y = m * dy / R³`):

- **Exterior (r > 0):** `∂_i Ξ^i = ∂_x(M x/r³) + ∂_y(M y/r³) + ∂_z(M z/r³) = 0`
  (standard result: divergence of x/r³ vanishes away from the origin)
- **At source (distributional):** `∂_i Ξ^i = 4π M δ³(x)`
  (Gauss's theorem: ∮ Ξ · dA = 4πM, same as for the Newtonian field)
- Therefore: `ρ_Ξ = (1/4π) × 4πM δ³(x) = M δ³(x)` — **recovers the point mass exactly**

This is a clean derivation with no softening calibration required. The quadratic invariant
`Θ_ij Θ^ij` is correctly relegated to the nonlinear energy-ledger sector per v14.4 Section
IV.A.d; it is not part of the geodesic_correspondence density step.

**Code status:** No new code is needed for the multi-source case. The divergence is linear, so
`∂_i Ξ^i_total = Σ_j ∂_i Ξ^i_j`. For the analytic point-source Poisson shortcut (Step 4), the
delta-function form means `∇²Φ = 4πG Σ M_j δ³(x − x_j)` → `Φ = −G Σ M_j / |x − x_j|`
directly. `derived_invariant_I_contracted()` is **not used** by this step.

### B.3 Step 4 — Poisson Solve

For point sources with masses M_i at positions x_i, the analytic solution to ∇²Φ = 4πG ρ with
ρ = Σ M_i δ³(x − x_i) is:

```
Φ(x) = −G Σ M_i / |x − x_i|
```

This is already computed implicitly by the Newtonian package
(`newtonian.cpp::NewtonianPackage::compute_accelerations()`, lines 12–21) through the sum over
sources. The Newtonian package does NOT store Φ explicitly; it directly forms −∂Φ/∂x_i.

**No numerical Poisson solver is needed.** The analytic solution for point sources means the
entire Step 4 is bypassed — Φ is never formed as an intermediate; its gradient is formed directly
in Step 5.

### B.4 Step 5 — Gradient of Φ

For point sources:

```
∂_i Φ = G Σ M_j (x_i − x_j) / |x − x_j|³    (with softening)
a^i = −∂^i Φ = −G Σ M_j (x_i − x_j) / |x − x_j|³
```

This is exactly what `NewtonianPackage::compute_accelerations()` computes (verified: lines 14–21
of `newtonian.cpp`). It uses `G_SI = 6.6743e-11`, BH at origin, optional star-star, and
`apply_plummer_softening()` for the softened form. The `v11_weak_field_truncation` path also
computes this, with a configurable `alpha_si` scalar (default −G).

**The `geodesic_correspondence` route is, at this level, numerically identical to the Newtonian
path.** The distinction is entirely in the theoretical chain documented.

### B.5 Minimal New Code Needed

The implementation is a **rename with parameter lock**, not a new parallel branch:

1. **Rename dispatch string** in `compute_accelerations()`: change the `if` condition from
   `"v11_weak_field_truncation"` to `"geodesic_correspondence"`; add a second `if` branch
   catching `"v11_weak_field_truncation"` that emits a deprecation warning and falls through to
   the renamed path (or runs in free-parameter mode if `alpha_si ≠ −G_SI`). ~10–15 lines.
2. **Lock `alpha_si`**: Remove it from the config surface for the `geodesic_correspondence` mode;
   hardcode `−TPF_G_SI` at the call site. ~5 lines.
3. **Config validation**: Update any allowlist of valid `tpf_dynamics_mode` strings to include
   `"geodesic_correspondence"` and retain `"v11_weak_field_truncation"` as a known alias.
4. **Documentation block**: Add a function-level comment to
   `compute_v11_weak_field_correspondence_accelerations()` (or its renamed wrapper) explaining
   the five-step chain and the "consistent with regularization" caveat for Step 3.

Optional but recommended:
5. **Diagnostic log**: At init, log which mode is active and confirm alpha_si = −G_SI for the
   geodesic path. ~5 lines.

---

## C. Effort Estimate

| Item | Lines | Files |
|------|-------|-------|
| Rename dispatch string + deprecation alias guard | 10–15 | `tpf_core_package.cpp` |
| Lock `alpha_si` to `−G_SI`; remove from config surface for this mode | 5–8 | `tpf_core_package.cpp`, `config.hpp` |
| Config validation allowlist update | 3–5 | `config.hpp`, `defaults.cfg` |
| Documentation block (five-step chain + regularization caveat) | 15–20 | `tpf_core_package.cpp` |
| Update unit test: routing smoke (new mode string) | 10–15 | `tests/unit/test_routing_semantics.cpp` |
| Update integration tests: replace v11 mode string with geodesic_correspondence | 5–10 | existing smoke tests |
| **Total** | **~50–75 lines** | **3–4 files touched** |

**Tests at risk:** `test_routing_semantics.cpp` tests the enumeration of valid mode strings; adding
a new mode will not break existing tests but the "rejects unknown mode" test may need updating if
it uses a positive allowlist. `tpfcore_vdsg_manifest_smoke.sh` and `tpfcore_legacy_vdsg_canonical_smoke.sh`
check manifest / mode names; they should be checked for hardcoded mode lists.

---

## D. Honest Assessment of Issues

### D.2 Divergence concern resolved by linear flux correspondence (Paper A v14.4 Appendix D)

The earlier concern — that `Θ_ij Θ^ij ~ (GM)²/r⁶` diverges at r→0, requiring softening
calibration — **no longer applies to the geodesic_correspondence route**. Paper A v14.4
Appendix D resolves this by replacing the quadratic density proxy with the linear flux divergence:

```
ρ_Ξ = (1/4π) ∂_i Ξ^i
```

For the point-source ansatz this gives `ρ_Ξ = M δ³(x)` exactly — no softening parameter
enters the density recovery, and no calibration is required (see B.2 for derivation).

**Where the quadratic invariant still belongs:** `Θ_ij Θ^ij` is correctly placed in the
nonlinear energy-ledger sector per v14.4 Section IV.A.d (and Paper B App. E). It quantifies
the configuration-gradient energy density, not the matter source density. The two quantities
serve different roles in the theory and should not be conflated.

**What this means for the implementation:** `tpf_kappa` and `derived_invariant_I_contracted()`
are not used by `geodesic_correspondence`. The density step is handled analytically via the
flux divergence identity; `tpf_kappa` remains an uncalibrated placeholder relevant only to the
legacy `direct_tpf` and radial-Poisson paths.

### D.3 Static Dust Sufficiency

For galaxy simulation at typical velocities v/c ≲ 10⁻³ (galactic orbital speeds ~200 km/s vs
c = 3×10⁸ m/s), the static dust approximation is excellent. Pressure corrections enter at
O(v²/c²) ~ 10⁻⁶ and are negligible compared to Newtonian-level accuracy. The static-dust
correspondence is sufficient for all current use cases in this simulator.

Momentum terms in the stress-energy do matter for the VDSG extension (which models a
velocity-dependent correction), but that is an additive modifier on top of the base acceleration
and would remain independent of the geodesic_correspondence route.

### D.4 Geodesic in Discretized Form

In the static weak-field linearized metric (`ds² = −(1 + 2Φ/c²)c² dt² + (1 − 2Φ/c²) δ_ij dx^i dx^j`),
the geodesic equation for a slow particle reduces to:

```
d²x^i/dt² = −∂^i Φ   (leading order, v/c → 0)
```

In the simulator's Verlet/leapfrog integrator, the acceleration vector is consumed directly as
`d²x/dt²`. Setting `ax = −∂_x Φ` is therefore the exact intended discretization. No additional
integration step is required; the geodesic equation IS the Newtonian force law in this limit.
The time-stepping scheme is unchanged.

---

## E. Recommendation

### E.1 Smallest Viable First Implementation

**Rename `v11_weak_field_truncation` to `geodesic_correspondence`; do not add a parallel route.**

The plan previously proposed a separate `geodesic_correspondence` dispatch branch delegating to the
v11 kernel. That is unnecessary — the routes would share an identical call graph, identical kernel,
and identical guards. The actual minimal change is:

1. **Rename** the dispatch string `"v11_weak_field_truncation"` → `"geodesic_correspondence"` in
   `tpf_core_package.cpp` and all config/manifest references.
2. **Lock `alpha_si` to `−G_SI`** in the renamed route — it is now a derived constant (the
   geodesic chain derives `a = −G Σ M/r²` from first principles), not a tunable calibration
   parameter. Remove `alpha_si` from the config surface for this mode.
3. **Keep `"v11_weak_field_truncation"` as a deprecation alias** in the dispatch: if the mode
   string is `"v11_weak_field_truncation"` and `alpha_si == −G_SI`, emit a deprecation warning and
   forward to `geodesic_correspondence`. If `alpha_si ≠ −G_SI`, run as-is (Paper B calibration
   use) with a notice that it is in free-parameter calibration mode, not the foundations-grounded
   path.
4. **Add a documentation block** at the top of the renamed method explaining the five-step chain
   (Ξ → Θ → ρ_Ξ → Poisson → geodesic) with the v14.4 Appendix D linear flux correspondence:
   `ρ_Ξ = (1/4π) ∂_i Ξ^i`. For point sources this recovers M δ³(x) exactly; no softening
   calibration is required. Note that Θ_ij Θ^ij is the nonlinear energy-ledger quantity
   (v14.4 §IV.A.d), not the density proxy used here.

This is approximately 15–20 lines of changes (dispatch rename, parameter lock, alias guard,
documentation). The key deliverable is the correct theoretical framing, not new numerical code.

Route names should describe physical content, not paper-version history. Severing the v11 label
from this route is the primary architectural gain.

### E.2 What to Defer

- **Full numerical 2D Poisson solve** for smooth density distributions: not needed for point
  sources; defer until smooth-matter (N-body mesh) support is added.
- **Higher-order metric terms** (post-Newtonian corrections at O(v²/c²)): defer; below noise
  for current galaxy simulations.
- **Non-static corrections** (retardation, radiation reaction): defer; these require the time-
  dependent sector of TPF, which is not yet in scope.
- **Θ_ij Θ^ij energy-ledger implementation**: the nonlinear ledger quantity (v14.4 §IV.A.d) is
  not part of the geodesic_correspondence density step; defer its implementation to when the
  energy-ledger sector is developed separately.
- **VDSG additive extension on this route**: defer until the base route is validated; can be
  added trivially later (same pattern as `direct_tpf` + VDSG).

### E.3 Handling Existing Routes

| Route | Recommendation | Reason |
|-------|---------------|--------|
| `xi_kernel_deformed` | Keep (default) | Active development path |
| `direct_tpf` | Keep with deprecation warning | Wrong readout step; useful for C_μν diagnostics |
| `v11_weak_field_truncation` | Rename → `geodesic_correspondence`; keep as deprecation alias | Route name should reflect physical content, not paper version. Free-parameter calibration variant retained via alias. |
| `legacy_readout` | Keep (existing deprecation) | Historical comparison |
| `geodesic_correspondence` | New canonical name for the renamed v11 route | Correct five-step theoretical chain; alpha_si locked to −G_SI |

No existing routes are removed. `v11_weak_field_truncation` as a config string remains functional
as a migration alias, emitting a deprecation notice.

### E.4 Implementation Chain Summary

The five steps of the route and their code grounding:

1. **Ξ** — computed by `provisional_point_source_field()` (`source_ansatz.cpp:12–37`), Paper A Eq. (19).
2. **Θ** — computed by same function, Paper A Eqs. (20)–(25).
3. **ρ_Ξ** — linear flux correspondence `ρ_Ξ = (1/4π) ∂_i Ξ^i` (Paper A v14.4 Appendix D);
   for point sources recovers M δ³(x) via Gauss's theorem; no new code required.
4. **Poisson / Φ** — analytic solution `Φ = −G Σ M_i / |x − x_i|` for point sources;
   gradient computed by `compute_v11_weak_field_correspondence_accelerations()` kernel.
5. **Geodesic / a** — static weak-field linearized metric gives `a^i = −∂^i Φ`;
   consumed directly by the Verlet integrator unchanged.

Behavior of the route in non-trivial test cases (smooth matter, multi-source, complex boundaries)
is to be determined empirically once implemented. The plan does not predict numerical results in advance.

---

## Appendix: Key File Locations

| File | Function | Lines | Purpose |
|------|----------|-------|---------|
| `engine/physics/TPFCore/source_ansatz.cpp` | `provisional_point_source_field()` | 12–37 | Ξ and Θ from single point source |
| `engine/physics/TPFCore/source_ansatz.cpp` | `compute_invariant_I()` | 74–80 | I = Θ_μν Θ^μν − λ Θ² |
| `engine/physics/TPFCore/source_ansatz.hpp` | `LAMBDA_4D` | 58 | λ = 0.25 (fixed) |
| `engine/physics/TPFCore/field_evaluation.cpp` | `evaluate_provisional_field_multi_source()` | 40–69 | Multi-source Ξ, Θ accumulation |
| `engine/physics/TPFCore/direct_tpf_baseline.cpp` | `compute_principal_Cij_from_eq10_baseline()` | 29–50 | C_μν field equation tensor |
| `engine/physics/TPFCore/direct_tpf_baseline.cpp` | `compute_xi_directed_tensor_readout()` | 52–66 | Wrong readout (C_μν → acceleration) |
| `engine/physics/TPFCore/derived_tpf_radial.cpp` | `derived_invariant_I_contracted()` | 48–52 | Pure Θ_ij Θ^ij (no λ term) |
| `engine/physics/TPFCore/derived_tpf_radial.cpp` | `build_tpf_gravity_profile()` | 102–182 | 1D radial Poisson ledger |
| `engine/physics/TPFCore/derived_tpf_radial.hpp` | `DerivedTpfPoissonConfig` | 21–32 | Poisson config struct |
| `engine/physics/TPFCore/tpf_core_package.cpp` | `compute_accelerations()` (dispatch) | 802–874 | Mode dispatch chain |
| `engine/physics/TPFCore/tpf_core_package.cpp` | `compute_xi_kernel_deformed_accelerations()` | 650–800 | Active default path (wrong chain) |
| `engine/physics/TPFCore/tpf_core_package.cpp` | `compute_direct_tpf_accelerations()` | 562–585 | direct_tpf path |
| `engine/physics/TPFCore/tpf_core_package.cpp` | `compute_v11_weak_field_truncation_accelerations()` | 587–595 | v11 path |
| `engine/physics/TPFCore/tpf_core_package.hpp` | `derived_poisson_cfg_` | 174 | Poisson config member |
| `engine/physics/TPFCore/v11_weak_field_correspondence.cpp` | `compute_v11_weak_field_correspondence_accelerations()` | 415–456 | Newtonian-equivalent acceleration (α·M/r²) |
| `engine/physics/Newtonian/newtonian.cpp` | `compute_accelerations()` | 12–21 | Ground-truth Newtonian: a = −G Σ M/r² |
| `engine/config.hpp` | `tpf_dynamics_mode` | 161 | Mode dispatch config key (default: xi_kernel_deformed) |
| `engine/config.hpp` | `tpf_kappa` | 215 | κ coupling (default 1e32, uncalibrated) |
