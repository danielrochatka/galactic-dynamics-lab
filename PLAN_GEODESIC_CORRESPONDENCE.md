# Plan: geodesic_correspondence Acceleration Route
## Date: 2026-05-02
## Status: Investigation — read-only, no code modified

## Executive Summary

The proposed `geodesic_correspondence` route implements the chain matter → Ξ → ρ_Ξ → Poisson →
Φ → a = −∇Φ per TPF_FOUNDATIONS.md §4. The density step uses the linear flux correspondence
`ρ_Ξ = (1/4π) ∂_i Ξ^i` (Paper A v14.4 Appendix D). In the singular (ε → 0) limit this recovers
M δ³(x) exactly; the softened implementation produces a regularized Plummer density that integrates
to M (see B.2). The runtime route uses the analytic closed-form shortcut of the correspondence
chain for point sources and does not numerically compute ρ_Ξ, solve Poisson, or differentiate Φ
on a grid. All necessary infrastructure either already exists in the codebase or can be borrowed
unchanged from the Newtonian package. Behavior in non-trivial test cases is to be determined
empirically once implemented.

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

**κ terminology note:** Two distinct things share the name "kappa" in this codebase and should
not be conflated:
- **κ = 8πG/c⁴**: the formal tensor correspondence coupling per Paper A v14.5; a fixed physical
  constant derived from the field equation structure.
- **`tpf_kappa`**: the legacy simulator config knob (default 1e32) used by `direct_tpf` and the
  radial-Poisson paths above. It is a numerical parameter specific to those legacy paths, not the
  physical κ.

`geodesic_correspondence` uses G directly through the analytic Poisson solution for point sources.
`tpf_kappa` is irrelevant to it.

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
catching the old string `"v11_weak_field_truncation"` as a deprecation alias using
**config-key-presence semantics** (not float equality):

- If the config does **not** explicitly set `alpha_si` (key absent or default): emit a deprecation
  notice and forward to `geodesic_correspondence`.
- If the config **explicitly sets** `alpha_si`: honor it with a "free-parameter legacy mode"
  warning; do not forward.

`geodesic_correspondence` ignores any `alpha_si` config value entirely; it uses `−G_SI` internally.
The presence check must use the config-loading layer's key-presence flag, not a numerical
comparison — float equality against `−G_SI` is fragile across config-loading paths.

No new private method is required; `compute_v11_weak_field_truncation_accelerations()` is renamed
to `compute_geodesic_correspondence_accelerations()` in place. Also update `"geodesic_correspondence"`
in any validation or manifest enumerations. The `init_from_config()` block at line 113 already
handles a similar alias rename for `"weak_field_correspondence"` and serves as the pattern.

### B.2 Step 3 — ρ from linear flux divergence of Ξ (Paper A v14.4 Appendix D)

**Supersedes the earlier Θ_ij Θ^ij approach.** Paper A v14.4 Appendix D establishes the
correct correspondence as the linear flux divergence:

```
ρ_Ξ = (1/4π) ∂_i Ξ^i
```

The simulator uses **softened** kernels (Plummer softening with parameter ε). The actual ansatz
implemented in `provisional_point_source_field()` (verified: `xi.x = m * dx / R³` where
`R² = dx² + dy² + eps²`) is:

```
Ξ^i = M x^i / (r² + ε²)^(3/2)
```

**Singular limit (ε → 0):** `Ξ^i → M x^i / r³`, and distributionally:
- Exterior (r > 0): `∂_i Ξ^i = 0`
- At source: `∂_i Ξ^i = 4π M δ³(x)` (Gauss's theorem)
- Therefore `ρ_Ξ = M δ³(x)` — recovers the point mass exactly

**Softened implementation (ε > 0):** Taking the divergence of the softened kernel:
```
∂_i Ξ^i = 3Mε² / (r² + ε²)^(5/2)
```
giving a regularized Plummer-like density:
```
ρ_ε(r) = (1/4π) ∂_i Ξ^i = 3Mε² / [4π (r² + ε²)^(5/2)]
```
This integrates to M over all space (standard Plummer integral) and produces the standard
softened Newtonian acceleration. It is **not** literally M δ³(x); it is a smooth regularization
that approaches M δ³(x) as ε → 0. No free calibration parameter is introduced — ε is already
fixed by the existing softening policy — but the density is ε-dependent.

The quadratic invariant `Θ_ij Θ^ij` is correctly relegated to the nonlinear energy-ledger sector
per v14.4 §IV.A.d; it is not part of the geodesic_correspondence density step.

**Code status:** No new code is needed for the multi-source case. The divergence is linear, so
`∂_i Ξ^i_total = Σ_j ∂_i Ξ^i_j`. The runtime route uses the analytic closed-form shortcut (see
B.5); `derived_invariant_I_contracted()` is **not used** by this step.

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
   `"v11_weak_field_truncation"` to `"geodesic_correspondence"`; add a deprecation alias branch
   using **config-key-presence** (not float equality) to forward or run in legacy mode (see B.1).
   ~10–15 lines.
2. **Lock `alpha_si`**: `geodesic_correspondence` ignores `alpha_si` from config entirely;
   hardcode `−TPF_G_SI` at the call site. ~5 lines.
3. **Config validation**: Update allowlist of valid `tpf_dynamics_mode` strings to include
   `"geodesic_correspondence"`; retain `"v11_weak_field_truncation"` as a known alias.
4. **Contamination guards**: Add rejection/warning for incompatible config keys (see B.6). ~10–15 lines.
5. **Analytic shortcut note**: Add a documentation block making clear that the runtime route uses
   the analytic closed-form shortcut for point sources. It does **not** numerically compute ρ_Ξ,
   solve Poisson, or differentiate Φ on a grid. The five-step chain (Ξ → ρ_Ξ → Poisson →
   geodesic) is the theoretical grounding; the runtime collapses it via the analytic Poisson
   identity for point sources.
6. **Manifest entry**: Record `acceleration_code_path = "geodesic_correspondence: rho_Xi -> Poisson analytic -> geodesic"` in run_info output. ~5 lines.

### B.6 Contamination Guards

`geodesic_correspondence` is a paper-facing route. The following config keys must be actively
rejected (throw) or reported as inactive with a loud warning to prevent accidental contamination
of paper-facing runs by leftover experimental settings:

| Config key | Condition that triggers guard | Behavior |
|---|---|---|
| `tpf_vdsg_coupling` | value ≠ 0 | throw: "geodesic_correspondence rejects nonzero tpf_vdsg_coupling" |
| `tpf_global_accel_shunt_enable` | true | throw: "geodesic_correspondence rejects tpf_global_accel_shunt_enable=true" |
| `tpf_cooling_fraction` | value > 0 | throw: "geodesic_correspondence rejects positive tpf_cooling_fraction" |
| `tpfcore_enable_provisional_readout` | true (if it alters acceleration) | throw: "geodesic_correspondence rejects tpfcore_enable_provisional_readout=true" |

Fail loudly (throw) rather than silently producing contaminated results. The pattern is already
established in the existing `v11_weak_field_truncation` dispatch block (lines ~848–873 of
`tpf_core_package.cpp`), which uses the same throw-on-incompatible-config approach.

---

## C. Effort Estimate

| Item | Lines | Files |
|------|-------|-------|
| Rename dispatch string + config-key-presence alias guard (B.1) | 15–20 | `tpf_core_package.cpp` |
| Lock `alpha_si` internally; ignore from config for this mode (B.5) | 5–8 | `tpf_core_package.cpp`, `config.hpp` |
| Config validation allowlist update | 3–5 | `config.hpp`, `defaults.cfg` |
| Contamination guards — four throw blocks (B.6) | 10–15 | `tpf_core_package.cpp` |
| Analytic shortcut documentation block + manifest entry (B.5) | 15–20 | `tpf_core_package.cpp` |
| Unit and integration tests (seven tests below) | 80–110 | test files (see below) |
| **Total** | **~130–180 lines** | **4–5 files touched** |

**Required tests:**

1. `geodesic_correspondence` acceleration equals `NewtonianPackage` acceleration within numerical
   tolerance for a one-body point-source setup.
2. `geodesic_correspondence` equals `NewtonianPackage` for a galaxy configuration with
   `star_star=false`.
3. `geodesic_correspondence` equals `NewtonianPackage` for a galaxy configuration with
   `star_star=true` (use small N to keep the test fast).
4. Setting `alpha_si` in config does not affect `geodesic_correspondence` acceleration (parameter
   is silently ignored for this mode).
5. `"v11_weak_field_truncation"` mode string still runs (when `alpha_si` is not explicitly set)
   and emits a deprecation notice in stderr.
6. Nonzero `tpf_vdsg_coupling`, `tpf_global_accel_shunt_enable=true`, or positive
   `tpf_cooling_fraction` is rejected with a thrown exception when `geodesic_correspondence`
   is selected.
7. `manifest` / `run_info` records
   `acceleration_code_path = "geodesic_correspondence: rho_Xi -> Poisson analytic -> geodesic"`.

**Tests at risk:** `test_routing_semantics.cpp` "rejects unknown mode" test may need updating if it
uses a positive allowlist. `tpfcore_vdsg_manifest_smoke.sh` and
`tpfcore_legacy_vdsg_canonical_smoke.sh` should be checked for hardcoded mode strings.

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
are not used by `geodesic_correspondence`. The density step is handled via the analytic
closed-form shortcut; `tpf_kappa` (the legacy config knob, not the physical κ = 8πG/c⁴) is
relevant only to the `direct_tpf` and radial-Poisson paths. See A.4 for the κ terminology
distinction.

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
3. **Keep `"v11_weak_field_truncation"` as a deprecation alias** using config-key-presence
   semantics (not float equality — float equality against `−G_SI` is fragile across config-loading
   paths): if the config does **not** explicitly set `alpha_si`, emit a deprecation warning and
   forward to `geodesic_correspondence`; if the config **explicitly sets** `alpha_si`, honor it
   with a "free-parameter legacy mode" warning.
4. **Add a documentation block** at the top of the renamed method stating: the theoretical chain
   is Ξ → ρ_Ξ → Poisson → geodesic (Paper A v14.4 Appendix D; TPF_FOUNDATIONS.md §4). The
   runtime uses the analytic closed-form shortcut for point sources and does **not** numerically
   compute ρ_Ξ, solve Poisson, or differentiate Φ on a grid. The five-step chain is theoretical
   grounding; the runtime collapses it via the analytic Poisson identity. In the softened
   implementation ρ_ε(r) = 3Mε²/[4π(r²+ε²)^(5/2)] integrates to M; Θ_ij Θ^ij is the nonlinear
   energy-ledger quantity (v14.4 §IV.A.d), not the density proxy.

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

**Default not changed in Stage 1.** Changing the simulator default from `xi_kernel_deformed` to
`geodesic_correspondence` is a separate future PR. Existing configs that explicitly set
`tpf_dynamics_mode = xi_kernel_deformed` continue to work with no migration required. Paper-facing
validation runs must explicitly set `tpf_dynamics_mode = geodesic_correspondence`; there is no
automatic migration.

### E.4 Implementation Chain Summary

The five steps of the route and their code grounding:

1. **Ξ** — computed by `provisional_point_source_field()` (`source_ansatz.cpp:12–37`), Paper A Eq. (19).
2. **Θ** — computed by same function, Paper A Eqs. (20)–(25).
3. **ρ_Ξ** — linear flux correspondence `ρ_Ξ = (1/4π) ∂_i Ξ^i` (Paper A v14.4 Appendix D).
   Singular limit: M δ³(x). Softened implementation: ρ_ε(r) = 3Mε²/[4π(r²+ε²)^(5/2)], integrates
   to M. **Runtime uses the analytic closed-form shortcut; ρ_Ξ is not numerically computed.**
4. **Poisson / Φ** — analytic identity `Φ = −G Σ M_i / |x − x_i|` (no grid solve);
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
| `engine/config.hpp` | `tpf_kappa` | 215 | Legacy config knob (default 1e32) for direct_tpf/radial-Poisson paths; distinct from physical κ = 8πG/c⁴; not used by geodesic_correspondence |
