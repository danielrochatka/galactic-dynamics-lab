# Theory–Implementation Correspondence Audit
## Transformational Physics Framework (TPF) Gravitational Simulator

**Audit date:** 2026-05-02  
**Auditor:** Claude (read-only; no source files modified)  
**Branch audited:** `experimental/vdsg-next`  
**Paper A:** "A Configuration-Gradient Formulation for Static Weak-Field Gravitation" (v14.3)  
**Paper B:** "Gravitational Energy Reframed: A Transformational Physics Framework Approach" (v11)

---

## Preamble

This audit compares the mathematical content of both papers against the simulator source code.
All findings are read-only observations; no files were modified.
Paper A is deliberately scoped to the static weak-field sector.
Paper B is the broader earlier version.
The simulator is demonstrably ahead of Paper A on several axes.

---

## Revision Notes (Second-Pass Revision)

This document is a second-pass revision of the original 2026-05-01 audit. The first pass contained
material errors traceable to trusting code comments and file-level documentation without
independently verifying function bodies. Code comments in this codebase are **not authoritative**:
some were written during a period when the author had been given incorrect guidance about which
paths were "theoretically grounded." Several first-pass findings were inverted (a correct path was
labelled blocking; an incorrect path was labelled a match).

**Ground rules applied in this revision:**
- Behavioral claims are either "verified by reading function body" or explicitly flagged as
  "stated in comment (not independently verified)."
- Code comments are treated as suspect unless corroborated by the function body.
- Function names are not trusted — the function body is read.

**Author corrections incorporated:**
1. Ξ_μ (Xi / Configuration Displacement Field) is the primary field in both papers. The correct
   acceleration readout is `a = −K * Ξ` (direct from the displacement field, Newtonian r⁻² falloff).
2. C_μν (Paper A Eq. 15) is the **field equation**, not an acceleration readout. Using C_μν
   directly as acceleration produces r⁻⁶ falloff — physically wrong for gravity.
3. Paper A Eq. (29), which claims `a = −(C_xx û_x + C_xy û_y)`, is incorrect. The author has
   identified this as an AI-injected error and will revise the paper.
4. `direct_tpf` uses C_μν as an acceleration readout and is the misinterpretation. It should be
   quarantined or renamed.
5. `xi_kernel_deformed` at `xi_kernel_coupling = 0` is the paper-grounded Ξ readout. At
   coupling > 0 it adds VDSG deformation — a legitimate research extension by the author.

---

## Finding 1: Lambda (λ) — correct value, correct locus

**Category:** match-A  
**Severity:** informational  
**Paper A reference:** Eq. (8)–(9), λ = 1/n = 1/4 in n = 4  
**Paper B reference:** Eq. (4), same  
**Code location:** `engine/physics/TPFCore/source_ansatz.hpp:58`, `v11_weak_field_correspondence.cpp:23`

### Description
`constexpr double LAMBDA_4D = 0.25` is declared once in `source_ansatz.hpp` and imported wherever λ appears.
`v11_weak_field_correspondence.cpp` independently hard-codes `kLambda4D = 0.25` with the comment "Paper v11: λ = 1/n in n = 4 dimensions (fixed; not tunable)".
The `direct_tpf_baseline.cpp` receives λ at call time via `tpfcore::LAMBDA_4D` (line 573 of `tpf_core_package.cpp`).

### Theory comparison
Both papers agree: λ = 1/4 in 4D spacetime dimensions.
Paper A Box (9), Paper B Box (4).
Code matches exactly. The constant is not tunable, which is correct per both papers.

### Recommendation
No action required. Add an explicit assertion `static_assert(LAMBDA_4D == 0.25)` if λ-consistency across translation units is ever a concern.

---

## Finding 2: Theta tensor — 3D Hessian of a softened potential (provisional ansatz)

**Category:** match-A (for what the paper calls the "current simulator realization")  
**Severity:** minor  
**Paper A reference:** Eqs. (18)–(25) (Sec. IV.A — "Current simulator realization of the weak-field sector")  
**Paper B reference:** Eq. (15), Θ_{ij} = ∂_i ∂_j φ for the curl-free weak-field case  
**Code location:** `engine/physics/TPFCore/source_ansatz.cpp:12–37`, `engine/physics/TPFCore/source_ansatz.hpp:34`

### Description
`provisional_point_source_field()` computes:
```
Xi_x  = m * dx / R^3,   Xi_y  = m * dy / R^3      (Eq. 19 Paper A)
Theta_xx = m*(r2 - 3*dx^2) / R^5                    (Eq. 20 Paper A)
Theta_xy = -3*m*dx*dy / R^5                          (Eq. 21 Paper A)
Theta_xz = 0                                         (Eq. 22 Paper A — by symmetry on z=0)
Theta_yy = m*(r2 - 3*dy^2) / R^5                    (Eq. 23 Paper A)
Theta_yz = 0                                         (Eq. 24 Paper A)
Theta_zz = m*r2 / R^5  = m/R^3                      (Eq. 25 Paper A)
```
where `R^2 = dx^2 + dy^2 + eps^2` (softened, Eq. 18 Paper A).

### Theory comparison
The code exactly reproduces Paper A Eqs. (20)–(25) including Theta_zz = m/R^3.
Paper A explicitly flags this as a "provisional source ansatz" rather than an on-shell solution.
Paper B Eq. (15) gives the same Θ_{ij} = ∂_i ∂_j φ in the weak-field curl-free case.
The softening ε breaks the exact vacuum harmonic condition; Paper A acknowledges this (Sec. IV.A paragraph 4).

### Recommendation
The code matches the paper's stated provisional ansatz correctly.
The label "provisional" in both paper and code is accurate and should be preserved in any publication description.

---

## Finding 3: Invariant I — 3D Euclidean sum including Theta_zz (NOT a true 4D sum)

**Category:** divergence  
**Severity:** concerning  
**Paper A reference:** Eq. (5) I = Θ_{μν} Θ^{μν} − λΘ² (4D covariant); Eq. (34) reduces to I ≃ Θ_{ij}Θ^{ij} − λΘ² in static limit  
**Paper B reference:** Eq. (3) I = Θ_{μν}Θ^{μν} − λΘ², same definition  
**Code location:** `engine/physics/TPFCore/source_ansatz.cpp:74–79`

### Description
```cpp
double compute_invariant_I(const Theta3D& theta) {
  double Theta_mn_Theta_mn =
      theta.xx*theta.xx + theta.yy*theta.yy + theta.zz*theta.zz +
      2.0*(theta.xy*theta.xy + theta.xz*theta.xz + theta.yz*theta.yz);
  double Theta_trace = theta.trace();  // xx + yy + zz
  return Theta_mn_Theta_mn - LAMBDA_4D * Theta_trace * Theta_trace;
}
```
The sum runs over all nine components of the 3×3 spatial Hessian using the 3D Frobenius norm (all positive signs).
It does NOT include temporal components (Θ_{tt}, Θ_{0i}).

Note: `compute_invariant_I` is called only within the `direct_tpf` path and for the static residual benchmark. The `xi_kernel_deformed` path (the paper-grounded dynamics path at coupling=0) does NOT call this function — verified by reading `compute_xi_kernel_deformed_accelerations` function body.

### Theory comparison
Paper A Eq. (34) explicitly says: "because the spatial sector dominates in the static approximation, I ≃ Θ_{αβ}Θ^{αβ} − λΘ² ≃ Θ_{ij}Θ^{ij} − λΘ²". The simulator's spatial Frobenius sum is consistent with this approximation.

Since `xi_kernel_deformed` (the correct dynamics path) does not use I at all, this concern applies only to `direct_tpf` (which should be quarantined per Finding 13) and the static residual benchmark.

### Recommendation
Document in code that `compute_invariant_I` computes the static spatial approximation I ≃ Θ_{ij}Θ^{ij} − λΘ² (Paper A Eq. 34), valid under static + weak-field assumptions. Note that it is not used in the paper-grounded `xi_kernel_deformed` path. The severity of this finding is moderated by the fact that the correct dynamics path is unaffected.

---

## Finding 4: Field equation C_{μν} — computation correct AS A FIELD EQUATION

**Category:** match-A (explicitly partial)  
**Severity:** informational  
**Paper A reference:** Eq. (15) — full C_{μν}; Eqs. (26)–(28) — planar principal-part components  
**Paper B reference:** Eq. (10) — full C_{μν} with ΔC_{μν}  
**Code location:** `engine/physics/TPFCore/direct_tpf_baseline.cpp:29–65`

### Description
`compute_principal_Cij_from_eq10_baseline()` implements exactly Paper A Eqs. (26)–(28) — verified by reading function body:
```
C_xx^(principal) = κ (Θ_xx² + Θ_xy² + Θ_xz² − λΘΘ_xx − ½I)
C_xy^(principal) = κ (Θ_xx Θ_xy + Θ_xy Θ_yy + Θ_xz Θ_yz − λΘΘ_xy)
C_yy^(principal) = κ (Θ_xy² + Θ_yy² + Θ_yz² − λΘΘ_yy − ½I)
C_zz^(principal) = κ (Θ_xz² + Θ_yz² + Θ_zz² − λΘΘ_zz − ½I)
```
`compute_deltaC_placeholder_zero()` returns all zeros — ΔC_{μν} = 0.

The computation of C_μν correctly implements Paper A Eqs. (26)–(28) **as a field equation**. This is fine — C_μν is indeed the field equation (Paper A Eq. 15). The computation is not wrong.

The ERROR is using C_μν as an acceleration readout (see Finding 13 and Finding 21). C_μν ~ κ*(m/r³)² gives r⁻⁶ falloff, not r⁻² gravity. C_zz is computed but unused in the 2D readout.

### Theory comparison
Paper A explicitly retains ΔC_{μν} symbolically but does not expand it (Sec. IV.A penultimate paragraph).
Paper B Eq. (10) also keeps ΔC_{μν} symbolic.
The C_μν computation is a correct implementation of the field equation principal-part. It is appropriate for residual/diagnostic purposes. It must not drive dynamics.

### Recommendation
The C_μν computation in `direct_tpf_baseline.cpp` can be retained as a diagnostic (field equation residual). C_zz is computed but unused in the 2D readout — document as dead code or remove from the acceleration path. Flag C_μν output clearly as a "field equation value, not an acceleration."

---

## Finding 5: Paper A Eq. (29) acceleration readout — PAPER ERROR

**Category:** divergence (paper error, not code error)  
**Severity:** blocking (paper must be revised before submission; the code implementing this is `direct_tpf`)  
**Paper A reference:** Eq. (29) — a_x = −(C_xx^(p) Ξ̂_x + C_xy^(p) Ξ̂_y)  
**Paper B reference:** No corresponding acceleration readout in Paper B  
**Code location:** `engine/physics/TPFCore/direct_tpf_baseline.cpp:52–65`

### Description
Paper A Eq. (29) states that the acceleration readout is obtained by projecting C_μν onto the normalized Ξ direction. The code in `compute_xi_directed_tensor_readout` implements this exactly — verified by reading function body:
```cpp
readout.ax = -(principal_cij.c_xx * u_x + principal_cij.c_xy * u_y);
readout.ay = -(principal_cij.c_xy * u_x + principal_cij.c_yy * u_y);
```
where u_x = ξ.x / |ξ|, u_y = ξ.y / |ξ|.

However, Paper A Eq. (29) is **incorrect**. The author has identified this as an AI-injected error
in the paper that swapped tensor math for Newtonian math. C_μν is the field equation (Paper A
Eq. 15), not an acceleration source. Projecting C_μν onto Ξ̂ gives r⁻⁶ falloff (see Finding 13
and Finding 21). The correct readout is `a = −K * Ξ`, which is what `xi_kernel_deformed` at
coupling=0 implements (see Finding 7).

**The code is correct on this point only insofar as it faithfully implements Eq. (29). The equation
itself is wrong. Paper A Eq. (29) requires revision.**

### Theory comparison
The first-pass audit labelled this finding "match-A / informational." That was wrong: it matched a
paper equation that is itself erroneous. The match to an erroneous equation does not constitute a
correct implementation of TPF dynamics.

### Recommendation
Revise Paper A Eq. (29) before submission. The `direct_tpf_baseline.cpp` implementation of this
equation should be quarantined as `delta_C_diagnostic` (see Finding 13). The paper-grounded
acceleration readout is `a = −K * Ξ`, implemented correctly by `xi_kernel_deformed` at coupling=0.

---

## Finding 6: 4D-to-2D Projection Audit — Θ_zz systematic (diagnostic concern only)

**Category:** ambiguity  
**Severity:** minor (diagnostic concern; does not affect the correct dynamics path)  
**Paper A reference:** Eq. (3) — covariant Θ_{μν} = ∇_μ Ξ_ν; Eq. (25) Θ_zz = m/R³; Sec. IV.A explicitly describes 2D reduction  
**Paper B reference:** No explicit 2D projection discussion  
**Code location:** `engine/physics/TPFCore/source_ansatz.cpp:12–36`, `engine/physics/TPFCore/direct_tpf_baseline.cpp:29–65`

### Description
The path from 4D theory to 2D simulation via Θ, I, and C_μν involves a systematic where Θ_zz ≠ 0 on the z = 0 plane (Θ_zz = m/R³, Paper A Eq. 25). This enters the invariant I (and hence C_xx, C_yy via the ½I term) but has no corresponding acceleration component along z.

**Key limitation of this finding's scope:** The `xi_kernel_deformed` path — the paper-grounded
dynamics path — does NOT compute Θ, I, or C_μν at all (verified by reading
`compute_xi_kernel_deformed_accelerations` function body). The Θ_zz systematic is therefore
relevant only to `direct_tpf` (which should be quarantined) and to C_μν residual diagnostics.

The ~17% Θ_zz contribution to I (Θ_zz² / I_total ~ 1/6 near the source) would shift the C_μν
force magnitude in `direct_tpf` by that fraction, but since `direct_tpf` uses C_μν as an
acceleration readout — which is itself wrong — the systematic is a second-order concern relative
to the primary error in that path.

### Theory comparison
Paper A explicitly describes and authorizes this 2D-plane implementation in Sec. IV.A, calling it a "planar provisional ansatz." The Θ_zz-in-I effect is a known, paper-acknowledged feature of the provisional ansatz, not an undocumented deviation.

### Recommendation
For C_μν residual diagnostic work using `direct_tpf_baseline`, document the ~17% Θ_zz contribution to I as a known systematic. For the paper-grounded `xi_kernel_deformed` dynamics path, this finding is not applicable.

---

## Finding 7: Active dynamics path — xi_kernel_deformed — PAPER-GROUNDED AT coupling=0

**Category:** match-A (at coupling=0) / extension-research (at coupling>0)  
**Severity:** informational (coupling=0 case) / minor (coupling>0 VDSG case)  
**Paper A reference:** Eq. (2) X^μ(x) = x^μ + Ξ^μ(x); Eq. (19) Ξ_i = m·d_i/R³ (the Xi field)  
**Paper B reference:** same Ξ field definition  
**Code location:** `engine/physics/TPFCore/tpf_core_package.cpp:650–799`, `engine/physics/TPFCore/tpf_core_package.hpp:202–224`

### Description

**At `xi_kernel_coupling = 0` (pure Ξ readout):**

Verified by reading `accumulate_source` lambda in `compute_xi_kernel_deformed_accelerations`
(lines ~698–760):
- `factor_raw = xi_kernel_coupling_ * (gamma_rel - 1.0)` — when `xi_kernel_coupling_ == 0.0`,
  `factor_raw = 0.0` regardless of velocity or kernel mode
- With `factor_raw = 0`: `1/(1 + 0) = 1.0`, so `metric_scale = 1.0` (within clamp bounds)
- `alpha = (metric_scale - 1.0) = 0.0`, therefore `gx = dx`, `gy = dy`
- `xi_sx = src.mass * dx * inv_r3 = src.mass * dx / r³`
- `dax = -xi_motion_readout_scale_ * xi_sx = -K * src.mass * dx / r³`

This is pure Ξ readout: `a = −K * Ξ`, giving Newtonian r⁻² force. This holds for ALL
`xi_kernel_mode_` values when `xi_kernel_coupling_ == 0`. This is the paper-grounded path per
Paper A Eq. (2) and Eq. (19): the acceleration reads directly from the Ξ field.

**At `xi_kernel_coupling > 0` (VDSG deformation active):**

The metric_scale becomes velocity-dependent. The displacement vector (dx, dy) is deformed into
(gx, gy) before computing Ξ. This adds velocity-dependent corrections — a legitimate research
extension by the author exploring VDSG dynamics beyond the static-scope Paper A.

The active rotation-curve config (`configs/my.local.cfg`) uses `xi_kernel_coupling = 5000` and
`xi_kernel_mode = metric_transverse_wake`. This is the VDSG research extension, not the baseline
coupling=0 TPF path.

### Theory comparison
At coupling=0: **confirmed paper-grounded Ξ readout** (verified by reading function body).
At coupling>0: velocity-dependent VDSG deformation, a legitimate research extension by the author.
Paper A is static-scoped so velocity terms are outside its stated scope, but this does not make
them physically wrong — it means Paper A does not yet describe them. Paper B's dynamical terms
differ in form (Θ_{tt}/Θ_{rr}) but do not prohibit this kind of VDSG extension.

The first-pass audit labelled this finding "blocking / not backed by either paper." That assessment
was wrong and stemmed from trusting misleading code comments without reading function bodies. The
coupling=0 path is the primary paper-grounded path.

### Recommendation
For paper-facing runs testing Newtonian/TPF weak-field correspondence: use
`xi_kernel_coupling = 0` (or confirm it is 0 in the config). This gives the clean
`a = −K * Ξ` readout that Paper A derivation leads to.
For VDSG research runs: keep coupling>0 but label results as "VDSG research extension."
See also Findings 8–10 for the specific sub-parameters of the VDSG deformation.

---

## Finding 8: xi_kernel_factor_mode — "gamma_minus_one" (special relativity, not in papers)

**Category:** code-only  
**Severity:** concerning  
**Paper A reference:** none (static scope; all velocity terms excluded)  
**Paper B reference:** none (no relativistic factor of this form derived)  
**Code location:** `engine/physics/TPFCore/tpf_core_package.cpp:696–700`

### Description
```cpp
const double gamma_rel = 1.0 / std::sqrt(1.0 - beta_for_gamma * beta_for_gamma);
if (xi_kernel_factor_mode_ == "beta_power")
  factor_raw = xi_kernel_coupling_ * std::pow(beta_effective, xi_kernel_beta_power_);
else   // "gamma_minus_one"
  factor_raw = xi_kernel_coupling_ * (gamma_rel - 1.0);
```
The `gamma_minus_one` factor is the standard Lorentz kinetic-energy factor (γ−1). It is zero at rest and grows as v²/2c² for non-relativistic speeds. For galaxy-scale speeds (v ~ 200 km/s, β ~ 7×10⁻⁴), γ−1 ~ 2.5×10⁻⁷. Multiplied by coupling 5000, factor_raw ~ 1.25×10⁻³ — small but nonzero.

This factor only activates when `xi_kernel_coupling > 0` (the VDSG extension). At coupling=0, factor_raw=0 and this mode has no effect (verified by reading function body).

### Theory comparison
This is a research ansatz. No TPF derivation in either paper produces a factor of (γ−1) as a kernel modifier. In GR, gravitomagnetic effects enter at order v/c through metric perturbation terms. Paper B's Earth-Moon example uses Θ_{rr} + Θ_{tt} (Eq. 48) — a trace condition, not a Lorentz boost.

### Recommendation
Document xi_kernel_factor_mode = "gamma_minus_one" as **exploratory VDSG extension, not derived from either paper**. If the intent is to approximate a weak gravitomagnetic correction, that needs a derivation from the TPF field equations. Only active when xi_kernel_coupling > 0.

---

## Finding 9: xi_kernel_coupling = 5000 — no physical derivation

**Category:** code-only  
**Severity:** concerning  
**Paper A reference:** none  
**Paper B reference:** none  
**Code location:** `configs/my.local.cfg:108`

### Description
The config comment says "At realistic SI speeds, beta_rel may be much larger than in toy tests, so start lower than the toy-scale 1e8 value." This is explicitly trial-and-error tuning language. The coupling 5000 is dimensionless when used with gamma_minus_one mode. No calibration procedure, no correspondence target, no physical argument appears anywhere in the code or comments.

### Theory comparison
κ = 16πG in Papers A and B (the actual theory coupling). The xi_kernel_coupling is an entirely separate dimensionless VDSG deformation strength with no paper analog. The two parameters are in different subsystems and serve different purposes.

### Recommendation
Either derive xi_kernel_coupling from first principles (e.g., from a desired magnitude of velocity-dependent correction at a fiducial radius and speed) or clearly label it as a fit parameter. Do not confuse it with κ from the papers.

---

## Finding 10: metric_clamp (0.5–2.0) — numerical regularization

**Category:** code-only (regularization)  
**Severity:** minor  
**Paper A reference:** none  
**Paper B reference:** none  
**Code location:** `configs/my.local.cfg:117–118`, `engine/physics/TPFCore/tpf_core_package.cpp:706–709`

### Description
```
tpf_4d_xi_kernel_metric_min = 0.5
tpf_4d_xi_kernel_metric_max = 2.0
metric_scale = clamp(metric_min, metric_max, 1.0 / (1.0 + factor_raw))
```
The clamp prevents the effective metric deformation from going below 0.5 or above 2.0 — i.e., forces can be at most doubled or halved relative to the undeformed case, regardless of velocity.

At coupling=0, factor_raw=0 so metric_scale=1.0, which is within the clamp range — the clamp has no effect on the paper-grounded path (verified by reading function body).

### Theory comparison
Neither paper describes any metric deformation clamping. This is a numerical stability measure for the VDSG extension. The config comment "Clamp metric deformation so the run does not instantly explode" confirms its regularization purpose.

### Recommendation
Classified as regularization. Document explicitly that the clamp has no theoretical justification and is a numerical safeguard for the VDSG extension. Verify that clamp events are rare in reported runs.

---

## Finding 11: velocity-dependent terms (VDSG) — legitimate research extension

**Category:** extension-research (exploratory, beyond Paper A static scope)  
**Severity:** concerning (for direct paper-attribution of results)  
**Paper A reference:** none (static scope explicitly)  
**Paper B reference:** none (Paper B's dynamical terms are Θ_{tt}/Θ_{rr}/Θ_{tr}, not a Doppler factor)  
**Code location:** `engine/physics/TPFCore/extensions_vdsg.cpp:43–163`

### Description
The VDSG (velocity-dependent structural gradient) extension computes an additive force modification:
```
factor = f(v_rel, beta_rad, mass_gate, weak_field_gate)
delta_a = softened_Newtonian_force * (factor - 1)
```
Modes include:
- `legacy_speed`: factor = 1 + λ_eff * |v_rel|/c  (linear in total speed)
- `radial_doppler_rational`: factor = 1/(1 − x) where x ∝ radial β
- `radial_doppler_exp`: factor = exp(x)
- `radial_doppler_bounded`: factor = 1 + A * tanh(x/A)

The VDSG subsystem is NOT active in the current rotation-curve config (my.local.cfg sets `tpf_vdsg_coupling = 0`).

VDSG is a legitimate research extension by the author — it explores velocity-dependent corrections in the TPF framework that Paper A does not yet describe (because Paper A is static-scoped). It is not an unauthorized invention; the author is aware of its status as a research extension beyond the current paper scope.

### Theory comparison
Paper B Section XI and Eq. (48) use Θ_{rr} and Θ_{tt} in a heuristic time-stepped collapse equation. These are tensor-component-based, not Doppler-factor-based. There is no TPF derivation for the Doppler-rational or exponential factors in the current papers. However, absent does not mean wrong — it means the derivation is future work.

### Recommendation
Classify VDSG explicitly as "exploratory research extension beyond Paper A static scope." Ensure that any published rotation-curve results explicitly state whether VDSG was active. The code's existing "VDSG: exploratory additive extension" labeling is appropriate.

---

## Finding 12: v11_weak_field_truncation — does "v11" refer to Paper B (v11)?

**Category:** match-B  
**Severity:** informational  
**Paper A reference:** Eqs. (36)–(39) — weak-field Poisson correspondence  
**Paper B reference:** Sec. II.B (Newtonian Limit), Eqs. (42)–(47) for Earth-Moon benchmark  
**Code location:** `engine/physics/TPFCore/v11_weak_field_correspondence.cpp:22`, `engine/physics/TPFCore/tpf_core_package.cpp:587–594`

### Description
The comment in `v11_weak_field_correspondence.cpp:22` reads: "Paper v11: λ = 1/n in n = 4 dimensions."
The function `compute_v11_weak_field_correspondence_accelerations` implements:
```
a_x = alpha_SI * m * dx / R³
```
where `alpha_SI = -G` by default. This is Newtonian gravity with coupling alpha instead of G.
The Earth-Moon benchmark (`run_earth_moon_line_correspondence_benchmark`) implements Paper B Eqs. (42)–(47) explicitly, including:
- φ(x) = α(M_E/x + M_M/(D−x))    [Paper B Eq. 44]
- a_TPF(x) = −dφ/dx + Ω²(x − x_b) [Paper B Eq. 45]

### Theory comparison
"v11" explicitly refers to Paper B (which is dated v11). The correspondence helper correctly references and implements Paper B's weak-field reduction.
Paper A's content (Eqs. 38–39) gives the same correspondence statement.
This path and `xi_kernel_deformed` at coupling=0 both implement `a = −G * m * dx / R³` — they converge to the same Ξ readout.

### Recommendation
This path correctly implements the paper-backed weak-field correspondence. It is NOT used for rotation-curve work (my.local.cfg uses xi_kernel_deformed). Both this path and xi_kernel_deformed at coupling=0 are paper-grounded for the weak-field static case.

---

## Finding 13: direct_tpf path — MISINTERPRETATION, QUARANTINE RECOMMENDED

**Category:** divergence  
**Severity:** blocking  
**Paper A reference:** Eq. (15) C_μν (field equation); Eqs. (26)–(28) principal-part computation; Eq. (29) (INCORRECT — author will revise)  
**Paper B reference:** Eq. (10) — full C_{μν}  
**Code location:** `engine/physics/TPFCore/tpf_core_package.cpp:562–584`, `engine/physics/TPFCore/direct_tpf_baseline.cpp`

### Description
Verified by reading `direct_tpf_baseline.cpp` function bodies:

1. `compute_principal_Cij_from_eq10_baseline`: computes C_{ij} = κ(Θ_{ik}Θ^{kj} − λ tr(Θ)Θ_{ij} − ½I δ_{ij})
   — this is the field equation (Paper A Eq. 15).
2. With Θ_{ij} ~ m/r³ (Hessian of 1/r potential), C_{ij} ~ κ*(m/r³)² = κm²/r⁶
3. `compute_xi_directed_tensor_readout`: `readout.ax = -(c_xx * u_x + c_xy * u_y)`
   — **C_μν used directly as acceleration**
4. Force ~ κm²/r⁶ — **r⁻⁶ falloff confirmed by reading function body**

This is physically wrong for gravity (should be r⁻²). The `direct_tpf` path is a misinterpretation
of the TPF field equation: it mistakes C_μν (the field equation) for an acceleration source.

The first-pass audit labelled this "match-A (dormant) / concerning" and recommended comparing it
against xi_kernel_deformed for rotation curves. That assessment was wrong. `direct_tpf` does not
produce TPF-consistent gravity; it produces r⁻⁶ force falloff.

Paper A Eq. (29), which appears to authorize using C_μν as an acceleration readout, is incorrect
(see Finding 5 and Finding 21). The author is aware and will revise the paper.

### Theory comparison
C_μν is the field equation (Paper A Eq. 15). It characterizes the configuration of the Ξ field.
It is not a force law. The correct acceleration comes from `a = −K * ∇Ξ = −K * Ξ` (for point
sources), which is what `xi_kernel_deformed` at coupling=0 computes.

### Recommendation
Quarantine `direct_tpf`:
1. Rename to `delta_C_diagnostic` or `delta_C_residual` to make clear it computes a field
   equation residual, not dynamics.
2. Remove it from valid `tpf_dynamics_mode` options for dynamics runs.
3. Add a compile-time static assertion or runtime error if `direct_tpf` is selected for any
   run that will be used to generate physical predictions.
4. The C_μν computation itself (`compute_principal_Cij_from_eq10_baseline`) is correct as a
   field-equation diagnostic and can be retained for that purpose.

---

## Finding 14: kappa (κ) — inconsistency is a consequence of direct_tpf being wrong

**Category:** divergence  
**Severity:** informational (relevant only to the incorrect direct_tpf path)  
**Paper A reference:** Eq. (38): κ = 16πG SI in the standard normalization; Box (39): Q_wf = ρ when κ = 16πG  
**Paper B reference:** Same calibration; Paper B Sec. II.B.c  
**Code location:** `engine/config.hpp:214` (default tpf_kappa = 1.0e32); `configs/my.local.cfg:84` (tpf_kappa = 1)

### Description
Paper A Eq. (38) sets κ = 16πG ≈ 2.09×10⁻⁹ SI to reproduce the Poisson sector.
The config default is `tpf_kappa = 1.0e32` — 41 orders of magnitude off.
The active rotation-curve config sets `tpf_kappa = 1`.

However, κ is used **only** by the `direct_tpf` path (verified by reading dispatch logic in
`compute_accelerations`). The paper-grounded path `xi_kernel_deformed` does not use κ at all.
Since `direct_tpf` should be quarantined (Finding 13), the κ inconsistency is a consequence of
that path being wrong, not an independent blocking issue.

For the correct path (`xi_kernel_deformed`), the analog of κ is `xi_motion_readout_scale` = G,
which is correctly set to 6.67430×10⁻¹¹ in my.local.cfg.

### Theory comparison
The κ = 1 mismatch in the active config is irrelevant to the `xi_kernel_deformed` dynamics.
If `direct_tpf` is retained as a diagnostic tool after quarantine, κ should be set to 16πG for
any residual-magnitude comparison.

### Recommendation
Downgraded from "blocking" to "informational" because the correct dynamics path does not use κ.
If `direct_tpf` is converted to a field-equation diagnostic, add a comment that κ = 16πG is
required for physically meaningful C_μν residual values.

---

## Finding 15: Temporal Xi mode — purely diagnostic, not active

**Category:** code-only (dormant)  
**Severity:** informational  
**Paper A reference:** none  
**Paper B reference:** none  
**Code location:** `engine/physics/TPFCore/tpf_core_package.hpp:218`, `configs/my.local.cfg:125–126`

### Description
`xi_temporal_mode` and `xi_temporal_coupling` appear in the config and class but are set to "off" / 0.0 in the active config. The code reads these parameters but the accumulate_source lambda only references them indirectly through wake kinematics. The "spacetime_metric" kernel mode would use temporal Xi components (Xi_t ≠ 0) but is not active.

### Theory comparison
Neither paper derives a temporal Xi contribution to the 2D planar force beyond the static approximation. Θ_{0μ} ≈ 0 is one of the stated static-limit assumptions in both papers.

### Recommendation
Keep as dormant exploratory code. Mark clearly as "not paper-grounded; exploratory extension for future dynamic regime testing."

---

## Finding 16: Energy ledger — Paper B Appendix E vs no code implementation

**Category:** extension-B (present in paper, absent in code)  
**Severity:** informational  
**Paper A reference:** none (energy ledger explicitly out of scope; Appendix A discusses dimensional conventions only)  
**Paper B reference:** Appendix E "Optional GR-Ledger Translation (Derived Bookkeeping)," Eqs. (6)–(11)  
**Code location:** none for energy ledger; `compute_potential_energy` always returns 0.0 (tpf_core_package.hpp:55)

### Description
Paper B Appendix E defines:
- ρ_ledger ≡ Θ_{ij}Θ^{ij}   [Eq. (7) Paper B App E]
- M = ∫ ρ_ledger d³x = ∫ Θ_{ij}Θ^{ij} d³x   [Eq. (8)]
- E = M c²   [Eq. (9)]
- E = c² ∫ Θ_{ij}Θ^{ij} d³x   [Eq. (10)]
- Noether charge E_Noether = ∫ J⁰ d³x   [Eq. (11)]

The code's `compute_potential_energy` returns 0.0 always (line 55 of `tpf_core_package.hpp`). No Θ-derived energy is logged anywhere.

### Theory comparison
Paper B makes the energy ledger explicit and derivable from the configuration gradients already computed. The required quantity Θ_{ij}Θ^{ij} is already computed as part of `compute_invariant_I` (without the −λΘ² term). The infrastructure to compute this exists; it just isn't exposed.

### Recommendation
Minimal energy ledger implementation:
1. In `compute_invariant_I`, additionally return the Frobenius norm squared Θ_{ij}Θ^{ij}.
2. In diagnostic output passes, sum Θ_{ij}Θ^{ij} over all sources and multiply by c² and an appropriate integration volume to estimate the configuration-gradient energy.
3. Log this as `E_ledger_J` alongside the standard kinetic energy.
This requires ~10 lines of additional output code and no changes to dynamics.

---

## Finding 17: Static residual benchmark (tpf_4d_static_residual) — not wired to dynamics

**Category:** match-A (for its stated scope — benchmark only)  
**Severity:** informational  
**Paper A reference:** Eq. (11) — field equation ∇_μ(Θ^μ_ν − λδ^μ_ν Θ) = 0  
**Paper B reference:** Eq. (9), same  
**Code location:** `engine/physics/TPFCore/tpf_4d_static_residual.cpp`

### Description
`evaluate_static_configuration_residual_4d` computes the full 4D spatial divergence of A^μ_ν = Θ^μ_ν − λδ^μ_ν Θ on a 3D grid, using finite differences. The residual measures how far the provisional ansatz is from satisfying the field equation. This runs as a benchmark harness but is explicitly "not wired into dynamics."

### Theory comparison
This correctly implements the field equation residual audit from Paper A Eq. (11) / Paper B Eq. (9) in 3D. The benchmark is honest: the softened Phi ansatz produces non-zero residuals near the source (as the code comments document), because a softened monopole is not an exact solution.

### Recommendation
No action required for dynamics. Consider running this benchmark with the actual Theta field from a production run to characterize how far off-shell the provisional ansatz is at different radii.

---

## Finding 18: Softening policy — numerical regularization only

**Category:** match-A (regularization acknowledged)  
**Severity:** informational  
**Paper A reference:** ε in Eqs. (18)–(25) — "source softening scale"; Paper A calls it "a regularized implementation scaffold"  
**Paper B reference:** none explicitly  
**Code location:** `engine/softening_policy.cpp`

### Description
Softening uses Plummer softening (`R² = r² + ε²`). The `resolve_softening` function implements collisionless softening (mean inter-particle separation / factor), with structural caps (inner_cap, outer_cap, min, max). The auto-softening is based on particle count and annulus geometry, not physics.

### Theory comparison
Paper A acknowledges softening as a numerical regularization. Neither paper requires a specific softening profile. The code's structural caps for galaxy mode are numerical stability measures with no theoretical counterpart.

### Recommendation
Document the softening as numerical regularization. For any precision weak-field test (e.g., xi_kernel_deformed at coupling=0 correspondence check), use analytical softening (fixed ε) rather than auto-softening to control the Ξ field profile systematically.

---

## Finding 19: Schwarzschild exterior and bounce metric — Paper B only, not in code

**Category:** extension-B (in paper, not in code)  
**Severity:** informational  
**Paper A reference:** none  
**Paper B reference:** Secs. VIII–X, Eqs. (27)–(39) — bounce metric, TPF horizon  
**Code location:** none

### Description
Paper B Secs. VIII–X construct an explicit regular-core ("bounce") metric:
f_TPF(r) = 1 − 2GMr² / (r³ + 2GMℓ²)
and a TPF horizon criterion based on divergence of I(r). Neither of these appears anywhere in the simulator code base.

### Theory comparison
These are static, spherically symmetric solutions to the full 4D field equation. The current simulator uses 2D planar dynamics with a provisional ansatz. The bounce metric and horizon construction are beyond the current simulator scope.

### Recommendation
Note for future work: the bounce metric provides a natural strong-field completion that could eventually replace the softened point-source ansatz as the near-source field profile. This is Paper B content deliberately cut from Paper A.

---

## Finding 20: Noether symmetry / conserved charge — Paper B only, not in code

**Category:** extension-B (in paper, not in code)  
**Severity:** informational  
**Paper A reference:** Eq. (16) — "structural conservation statement" from diffeomorphism invariance  
**Paper B reference:** Sec. V, Eqs. (8)–(11) in App. E — Noether charge Q, conserved current J^μ  
**Code location:** none

### Description
Paper B Section V discusses the conserved Noether charge Q = ∫ J⁰ d³x associated with time-translation invariance of S[g, Ξ]. This is the TPF analog of energy conservation. Paper A mentions a divergence identity ∇^μ C_{μν} = 0 (Eq. 16) but does not compute the Noether charge.

Neither paper's energy conservation/Noether charge is monitored in the simulator (see also Finding 16).

### Recommendation
Implementing the Noether charge monitor would both serve as a simulation sanity check and demonstrate Paper B's claim that energy conservation emerges from configuration-gradient invariance without introducing gravitational energy as a primitive.

---

## Finding 21: Paper A Eq. (29) contains an AI-injected error

**Category:** divergence  
**Severity:** blocking (paper must be revised before submission)  
**Paper A reference:** Eq. (29) — stated as `a_x = −(C_xx^(p) Ξ̂_x + C_xy^(p) Ξ̂_y)`  
**Code location:** `engine/physics/TPFCore/direct_tpf_baseline.cpp:52–65`

### Description
Paper A Eq. (29) states that the acceleration readout is obtained by projecting C_μν (the principal-part field equation tensor) onto the normalized Ξ direction. The author has identified this equation as incorrect — it was inserted into the paper through AI deception at an earlier time. The author is aware and will revise the paper.

**Why it is wrong (verified by reading function body):**
- C_μν is the field equation (Paper A Eq. 15). It characterizes the configuration gradient structure.
- Θ_{ij} ~ m/r³ (Hessian of the 1/r potential), so C_{ij} ~ κ*(m/r³)² = κm²/r⁶
- Reading `direct_tpf_baseline.cpp:63–64`: `readout.ax = -(c_xx * u_x + c_xy * u_y)` — C used directly as acceleration
- Therefore force ~ κm²/r⁶ — **r⁻⁶ falloff, not r⁻² gravity**

This is confirmed by the code: `direct_tpf` faithfully implements Eq. (29) and therefore
also produces the wrong r⁻⁶ falloff. The equation itself is the error.

**The correct readout (Paper A Eq. 19, Paper A Eq. 2):**
The Ξ field for a point source is Ξ_i = m·d_i/R³. The acceleration reads directly:
`a = −K * Ξ = −K * m * d/R³`, giving r⁻² Newtonian force.
This is implemented correctly in `xi_kernel_deformed` at `xi_kernel_coupling = 0`
(verified by reading `compute_xi_kernel_deformed_accelerations` function body).

### Recommendation
1. **Revise Paper A Eq. (29) before submission.** Replace C_μν projection readout with
   `a_i = −K * Ξ_i`, consistent with the paper's own Eq. (2) and Eq. (19).
2. **Quarantine `direct_tpf_baseline.cpp` implementation.** Rename the path to
   `delta_C_diagnostic` or `delta_C_residual` to make clear it computes a field equation
   residual value, not a dynamics acceleration.
3. The `xi_kernel_deformed` path at coupling=0 correctly implements the actual paper-grounded
   acceleration readout and should be identified as such in Paper A's "Simulation" section.

---

## Summary Tables

---

### Table 1: Paper A Correspondence

| Eq./Concept | Code Location | Status | Notes |
|---|---|---|---|
| λ = 1/4 (Eq. 9) | `source_ansatz.hpp:58` | MATCH | Exact; verified |
| X^μ(x) = x^μ + Ξ^μ(x) (Eq. 2) | `tpf_core_package.cpp:717–718` | MATCH | xi_kernel_deformed at coupling=0 reads directly from Ξ — verified by reading function body |
| Θ_{μν} = ∇_μΞ_ν (Eq. 1/3) | `source_ansatz.cpp:23–35` | MATCH | 3D Hessian on z=0, provisional ansatz |
| Θ_{xx–zz} components (Eqs. 20–25) | `source_ansatz.cpp:23–35` | MATCH | Including Θ_zz = m/R³ |
| I = Θ_{μν}Θ^{μν} − λΘ² (Eq. 5) | `source_ansatz.cpp:74–79` | MATCH (static approx, diagnostic only) | Spatial sum; consistent with Eq. (34); not used in correct dynamics path |
| C_{μν}^(principal) (Eqs. 26–28) | `direct_tpf_baseline.cpp:29–50` | MATCH as field equation | ΔC = 0 deliberately; C_μν is correct as field eq, WRONG as acceleration source |
| Acceleration readout (Eq. 29) | `direct_tpf_baseline.cpp:52–65` | DIVERGENCE — paper error | Eq. (29) is incorrect (AI-injected); C_μν gives r⁻⁶ falloff; author will revise |
| `a = −K * Ξ` (correct readout per Eqs. 2, 19) | `tpf_core_package.cpp:757–758` | MATCH — xi_kernel_deformed at coupling=0 | Verified by reading function body; Newtonian r⁻² force |
| Poisson correspondence κ = 16πG (Eq. 38) | `v11_weak_field_correspondence.cpp:26` | MATCH (in v11 path) | κ unused in correct xi_kernel_deformed path |
| Q_wf = ρ (Eq. 39) | `v11_weak_field_correspondence.cpp` | MATCH (correspondence only) | |
| ΔC_{μν} symbolic (Eq. 14–15) | `direct_tpf_baseline.cpp:20–27` | MATCH (zeroed) | Paper authorizes omission in weak-field |
| Static, weak-field, non-relativistic scope | xi_kernel_deformed at coupling=0 | MATCH | coupling=0 path is static/Newtonian; coupling>0 is VDSG research extension |
| Planar ansatz on z=0 (Sec. IV.A) | `source_ansatz.hpp:6–23` | MATCH | |
| ∇^μC_{μν} = 0 (Eq. 16) | Not monitored | MISSING from code | |

---

### Table 2: Paper B Content in Code but Cut from Paper A

| Paper B Content | Code Location | Status |
|---|---|---|
| v11 Earth-Moon weak-field benchmark (Eqs. 42–47) | `v11_weak_field_correspondence.cpp:265–403` | Present in code, absent from Paper A |
| v11 reference label (v11_weak_field_truncation) | `tpf_core_package.cpp:8,587` | Code explicitly names Paper B version |
| Orbit coherence condition Θ_{tt} + Θ_{rr} = 0 (Paper B Case 1) | Not in code | In Paper B Sec. XI.E; heuristic only |
| Noether charge / conserved energy (Paper B Sec. V, App. E) | Not implemented | In Paper B; explicitly cut from Paper A |
| Energy-Theta ledger E = c² ∫ Θ_{ij}Θ^{ij} d³x (App. E, Eq. 10) | Not implemented | Infrastructure exists (Frobenius norm in source_ansatz.cpp:75) |
| Schwarzschild exterior consistency (Paper B Secs. VIII) | Not in code | Paper B only; not in simulator |
| Bounce metric / TPF horizon (Paper B Secs. IX–X) | Not in code | Paper B only; strong-field territory |
| Emergence of G from sub-atomic structure (Paper B Sec. XII) | Not in code | Paper B Sec. XII; very speculative |
| Nonlinearity index N (Paper B Eq. 12) | Not in code | Paper B only |

---

### Table 3: Code-Only Territory (Beyond Both Papers)

| Parameter / Feature | Config Key | Default | Active Value | Classification | Notes |
|---|---|---|---|---|---|
| xi_kernel_coupling = 0 | tpf_4d_xi_kernel_coupling | 0.0 | 0.0 (baseline) | paper-grounded | coupling=0 is Ξ readout = Newtonian, verified |
| xi_kernel_coupling > 0 (VDSG) | tpf_4d_xi_kernel_coupling | 0.0 | 5000 (active cfg) | research extension | VDSG deformation; no paper derivation yet |
| Xi readout scale | tpf_4d_xi_motion_readout_scale | 1e-12 | 6.674e-11 (≈G) | calibration | Correct: K = G for Newtonian correspondence |
| Kernel factor mode: gamma_minus_one | tpf_4d_xi_kernel_factor_mode | beta_power | gamma_minus_one | exploratory | SR γ−1 factor; not in TPF papers; only active when coupling>0 |
| Kernel factor mode: beta_power | tpf_4d_xi_kernel_factor_mode | beta_power | — | exploratory | β^p factor; not in TPF papers |
| beta_power exponent | tpf_4d_xi_kernel_beta_power | 1.0 | 1.0 | fit parameter | |
| metric_min clamp | tpf_4d_xi_kernel_metric_min | 0.1 | 0.5 | regularization | No effect at coupling=0 (verified) |
| metric_max clamp | tpf_4d_xi_kernel_metric_max | 10.0 | 2.0 | regularization | Numerical stability for VDSG only |
| xi_temporal_mode | tpf_4d_xi_temporal_mode | off | off | dormant | Not paper-grounded |
| xi_temporal_coupling | tpf_4d_xi_temporal_coupling | 0.0 | 0.0 | dormant | |
| xi_kernel_mode modes (scalar_beta, metric_radial, etc.) | tpf_4d_xi_kernel_mode | off | metric_transverse_wake | exploratory | No effect at coupling=0 (verified); VDSG extension |
| VDSG coupling | tpf_vdsg_coupling | 1e-20 | 0.0 (inactive) | exploratory | Velocity-dependent; no paper backing |
| VDSG mode (radial_doppler_rational etc.) | tpf_vdsg_mode | legacy_speed | — | exploratory | Doppler factors not in TPF |
| direct_tpf mode | tpf_dynamics_mode | xi_kernel_deformed | (not active) | QUARANTINE | r⁻⁶ falloff — misinterpretation; see Findings 13, 21 |
| kappa in active config | tpf_kappa | 1e32 | 1 | informational | Unused by xi_kernel_deformed; only relevant if direct_tpf retained as diagnostic |
| Global accel shunt | tpf_global_accel_shunt_enable | false | false | regularization | Numerical cap |
| Cooling fraction | tpf_cooling_fraction | 0.2 | 0 (inactive) | regularization | Velocity damping; no physics basis |
| Auto-softening caps | auto_softening_inner_cap | — | 0.05 | regularization | Structural stability only |
| provisional_readout closures | tpfcore_readout_mode | tensor_radial_projection | — | legacy | Not paper-grounded |

---

### Table 4: Severity-Sorted Action Items

| Priority | Severity | Finding | Short Description | Action |
|---|---|---|---|---|
| 1 | blocking | Finding 21 | Paper A Eq. (29) is incorrect (AI-injected error) | Revise Paper A Eq. (29) before submission; replace C_μν projection with `a = −K * Ξ` |
| 2 | blocking | Finding 5 | Paper A Eq. (29) used as code readout in direct_tpf | Quarantine/rename direct_tpf to delta_C_diagnostic; remove from acceleration paths |
| 3 | blocking | Finding 13 | direct_tpf uses C_μν as acceleration → r⁻⁶ falloff | Same as above; add runtime/compile error if direct_tpf selected for dynamics |
| 4 | concerning | Finding 8 | gamma_minus_one factor not derived from TPF | Document as exploratory VDSG extension; note only active when coupling>0 |
| 5 | concerning | Finding 9 | xi_kernel_coupling = 5000 has no derivation | Derive from first principles or label as fit parameter for VDSG research runs |
| 6 | concerning | Finding 11 | VDSG velocity-dependent forces beyond Paper A scope | Classify as research extension in all outputs; ensure off for pure-TPF paper-facing runs |
| 7 | minor | Finding 3 | I computed as static spatial sum (not used in correct path) | Document clearly; note function not called by xi_kernel_deformed |
| 8 | minor | Finding 4 | C_zz computed but unused in 2D readout | Remove from acceleration path or document as dead code |
| 9 | minor | Finding 6 | Θ_zz in I causes ~17% systematic in direct_tpf diagnostic | Document systematic; applies only to C_μν diagnostic, not dynamics |
| 10 | minor | Finding 10 | metric_clamp has no theoretical justification | Document as regularization for VDSG; verify no effect at coupling=0 (confirmed) |
| 11 | informational | Finding 14 | κ = 1 in active config | Irrelevant to xi_kernel_deformed; only relevant if direct_tpf retained as diagnostic |
| 12 | informational | Finding 16 | Energy ledger (Paper B App. E) not implemented | Implement minimal E_ledger = c² Σ Θ_{ij}² as optional diagnostic |
| 13 | informational | Finding 20 | Noether charge not monitored | Implement as simulation sanity check |
| 14 | informational | Finding 1 | λ = 1/4 correct | No action |
| 15 | informational | Finding 2 | Theta ansatz correct per Paper A | No action |
| 16 | informational | Finding 7 | xi_kernel_deformed at coupling=0 is paper-grounded Ξ readout | Identify explicitly as such in Paper A's simulation section |
| 17 | informational | Finding 12 | v11_weak_field_truncation correctly implements Paper B | No action |
| 18 | informational | Finding 15 | Temporal Xi mode dormant | Keep as exploratory with clear label |
| 19 | informational | Finding 17 | 4D static residual benchmark correct but not wired to dynamics | No action; clarify scope |
| 20 | informational | Finding 18 | Softening is numerical regularization | No action; label correctly |
| 21 | informational | Finding 19 | Bounce metric / TPF horizon (Paper B) not in code | Note for future work |

---

*End of audit. No source files were modified. This is a second-pass revision; see Revision Notes section.*
