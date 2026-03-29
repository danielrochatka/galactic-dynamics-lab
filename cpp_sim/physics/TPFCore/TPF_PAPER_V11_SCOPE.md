# TPF manuscript (v11) vs this codebase

**Reference:** *Gravitational Energy Reframed: A Transformational Physics Framework Approach* (author’s v11 PDF, dated March 6, 2026 in the document). This file is the **scope contract** for what is paper-backed, what is an explicit placeholder in the paper, and what is **simulator-only experiment**.

Re-verify after any edit to `provisional_readout.cpp`, `derived_tpf_radial.*`, or `tpf_core_package.cpp`.

---

## Tier 1 — Paper core (definitions and static / quasi-static claims)

| Manuscript object | Equation / text | Role in code |
|-------------------|-----------------|--------------|
| Configuration map | X^μ = x^μ + Ξ^μ | Not evolved as an independent field; **Ξ** enters via **Φ = −M/R** ansatz → **∇Ξ-like** structure in `source_ansatz` (provisional). |
| **Θ_μν = ∇_μ Ξ_ν** | (2) | **Theta3D** = 3D Hessian of **Φ = −M/R** on **z = 0** (`provisional_point_source_field`, `evaluate_derived_theta` for SI ledger). |
| **I = Θ_μν Θ^μν − λ Θ²**, **λ = 1/4** | (3)–(4) | **`compute_invariant_I`** in `source_ansatz.cpp` (uses **LAMBDA_4D**). **κ–ledger** density in `build_tpf_gravity_profile` uses **`std::abs(compute_invariant_I(θ_tot))`** so negative **I** does not flip **ρ** sign. |
| Configuration EOM | **∇_μ(Θ^μν − λ δ^μν Θ) = 0** (9) | **Single-source residual** in `provisional_point_source_residual` (inspection). **Not** enforced as a dynamical update for **Ξ** during N-body. |
| **C_μν** / ledger | (5)–(6), (10) | **Not** implemented as a coupled metric solver. Simulator uses **readout closures** and/or **effective radial force**, not full **C_μν = κ(...)** time evolution. |
| Weak-field Poisson calibration | Paper Sec. II | **Diagnostic**: `run_weak_field_calibration` compares TPF readout vs **NewtonianPackage** on an axis. |
| Bounce / regular core | Sec. IX (paper) | **`get_tpf_mass_at_r`**: R⁶/(R⁶+R_s⁶) on **enclosed baryonic mass**; separate from **κ** ledger. |

---

## Tier 2 — Paper placeholder (**∆C_μν**)

The paper writes **C_μν** with an explicit **∆C_μν** term: contributions from varying the connection in **Θ_μν = ∇_μ Ξ_ν**, kept **symbolic** and deferred to **future work** (not expanded for weak-field applications).

**There is no single function named ∆C in the code.** Stand-ins and related experiment hooks:

| Idea | Code location | Notes |
|------|---------------|--------|
| Higher-order / non-Newtonian readout | `apply_tensor_radial_closure`, `apply_experimental_radial_r_scaling_closure` | Map **Θ → a** without the hybrid **M_eff** shell model. |
| **κ · I → ρ → M_eff** shell Poisson | `build_tpf_gravity_profile`, `radial_acceleration_scalar_derived` | **Closure**, not derived from expanded **∆C_μν**. Uses **manuscript I** (via **`compute_invariant_I`**) for **ρ_raw**. |
| Velocity-dependent rescaling | `accumulate_velocity_deformed_centripetal_gravity` in `tpf_core_package.cpp` | When **`tpf_gdd_coupling ≠ 0`**, **tensor readout is skipped**; SI Newtonian-style **G M / r²** with **doppler_scale** and **|a| shunt**. |
| Numerical stabilization | ρ shunts, `apply_global_accel_magnitude_shunt` | Engineering, not theory. |

---

## Tier 3 — Code experimental (must be labeled in papers)

- **Provisional readout** (`provisional_readout.cpp`): all **`apply_*_closure`** paths.
- **Radial cooling** (`simulation.cpp`): non-Hamiltonian **v_radial** damping for part of the run.
- **Binning**, **`tpf_kappa`**, **`tpf_poisson_*`**, softening: discretization / tuning.
- **Regime diagnostics** (`regime_diagnostics.hpp`): **reporting only**.

---

## Readout modes (current routing)

| Config string | Entry | **Particle acceleration** |
|---------------|--------|---------------------------|
| `derived_tpf_radial_readout` | `is_derived_tpf_radial_readout_mode` → `apply_derived_tpf_radial_readout_closure` | **Radial only**: **a = a_s r̂**, **a_s** from **`radial_acceleration_scalar_derived`** (bounced baryons + **M_eff**). |
| `tr_coherence_readout` | **Same** as above | **Same closure** as `derived_tpf_radial_readout`. **`theta_tt` / `theta_tr` / `provisional_tangential_readout`** are **diagnostics only**; they are **not** added to **a_x, a_y**. |
| `tensor_radial_projection` (+ `_negated`) | `apply_tensor_radial_closure` | **Θ·r̂**-style superposed contribution. |
| `experimental_radial_r_scaling` | `apply_experimental_radial_r_scaling_closure` | Radial closure from **−θ_rr × r** scaling. |

When **`tpf_gdd_coupling ≠ 0`**, **`TPFCorePackage::compute_accelerations`** uses **only** the GDD centripetal path (stderr notice once).

---

## Equation-to-function map (quick)

| Paper | Function / file |
|-------|-----------------|
| **I** (scalar invariant) | `compute_invariant_I` — `source_ansatz.cpp` |
| **Θ** from −M/R (code units) | `provisional_point_source_theta`, `evaluate_provisional_field_*` — `source_ansatz` / `field_evaluation` |
| **Θ** from −GM/r (SI, ledger) | `evaluate_derived_theta`, `sum_derived_theta_at_point` — `derived_tpf_radial.cpp` |
| **ρ ∝ κ \|I\|** (ledger) | `build_tpf_gravity_profile` — `derived_tpf_radial.cpp` |
| **a_r = −G (M_bary + M_eff) / r²** | `radial_acceleration_scalar_derived` — `derived_tpf_radial.cpp` |
| Residual (weak-field) | `provisional_point_source_residual` — `source_ansatz.cpp` |
| Integrator | `velocity_verlet` — `integrator.cpp` (unchanged) |

**Legacy diagnostic:** `derived_invariant_I_contracted` = **Θ_xx² + …** only (no **−λ Θ²**). **Do not** treat it as manuscript **I**; κ **ledger** uses **`compute_invariant_I`**.

---

## Units

- **NewtonianPackage**: **G ≡ 1** (simulation mass–length–time convention).
- **TPF hybrid / GDD / derived θ**: **SI** via **`TPF_G_SI`** and **`c`** from `config.hpp`.

---

## What you can say in a methods section (safe)

- The simulator evaluates **Θ** and **I** from a **stated provisional potential ansatz**, applies **optional readout closures** to obtain **particle acceleration**, and integrates with **Verlet**.
- **Manuscript Eq. (3)** for **I** is used in the **κ–ledger** shell model (**Tier 2 closure**).
- **∆C_μν** and full **C_μν** dynamics are **not** implemented as written in v11.

## Overclaiming (avoid)

- That **`tr_coherence_readout`** applies **tangential coherence acceleration** (it does **not** in the shared derived closure).
- That **GDD-on** runs are “tensor TPF readout” dynamics (readout is **bypassed**).
- That galaxy rotation results are **claimed v11 results** (paper scopes dynamics / DM as **future work**).
