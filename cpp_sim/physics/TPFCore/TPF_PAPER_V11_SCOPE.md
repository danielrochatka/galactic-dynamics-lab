# TPF manuscript (v11) vs this codebase

**Reference:** *Gravitational Energy Reframed: A Transformational Physics Framework Approach* (author’s v11 PDF, dated March 6, 2026 in the document). This file is a **scope contract** for what this repo implements vs what is **explicit** manuscript structure vs **closure / ledger / simulator** layers. It is **not** a substitute for the PDF: **narrative claims about what the manuscript “says” or “defers” must be checked against the v11 PDF** unless they are grounded in this file or in code comments that already cite behavior.

**Ground truth:** **current C++ behavior** (see `TPFCorePackage::compute_accelerations`, `provisional_readout.cpp`, `derived_tpf_radial.*`). Where this document summarizes a manuscript equation, **verify numbering and text in the PDF** if you need a publication citation.

Re-verify after any edit to `provisional_readout.cpp`, `derived_tpf_radial.*`, or `tpf_core_package.cpp`.

---

## Labels used here

| Label | Meaning |
|-------|---------|
| **Tier 1** | Manuscript-aligned **definitions / static** structure as implemented in this repo (Ξ/Θ/I, ansatz, residual inspection). |
| **Tier 2** | **Closure / ledger** layers that map **Θ, I** into effective forces or readout — **not** full **C_μν** dynamics as expanded in the manuscript. |
| **Tier 3** | **Simulator-only** experimental: provisional readouts, cooling, binning, regime reporting. |
| **Engineering** | Numerical stabilization (shunts, caps) — **not** theory. |

---

## Tier 1 — Paper core (definitions and static / quasi-static claims)

| Manuscript object | Equation / text | Role in code |
|-------------------|-----------------|--------------|
| Configuration map | X^μ = x^μ + Ξ^μ | Not evolved as an independent field; **Ξ** enters via **Φ = −M/R** ansatz → **∇Ξ-like** structure in `source_ansatz` (provisional). |
| **Θ_μν = ∇_μ Ξ_ν** | (2) | **Theta3D** = 3D Hessian of **Φ = −M/R** on **z = 0** (`provisional_point_source_field`, `evaluate_derived_theta` for SI ledger). |
| **I = Θ_μν Θ^μν − λ Θ²**, **λ = 1/4** | (3)–(4) | **`compute_invariant_I`** in `source_ansatz.cpp` (uses **LAMBDA_4D**). **κ–ledger** density in `build_tpf_gravity_profile` uses **`std::abs(compute_invariant_I(θ_tot))`** so negative **I** does not flip **ρ** sign. |
| Configuration EOM | **∇_μ(Θ^μν − λ δ^μν Θ) = 0** (9) | **Single-source residual** in `provisional_point_source_residual` (inspection). **Not** enforced as a dynamical update for **Ξ** during N-body. |
| **C_μν** / ledger | (5)–(6), (10) | **Not** implemented as a coupled metric solver. Simulator uses **readout closures** and/or **effective radial force**, not full **C_μν = κ(...)** time evolution. |
| Weak-field Poisson calibration | Sec. II (verify in PDF) | **Diagnostic**: `run_weak_field_calibration` compares TPF readout vs **NewtonianPackage** on an axis. |
| Bounce / regular core | Sec. IX (verify in PDF) | **`get_tpf_mass_at_r`**: R⁶/(R⁶+R_s⁶) on **enclosed baryonic mass**; separate from **κ** ledger. |

---

## Tier 2 — Paper placeholder (**∆C_μν**)

**Repo summary (verify against v11 PDF):** the manuscript introduces **C_μν** with an explicit **∆C_μν** term tied to varying the connection in **Θ_μν = ∇_μ Ξ_ν**; the full expansion is **not** implemented here as weak-field dynamics.

**There is no single function named ∆C in the code.** Stand-ins and related experiment hooks:

| Idea | Code location | Notes |
|------|---------------|--------|
| Higher-order / non-Newtonian readout | `apply_tensor_radial_closure`, `apply_experimental_radial_r_scaling_closure` | Map **Θ → a** without the hybrid **M_eff** shell model. |
| **κ · I → ρ → M_eff** shell Poisson | `build_tpf_gravity_profile`, `radial_acceleration_scalar_derived` | **Closure**, not derived from expanded **∆C_μν**. Uses **manuscript I** (via **`compute_invariant_I`**) for **ρ_raw**. |
| Velocity-dependent rescaling | `accumulate_velocity_deformed_centripetal_gravity` in `tpf_core_package.cpp` | When **`tpf_vdsg_coupling ≠ 0`**, **tensor readout is skipped**; SI Newtonian-style **G M / r²** with **doppler_scale** and **|a| shunt**. |
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

When **`tpf_vdsg_coupling ≠ 0`**, **`TPFCorePackage::compute_accelerations`** uses **only** the VDSG centripetal path (stderr notice once).

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

## Unresolved (requires v11 PDF or author confirmation)

This repo does **not** embed the full manuscript. The following are **not** settled here:

- **Exact** manuscript wording or section references for what is “future work” vs asserted for **galaxy-scale** dynamics, **dark matter**, or **rotation curves** (do **not** cite this file as if it quoted the PDF).
- Equation **numbering** in the PDF vs the labels used in the tables above (cross-check when publishing).

---

## Units

- **NewtonianPackage**: **G ≡ 1** (simulation mass–length–time convention).
- **TPF hybrid / VDSG / derived θ**: **SI** via **`TPF_G_SI`** and **`c`** from `config.hpp`.

---

## What you can say in a methods section (safe, code-grounded)

- The simulator evaluates **Θ** and **I** from a **stated provisional potential ansatz**, applies **optional readout closures** or **VDSG** per **`TPFCorePackage::compute_accelerations`**, and integrates with **Verlet**.
- **Invariant I** from **`compute_invariant_I`** (same combination as manuscript **Eq. (3)** in this repo’s implementation) feeds the **κ–ledger** shell model (**Tier 2 closure**). **Verify equation numbers against the PDF** in a paper.
- **∆C_μν** and full **C_μν** dynamics are **not** implemented as expanded weak-field physics in this codebase.

---

## Overclaiming (avoid)

- That **`tr_coherence_readout`** applies **tangential coherence acceleration** in **ax, ay** (it does **not**: **current code** uses the same **radial-only** derived closure as **`derived_tpf_radial_readout`**; tangential quantities are **diagnostics-only** on that path).
- That **VDSG-on** runs use **tensor readout** for **ax, ay** (readout **does not** supply those accelerations when **`tpf_vdsg_coupling ≠ 0`**).
- That **galaxy** or **rotation-curve** results from this simulator are **automatically** “v11 manuscript results” without checking **both** this code and the **v11 PDF** (see **Unresolved** above).
