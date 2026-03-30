# TPFCore physics package

**TPFCore** is a **compiled physics package** inside the **Galactic Dynamics Lab** repository (`cpp_sim/physics/TPFCore/`). It implements the **primitive TPF field structure** from the manuscript (configuration displacement **Ξ**, gradient tensor **Θ**, scalar invariant **I**) on the simulation plane, plus **labeled experimental layers** for motion and diagnostics.

It is **not** the old removed “weak-field Newtonian-like TPF” package. For **lab-wide overview** and **engine** behavior, see **[../../../README.md](../../../README.md)** and **[../../README.md](../../README.md)** (`cpp_sim`). Here: **what this package is**, **what is paper-aligned vs exploratory**, and **how to read branch labels** in manifests.

**Manuscript v11 vs simulator tiers:** **[TPF_PAPER_V11_SCOPE.md](TPF_PAPER_V11_SCOPE.md)**.

---

## Paper-backed core (honest scope)

- **Ξ** — displacement field from the potential ansatz (see `source_ansatz.*`).
- **Θ_{μν} = ∇_μ Ξ_ν** — symmetric tensor; evaluated as **3D Hessian** of a softened Coulomb **Φ** on **z = 0** (field and sources in the plane; `Theta_xz = Theta_yz = 0` on slice).
- **I = Θ_{μν}Θ^{μν} − λ Θ²** with **λ = 1/4** in 4D (**fixed**, not a tunable “theory knob” in the sense of fitting data).
- **Field-equation residual** (inspection modes): weak-field configuration equation in the form used in code; analytic derivatives; **structural sanity check**, not proof of full nonlinear dynamics.

---

## Ansatz vs closure vs diagnostics (distinctions)

| Layer | Role |
|-------|------|
| **Ansatz** | **Φ = −M/R**, **R² = dx²+dy²+eps²**; **Ξ**, **Θ** from closed-form derivatives (`source_ansatz.*`). Provisional where the manuscript leaves the full source unspecified. |
| **Closure (acceleration)** | **Current code:** `TPFCorePackage::compute_accelerations` always builds **TPF readout baseline** **`ax, ay`** from **`compute_provisional_readout_acceleration`** when provisional readout is on. If **`tpf_vdsg_coupling ≠ 0`**, it **adds** the SI velocity modifier from **`accumulate_vdsg_velocity_modifier`** (Newtonian magnitude × (**doppler_scale − 1**) per interaction). Paths are **exploratory** unless stated otherwise. |
| **Diagnostics** | CSVs, debug columns, and **`ReadoutDiagnostics`**: on **derived-radial** readout modes, **theta_tt** / **theta_tr** / **provisional_tangential_readout** are **not** added to **ax, ay** (only radial **a_s** is). **VDSG** contributes an additive SI excess on top of that baseline, not a replacement readout. |

---

## Dynamics branch vs metrics branch (audit language)

The simulator exposes **resolved strings** in **`run_info.txt`** and **`render_manifest.json`** (computed in **`render_audit.cpp`**):

- **`active_dynamics_branch`** — **Readout identity** for dynamics: **`TPF_readout_acceleration:<mode>`** when **`tpfcore_enable_provisional_readout`** is on (independent of **`tpf_vdsg_coupling`**). If provisional readout is off, dynamics are disabled for TPFCore acceleration.
- **`active_metrics_branch`** — Identity of the **configured readout** label (e.g. `tpfcore_readout:derived_tpf_radial_readout`).

**Integrator accelerations** are **baseline readout + optional VDSG modifier** when coupling is nonzero. **`acceleration_code_path`** appends **`+ accumulate_vdsg_velocity_modifier`** when **`tpf_vdsg_coupling ≠ 0`**.

---

## VDSG (Velocity-Deformed Spacetime Gradient)

**VDSG** (**Velocity-Deformed Spacetime Gradient**) is an **exploratory**, **velocity-dependent** **additive** correction: per interaction, Newtonian magnitude **a_N = G M / r²** is scaled by **doppler_scale = 1 + λ_eff (v·r̂)/c**; the code adds only the **excess** **a_N (doppler_scale − 1)** on top of the **TPF readout baseline** (see **`accumulate_vdsg_velocity_modifier`** in `tpf_core_package.cpp`). **`active_dynamics_branch`** stays **`TPF_readout_acceleration:<mode>`**; manifests record the modifier via **`acceleration_code_path`**.

**Legacy alias (once):** the parser accepts **`tpf_gdd_coupling`** as an alias for **`tpf_vdsg_coupling`** (historical name). **Canonical key:** **`tpf_vdsg_coupling`**. Manifests note the rename for audit.

---

## Readout modes (summary)

Details and column semantics: **`provisional_readout.cpp`**, **`TPF_PAPER_V11_SCOPE.md`**.

- **`tensor_radial_projection`** (and negated variant) — Exploratory **Θ·r̂**-style superposed acceleration from **`apply_tensor_radial_closure`**.
- **`derived_tpf_radial_readout`**, **`tr_coherence_readout`** — **Current code:** both match **`is_derived_tpf_radial_readout_mode`** (`derived_tpf_radial.hpp`) and call **`apply_derived_tpf_radial_readout_closure`**. **Particle accelerations** are **purely radial**: `ax = a_s (x/r)`, `ay = a_s (y/r)` with **`a_s`** from **`radial_acceleration_scalar_derived`**. **theta_tt**, **theta_tr**, and **provisional_tangential_readout** are computed **only** into **`ReadoutDiagnostics`** (and related diagnostics); they are **not** added to **ax, ay** on this path.
- **`experimental_radial_r_scaling`** — Separate closure (**`apply_experimental_radial_r_scaling_closure`**); see scope doc.

When **`tpf_vdsg_coupling ≠ 0`**, **`compute_accelerations`** adds the VDSG velocity modifier **after** the readout baseline; **`apply_global_accel_magnitude_shunt`** applies to the **sum** when the modifier is active.

---

## Config keys (`defaults.cfg`)

Package defaults live in **`defaults.cfg`** in this directory. Important keys (non-exhaustive; see file and `config.hpp`):

- **`tpfcore_enable_provisional_readout`** — Required for dynamical modes with TPFCore.
- **`tpfcore_readout_mode`**, **`tpfcore_readout_scale`**, **`tpfcore_theta_tt_scale`**, **`tpfcore_theta_tr_scale`**
- **`tpf_kappa`**, **`tpf_poisson_bins`**, **`tpf_poisson_max_radius`**, **`tpf_cooling_fraction`**
- **`tpf_vdsg_coupling`**, **`tpf_vdsg_mass_baseline_kg`**
- Inspection: **`tpfcore_probe_radius_*`**, **`tpfcore_dump_*`**, **`tpfcore_source_softening`**

**Simulation-wide** keys (galaxy ICs, softening, `n_stars`, etc.) belong to the **application** config, not this file alone — see **[../../README.md](../../README.md)**.

---

## CSV outputs (inspection)

**`theta_profile.csv`**, **`invariant_profile.csv`**, **`field_summary.txt`** — produced in **inspection** modes (`tpf_single_source_inspect`, …). Column tables and symmetry expectations were in previous versions of this README; the **authoritative** column definitions are in code comments and CSV headers. Symmetry: e.g. **residual_y ≈ 0** on the +x axis for single-source at origin.

**`tpf_readout_debug.csv`** — Dynamical runs when **`tpfcore_dump_readout_debug`**: mode-dependent columns for diagnosing radial vs tangential acceleration components.

---

## What TPFCore claims vs does not claim

**Claims:**

- Implements the **primitive** Ξ / Θ / **I** structure and **λ = 1/4** as stated.
- Exposes **explicit** exploratory closures for motion and **VDSG** with **stderr / run_info / manifest** branch labeling.

**Does not claim:**

- Full **C_μν** / metric evolution or **∆C_μν** expansion as in the manuscript’s deferred sections.
- That **galaxy** runs with cooling, readout, or VDSG constitute **final** TPF predictions — they are **experiments** with **documented** stand-ins and **audit** metadata.

**Honesty:** Higher-order / galaxy-scale behavior that depends on **provisional** readout, **artificial cooling**, or **heuristic** VDSG normalization must be interpreted as **experimental**, not as proof of the closed theory. Use **`render_manifest.json`**, **`run_info`**, and this README’s branch distinction to avoid conflating **measurement** with **dynamics**.
