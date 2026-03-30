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
| **Ansatz** | **Φ = −M/R**, **R² = dx²+dy²+eps²**; **Ξ**, **Θ** from closed-form derivatives (`source_ansatz.*`). Provisional where the paper leaves the full source unspecified. |
| **Closure** | Maps **Θ** (and ledgers) into **acceleration** or into **SI-style VDSG** — **exploratory** unless stated otherwise; isolated in `provisional_readout.*` / `tpf_core_package.*`. |
| **Diagnostics** | CSVs and stderr ledgers that record tensor components, regime labels, etc.; may **not** be the same as the integrator’s **ax, ay** when **VDSG** drives dynamics. |

---

## Dynamics branch vs metrics branch (audit language)

The simulator exposes **resolved strings** in **`run_info.txt`** and **`render_manifest.json`**:

- **`active_dynamics_branch`** — What **actually integrates** particle motion (e.g. **VDSG centripetal** when `tpf_vdsg_coupling ≠ 0`, else readout-based acceleration from `tpfcore_readout_mode` when provisional readout is on).
- **`active_metrics_branch`** — Identity of the **provisional readout / tensor diagnostic** path (e.g. `tpfcore_readout:derived_tpf_radial_readout`), which may still run for profiles and CSVs **even when VDSG supersedes accelerations**.

Do not equate **configured** `tpfcore_readout_mode` with **actual** acceleration routing: when VDSG is active, tensor readout is **not** used for **ax, ay** (see `TPFCorePackage::compute_accelerations`). **`acceleration_code_path`** in the manifest names the C++ symbol-level path for audit.

---

## VDSG (Velocity-Deformed Spacetime Gradient)

**VDSG** stands for **Velocity-Deformed Spacetime Gradient**: the **exploratory** closure path in this package that applies a **velocity-dependent rescaling** (e.g. **`doppler_scale`** tied to **`tpf_vdsg_coupling`**) on top of a **centripetal / SI-style** acceleration route, implemented in **`accumulate_velocity_deformed_centripetal_gravity`**. In manifests and **`run_info`**, labels such as **`VDSG_centripetal_SI`** refer to this dynamics branch.

### Naming and legacy alias

- **Canonical config key:** **`tpf_vdsg_coupling`** (λ in **doppler_scale** for the VDSG SI path; effective λ may be mass-normalized per interaction — see code and run_info).
- **Legacy alias:** **`tpf_gdd_coupling`** is accepted in the config parser and maps to **`tpf_vdsg_coupling`**. Manifests record the alias for audit.

---

## Readout modes (summary)

Details and column semantics: **`provisional_readout.cpp`**, **`TPF_PAPER_V11_SCOPE.md`**.

- **`tensor_radial_projection`** (and negated variant) — Exploratory spatial projection; not the paper’s preferred t–r coherence story.
- **`derived_tpf_radial_readout`**, **`tr_coherence_readout`** — Share the **hybrid radial** closure where `is_derived_tpf_radial_readout_mode` holds; **particle acceleration** is **radial** from **`radial_acceleration_scalar_derived`** and the κ–**I** ledger. **Theta_tt / Theta_tr** (and similar) can appear in **diagnostics**, not necessarily in **a** for that path.
- **`experimental_radial_r_scaling`** — Experimental (see code / scope doc).

When **`tpf_vdsg_coupling ≠ 0`**, **`accumulate_velocity_deformed_centripetal_gravity`** defines dynamics; readout remains available for **measurement / CSV**, not for substituting Newtonian −∇Φ in that regime.

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
