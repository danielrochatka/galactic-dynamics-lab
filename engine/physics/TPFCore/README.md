# TPFCore physics package

**TPFCore** is a **compiled physics package** inside the **Galactic Dynamics Lab** repository (`engine/physics/TPFCore/`). It contains the runtime TPF dynamics routes used by the simulator plus staged field/residual benchmark paths used for isolated static diagnostics.

It is **not** the old removed “weak-field Newtonian-like TPF” package. For **lab-wide overview** and **engine** behavior, see **[../../../README.md](../../../README.md)** and **[../../README.md](../../README.md)** (`engine`). Here: **what this package is**, **what is paper-aligned vs exploratory**, and **how to read branch labels** in manifests.

Current spike branch status: runtime dynamics still use the existing projected-vector / spatial-tensor runtime path (`legacy_readout`, `direct_tpf`, `v11_weak_field_truncation`). Stage 1–4 added an isolated `tpf_4d_static_residual_benchmark` path with static Xi4/ordered-Theta4 field objects and full spatial-support residual diagnostics over an x/y/z stencil. Stage 5 adds view-plane export/plot outputs rendered from that full spatial-support static 4D residual evaluation. Stage 7A adds `tpf_4d_static_motion_readout_benchmark`, a candidate probe-motion readout benchmark from frozen static 4D field data using `GravityStaticMotionReadout_v1` (principal tensor readout with DeltaC omitted for this named benchmark).

**Manuscript v11 vs simulator tiers:** **[TPF_PAPER_V11_SCOPE.md](TPF_PAPER_V11_SCOPE.md)**.

**Audit evidence records (contamination/sanity):** **[audits/](audits/README.md)**.

---

## Runtime dynamics vs static 4D benchmark scope

- **Ξ** — displacement field from the potential ansatz (see `source_ansatz.*`).
- **Runtime dynamics path (current)**: uses existing projected-vector / spatial-tensor structures and route-specific acceleration closures.
- **Static 4D benchmark path (Stage 1–4, spike branch)**: diagnostic-only path that evaluates static Xi4 / ordered Theta4 and a full spatial-support residual diagnostic (x/y/z stencil), with no particle integration and no acceleration path.
- **Static 4D motion-readout benchmark (Stage 7A, spike branch)**: simulator-exploratory path that evaluates static Xi4 / ordered Theta4 field data and computes candidate probe-motion readout accelerations via `GravityStaticMotionReadout_v1`; this is a candidate future paper-expansion path if benchmark behavior continues to pass regression tests.
- **I = Θ_{μν}Θ^{μν} − λ Θ²** with **λ = 1/4** in 4D (**fixed**, not a tunable “theory knob” in the sense of fitting data).
- **Static residual benchmark outputs**: `tpf_4d_static_residual_summary.txt`, `tpf_4d_static_residual_slice.csv` (legacy alias for central xy view-plane), `tpf_4d_static_residual_slice_xy.csv`, `tpf_4d_static_residual_slice_xz.csv`, `tpf_4d_static_residual_slice_yz.csv`, `tpf_4d_static_residual_sources.csv`, `tpf_4d_static_residual_bins_nearest_source.csv`, and `tpf_4d_static_residual_bins_origin.csv`; slice CSVs are central view-plane diagnostic renderings from full spatial-support static 4D residual evaluation, and binned CSVs provide quantitative static-field residual evidence for cross-run regression comparison.
- **Optional plot outputs**: when `--plot` is supplied and plotting dependencies are available, `plot_tpf_4d_static_residual.py` emits PNG view-plane diagnostic renderings from full spatial-support static 4D residual evaluation.
- **Exercised and unexercised scope**: plots visualize the static 4D residual benchmark. Dynamics, moving sources, source velocities, time evolution, source worldlines, physical coupling, orbital behavior, and DeltaC closure are outside this benchmark’s exercised scope.

---

## Ansatz vs closure vs diagnostics (distinctions)

| Layer | Role |
|-------|------|
| **Ansatz** | **Φ = −M/R**, **R² = dx²+dy²+eps²**; **Ξ**, **Θ** from closed-form derivatives (`source_ansatz.*`). Provisional where the manuscript leaves the full source unspecified. |
| **Closure (acceleration)** | **Current code is route-dependent:** **`direct_tpf`** is the tensor principal-part route (field_evaluation → legacy spatial tensor objects → principal_Cij → tensor_projection; Theta/I/kappa baseline; DeltaC omitted in current implementation scope) with optional additive VDSG extension. **`v11_weak_field_truncation`** is the explicit weak-field correspondence helper (alpha_si path; legacy/benchmark compatibility). **`legacy_readout`** uses readout baseline from **`compute_provisional_readout_acceleration`**, then **`accumulate_vdsg_velocity_modifier`** (no-op when λ = 0), then optional **`apply_global_accel_magnitude_shunt`**. Modifier uses **doppler_scale = 1 + λ_eff |v_rel|/c** per interaction. |
| **Diagnostics** | CSVs, debug columns, and **`ReadoutDiagnostics`**: on **derived-radial** readout modes, **theta_tt** / **theta_tr** / **provisional_tangential_readout** are **not** added to **ax, ay** (only radial **a_s** is). **VDSG** contributes an additive SI excess on top of that baseline, not a replacement readout. |

---

## Dynamics branch vs metrics branch (audit language)

The simulator exposes **resolved strings** in **`run_info.txt`** and **`render_manifest.json`** (computed in **`render_audit.cpp`**):

- **`active_dynamics_branch`** — runtime branch identity (`direct_tpf` tensor principal-part route vs `v11_weak_field_truncation` correspondence helper vs `legacy_readout` provisional path).
- **`active_metrics_branch`** — matching metrics branch identity for that runtime route.

**Integrator accelerations** depend on routing: **`direct_tpf`** uses the tensor principal-part Theta/I/kappa baseline (DeltaC omitted in current implementation scope; VDSG optional additive extension; readout/shunt/cooling rejected), **`v11_weak_field_truncation`** is the weak-field correspondence helper (alpha_si path, legacy/benchmark compatibility), while **`legacy_readout`** uses baseline readout + optional VDSG + optional global shunt.

---

## VDSG (Velocity-Deformed Spacetime Gradient)

**VDSG** (**Velocity-Deformed Spacetime Gradient**) is an **exploratory**, **velocity-dependent** **additive** correction: per interaction, **doppler_scale = 1 + λ_eff |v_rel| / c** (relative speed of the interaction, not **v·r̂** — so tangential / circular motion still couples); excess **a_N (doppler_scale − 1)** is added along the Newtonian line on top of the **TPF readout baseline** (see **`accumulate_vdsg_velocity_modifier`** in `tpf_core_package.cpp`). **`apply_global_accel_magnitude_shunt`** runs **after** every TPFCore acceleration evaluation (**same for λ = 0 and λ ≠ 0**). **`active_dynamics_branch`** stays **`TPF_readout_acceleration:<mode>`**; **`acceleration_code_path`** lists the full pipeline.

**Legacy alias (once):** the parser accepts **`tpf_gdd_coupling`** as an alias for **`tpf_vdsg_coupling`** (historical name). **Canonical key:** **`tpf_vdsg_coupling`**. Manifests note the rename for audit.

---

## Readout modes (summary)

Details and column semantics: **`provisional_readout.cpp`**, **`TPF_PAPER_V11_SCOPE.md`**.

- **`tensor_radial_projection`** (and negated variant) — Exploratory **Θ·r̂**-style superposed acceleration from **`apply_tensor_radial_closure`**.
- **`derived_tpf_radial_readout`**, **`tr_coherence_readout`** — **Current code:** both match **`is_derived_tpf_radial_readout_mode`** (`derived_tpf_radial.hpp`) and call **`apply_derived_tpf_radial_readout_closure`**. **Particle accelerations** are **purely radial**: `ax = a_s (x/r)`, `ay = a_s (y/r)` with **`a_s`** from **`radial_acceleration_scalar_derived`**. **theta_tt**, **theta_tr**, and **provisional_tangential_readout** are computed **only** into **`ReadoutDiagnostics`** (and related diagnostics); they are **not** added to **ax, ay** on this path.
- **`experimental_radial_r_scaling`** — Separate closure (**`apply_experimental_radial_r_scaling_closure`**); see scope doc.

In **`legacy_readout`**, **`compute_accelerations`** adds the VDSG modifier vector (zeros if λ = 0) and can apply optional global **`|a|` shunt** (off by default). In canonical **`direct_tpf`**, those exploratory/stabilizer knobs are rejected.

---

## Config keys (`defaults.cfg`)

Package defaults live in **`defaults.cfg`** in this directory. Important keys (non-exhaustive; see file and `config.hpp`):

- **`tpfcore_enable_provisional_readout`** — Required for dynamical modes with TPFCore.
- **`tpfcore_readout_mode`**, **`tpfcore_readout_scale`**, **`tpfcore_theta_tt_scale`**, **`tpfcore_theta_tr_scale`**
- **`tpf_kappa`**, **`tpf_poisson_bins`**, **`tpf_poisson_max_radius`**, **`tpf_cooling_fraction`**
- **`tpf_vdsg_coupling`**, **`tpf_vdsg_mass_baseline_kg`**
- **`tpf_global_accel_shunt_enable`**, **`tpf_global_accel_shunt_fraction`**, **`tpf_accel_pipeline_diagnostics_csv`**
- Inspection: **`tpfcore_probe_radius_*`**, **`tpfcore_dump_*`**, **`tpfcore_source_softening`**

**Simulation-wide** keys (galaxy ICs, validation scenario ICs, softening, `n_stars`, etc.) belong to the **application** config + scenario resolver (`engine/scenario_defaults.cpp`, `engine/resolved_scenario.cpp`), not this file alone — see **[../../README.md](../../README.md)**.

---

## CSV outputs (inspection)

**`theta_profile.csv`**, **`invariant_profile.csv`**, **`field_summary.txt`** — produced in **inspection** modes (`tpf_single_source_inspect`, …). Column tables and symmetry expectations were in previous versions of this README; the **authoritative** column definitions are in code comments and CSV headers. Symmetry: e.g. **residual_y ≈ 0** on the +x axis for single-source at origin.

**`tpf_readout_debug.csv`** — Dynamical runs when **`tpfcore_dump_readout_debug`**: mode-dependent columns for diagnosing radial vs tangential acceleration components.

**`tpf_4d_static_residual_summary.txt`**, **`tpf_4d_static_residual_slice.csv`**, **`tpf_4d_static_residual_slice_xy.csv`**, **`tpf_4d_static_residual_slice_xz.csv`**, **`tpf_4d_static_residual_slice_yz.csv`**, **`tpf_4d_static_residual_sources.csv`**, **`tpf_4d_static_residual_bins_nearest_source.csv`**, **`tpf_4d_static_residual_bins_origin.csv`** — benchmark artifacts from `simulation_mode=tpf_4d_static_residual_benchmark`. The static residual computation evaluates full 4D tensor objects over 3D spatial support; slice CSVs are central view-plane diagnostic renderings from that full spatial-support evaluation, and bin CSVs are derived quantitative summaries of the static 4D residual benchmark.

**`tpf_4d_static_motion_readout_summary.txt`**, **`tpf_4d_static_motion_readout_probe_grid.csv`**, **`tpf_4d_static_motion_readout_bins_origin.csv`** — benchmark artifacts from `simulation_mode=tpf_4d_static_motion_readout_benchmark`. This path evaluates frozen static 4D field objects and computes candidate probe-motion readout accelerations from the principal spatial tensor construction (`GravityStaticMotionReadout_v1`).

**`tpf_4d_xi_motion_probe_summary.txt`**, **`tpf_4d_xi_motion_probe_trajectories.csv`**, **`tpf_4d_xi_motion_probe_initial_readout.csv`** — benchmark artifacts from `simulation_mode=tpf_4d_xi_motion_probe_benchmark`. This path advances dynamic probes with `GravityXiMotionReadout_v1` (`a=-K_xi*Xi_spatial`) using fixed-source Stage 7B field evaluation and writes trajectory readout samples for each probe over time.

### Stage 8A Xi-kernel deformation (benchmark-only)

`GravityXiKernelDeformation_v1` adds an isolated, configurable kernel deformation stage inside `tpf_4d_xi_motion_probe_benchmark` before Xi acceleration readout. The benchmark now computes per-source Xi contributions, applies a configured deformation mode (`off`, `scalar_beta`, `metric_radial`, `metric_velocity`, `spacetime_metric`), sums `Xi_eff`, then reads acceleration strictly as `a=-K_xi*Xi_eff_spatial`.

- `off` preserves Stage 7B behavior exactly (same Xi kernel and same readout equation).
- `spacetime_metric` can emit an `Xi_t`/Xi0 diagnostic (`tpf_4d_xi_temporal_mode=norm_scaled`) but does not feed Xi0 into acceleration in Stage 8A.
- Active velocity-dependent kernel deformation (`tpf_4d_xi_kernel_mode!=off` with nonzero `tpf_4d_xi_kernel_coupling`) requires `tpf_4d_xi_motion_integrator=semi_implicit_euler`; `velocity_verlet` remains valid for off/identity (Stage 7B-equivalent) behavior.
- This is distinct from older additive acceleration VDSG paths; Stage 8A deforms Xi kernel evaluation before readout and does not append acceleration terms.

**Optional trajectory visualization PNGs** — `plot_tpf_4d_xi_motion_probe.py` reads the already-computed Xi-direct trajectory CSV and can emit:
`tpf_4d_xi_motion_probe_xy_trajectories.png`,
`tpf_4d_xi_motion_probe_xz_trajectories.png`,
`tpf_4d_xi_motion_probe_yz_trajectories.png`,
`tpf_4d_xi_motion_probe_radius_vs_time.png`,
`tpf_4d_xi_motion_probe_acceleration_norm_vs_time.png`.

---

## What TPFCore claims vs does not claim

**Claims:**

- Implements the **primitive** Ξ / Θ / **I** structure and **λ = 1/4** as stated.
- Exposes **explicit** exploratory closures for motion and **VDSG** with **stderr / run_info / manifest** branch labeling.

**Does not claim:**

- Full **C_μν** / metric evolution or **∆C_μν** expansion as in the manuscript’s deferred sections.
- That **galaxy** runs with cooling, readout, or VDSG constitute **final** TPF predictions — they are **experiments** with **documented** stand-ins and **audit** metadata.

**Honesty:** Higher-order / galaxy-scale behavior that depends on **provisional** readout, **artificial cooling**, or **heuristic** VDSG normalization must be interpreted as **experimental**, not as proof of the closed theory. Use **`render_manifest.json`**, **`run_info`**, and this README’s branch distinction to avoid conflating **measurement** with **dynamics**.
