# TPFCore Calculation Audit — 2026-04-29

## Executive summary
- `xi_kernel_deformed` runtime routing is isolated to `compute_xi_kernel_deformed_accelerations` and returns before legacy/direct/v11 paths.
- Xi sign convention and readout sign match declared contract: `d = target-source`, `Xi_base = m*d/(|d|^2+eps^2)^(3/2)`, `a=-K_xi*Xi_eff`.
- Wake/continuous transverse gates match documented formulas in code.
- **Definite mismatch found:** derived radial path always includes enclosed stellar mass in `radial_acceleration_scalar_derived`, independent of `enable_star_star_gravity`; this requires metadata truth fix and theory decision on physics behavior.

## Route table
| Route | Runtime entry | Status |
|---|---|---|
| xi_kernel_deformed | `TPFCorePackage::compute_xi_kernel_deformed_accelerations` | Active and isolated |
| direct_tpf | `compute_direct_tpf_accelerations` (+optional VDSG) | Separate branch |
| legacy_readout/provisional | `eval_accel_pipeline` | Separate branch |
| v11_weak_field_truncation | `compute_v11_weak_field_truncation_accelerations` | Separate branch |
| derived_tpf_radial readout | via provisional readout mode closures | Diagnostic/readout branch |
| additive VDSG | direct/legacy branches only | not used in xi branch |

## Formula table (implemented)
- off: `Xi_eff = Xi_base`.
- scalar_beta: scales Xi magnitude by `1 + factor_raw`.
- metric_radial: anisotropic metric along radial axis `n=r_hat`.
- metric_velocity: anisotropic metric along relative velocity axis.
- metric_transverse_wake: `beta_pass = |v_trans|/c`; gate `0.5*(1+tanh((v_rad/v_trans-0.10)/0.05))`.
- metric_transverse_continuous: gate `0.5*(1+tanh(v_rad/v_trans))`.
- spacetime_metric: same kernel metric machinery plus temporal coupling branch.

## Config truth table (condensed)
- Active in xi runtime: `tpf_dynamics_mode`, `tpf_4d_xi_motion_readout_scale`, `tpf_4d_xi_kernel_mode`, `tpf_4d_xi_kernel_coupling`, `tpf_4d_xi_kernel_factor_mode`, `tpf_4d_xi_kernel_beta_power`, `tpf_4d_xi_kernel_metric_min/max`, `tpf_4d_xi_temporal_mode`, `tpf_4d_xi_temporal_coupling`.
- Legacy/direct only: `tpf_vdsg_coupling`, `tpfcore_enable_provisional_readout`, `tpfcore_readout_mode`, `tpfcore_readout_scale`, `tpfcore_theta_tt_scale`, `tpfcore_theta_tr_scale`, `tpf_kappa`, `tpf_cooling_fraction`, `tpf_global_accel_shunt_enable/fraction`.
- Diagnostic-only for xi route: `tpf_xi_kernel_dump_field_diagnostics`.

## Findings by severity
### P1
1. **[B/C/F] Derived radial star-star semantics mismatch**
   - File/function: `engine/physics/TPFCore/derived_tpf_radial.cpp`, `radial_acceleration_scalar_derived`.
   - Quote: `double M_stars_enc = enclosed_stellar_mass_cyl(state, r_cyl);`
   - Current behavior: enclosed stars are always included in baryonic bounced mass.
   - Expected/declared behavior: should either respect `enable_star_star_gravity` or metadata must clearly state it does not.
   - Safe fix?: metadata/reporting truth fix is safe; formula change needs theory decision.

### P2
2. **[D] Wake gate constants hard-coded (0.10, 0.05)**
   - File/function: `runtime_package_helpers.cpp`, `compute_xi_wake_kinematics`.
   - Behavior: constants are implementation settings, not config-driven.
   - Expected: acceptable if intentional; document as fixed contract.
   - Safe fix?: no formula change; documentation-only.

### P3
3. **[A] Xi routing/readout sign contract verified**
   - File/function: `tpf_core_package.cpp`, `compute_xi_kernel_deformed_accelerations`.
   - Behavior: inward acceleration for positive `K_xi` verified by explicit runtime sign audit and tests.

## Tests added
- `test_xi_kernel_runtime_audit.cpp`
  1. Xi-off sign/direction with `K_xi=1` matches Newtonian shape/sign with identical softening.
  2. `star_star=false` ignores stellar source-mass changes.
  3. `star_star=true` responds to stellar source-mass changes.

## Tests still missing (from requested matrix)
- explicit monotonic strengthening vs coupling with clamp-saturation checks on runtime integrator path.
- explicit metric_min saturation sweep assertion in a single deterministic unit test.

## Theory-decision questions for Dan
1. Should `derived_tpf_radial` continue to include enclosed stellar mass when `enable_star_star_gravity=false`, or should this path honor the global star-star toggle?
2. Should wake gate constants remain fixed (`0.10`, `0.05`) or be elevated to explicit config keys with provenance labels?
