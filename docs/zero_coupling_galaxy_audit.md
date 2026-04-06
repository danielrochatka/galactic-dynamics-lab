# Zero-coupling galaxy-mode audit (falsification-first)

Scope: `tpf_dynamics_mode=direct_tpf`, `tpf_kappa=0`, `tpf_vdsg_coupling=0` in galaxy runs.

Findings:

1. **Force write path in direct_tpf with both couplings zero**
   - `compute_direct_tpf_accelerations` computes `ax_raw/ay_raw` from the principal-part tensor terms and then writes
     `ax[i] = kappa * ax_raw; ay[i] = kappa * ay_raw`.
   - With `kappa=0`, this write makes baseline `ax/ay` zero per particle.
   - The additive VDSG stage calls `accumulate_vdsg_velocity_modifier`, which immediately returns zeros when coupling is zero.
   - Therefore final accelerations are zero (barring NaN propagation from pathological inputs).

2. **No Newtonian fallback in direct_tpf dynamics path**
   - `TPFCorePackage::compute_accelerations` routes direct_tpf runs only through
     `compute_direct_tpf_accelerations` + `apply_vdsg_additive_extension` and returns.
   - Newtonian dynamics are not called in this branch.

3. **Galaxy initialization seeds nonzero circular velocity by design**
   - `initialize_galaxy_disk` sets `vx/vy` from `v_circ` in both init branches:
     - derived-radial init: `v_circ = sqrt(abs(a_s) * r)`
     - default init: `v_circ = sqrt(G * enclosed_mass / r)`
   - Both then assign `state.vx/vy = initial_velocity_scale * ...`.

4. **Acceleration buffers are not stale-reused as a source of force**
   - Integrator calls `physics->compute_accelerations(...)` each step before position update and again after drift for velocity update.
   - Direct TPF path explicitly reassigns `ax/ay` with zeros at entry before computing.

5. **No BH additive force branch when couplings are zero in direct_tpf**
   - BH enters field construction for direct_tpf baseline but is multiplied out by `kappa=0` at final write.
   - VDSG BH term is gated by nonzero `tpf_vdsg_coupling`; with zero coupling, contribution is zero.

6. **Integrator behavior with zero acceleration**
   - Velocity-Verlet becomes constant-velocity drift:
     - `x += vx*dt + 0.5*ax*dt^2` (ax=0)
     - `vx += 0.5*(ax + ax_new)*dt` (both zero)
   - This preserves initial velocities and produces outward motion for outward radial components present at init/noise.

First explanatory point for nonzero motion in this scenario:
- `initialize_galaxy_disk` nonzero circular-velocity assignment (`v_circ` -> `state.vx/vy`) prior to integration.
- With zero acceleration, this initial velocity field persists and advects stars.
