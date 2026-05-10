# Pushback Artifact — Architecture Reset Phase A (Config/Mode Isolation)

Date: 2026-05-10
Branch: `tpf-xi-theta-v1`

## θ_actual measurement (required pre-change scan)
Command run:

```bash
rg -n "TPF|tpf_|TPFCore|Xi|xi_|Theta|theta_|VDSG|direct_tpf|legacy_readout|geodesic_correspondence|v11_weak_field_truncation|C_mu_nu|C_μν" \
  engine/config.hpp engine/config.cpp
```

Measured counts:
- Total matches: **432**
- `engine/config.hpp`: **176**
- `engine/config.cpp`: **256**

Top coupling categories observed:
1. **Typed TPF fields in generic `Config`** (`tpf_*`, Xi/Theta/VDSG/readout/runtime knobs).
2. **TPF-only utility/benchmark mode ownership in generic `SimulationMode`** enum.
3. **TPF key parsing/validation in generic parser** (`apply_config_kv` branches).
4. **TPF key serialization in generic serializer** (`serialize_config_kv` branches).
5. **TPF helper naming/dispatch logic in generic config layer**.

## Exact blocker
A safe Phase A implementation cannot be completed *within the allowed file scope* without touching additional engine files that currently depend directly on TPF-typed members and TPF enum modes from generic `Config`.

The required compatibility-preserving migration is cross-cutting: moving TPF fields/modes out of `config.hpp/cpp` requires all existing consumers to switch from direct typed-field access to package-owned parsed config and package-owned mode routing.

Those consumers include (non-exhaustive):
- `engine/main.cpp`
- `engine/simulation.cpp`
- `engine/output.cpp`
- `engine/render_audit.cpp`
- `engine/force_compare.cpp/.hpp`
- `engine/galaxy_init.cpp/.hpp`
- `engine/resolved_scenario.cpp`
- generic unit/integration tests asserting TPF config/mode behavior in `engine/tests/**`

Because the instruction explicitly restricts touching those files except minimally for compile, fully removing TPF ownership from `config.hpp/cpp` would either break compile or silently break runtime compatibility.

## Exact config/mode compatibility risk
If TPF fields/modes are removed from generic `Config`/`SimulationMode` in this PR without coordinated consumer migration:
- Existing TPF config files would parse into storage that runtime code does not read.
- Existing user-facing TPF mode strings would no longer map to executable paths correctly.
- Run-info/render/audit paths relying on current mode and typed-field values would drift or fail.
- Unit/integration behavior coverage would become invalid (false pass or compile failures).

This violates the hard requirement to preserve legacy config/mode behavior while relocating ownership.

## Smallest coherent implementation slice that reduces θ safely
A safe, coherent slice that materially reduces θ while preserving behavior requires:
1. Add generic package option storage + generic mode token in `Config`.
2. Add package config/mode hook surfaces in `physics_package.hpp`.
3. Implement TPFCore-owned typed config adapter/parser/serializer under `engine/physics/TPFCore/**`.
4. Migrate at least the direct config/mode consumers (`main/simulation/output/resolved_scenario` minimum) to package hooks.
5. Move TPF config/mode tests from generic `engine/tests/unit/test_config.cpp` into `engine/physics/TPFCore/tests/**` and keep generic tests package-agnostic.

Without step 4, removing TPF ownership from `config.hpp/cpp` cannot preserve behavior.

## Expected θ reduction from that coherent slice
- `engine/config.hpp` + `engine/config.cpp`: from **432** toward **0** TPF matches (target), with any temporary residual references limited to generic package dispatch compatibility.
- Full engine boundary scan (excluding `TPFCore`) should drop significantly once consumer migration lands in same slice.

## Why this is pushback (and not scaffold-only)
This artifact records measured θ_actual, concrete blockers, compatibility risk, and the minimum coherent implementation sequence required to avoid breaking existing behavior. It is provided because implementing only `config.hpp/cpp` changes under current scope constraints would be unsafe and non-compliant with stated compatibility requirements.
