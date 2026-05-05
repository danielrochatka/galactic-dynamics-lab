# TPF Conceptual Foundations

## Author: Daniel Rochatka

## Status: Canonical reference document — anchor for paper, code, and AI work

This document captures the conceptual framework of TPF in the author's own conception, before
any later formalization or simulator shortcut is allowed to drift the theory. Any paper claim,
simulation route, or AI-assisted rewrite must trace back to these foundations. Drift from this
document is a bug unless the foundations are intentionally revised.

TPF is built on three core commitments:

1. Recursive geometric state evolution is fundamental.
2. Spacetime is not treated as fundamental.
3. Energy is a ledger/readout of configuration change, not a stored substance.

General relativity remains the effective macroscopic correspondence language in the regimes
where its tested predictions apply. TPF does not begin with spacetime curvature as a primitive
substance. It begins with recursive configuration states and the geometric-gradient structure
derived from them.

---

## 1. Primary Objects

### Ξ_μ / Ξ^μ — Configuration-State Field

Ξ is the primitive configuration-state field.

At a given state interval `n`, the system is represented by a configuration state:

```text
Ξ_n
```

Transitions occur between states:

```text
Ξ_n -> Ξ_{n+1}
```

Ξ is therefore a state object, not itself the transition. The transition is the recursive
update from one Ξ state to the next.

In continuum notation, Ξ may be represented as a vector field:

```text
Ξ^μ(x)
```

This continuum representation should not be mistaken for a commitment that spacetime is
fundamental. Coordinates are the representation language used to express the state and compare
with standard physics.

Ξ is application-dependent:

- For gravity: Ξ encodes gravitational configuration constraints.
- For EM: Ξ would encode electromagnetic configuration constraints.
- For other fields: Ξ encodes domain-specific configuration structure.

In simulation, Ξ may be represented continuously, on a discrete grid, or at lattice sites
depending on scope. These are representations of the same configuration-state object, not
different ontological objects.

For n-body intuition, Ξ is the summed configuration/influence field at each sampled point.
At a spatial point `x`, each source contributes a vector-valued field contribution, and those
contributions sum to the total Ξ at that point. Thus Ξ tells the simulator the net local
direction/strength of the configuration field in the selected sector.

A useful visualization is `(x, y, z, field strength)`, but mathematically this should be read
as a field over three-dimensional space, not as a literal fourth spatial coordinate. More
precisely, the field is:

```text
Ξ_i(x, y, z)
```

or, at a state interval:

```text
Ξ_i^(n)(x, y, z)
```

The field strength or magnitude `|Ξ|` is a value assigned to each spatial point. It can be
visualized as a landscape or spacetime-like map of curves, but the value is not itself an
additional spatial dimension.

### Θ_μν — Configuration-Gradient Tensor

Θ is the covariant derivative of Ξ:

```text
Θ_μν = ∇_μ Ξ_ν
```

or, in a spatial weak-field implementation:

```text
Θ_ij = ∂_i Ξ_j
```

Θ is derived from Ξ. It is not a separate primitive field.

Θ describes how the configuration state varies locally. It is the core geometric-gradient
object of TPF. It reflects how matter/source structure limits or shapes possible changes
around it.

Θ is the natural place to describe:

- local gradient structure,
- tidal structure,
- curvature-like behavior,
- extended-body response candidates,
- nonlinear configuration-intensity diagnostics,
- sector-specific geometric response.

In regimes where both GR and TPF descriptions apply, GR spacetime curvature and TPF
configuration-gradient structure describe the same observable response behavior in different
languages. TPF treats Θ as the lower-level geometric-gradient object. GR expresses the
corresponding macroscopic behavior as spacetime curvature.

Important distinction:

```text
Ξ = configuration-state field
Θ = derived configuration-gradient structure
Ξ_n -> Ξ_{n+1} = recursive state transition
```

Do not collapse these three ideas into one.

Operationally:

```text
Ξ tells the net field value at a point.
Θ tells how that net field changes around that point.
```

This distinction matters in cancellation regions. Two equal source weights can produce a
midpoint where the summed Ξ field vanishes:

```text
Ξ_total = 0
```

but the gradient structure need not vanish:

```text
Θ = ∇Ξ ≠ 0
```

The midpoint is balanced, but it is not geometrically featureless. A small displacement from
the midpoint breaks the balance. Θ captures that local topology: the tidal, curvature-like,
and instability structure of the summed field.

### C_μν — TPF Field-Response Tensor

C_μν is the TPF field-response tensor used in the full geometric route. It is analogous in role
to Einstein's G_μν, but it is built from TPF configuration-gradient machinery rather than from
energy-momentum as a primitive source.

The intended full field-equation route is:

```text
source configuration
  -> Ξ
  -> Θ = ∇Ξ
  -> C_μν[Θ machinery]
  -> effective geometry / GR correspondence structure
  -> geodesic-like motion
```

C_μν is a field-equation object, not an acceleration formula.

A schematic field relation is:

```text
C_μν = κ[Θ machinery]
```

with an on-shell divergence identity analogous in role to the Bianchi identity:

```text
∇^μ C_μν = 0  (on shell)
```

This is the publishable geometric route of TPF. It should not be confused with every current
simulator closure or weak-field benchmark shortcut.

---

## 2. Primitive Inputs

In the present gravitational sector, the primitive inputs are:

- configuration state Ξ,
- source weight / configuration weight,
- source position,
- source type,
- gravitational coupling scale G.

The local near-surface scale `g` is not an independent primitive coupling. It is a derived
local approximation, for example:

```text
g ≈ G M_Earth / R_Earth^2
```

Mass/source weight is a primitive input in the current gravitational sector. TPF does not need
to derive G or source weight from deeper principles in the present paper.

Energy is not primitive. Density, stress-energy, potential energy, and momentum are also not
primitive drivers of the state update. They are conventional ledger/readout quantities used in
appropriate correspondence regimes.

---

## 3. How Particles Move

The full geometric route is:

```text
matter/source configuration
  -> Ξ
  -> Θ
  -> C_μν / effective geometry
  -> geodesic-like motion
```

In this route, particles do not move because energy acts as a substance or because Newtonian
force is primitive. They move according to the effective geometry generated by recursive
configuration-gradient structure.

This is the route most directly comparable to GR.

However, current simulator routes may use weaker correspondence closures. A weak-field
Xi-direct route may read point-matter acceleration from Ξ:

```text
a_i(x) = -G Ξ_i(x)
```

or, in a route with a configured readout scale:

```text
a = -K_xi Ξ_eff,spatial
```

This is a calibrated weak-field correspondence readout, not the primitive ontology of TPF.

Therefore:

```text
Full theory route:
  Ξ -> Θ -> C_μν / geometry -> geodesic-like motion

Weak-field simulator closure:
  Ξ -> acceleration readout -> position update
```

The second route is useful for correspondence testing and simulation, but it must not be
presented as the complete TPF field-equation route.

---

## 4. Source Specification

Ξ is defined by source configuration, source weights, sector rules, and boundary conditions.

A closed matter-to-Ξ map for all regimes is currently open. The present weak-field
implementation may use correspondence-calibrated source forms, but those forms should not be
mistaken for the complete general definition of Ξ.

For a weak-field point-source gravitational correspondence, a source contribution may be
represented schematically as:

```text
r_a = x - x_a
Ξ_a,i(x) = W_a r_a,i / |r_a|^3
```

where:

- `x` is the field/sample position,
- `x_a` is the source position,
- `W_a` is the source weight,
- softening may be applied in code for numerical regularization.

This weak-field expression is a correspondence-sector source rule. It should not be presented
as the final general TPF source law.

The Newtonian gradient is what the weak-field correspondence may reproduce. It is not the full
definition of Ξ in all regimes.

---

## 5. Field Aggregation and Multi-Source Systems

In the current simulator architecture, each sector constructs its field representation by
summing the per-source contributions assigned by that sector:

```text
Ξ_total^(s)(x) = Σ_a Ξ_a^(s)(x)
```

where `s` labels the sector.

This aggregation rule is an implementation contract for current TPFCore field construction.
It should not be overread as a final universal claim that every physical sector in the full
theory is linearly superposed in every regime.

Nonlinear behavior may enter through:

- the per-source contribution function,
- boundary conditions,
- sector coupling,
- kernel deformation,
- the motion/update closure,
- the field-response tensor,
- the ledger/readout map.

The summation layer should remain explicit and auditable.

### Cancellation and topology

In a multi-source system, field contributions may cancel at specific points. For example, two
equal source weights placed symmetrically can create a midpoint where the net Ξ field is zero.

For a weak-field point-motion readout, this means:

```text
a_i = -G Ξ_i = 0
```

at the exact balance point.

But this does not mean the field geometry has disappeared. The surrounding field can still
have nonzero gradient structure:

```text
Θ_ij = ∂_i Ξ_j
```

A cancellation point is therefore a useful example of the difference between field value and
field topology:

```text
Ξ = net local value / direction / strength
Θ = local variation / topology / tidal structure of Ξ
```

This is why Θ is not redundant with Ξ. Θ records how the summed field changes across nearby
points, even when the summed field itself vanishes at a point.

---

## 6. Recursion and State Evolution

TPF is a recursive geometric state framework. The primitive conceptual sequence is:

```text
Ξ_1, Ξ_2, ..., Ξ_n, Ξ_{n+1}, ...
```

The central physical idea is not that an object moves through a pre-existing spacetime
substance. The central idea is that the system evolves through recursively related geometric
configuration states.

The full theoretical transition is:

```text
Ξ_n -> Ξ_{n+1}
```

A complete primitive transition law remains part of the deeper TPF program. Current simulator
routes may implement specific closures or correspondence shortcuts rather than the complete
recursive state law.

### Current quasi-static field-rebuild closure

The current quasi-static n-body simulator may use the following operational loop:

1. Read current source positions and properties.
2. Each source writes its sector contribution.
3. Sum contributions to obtain the current Ξ field.
4. Derive Θ = ∇Ξ when needed.
5. Objects read the appropriate sector field/readout.
6. Objects move according to the selected dynamics closure.
7. Updated source positions become input for the next step.
8. Rebuild Ξ from the updated source configuration.

In this simulator closure:

```text
Ξ_{n+1}
```

is rebuilt from the updated source configuration after the step.

This is a quasi-static / correspondence-sector implementation. It should not be confused with
the full primitive Ξ_n -> Ξ_{n+1} transition law.

---

## 7. Energy

Energy is not a primitive in TPF.

Energy is a scalar ledger/readout associated with weighted configuration change. It is
bookkeeping over configuration transitions, not a stored substance.

The ledger is the accounting/readout structure associated with weighted intervals of Ξ-state
change.

Energy is one scalar expression of that ledger. Density, stress-energy, potential energy, and
momentum are also conventional readouts or representations of the ledger in specified sectors
or correspondence regimes.

### Bucket of water example

A bucket at height `h` does not "contain" energy `Mgh`.

It is a configuration at one state interval relative to a boundary/reference configuration.
When dropped, the system transitions through new configurations.

The quantity `Mgh` is the scalar bookkeeping label assigned to the weighted geometric relation
between the suspended configuration and the boundary/reference configuration.

Neither Ξ_n nor Ξ_{n+1} contains energy as a substance. The scalar energy value is a ledger
readout of the weighted relation/transition.

In this example:

```text
configuration relation:
  separation h from boundary/reference state

source weight:
  M

local field scale:
  g ≈ G M_Earth / R_Earth^2

ledger scalar:
  Mgh
```

### Battery example

A battery is a chemical configuration with available electron states. Connecting a circuit
creates a vector space for electrons to traverse. The electron traversal is a configuration
change.

"Energy loss" is the bookkeeping of this configuration change. The battery never stored energy
as a substance; it had a configuration that permitted certain transitions.

### Ledger readout examples

Different ledger readouts answer different physical questions:

```text
Mgh-type ledger:
  scalar energy value from weighted geometric separation relative to a boundary/reference state

ρ_Ξ:
  ρ_Ξ = (1 / 4π) ∂_i Ξ^i
  linear flux density/source readout in the static weak-field sector

Θ_ij Θ^ij:
  nonlinear configuration-intensity readout

stress-energy:
  macroscopic tensor representation of the ledger in the relativistic / GR correspondence regime
```

No single ledger projection exhausts the full accounting structure. Each readout answers a
specific physical question in a specified sector or correspondence regime.

The ledger does not drive the state transition. It is an accounting/readout layer.

Important correction:

```text
ρ_Ξ = (1 / 4π) ∂_i Ξ^i
```

is the leading scalar density/source readout in the static weak-field sector.

```text
Θ_ij Θ^ij
```

is a nonlinear configuration-intensity readout. It should not be used as the leading Newtonian
density/source readout.

---

## 8. Source, Test Object, and Motion Roles

Sources and test objects are practical role labels, not separate ontological categories.

A massive object may both:

- write to a sector field as a source,
- read from the resulting total field as a responding configuration.

The distinction is operational. It depends on what is being modeled and what approximation is
being used.

For weak-field point matter:

```text
motion may read Ξ
```

For tidal, extended-body, nonlinear, or curvature-like diagnostics:

```text
structure reads Θ = ∇Ξ
```

This is the practical simulator distinction:

```text
Ξ determines the net local point-motion readout in the weak-field closure.
Θ determines the local shape/topology of the field around that point.
```

For light or EM propagation, TPFCore needs a sector-specific map from Ξ / Θ into EM propagation
behavior. Do not silently reuse the point-matter acceleration readout for light.

If a route reads acceleration directly from Θ, it must define the closure or contraction that
maps Θ to acceleration. Otherwise, Θ should be treated as gradient/tidal/curvature-like
structure, not as ordinary point acceleration.

---

## 9. Multi-Field Vision

TPF is intended to apply beyond gravity.

Candidate sectors include:

- gravity,
- electromagnetism,
- exploratory velocity-dependent geometric sectors,
- future lattice/material/configuration sectors.

The shared structural pattern is:

```text
Ξ^(s) -> Θ^(s) = ∇Ξ^(s)
```

where `s` is the sector.

The statement "TPF can include gradients beyond spacetime" means that different sectors may
define different configuration fields and corresponding gradient structures.

It does not mean that arbitrary gradients with incompatible units or meanings can be added
without a coupling rule.

Any combined response structure must define:

- which sector fields are included,
- how their gradients are represented,
- what coupling maps are used,
- what readout or transition rule follows.

Θ, intensity objects, and C_μν-type response structures are the universal structural ideas.
Their sector-specific meanings must be defined explicitly.

---

## 10. Correspondence with Standard Physics

In the static weak-field sector with point-source contributions and inverse-square behavior,
TPFCore should reproduce the Newtonian point-motion correspondence through:

```text
Ξ_i(x) = Σ_a W_a r_a,i / |r_a|^3
a_i(x) = -G Ξ_i(x)
```

This is the weak-field correspondence endpoint. It is not the full TPF primitive engine.

GR remains the effective macroscopic correspondence theory in regimes where its tested
predictions apply. TPF differs in primitive interpretation:

```text
GR:
  spacetime geometry + stress-energy

TPF:
  recursive configuration states + derived gradient structure + ledger readouts
```

In shared tested regimes, TPF and GR should agree on observables. In extensions beyond
established correspondence regimes, TPF may allow different predictions, but those require
explicit sector definitions, mathematical closure, and tests.

TPF and GR are therefore not in competition as numerical descriptions inside GR's tested
regime. They differ in what they treat as primitive.

---

## 11. What TPF Is NOT

- TPF is not Newton's force law with different notation.
- TPF is not energy conservation with relabeled variables.
- TPF is not a claim that energy is a stored substance.
- TPF is not a claim that spacetime curvature is fundamental.
- Ξ is not a transition by itself; transitions occur between Ξ_n and Ξ_{n+1}.
- Θ is not a separate primitive field; Θ is derived from Ξ.
- C_μν is not an acceleration formula.
- The weak-field Xi-direct acceleration readout is not the full primitive theory.
- The current simulator is not proof of the full recursive TPF engine.
- VDSG or velocity-dependent kernel deformation is exploratory unless explicitly formalized and tested.

---

## 12. Implementation Implications

The full geometric simulator route should:

1. Specify Ξ from source configuration, source rules, and boundary conditions.
2. Compute Θ_μν = ∇_μ Ξ_ν.
3. Compute intensity / response diagnostics as needed.
4. Compute C_μν or an effective geometric correspondence structure.
5. Move particles via geodesic-like motion in the resulting effective geometry.

The weak-field correspondence simulator route may:

1. Build Ξ from source positions and source weights.
2. Sum source contributions.
3. Read point-matter acceleration from Ξ.
4. Move objects.
5. Rebuild Ξ from the updated source configuration.

These are different routes. They must be labeled separately.

The simulator must not:

1. Set Ξ = ∇Φ_Newton and call that the complete TPF theory.
2. Use C_μν as a direct acceleration formula.
3. Treat Xi-direct acceleration as the full primitive ontology.
4. Treat Θ_ij Θ^ij as the leading weak-field density/source readout.
5. Add velocity-dependent kernel modifications to Ξ and call that confirmed TPF physics.
6. Allow a ledger readout to become a hidden transition law.

Before adding or changing a route, answer:

```text
What is the state variable?
What is derived from the state?
How is the field constructed?
What determines motion or the next state?
What is merely a diagnostic or ledger readout?
Which sector does this route claim to represent?
```

If the answer is unclear, the route should be marked exploratory or deferred.

---

## 13. Current TPFCore Route Alignment

Current runtime and benchmark routes must be interpreted according to their closure.

### xi_kernel_deformed

Xi-direct runtime route.

- Computes per-source Ξ.
- Optionally applies Xi-kernel deformation.
- Sums Ξ_eff.
- Reads acceleration as:

```text
a = -K_xi * Ξ_eff,spatial
```

This is a weak-field/correspondence-style simulator route. It should not be described as the
full primitive TPF field-equation route.

### direct_tpf

Tensor principal-part route.

This route consumes field-evaluation objects and may use Θ / intensity / principal-tensor
structure for a paper-facing or benchmark-facing closure. It should remain clearly separated
from the Xi-direct acceleration route.

### v11_weak_field_truncation

Correspondence-control route.

This route is useful for verifying weak-field behavior but should not be treated as validation
of the full primitive recursive TPF engine.

### legacy_readout

Deprecated / compatibility route.

Legacy provisional readout behavior should remain quarantined and explicitly labeled when used.

---

## 14. Naming Discipline

Use these labels consistently:

```text
Ξ:
  configuration-state field

Θ:
  derived configuration-gradient tensor

Ξ_n -> Ξ_{n+1}:
  recursive state transition

field rebuild:
  quasi-static simulator closure from updated source positions

ledger readout:
  derived accounting/readout quantity

source weight:
  configuration/source role; use "mass" when speaking in correspondence language

G:
  gravitational coupling scale

g:
  derived local field scale, not independent primitive

C_μν:
  field-response tensor, not acceleration formula
```

Avoid saying:

```text
Ξ is a transition field.
Θ stores forces.
Energy drives the transition.
Test objects always read Θ.
Ξ-direct acceleration is the full primitive theory.
Spacetime curvature is fundamental in TPF.
```

Prefer saying:

```text
Ξ is the configuration-state field.
Θ is the derived configuration-gradient tensor.
Point-matter acceleration may read Ξ in the weak-field correspondence sector.
Θ supplies gradient, tidal, curvature-like, and nonlinear response structure.
Energy is a scalar ledger readout.
The current n-body simulator may rebuild Ξ each step as a quasi-static closure.
The full geometric route proceeds through Θ/C/effective geometry/geodesic-like motion.
```

---

## 15. Drift Detection

If any AI session, paper revision, or code refactor produces output inconsistent with this
document, that is a bug to be investigated, not a correction to be accepted.

The most common drift errors are:

1. Turning Ξ back into a transition instead of a state.
2. Treating energy as primitive.
3. Treating the ledger as a driver of dynamics.
4. Treating Θ_ij Θ^ij as the leading weak-field density source.
5. Treating Xi-direct acceleration as the full primitive theory.
6. Treating the full geometric route as already implemented by every simulator mode.
7. Using "spacetime curvature" as if it were TPF's primitive ontology.
8. Allowing simulator shortcuts to rewrite the foundations.
9. Treating exploratory VDSG/kernel deformation as validated physics.
10. Letting GR correspondence language erase the lower-level TPF ontology.
11. Collapsing Ξ and Θ into the same object, especially in cancellation regions where Ξ may vanish while Θ remains nonzero.

The framework should remain layered:

```text
Core ontology:
  recursive Ξ-state mechanics and derived Θ-gradient structure

Full geometric route:
  Ξ -> Θ -> C_μν -> effective geometry -> geodesic-like motion

Weak-field simulator closure:
  source positions -> Ξ field -> Xi-direct acceleration readout -> position update -> field rebuild

Ledger layer:
  energy, density, stress-energy, potential, and momentum as derived readouts
```