# TPF Conceptual Foundations
## Author: Daniel Rochatka
## Status: Reference document — anchor for all paper, code, and AI work

This document captures the conceptual framework of TPF in the
author's own conception, before any AI formalization. Any 
formalization, simulation, or paper claim must trace back to 
these foundations. Drift from this document is a bug.

## 1. Primary Objects

### Ξ_μ (Configuration Displacement Field)
A vector field representing the infinitesimal transformation 
between nearby configurations. Ξ is fundamental, not derived 
from a potential. It encodes the configuration of constraints 
acting at a point.

Ξ is application-dependent:
- For gravity: Ξ encodes gravitational configuration constraints
- For EM: Ξ would encode electromagnetic configuration constraints
- For other fields: domain-specific configurations

### Θ_μν (Configuration Gradient Tensor)
The covariant derivative of Ξ:
    Θ_μν = ∇_μ Ξ_ν

Θ describes how nearby configurations shift with respect to 
the background geometry. It reflects how matter limits 
changes around it. 

**Θ_μν REPLACES T_μν as the primitive source term in the 
field relation.** This is the essence of TPF: source structure 
without invoking energy as primitive.

### C_μν (TPF Field Equation Tensor)
The TPF analog of Einstein's G_μν. Renamed from G to avoid 
GR confusion. Obeys an on-shell divergence identity 
analogous to the Bianchi identity:

    ∇^μ C_μν = 0  (on shell)

The field equation:
    C_μν = κ[Θ machinery]

This is a field equation, NOT an acceleration formula.

## 2. How Particles Move

Particles travel in **straight lines through 4D spacetime** — 
geodesic motion in the metric that solves the field equation.

Acceleration is NOT directly computed from Ξ. The chain is:

    matter → Ξ → Θ → I/C_μν → geometry → geodesic motion

The discreteness of Ξ in any computation is a resolution 
limitation, not a fundamental feature. Higher resolution 
Ξ → finer trajectories. With infinite-resolution Ξ, problems 
like three-body would be solvable. We are precision-limited, 
not theory-limited.

## 3. Energy

Energy is **not** a primitive in TPF. Energy is the measurement 
of the change in states of matter. It is bookkeeping over 
configuration transitions, not a stored substance.

### Bucket of water example
A bucket at height h does not "contain" energy mgh. It is a 
configuration Ξ(n). When dropped, it becomes Ξ(n+1). The 
quantity mgh is the bookkeeping label for the transition 
Ξ(n) → Ξ(n+1). Neither Ξ(n) nor Ξ(n+1) contains energy. 
Only the change does.

### Battery example
A battery is a chemical configuration with available electron 
states. Connecting a circuit creates a vector space for 
electrons to traverse. The electron traversal is a configuration 
change — the displacement field changes. "Energy loss" is the 
bookkeeping of this configuration change. The battery never 
"stored" energy; it had a configuration that permitted certain 
transitions.

### Mass-energy ledger (Paper B Appendix E)
In the static weak-field correspondence, the rest-mass density 
emerges as:
    ρ_ledger ≡ Θ_ij Θ^ij

And:
    E = c² ∫ Θ_ij Θ^ij d³x
    
This recovers E = Mc² as derived bookkeeping, not as a 
fundamental source.

## 4. Source Specification

Ξ is defined by the matter configuration and boundary 
conditions. Ξ is NOT defined as the gradient of a Newtonian 
potential. The Newtonian gradient is what we GET in the 
weak-field limit, not what Ξ IS by definition.

A closed matter-to-Ξ map for general regimes is currently 
open. Paper A explicitly acknowledges this. The current 
weak-field implementation uses Ξ_i = ∂_i Φ as a calibrated 
correspondence, not as the general definition of Ξ.

## 5. Multi-Field Vision

TPF is intended to apply beyond gravity:
- Gravity (current focus)
- Electromagnetism (planned next)
- Other fields where configuration changes can be defined

Ξ takes different forms in each domain. Θ, I, and C_μν are 
the universal structural objects.

## 6. What TPF Is NOT

- TPF is NOT Newton's force law with different notation
- TPF is NOT energy conservation with relabeled variables
- Ξ is NOT the acceleration field
- C_μν is NOT the acceleration formula
- TPF is NOT a numerical approximation scheme for Newtonian gravity

## 7. Implementation Implications

The simulator must:
1. Specify Ξ from matter configuration (not from ∇Φ_Newton)
2. Compute Θ_μν = ∇_μ Ξ_ν
3. Compute I and C_μν as field-equation diagnostics
4. Solve for the geometry (or its weak-field metric perturbation)
5. Move particles via geodesics of the resulting geometry

The simulator must NOT:
1. Set Ξ = ∇Φ_Newton and call it TPF
2. Use C_μν as a direct acceleration readout (this is wrong)
3. Use Ξ as a direct acceleration readout (this is also wrong)
4. Add velocity-dependent kernel modifications to Ξ as if 
   that produces TPF physics (this is phenomenology, not TPF)

## 8. Drift Detection

If any AI session, paper revision, or code refactor produces 
output inconsistent with this document, that is a bug to be 
investigated, not a correction to be accepted. The framework 
itself is more stable than its formalizations.