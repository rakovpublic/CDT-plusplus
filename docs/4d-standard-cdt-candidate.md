# 4D Standard-CDT Candidate Validation

This implementation treats `cdt -d4` as a standard-CDT candidate only when the
state passes the validator in `Foliated_triangulation_4.hpp`. The validator
requires an explicit persistent vertex/simplex incidence complex, periodic
time, connected closed `S3` spatial slices, non-negative simplex counts,
non-negative proposal multiplicities derived from local sites, reciprocal
gluing for combinatorial simplices, and exact agreement between `N4` and the
sum of `N41`, `N32`, `N23`, and `N14`.

Count-only 4D states are deliberately rejected as standard CDT candidates. They
remain useful for validation fixtures and checkpoint compatibility, but the
production proposal kernel does not infer local moves from aggregate counts.

The action convention follows the usual 4D CDT bare-coupling form

```text
S_E = -(kappa_0 + 6 Delta) N0
      + kappa_4 N4
      + Delta (2 N41_total + N32_total)
      + epsilon (N4 - target_N4)^2
```

where `N41_total = N41 + N14` and `N32_total = N32 + N23`.
The numerical setup follows the standard CDT requirements that Monte Carlo
moves preserve fixed topology, causality and detailed balance; see the CDT
transfer-matrix lecture notes for the seven-move 4D setup and for the
`cos^3` de Sitter spatial-volume profile and `N4^(1/4)`, `N4^(3/4)` scaling.

Sources:

- https://indico.tpi.uni-jena.de/event/76/attachments/38/132/Goerlich_-_Causal_Dynamical_Triangulations.pdf
- https://www.scholarpedia.org/article/Causal_Dynamical_Triangulation

## Move Catalogue

The source of truth is `Move_catalog_4.hpp`. At present the persistent
incidence kernel production-enables only the verified `2<->4` pair. The other
standard 4D Pachner-like CDT move names remain in the public enum for API and
checkpoint compatibility, but they enumerate zero proposal sites until their
local subcomplex replacement rules are independently implemented and validated.

| Move | Inverse | Proposal multiplicity | Invariant delta `(N0,N1,N2,N3,N4)` |
| --- | --- | --- | --- |
| `TWO_FOUR` | `FOUR_TWO` | legal internal tetrahedral facets | `(0,+1,+4,+5,+2)` |
| `FOUR_TWO` | `TWO_FOUR` | removable order-four causal edges | `(0,-1,-4,-5,-2)` |
| `THREE_THREE` | `THREE_THREE` | not implemented | `(0,0,0,0,0)` |
| `FOUR_SIX` | `SIX_FOUR` | not implemented | `(0,0,0,0,0)` |
| `SIX_FOUR` | `FOUR_SIX` | not implemented | `(0,0,0,0,0)` |
| `TWO_EIGHT` | `EIGHT_TWO` | not implemented | `(0,0,0,0,0)` |
| `EIGHT_TWO` | `TWO_EIGHT` | not implemented | `(0,0,0,0,0)` |

The class-resolved simplex-type deltas `N41/N32/N23/N14` are not fixed by the
move name alone; they are derived from the selected local site and recorded in
the accepted `MoveApplication`.

## Detailed Balance

For a proposed transition `A -> B`, the sampler uses

```text
acceptance = min(1, exp(-(S(B)-S(A))) * q(B -> A) / q(A -> B))
```

where `q` is computed from the current proposal inventory, not from historical
move frequencies. The enumerable detailed-balance tests build a small reachable
state graph, verify reverse transitions, and compare both sides of

```text
pi(A) q(A -> B) A(A -> B) = pi(B) q(B -> A) A(B -> A)
```

with `pi(T) = exp(-S(T))`.

## Phase Diagnostics

The single-run summary reports the conservative profile verdict and the
measured `cos^3` profile correlation for the measurements available in that
run. It does not report likelihood, AIC, or BIC values. A full
`c_ds_supported` finite-size claim uses `diagnose_c_ds_finite_size()` and
requires all of the following:

- centered ensemble profile is better fit by `cos^3` than collapsed or
  alternating-slice alternatives;
- finite-size width scales as `N4^(1/4)`;
- finite-size peak volume scales as `N4^(3/4)`;
- covariance/effective-action extraction has enough decorrelated samples;
- C_b diagnostics are not triggered.

Synthetic tests exercise the analysis layer. Production claims still require
real independent chains and decorrelated measurements written under
`results/<run-id>/`.
