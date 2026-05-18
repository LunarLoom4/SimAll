# Theory — Combustion

Combustion models live under [src/combustion/](../../src/combustion)
and implement `combustion::ICombustionModel`.

## Reaction modelling

| Model | Approach | File |
|---|---|---|
| Finite-rate (Arrhenius)        | direct evaluation of $k_f, k_r$ | `FiniteRateChemistry` |
| Eddy-Dissipation               | $\omega = \min(A\rho k/\varepsilon \cdot Y_F, B\rho k/\varepsilon \cdot Y_O/s)$ | `EddyDissipation` |
| Eddy-Dissipation Concept       | fine-structure reactor | `EddyDissipationConcept` |
| Flamelet (steady)              | $Z$ + scalar dissipation $\chi$ | `SteadyFlamelet` |
| Flamelet-Progress-Variable     | $Z, c$ + variance | `FpvModel` |
| Partially-Premixed             | $Z, c$ with $G$-equation | `PartiallyPremixed` |
| Conditional Moment Closure     | $\langle Y_k \mid Z\rangle$ | `CMC` |

## Species transport

Species mass fractions $Y_k$ are stored in `FieldRegistry` as a
vector field; one `solver::ScalarTransport` per species writes its
reaction source channel `S_Y_k`. The energy equation aggregates
heat-of-reaction into `S_T`.

## Mechanism reader

`combustion::ChemkinReader` parses Chemkin-II mechanisms and
thermodynamic/transport data, building a `KineticsLibrary` consumed
by all rate-based models.
