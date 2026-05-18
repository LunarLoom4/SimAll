# Boundary Conditions

| Type | Fields needed | Notes |
|---|---|---|
| velocity-inlet      | $\mathbf{u}$, $T$, turbulence quantities, $Y_k$ | scalar fields default to upstream value |
| mass-flow-inlet     | $\dot m$, $T$ | inlet normal velocity solved for |
| pressure-inlet      | $p_t$, $T_t$ | total quantities |
| pressure-outlet     | $p_s$ | back-flow values for scalars |
| outflow             | — | zero-gradient on all transported scalars |
| wall                | thermal BC (adiabatic / fixed T / fixed q) | slip / no-slip / moving |
| symmetry            | — | zero-flux for all scalars |
| periodic            | translation / rotation | matched pair |
| interface           | shared between zones | conjugate heat transfer |
| inlet-vent          | loss coefficient | combination of pressure inlet + porous jump |

Each BC is owned by a `solver::IBoundaryCondition` subclass and is
attached to a `cad::PersistentId` (not a face index) so it survives
geometry edits.
