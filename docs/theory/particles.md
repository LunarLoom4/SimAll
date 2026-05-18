# Theory — Particles

`particles::LagrangianTracker` evolves parcels under

$$
m_p \frac{d\mathbf{u}_p}{dt} =
\mathbf{F}_{drag} + \mathbf{F}_{lift} + \mathbf{F}_{pressure}
+ \mathbf{F}_{vm} + m_p \mathbf{g}
$$

with sub-models for drag (Schiller–Naumann, Morsi–Alexander),
lift (Saffman, Magnus), virtual mass, Brownian motion, and
turbulent dispersion (Discrete Random Walk).

## Sub-models

| Sub-model | File |
|---|---|
| Drag library              | `DragModels` |
| Spray break-up (KH-RT, TAB) | `SprayBreakup` |
| Droplet evaporation       | `Evaporation` |
| Wall-film deposition      | `WallFilm` |
| Particle-particle (DEM)   | `Dem` |

## Two-way coupling

The tracker accumulates per-cell momentum source $-\sum_p
\mathbf{F}_{drag}/V_P$ on `S_U`, energy source on `S_T`, and species
source on `S_Y_k` (evaporated mass). The solver picks them up on the
next outer iteration.
