# Theory — Multiphase Flow

Models live under [src/multiphase/](../../src/multiphase).

## Approaches

| Model | Description | File |
|---|---|---|
| Volume-of-Fluid (VOF)         | compressive interface (CICSAM / HRIC) | `VofSolver` |
| Geometric VOF (PLIC)          | piecewise-linear reconstruction | `PlicVof` |
| Level-Set                     | signed distance + reinitialisation | `LevelSet` |
| Coupled VOF + Level-Set       | sharp interface, conservative | `ClsVof` |
| Eulerian-Eulerian (n-fluid)   | per-phase momentum + drag closures | `EulerianMultifluid` |
| Mixture model                 | algebraic slip velocity | `MixtureModel` |
| Cavitation (Schnerr–Sauer)    | mass-transfer source on $\alpha$ | `Cavitation` |
| Phase-change (Lee)            | evaporation/condensation source | `LeeMassTransfer` |

## Surface tension

Continuum Surface Force (CSF) with curvature from a smoothed VOF
gradient; height-function curvature available for high accuracy.

## Coupling

Mass-transfer models push to the species channel `S_Y_k` and to the
energy channel `S_T` (latent heat). The momentum equation receives
$\sigma\kappa\nabla\alpha$ from the surface-tension term.
