# Monitors

Monitors are point, surface, or volume probes registered with
`solver::MonitorRegistry`. They evaluate every $N$ iterations and
broadcast a `MonitorSample` event over the EventBus.

## Built-in monitors

| Monitor | Quantity |
|---|---|
| ResidualMonitor       | scaled residual per equation |
| ForceMonitor          | pressure + viscous force on a wall set |
| MomentMonitor         | torque about a user point |
| MassFlowMonitor       | $\int \rho\mathbf{u}\cdot\mathbf{n}\,dS$ on a face zone |
| HeatFluxMonitor       | $\int q\,dS$ on a wall set |
| PointProbe            | any field at a point (trilinear interp) |
| SurfaceAverageMonitor | $\bar\phi = \int\phi\,dS / \int dS$ |
| VolumeIntegralMonitor | $\int\phi\,dV$ |

Monitors are persistent across runs and re-evaluated when the project
is reloaded.
