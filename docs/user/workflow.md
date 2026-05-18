# Workflow Tree

The left-hand dock of `simall_beta` shows a workflow tree mirroring
the project graph:

```text
Project
 ├─ Geometry
 │   ├─ Imported parts
 │   └─ Named selections
 ├─ Mesh
 │   ├─ Sizing fields
 │   ├─ Volume mesh
 │   └─ Boundary layer
 ├─ Physics
 │   ├─ Turbulence:       k-ω SST
 │   ├─ Combustion:       (none)
 │   ├─ Multiphase:       (none)
 │   ├─ Radiation:        (none)
 │   └─ Particles:        (none)
 ├─ Materials
 │   └─ Air @ 25 °C
 ├─ Boundary conditions
 │   ├─ inlet:    velocity-inlet 10 m/s
 │   ├─ outlet:   pressure-outlet 0 Pa
 │   └─ wall:     no-slip, T = 300 K
 ├─ Solver
 │   ├─ Scheme:   SIMPLE
 │   ├─ Linear:   AMG-PCG
 │   └─ Time:     steady
 ├─ Monitors
 └─ Reports
```

Each node maps to a `core::Command`; double-click to edit, right-click
for context actions. Changes propagate downstream via the project
graph — out-of-date nodes are shown with an orange badge.
