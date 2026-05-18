# Solver Pipeline

```text
meshing::Mesh + materials::MaterialDatabase
   │
   ▼
solver::FieldRegistry::initialize()      (allocate SoA arrays per equation)
   │
   ▼
solver::BoundaryConditionManager::apply()
   │
   ▼
┌───────────────────────────────────────────────────────────┐
│ Outer iteration (SIMPLE / PISO / Coupled)                 │
│   for each equation E in {U, p, T, k, ω, Y_k, ...}        │
│      • solver::Gradient::compute(E)                       │
│      • solver::Convection::assemble(E)                    │
│      • solver::Diffusion::assemble(E)                     │
│      • physics modules add sources via FieldRegistry      │
│      • solver::LinearSolver::solve(A·x = b)               │
│   end                                                     │
│   monitor residuals; check convergence                    │
└───────────────────────────────────────────────────────────┘
   │
   ▼
io::SolutionWriter (every N steps)   +   monitors broadcast on EventBus
```

## Coupling protocol

Each physics module owns one or more **source channels** registered
with `FieldRegistry`:

| Channel | Owner | Consumed by |
|---|---|---|
| `S_U`, `S_p`            | combustion, multiphase, particles, IBM | momentum equation |
| `S_T`                   | combustion, radiation, particles | energy equation |
| `S_k`, `S_omega`        | particles (turbulence modulation) | turbulence model |
| `S_Y_k`                 | combustion, multiphase mass transfer | species equation |

The solver never references a physics module directly. Modules push
to source channels; the solver pulls during assembly.
