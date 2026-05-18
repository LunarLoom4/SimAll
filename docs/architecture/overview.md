# Architecture Overview

SimAll Beta is built as a strict layered system. Every subsystem listed
in [the master blueprint](../../SimAll_Beta_Master_Architecture_Blueprint_Prompt.md)
lives in its own static library under `src/`, and its inter-module
edges are checked by CMake (no cycles permitted).

![Layers](architecture.puml)

## Subsystem inventory

| Layer | Module(s) | Owns |
|---|---|---|
| Presentation | `gui`, `visualization` | windowing, ribbon, panels, VTK actors |
| Workflow | `core` (Command, EventBus, ServiceLocator), `io::Project` | project graph, undo/redo, dependency invalidation |
| Physics | `turbulence`, `combustion`, `multiphase`, `radiation`, `heat_transfer`, `particles`, `acoustics`, `emag`, `adjoint`, `rom` | constitutive models, transport-equation closures |
| Numerics | `solver` | FV assembly, SIMPLE/PISO/Coupled, linear solvers, gradients, limiters |
| Discretisation | `meshing`, `amr`, `ibm`, `morphing`, `rotating`, `dynamics` | mesh data, mesh generators, IBM, mesh motion |
| CAD | `cad` | exact OpenCASCADE B-Rep, topology graph, healing |
| Foundation | `core`, `utilities`, `parallel`, `gpu`, `io`, `scripting`, `materials` | runtime services, math primitives, MPI/OMP/TBB, CUDA, file I/O, Python, materials DB |

## Non-negotiable contracts

1. **GUI mutates state only via `core::ICommand`.** Every change is
   undo/redo capable.
2. **Cross-module communication goes through `core::EventBus`.** No
   direct includes of GUI code from physics, and no direct includes of
   solver code from GUI.
3. **CAD ≠ visualization ≠ mesh.** Three independent representations
   live side by side; the picking bridge maps between them.
4. **Solver is physics-agnostic.** Models register source terms via
   `solver::FieldRegistry` channels — never patch the solver source.
5. **Storage is SoA + 64-byte aligned.** No AoS in hot loops.
6. **GUI thread never blocks.** All long work runs on
   `core::ThreadPool` / `QtConcurrent`.

See [Layers](layers.md) for the precise dependency graph and
[Data ownership](data-ownership.md) for ownership boundaries.
