# SimAll Beta

Industrial multiphysics CFD/CAE platform. C++20 / Qt6 / VTK / OpenCASCADE.

## Architecture

Strictly layered. Each subsystem is its own static library; cyclic deps are
forbidden by the build system:

```
utilities → core → io → parallel → gpu → materials
                                ↓
                              cad → meshing → solver
                                                ↓
                       turbulence | heat_transfer | radiation
                       combustion | multiphase    | scripting
                                                ↓
                                       visualization → gui
```

Ownership rules (enforced by header boundaries):

| Subsystem      | Owns                                           |
|----------------|------------------------------------------------|
| `cad`          | exact `TopoDS_Shape`, persistent topology IDs |
| `meshing`      | SoA mesh storage (`Mesh::nodes/faces/cells`)  |
| `solver`       | CSR matrices, field registry, BCs             |
| `visualization`| `vtkPolyData` actors only — NEVER exact geom  |
| `gui`          | Qt widgets only — never touches solver memory |

## Build

Requirements: CMake ≥ 3.22, MSVC 19.36+ or GCC 13+ / Clang 17+, Qt 6.5+,
VTK 9.2+, OpenCASCADE 7.6+, Eigen 3.4+.

```powershell
cmake --preset default
cmake --build build --config Release
build\bin\Release\simall_beta.exe
```

## Phased implementation status

See [`SimAll_Beta_Master_Architecture_Blueprint_Prompt.md`](SimAll_Beta_Master_Architecture_Blueprint_Prompt.md).
This Phase-1 commit delivers:

- Core (logger, event bus, service locator, undo/redo, application)
- Project graph + .simall serializer scaffold
- Qt6 ribbon + workflow tree + dynamic property editor + residual plot
- VTK viewport with picking → topology id translation
- OpenCASCADE STEP/IGES/STL import + healing + tessellation
- Mesh SoA storage (Section 6 layout) + surface-mesh ingest from CAD
- Finite-volume framework + CG linear solver
- Turbulence / heat / radiation / combustion / multiphase plugin contracts
- Parallel (MPI/OpenMP/TBB) and GPU façades
- Plugin ABI (`extern "C" CreatePlugin`)
- Catch2 unit tests
