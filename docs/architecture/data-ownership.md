# Data Ownership

Per spec §"DATA OWNERSHIP MODEL".

| Owner | Owns | Forbidden access |
|---|---|---|
| `cad::CadKernel` | exact `TopoDS_Shape`, persistent topology IDs | mesh nodes/cells, solver matrices |
| `visualization::SceneGraph` | `vtkActor`, `vtkPolyData`, scalar maps | exact geometry, mesh storage |
| `meshing::Mesh` | SoA node/face/cell arrays, zone table | solver matrices, CAD shapes |
| `solver::FieldRegistry` | per-cell scalar/vector fields | mesh topology, GUI widgets |
| `solver::CSRMatrix` | one matrix per equation per outer iteration | field arrays it does not own |
| `core::Project` | dependency graph, undo stack | physics arrays directly |

Cross-references travel as **handles** (`meshing::ZoneId`,
`cad::PersistentId`, `solver::FieldHandle`) — never raw pointers
between layers.
