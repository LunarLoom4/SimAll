# CAD Pipeline

```text
file (STEP/IGES/BREP)
   │
   ▼
io::CadReader  ──►  OCCT XCAFDoc  (TopoDS_Compound + colours + names)
   │
   ▼
cad::CadKernel::import()
   ├─ ShapeFix_Shape    (heal small gaps)
   ├─ ShapeUpgrade_UnifySameDomain (merge co-planar faces)
   └─ assign cad::PersistentId per face / edge / vertex
   │
   ▼
cad::TopologyGraph        (immutable read-only view used by GUI)
   │                       └── visualization tessellation (BRepMesh_IncrementalMesh)
   │
   ▼
meshing::BoundaryRecovery  (snaps surface mesh to exact OCCT geometry)
```

## Healing strategy

1. `ShapeFix_Shape` — close hairline gaps, fix wire orientation.
2. `ShapeFix_FixSmallEdges` — collapse edges < tolerance.
3. `ShapeUpgrade_UnifySameDomain` — merge co-planar/co-cylindrical
   faces produced by re-importing tessellated geometry.
4. `BRepCheck_Analyzer` — final quality report; warnings surface in
   the GUI's CAD pane.

## Persistent IDs

Every B-Rep entity receives a `cad::PersistentId` on first import.
The ID survives healing, naming, and topology edits — boundary
conditions and named selections attach to IDs rather than indices.
