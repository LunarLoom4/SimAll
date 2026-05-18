# Mesh Pipeline

```text
cad::TopologyGraph
   │
   ▼
meshing::SizingField           (curvature + proximity + user sources)
   │
   ▼
┌──────────────────────────────────────────────────────────┐
│ Volume generators (pick one)                             │
│  • Cartesian + cut-cell  (CutCellHexMesher)              │
│  • Tetrahedral Delaunay  (DelaunayTetMesher)             │
│  • Poly-hexcore          (PolyHexCoreMesher)             │
│  • Structured hex sweep  (SweepHexMesher)                │
└──────────────────────────────────────────────────────────┘
   │
   ▼
meshing::BoundaryLayer         (prism inflation, smoothing, collapse)
   │
   ▼
meshing::QualityImprover       (Laplace, optimisation-based smoothing)
   │
   ▼
meshing::Mesh (SoA, zone table, persistent boundary IDs)
   │
   ▼
parallel::MeshPartitioner      (METIS / ParMETIS)
   │
   ▼
solver::FieldRegistry          (per-cell SoA fields)
```

## AMR / IBM / Overset hooks

- `amr::AdaptiveRefiner` re-meshes between time steps based on
  user-selected error indicators.
- `ibm::ImmersedBoundary` overlays solid surfaces onto a stationary
  background mesh; no remeshing required for moving bodies.
- `meshing::OversetInterpolation` couples multiple background meshes.

All three operate on `meshing::Mesh` in place and notify subscribers
via `EventBus::MeshChanged`.
