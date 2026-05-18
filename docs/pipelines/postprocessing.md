# Post-processing Pipeline

Visualization is **read-only** on the solver's `FieldRegistry` and
the mesh.

```text
FieldRegistry  ──►  vis::FieldProbe (zero-copy view)
                       │
                       ▼
              vis::FilterGraph
                ├─ contour
                ├─ slice / clip / threshold
                ├─ streamlines / pathlines / LIC
                ├─ vector glyphs
                ├─ iso-surface
                ├─ volume rendering
                └─ statistics (avg, integral, flux)
                       │
                       ▼
                vtkActor (rendered in vis::Viewport)
```

## Update protocol

When the solver finishes an outer iteration it publishes
`EventBus::SolutionUpdated`. The visualization layer subscribes and
re-evaluates only the filters whose inputs changed. The Qt GUI
thread paints at the screen refresh rate; the filter pipeline runs
on `core::ThreadPool`.
