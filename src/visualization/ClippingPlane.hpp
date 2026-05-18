// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ClippingPlane.hpp
// Phase  : 3 / Week 14 (filters: half-space clip)
//
// Clip a triangle SurfaceMesh against an oriented plane.  Triangles whose
// vertices are entirely on the "above" side of the plane (n·(p-p0) >= 0)
// are kept verbatim; partially-straddling triangles are re-triangulated
// against the plane with linear scalar interpolation.  Vertex provenance
// is preserved so the renderer can keep using the original lookup table.
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

namespace simall::visualization {

struct Plane {
    util::Vec3d point;        // any point on the plane
    util::Vec3d normal;       // points to the "above" half-space
};

class ClippingPlane {
public:
    /// Clip `surface`.  When `keepBelow` is true the half-space test is
    /// inverted (handy for showing both halves side-by-side).  pointScalars
    /// are interpolated when present.
    static SurfaceMesh clip(const SurfaceMesh& surface,
                            const Plane& plane,
                            bool keepBelow = false);
};

}  // namespace simall::visualization
