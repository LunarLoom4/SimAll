// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/SectionCut.hpp
// Phase  : 3 / Week 14 (filters: planar slice of a tetrahedral volume)
//
// Compute the polygonal intersection of an oriented plane with each tet of
// a VolumeMesh.  Per-vertex scalars are interpolated onto the slice so the
// downstream renderer can colour the cut surface with the same lookup
// table as the volume.
// =============================================================================
#pragma once

#include "visualization/ClippingPlane.hpp" // brings in Plane
#include "visualization/VisualizationTypes.hpp"

namespace simall::visualization
{

class SectionCut
{
public:
    /// Slice `volume` with `plane`.  Returns a triangulated SurfaceMesh
    /// embedded in the cut plane, with pointScalars populated when the
    /// volume carries them.
    static SurfaceMesh slice(const VolumeMesh& volume, const Plane& plane);
};

} // namespace simall::visualization
