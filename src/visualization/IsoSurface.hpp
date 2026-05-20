// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/IsoSurface.hpp
// Phase  : 3 / Week 14 (filters: 3D iso-surface via marching tetrahedra)
//
// Per-tetrahedron Marching Tetrahedra (16 cases reducible to 2 unique
// configurations: 1 triangle when one vertex is on one side, 1 quad/2
// triangles when two vertices straddle the iso).  No ambiguity tables;
// deterministic on degenerate cases (vertex sits exactly on the iso).
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

namespace simall::visualization
{

class IsoSurface
{
public:
    /// Extract a triangulated iso-surface at `isovalue` from a tetrahedral
    /// volume mesh with per-vertex scalars.  The output is a SurfaceMesh
    /// whose pointScalars are filled with `isovalue` (constant) so that
    /// downstream coloring stays consistent.
    static SurfaceMesh extract(const VolumeMesh& volume, double isovalue);
};

} // namespace simall::visualization
