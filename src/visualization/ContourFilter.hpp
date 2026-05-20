// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ContourFilter.hpp
// Phase  : 3 / Week 14 (filters: surface isolines via marching triangles)
//
// Given a triangulated surface with per-vertex scalars, extract one or more
// isovalue curves as LineSet primitives.  Uses the 8-case marching-triangle
// table.  Linear interpolation along edges so contour vertices land exactly
// on the scalar level set within the precision of the input mesh.
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <vector>

namespace simall::visualization
{

class ContourFilter
{
public:
    /// Extract isolines at the supplied scalar levels.  The output LineSet
    /// stores 2-vertex polylines (one CSR group per intersection segment).
    /// pointScalars MUST have size == points.size().
    static LineSet extract(const SurfaceMesh& surface, const std::vector<double>& isovalues);

    /// Convenience: extract N uniformly-spaced isovalues between min..max.
    static LineSet extract_uniform(const SurfaceMesh& surface,
                                   double sMin,
                                   double sMax,
                                   std::size_t count);
};

} // namespace simall::visualization
