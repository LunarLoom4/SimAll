// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/FeatureEdges.hpp
// Phase  : 4.5 (Feature edge extraction)
//
// Extracts the edges of a model that should be preserved during meshing:
//   * Sharp edges:    dihedral angle ≥ threshold
//   * Boundary edges: edges adjacent to a single face (free edges)
//   * Material edges: edges shared by faces with different attributes
// Returns polylines suitable for surface-mesher input (mark edges that
// MUST appear in the mesh, no matter how coarse).
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"
#include "utilities/MathTypes.hpp"

#include <vector>

namespace simall::cad
{

struct FeatureEdgeOptions
{
    double dihedralAngleDeg = 30.0;
    bool includeBoundary = true;
    bool includeNonManifold = true;
    bool includeSharp = true;
};

struct FeatureEdgePolyline
{
    std::vector<util::Vec3d> points;
    util::PersistentId edgeId = 0;
    bool isSharp = false;
    bool isBoundary = false;
    bool isNonManifold = false;
};

class FeatureEdges
{
public:
    std::vector<FeatureEdgePolyline> extract(const ShapeHandle& shape,
                                             const FeatureEdgeOptions& opts = {});
};

} // namespace simall::cad
