// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/OctreeMesher.hpp
// Phase  : 5.6 — Octree-based hex volume mesher from a closed STL surface.
//
// Algorithm:
//   1. Build a bbox-aligned octree, refine to maxDepth where any STL
//      triangle intersects the cell AABB (boundary refinement).
//   2. For every leaf cell, classify as INSIDE / OUTSIDE / BOUNDARY using a
//      stochastic ray-cast inside test against the surface.
//   3. Emit INSIDE leaves as hex cells. BOUNDARY leaves are kept as
//      uniformly-refined hexes (true cut-cells require ANSA-grade clipping
//      and arrive in a later pass; this pass produces a watertight stepped
//      approximation that is still general 3-D and exact at the bulk).
//
// Output: standard polyhedral Mesh, ready for the SIMPLE solver.
// =============================================================================
#pragma once

#include "MeshStorage.hpp"
#include "StlImporter.hpp"

#include <array>

namespace simall::meshing {

struct OctreeMeshOptions {
    int    maxDepth          = 6;
    int    minDepthGlobal    = 3;
    double targetEdgeLength  = 0.0;   // 0 = auto from bbox / 2^maxDepth
};

class OctreeMesher {
public:
    /// Generate volumetric hex mesh from a closed STL surface.
    void mesh(const StlSurface& surface,
              OctreeMeshOptions opt,
              Mesh& out);
};

}  // namespace simall::meshing
