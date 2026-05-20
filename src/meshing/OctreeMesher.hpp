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

namespace simall::meshing
{

struct OctreeMeshOptions
{
    int maxDepth = 6;
    int minDepthGlobal = 3;
    double targetEdgeLength = 0.0; // 0 = auto from bbox / 2^maxDepth

    /// (Pass 13) True cut-cell extraction.  When false (legacy default),
    /// every boundary leaf is emitted as a full hex (stepped staircase
    /// approximation).  When true, each boundary leaf is clipped against
    /// the STL surface by:
    ///   1. Classifying its 8 corners (inside/outside) using the same
    ///      ray-cast test as the existing point-in-surface helper.
    ///   2. Locating the surface intersection on every hex edge with a
    ///      sign change via bisection (controlled by `edgeBisectIters`).
    ///   3. Per hex face, clipping the 4-vertex cycle to the inside half
    ///      (marching-squares).  Each cut face contributes a closed loop
    ///      segment that joins the cap polygon.
    ///   4. Emitting one general polyhedron per cut leaf whose faces are
    ///      the 6 clipped face polygons (zero-or-more vertices) plus the
    ///      cap polygon glued around the cut.
    /// All-outside leaves are dropped; all-inside boundary leaves degrade
    /// to a plain hex (no extra cost).
    bool enableCutCells = false;
    int edgeBisectIters = 12; ///< Edge-intersection precision (>=4).
};

class OctreeMesher
{
public:
    /// Generate volumetric hex mesh from a closed STL surface.
    void mesh(const StlSurface& surface, OctreeMeshOptions opt, Mesh& out);
};

} // namespace simall::meshing
