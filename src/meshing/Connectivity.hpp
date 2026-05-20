// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/Connectivity.hpp
// Phase  : 5 — face-list construction & owner/neighbour assignment.
//
// Builds full FV-style topology from a cell-vertex list:
//   - unique faces (no duplicates across cell boundaries)
//   - owner / neighbour cell per face (boundary faces get kBoundaryCell)
//   - per-cell face index list (CSR)
//   - per-face node index list (CSR)
//   - geometric metrics (area vector, centroid, volume) computed afterwards
//
// The builder is hash-based: each face is keyed by its sorted vertex tuple.
// Complexity: O(C * f) where f = faces/cell (≤ 6 for hex, 4 for tet).
// =============================================================================
#pragma once

#include "MeshStorage.hpp"

#include <cstdint>
#include <vector>

namespace simall::meshing
{

/// Description of a single cell prior to global face stitching.
struct CellDescriptor
{
    /// One inner vector per face, listing vertex indices in CCW order
    /// (outward normal convention). For a hex this contains 6 vectors of 4.
    std::vector<std::vector<NodeId>> faces;
    ZoneId zone = 0;
};

/// Constructs full Mesh connectivity from per-cell face descriptions, then
/// invokes Mesh::compute_geometry() to populate areas / centroids / volumes.
/// `boundary_face_zones` maps `(sorted vertex tuple) → zone id` for any face
/// that lies on a named boundary; faces not present default to zone 1.
class ConnectivityBuilder
{
public:
    static void build(Mesh& out,
                      const NodeStorage& nodes,
                      const std::vector<CellDescriptor>& cells);
};

} // namespace simall::meshing
