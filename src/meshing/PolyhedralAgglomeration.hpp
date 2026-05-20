// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/PolyhedralAgglomeration.hpp
// Phase  : 6.10 — Convert tetrahedral mesh to a polyhedral mesh via dual-
//          cell agglomeration (Ansys Fluent / STAR-CCM+ approach).
//
// Algorithm (Mathur-Murthy / Mavriplis dual mesh):
//   1. For each interior node of the input tet mesh, gather all tets that
//      share that node.
//   2. The polyhedral cell associated with the node is the union of the
//      barycentre-bounded sub-volumes (so-called dual cell or "median
//      dual").
//   3. Polyhedral cell faces lie on:
//        – edge mid-segments (tet edge midpoint <-> incident tet centroid)
//        – tet face barycentres
//   4. Resulting mesh is conforming, valid, and typically reduces cell
//      count by ~3-5× while improving solver convergence (Mavriplis,
//      AIAA J. 1994).
//
// In this scaffold we produce the *agglomeration map*: a list mapping
// each input tet to its target polyhedral cell.  Construction of explicit
// polyhedral face geometry is delegated to a host meshing driver.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"

#include <cstdint>
#include <vector>

namespace simall::meshing
{

struct AgglomerationProps
{
    std::size_t maxClusterSize = 8; // hard cap on tets per polyhedral
    double qualityThreshold = 0.2;
    bool seedByNode = true; // false ⇒ K-way agglom by metric
};

class PolyhedralAgglomeration
{
public:
    void initialize(AgglomerationProps props);

    /// Compute a tet→cluster map and write into `clusterId` (size = nCells).
    /// Returns the number of polyhedral clusters generated.
    std::size_t agglomerate(const Mesh& tetMesh, std::vector<std::int32_t>& clusterId);

    /// Builds a coarse mesh from the cluster map: each cluster becomes one
    /// cell, its faces are the boundary faces of the original cluster, and
    /// internal faces of the cluster are dropped.
    std::size_t emit_coarse_mesh(const Mesh& tetMesh,
                                 const std::vector<std::int32_t>& clusterId,
                                 Mesh& outMesh);

    const AgglomerationProps& props() const noexcept { return p_; }

private:
    AgglomerationProps p_{};
};

} // namespace simall::meshing
