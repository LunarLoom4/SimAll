// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/DomainPartition.hpp
// Phase  : 17.1 — global-mesh → per-rank subdomain partitioning.
//
// Glues the W11 graph partitioners (METIS / Scotch) to the W12 ghost
// machinery.  Workflow:
//
//   1. Compute cell-to-rank map from the cell-adjacency graph of the
//      global mesh.
//   2. Tag local-owned cells (rank == myRank) and gather ghost cells
//      from neighbour-rank cells across each cut face.
//   3. Build NeighbourComm specs suitable for GhostExchange::add_neighbour.
//
// The local subdomain mesh extraction itself is large and would touch
// every subsystem; here we expose the *partition map* and the *ghost
// communication plan* (the two artefacts the rest of the solver
// actually needs) without rewriting MeshStorage.
// =============================================================================
#pragma once

#include "meshing/MeshOps.hpp"
#include "meshing/MeshStorage.hpp"
#include "parallel/GhostExchange.hpp"
#include "parallel/MpiContext.hpp"

#include <cstdint>
#include <vector>

namespace simall::parallel
{

enum class PartitionerKind
{
    Metis,
    Scotch
};

struct DomainPartitionProps
{
    PartitionerKind kind = PartitionerKind::Metis;
    int nParts = 0; // 0 → use ctx.size()
    double imbalanceTol = 1.05;
};

struct DomainPlan
{
    std::vector<std::int32_t> cellRank;    // size = nGlobalCells
    std::vector<std::int32_t> localCells;  // global cell ids owned by myRank
    std::vector<std::int32_t> ghostCells;  // global cell ids in halo
    std::vector<NeighbourComm> neighbours; // communication plan
    std::size_t edgeCut = 0;
};

class DomainPartition
{
public:
    explicit DomainPartition(MpiContext& ctx);

    void initialize(DomainPartitionProps props);

    /// Partition the global mesh and build the rank-local domain plan.
    /// Returns the number of local owned cells.
    std::size_t partition(const meshing::Mesh& globalMesh, DomainPlan& outPlan);

    /// Convenience: configure a GhostExchange from the produced plan.
    void install_into(GhostExchange& gx, const DomainPlan& plan) const;

    /// Materialise the rank-local subdomain mesh (owned cells + one-deep
    /// ghost halo) from a previously built `plan`.  Thin wrapper around
    /// `meshing::ops::extract_subdomain` using `ctx.rank()` as the local
    /// rank.  Returns the standard SubdomainStats reporting owned/ghost
    /// counts plus the wrapped split_mesh result.
    meshing::ops::SubdomainStats extract_local_mesh(const meshing::Mesh& globalMesh,
                                                    const DomainPlan& plan,
                                                    meshing::Mesh& outMesh) const;

    const DomainPartitionProps& props() const noexcept { return p_; }

private:
    MpiContext& ctx_;
    DomainPartitionProps p_{};
};

} // namespace simall::parallel
