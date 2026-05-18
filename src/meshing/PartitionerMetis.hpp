// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/PartitionerMetis.hpp
// Phase  : 6.13 — METIS-style multilevel k-way graph partitioner.
//
// Re-implements the Karypis-Kumar multilevel graph partitioning algorithm
// (METIS, 1998) for partitioning the cell-adjacency graph of a CFD mesh.
//
// The three phases of the multilevel scheme:
//   1. **Coarsening (heavy-edge matching):** repeatedly merge pairs of
//      vertices joined by the heaviest edge until the graph is small
//      enough (≤ 2·nparts vertices remain).
//   2. **Initial partitioning:** greedy growth-region partitioning on the
//      coarsest graph; computed deterministically by BFS from k seed
//      vertices spaced by graph diameter.
//   3. **Uncoarsening + refinement:** project the partition back through
//      the matching tree; at each level apply Kernighan-Lin / Fiduccia-
//      Mattheyses passes to reduce edge-cut.
//
// API mirrors the METIS_PartGraphKway signature:
//   partition(adjacency, nParts, partOut)  →  edge_cut
//
// Does NOT link against the external METIS library to keep the build
// dependency-free; the algorithm is a complete in-tree re-implementation.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"

#include <cstdint>
#include <vector>

namespace simall::meshing {

struct MetisProps {
    int    nParts          = 4;
    double imbalanceTol    = 1.03;
    int    refineSweeps    = 4;
    std::uint64_t rngSeed  = 0xCAFEDEED'1357ULL;
};

class PartitionerMetis {
public:
    void initialize(MetisProps props);

    /// Build the cell-adjacency graph from a Mesh and partition it.
    /// `partOut[c]` ∈ [0, nParts).  Returns the edge-cut count.
    std::size_t partition(const Mesh& mesh, std::vector<std::int32_t>& partOut);

    /// Partition an arbitrary CSR graph directly.  `xadj`/`adjncy` follow
    /// METIS conventions (xadj.size() = nVtx + 1).
    std::size_t partition(const std::vector<std::int32_t>& xadj,
                          const std::vector<std::int32_t>& adjncy,
                          std::vector<std::int32_t>& partOut);

    const MetisProps& props() const noexcept { return p_; }

private:
    MetisProps p_{};
};

}  // namespace simall::meshing
