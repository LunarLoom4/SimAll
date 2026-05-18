// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/GhostExchange.hpp
// Phase  : 17.2 — halo-cell exchange between MPI ranks.
//
// After DomainPartition splits a global mesh into rank-local subdomains,
// each rank's flux/gradient assembly needs read-access to the field
// values of cells that live on neighbour ranks (the so-called "halo" or
// "ghost" layer).  GhostExchange:
//
//   1. Holds a per-neighbour-rank send list (local owned cells whose
//      values are needed remotely) and a per-neighbour-rank recv list
//      (ghost cell indices in the local array that receive remote data).
//   2. Implements pack → Isend/Irecv → Wait → unpack for scalar and
//      3-component vector fields.
//
// Cell-array layout convention: owned cells occupy [0, nOwned); ghost
// cells occupy [nOwned, nOwned + nGhost).  This keeps the inner
// assembly loops contiguous over owned cells.
// =============================================================================
#pragma once

#include "parallel/IFieldSynchronizer.hpp"
#include "parallel/MpiContext.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <cstdint>
#include <vector>

namespace simall::parallel {

struct NeighbourComm {
    int rank;                            // remote rank index
    std::vector<std::int32_t> sendCells; // local owned cell ids to send
    std::vector<std::int32_t> recvCells; // local ghost cell ids to receive into
};

class GhostExchange final : public IFieldSynchronizer {
public:
    explicit GhostExchange(MpiContext& ctx);

    void set_layout(std::size_t nOwned, std::size_t nGhost);
    void add_neighbour(NeighbourComm spec);
    void clear_neighbours();

    std::size_t neighbour_count() const noexcept { return nbrs_.size(); }
    const std::vector<NeighbourComm>& neighbours() const noexcept { return nbrs_; }

    // ------------- IFieldSynchronizer ----------------
    void sync_scalar(util::aligned_vector<double>& v) override;
    void sync_vector(util::aligned_vector<double>& cx,
                     util::aligned_vector<double>& cy,
                     util::aligned_vector<double>& cz) override;
    std::size_t owned_count() const override { return nOwned_; }
    std::size_t ghost_count() const override { return nGhost_; }

private:
    MpiContext&                  ctx_;
    std::size_t                  nOwned_ = 0;
    std::size_t                  nGhost_ = 0;
    std::vector<NeighbourComm>   nbrs_;
};

}  // namespace simall::parallel
