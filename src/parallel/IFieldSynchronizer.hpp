// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/IFieldSynchronizer.hpp
// Phase  : 17 — abstract halo-update interface.
//
// Pure abstract base shared between the solver (which has a non-owning
// pointer to one) and the parallel subsystem (which provides the
// concrete GhostExchange implementation).  Decoupling via this interface
// avoids a Solver → Parallel circular include chain.
// =============================================================================
#pragma once

#include "utilities/AlignedAllocator.hpp"

#include <cstddef>

namespace simall::parallel {

class IFieldSynchronizer {
public:
    virtual ~IFieldSynchronizer() = default;

    /// Exchange ghost values of a scalar cell-centred field.
    virtual void sync_scalar(util::aligned_vector<double>& cellData) = 0;

    /// Exchange ghost values of a three-component vector cell-centred field.
    virtual void sync_vector(util::aligned_vector<double>& cx,
                             util::aligned_vector<double>& cy,
                             util::aligned_vector<double>& cz) = 0;

    /// Number of local owned cells (sizes leading partition of the array).
    virtual std::size_t owned_count() const = 0;

    /// Number of ghost cells appended after the owned block.
    virtual std::size_t ghost_count() const = 0;
};

}  // namespace simall::parallel
