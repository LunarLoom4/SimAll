// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/PartitionerScotch.hpp
// Phase  : 6.14 — Scotch-style recursive bisection partitioner.
//
// Implements the recursive bisection scheme championed by Pellegrini's
// Scotch (1996) library:  partition the graph into 2 parts, then recurse
// on each half until log₂(nParts) levels are produced.  Bisection uses
// a Fiduccia-Mattheyses min-cut sweep.
//
// Distinct from PartitionerMetis (which is k-way at the coarsest level
// followed by direct k-way refinement), Scotch's recursive scheme often
// produces lower-aspect-ratio sub-domains and is the historical default
// for finite-element codes.  Both partitioners share the FieldRegistry-
// free, dependency-free in-tree implementation.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"

#include <cstdint>
#include <vector>

namespace simall::meshing {

struct ScotchProps {
    int    nParts        = 4;
    int    fmPasses      = 4;
    double imbalanceTol  = 1.05;
    std::uint64_t rngSeed= 0x5C07'C0DEULL;
};

class PartitionerScotch {
public:
    void initialize(ScotchProps props);

    std::size_t partition(const Mesh& mesh, std::vector<std::int32_t>& partOut);

    std::size_t partition(const std::vector<std::int32_t>& xadj,
                          const std::vector<std::int32_t>& adjncy,
                          std::vector<std::int32_t>& partOut);

    const ScotchProps& props() const noexcept { return p_; }

private:
    void bisect(const std::vector<std::int32_t>& xadj,
                const std::vector<std::int32_t>& adjncy,
                const std::vector<std::int32_t>& vmap,    // sub-vertex → global
                std::vector<std::int32_t>& partOut,
                std::int32_t partA, std::int32_t partB,
                std::uint64_t seed);

    ScotchProps p_{};
};

}  // namespace simall::meshing
