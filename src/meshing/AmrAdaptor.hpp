// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/AmrAdaptor.hpp
// Phase  : 6.15 — Adaptive Mesh Refinement driver (tag-based).
//
// Reads per-cell refinement tags (e.g. error indicators) from a
// FieldRegistry scalar field and refines or coarsens cells subject to a
// 2:1 face balance constraint.  Hex cells are split 1→8 (octree), tet
// cells 1→8 (octasection à la Bey 1995), prism cells 1→8 (Bey-Schwarz
// 1991).  Coarsening recombines 8 sibling cells when all 8 carry a
// negative tag and have the same parent.
//
// Parentage is stored in a side-table indexed by CellId so the Mesh
// structure itself remains untouched between adaptation cycles.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::meshing {

struct AmrProps {
    int    maxLevel          = 4;
    double refineThreshold   = 0.7;
    double coarsenThreshold  = 0.2;
    bool   balance21         = true;
};

struct AmrReport {
    std::size_t cellsRefined   = 0;
    std::size_t cellsCoarsened = 0;
    std::size_t cellsBalanced  = 0;
};

class AmrAdaptor {
public:
    void initialize(Mesh& mesh, AmrProps props);

    /// Run one adaptation cycle.  `tagFieldName` is the name of a scalar
    /// field in `fields` whose magnitude is compared against refine /
    /// coarsen thresholds.
    AmrReport adapt(solver::FieldRegistry& fields,
                    const std::string& tagFieldName);

    int      cell_level(CellId c) const;
    const AmrProps& props() const noexcept { return p_; }

private:
    void enforce_balance(std::vector<int>& flag, AmrReport& rep);

    Mesh*                                       mesh_ = nullptr;
    AmrProps                                    p_{};
    std::unordered_map<CellId, int>             level_;     // per-cell refinement level
    std::unordered_map<CellId, CellId>          parent_;    // child → parent
};

}  // namespace simall::meshing
