// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/OversetBc.cpp
// =============================================================================
#include "solver/bc/OversetBc.hpp"

namespace simall::solver::bc {

std::size_t OversetBc::apply(BcContext& ctx) {
    if (!ctx.matrix || !ctx.rhs) return 0;
    auto& A = *ctx.matrix;
    auto& b = *ctx.rhs;
    for (auto const& [cell, val] : donor_)
        applyDirichlet(A, b, cell, val);
    return donor_.size();
}

}  // namespace simall::solver::bc
