// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/PeriodicBc.cpp
// =============================================================================
#include "solver/bc/PeriodicBc.hpp"

namespace simall::solver::bc
{

void PeriodicBc::initialize(const meshing::Mesh& mesh)
{
    pairs_ = build_periodic_pairs(mesh, {xform_});
}

std::size_t PeriodicBc::apply(BcContext& ctx)
{
    if (!ctx.mesh || !ctx.matrix || !ctx.rhs || pairs_.empty())
        return 0;
    const auto& F = ctx.mesh->faces();
    auto& A = *ctx.matrix;
    auto& b = *ctx.rhs;

    std::size_t hits = 0;
    return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
        const int tw = (f < pairs_.size()) ? pairs_[f] : -1;
        if (tw < 0)
            return;
        const meshing::CellId c1 = F.owner[f];
        const meshing::CellId c2 = F.owner[static_cast<std::size_t>(tw)];
        // Symmetric off-diagonal coupling: A_{c1,c1} += K, A_{c1,c2} -= K
        // (we cannot insert new sparsity here; use a penalty-Dirichlet on the
        // current field value at the twin cell, valid for SIMPLE/PISO outer
        // iterations).
        const auto& field = ctx.fields ? ctx.fields->find_scalar(ctx.variable) : nullptr;
        const double twinVal = (field && c2 < field->size()) ? (*field)[c2] : 0.0;
        applyDirichlet(A, b, c1, twinVal);
        ++hits;
        (void) hits;
    });
}

} // namespace simall::solver::bc
