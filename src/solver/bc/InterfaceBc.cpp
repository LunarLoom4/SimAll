// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/InterfaceBc.cpp
// =============================================================================
#include "solver/bc/InterfaceBc.hpp"

#include <cmath>

namespace simall::solver::bc {

std::size_t InterfaceBc::apply(BcContext& ctx) {
    if (!ctx.mesh || !ctx.matrix || !ctx.rhs || !trace_) return 0;
    const auto& F = ctx.mesh->faces();
    auto& A = *ctx.matrix;
    auto& b = *ctx.rhs;
    const std::string& v = ctx.variable;

    if (v != "T") return 0;

    std::size_t localIdx = 0;
    return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
        const meshing::CellId c = F.owner[f];
        const double T_other = (localIdx < trace_->T_other.size())
                             ? trace_->T_other[localIdx] : 0.0;
        const double q_other = (localIdx < trace_->q_other.size())
                             ? trace_->q_other[localIdx] : 0.0;
        ++localIdx;
        // Dirichlet on the matched temperature; q is enforced by the other side.
        applyDirichlet(A, b, c, T_other);
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
        const double area = std::sqrt(Ax*Ax + Ay*Ay + Az*Az);
        b[c] += q_other * area;
    });
}

}  // namespace simall::solver::bc
