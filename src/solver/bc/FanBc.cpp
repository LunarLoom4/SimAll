// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/FanBc.cpp
// =============================================================================
#include "solver/bc/FanBc.hpp"

#include <cmath>

namespace simall::solver::bc
{

std::size_t FanBc::apply(BcContext& ctx)
{
    if (!ctx.mesh || !ctx.matrix || !ctx.rhs)
        return 0;
    const auto& F = ctx.mesh->faces();
    auto& b = *ctx.rhs;
    const std::string& v = ctx.variable;
    if (v != "p" && v != "p_rgh")
        return 0;

    auto* Ux = ctx.fields ? ctx.fields->find_vector("U") : nullptr;

    return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
        const meshing::CellId c = F.owner[f];
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
        const double area = std::sqrt(Ax * Ax + Ay * Ay + Az * Az);
        double Q = 0.0;
        if (Ux && c < Ux->size())
            Q = std::abs(Ux->x[c] * Ax + Ux->y[c] * Ay + Ux->z[c] * Az);
        const double dp = curve_.evaluate(Q);
        b[c] += dp * area;
    });
}

} // namespace simall::solver::bc
