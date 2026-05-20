// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/PorousJumpBc.cpp
// =============================================================================
#include "solver/bc/PorousJumpBc.hpp"

#include <cmath>

namespace simall::solver::bc
{

std::size_t PorousJumpBc::apply(BcContext& ctx)
{
    if (!ctx.mesh || !ctx.matrix || !ctx.rhs)
        return 0;
    const auto& F = ctx.mesh->faces();
    auto& A = *ctx.matrix;
    auto& b = *ctx.rhs;
    const std::string& v = ctx.variable;
    if (v != "p" && v != "p_rgh")
        return 0;

    auto* rho = ctx.fields ? ctx.fields->find_scalar("rho") : nullptr;
    auto* mu = ctx.fields ? ctx.fields->find_scalar("mu") : nullptr;
    auto* Ux = ctx.fields ? ctx.fields->find_vector("U") : nullptr;

    return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
        const meshing::CellId c = F.owner[f];
        const double r = (rho && c < rho->size()) ? (*rho)[c] : 1.225;
        const double mv = (mu && c < mu->size()) ? (*mu)[c] : 1.8e-5;
        double Umag = 0.0;
        if (Ux && c < Ux->size()) {
            const double ux = Ux->x[c], uy = Ux->y[c], uz = Ux->z[c];
            Umag = std::sqrt(ux * ux + uy * uy + uz * uz);
        }
        const double dp =
            (mv / p_.permeability + 0.5 * p_.inertialCoefficient * r * Umag) * Umag * p_.thickness;
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
        const double area = std::sqrt(Ax * Ax + Ay * Ay + Az * Az);
        // Add a pressure-jump source: increases RHS on the upstream cell.
        b[c] -= dp * area;
    });
}

} // namespace simall::solver::bc
