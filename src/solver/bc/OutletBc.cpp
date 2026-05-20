// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/OutletBc.cpp
// =============================================================================
#include "solver/bc/OutletBc.hpp"

namespace simall::solver::bc
{

std::size_t OutletBc::apply(BcContext& ctx)
{
    if (!ctx.mesh || !ctx.matrix || !ctx.rhs)
        return 0;
    const auto& F = ctx.mesh->faces();
    auto& A = *ctx.matrix;
    auto& b = *ctx.rhs;
    const std::string& v = ctx.variable;

    if (v == "p" && p_.kind == OutletKind::PressureStatic) {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            applyDirichlet(A, b, F.owner[f], p_.staticPressure);
        });
    }
    // Outflow: zero-gradient for everything; matrix unchanged.
    // The solver's velocity coupling enforces ∂U/∂n = 0 explicitly.
    if (v.size() >= 3 && v[0] == 'U' && v[1] == '.')
        return 0;
    if (v == "T")
        return 0;
    if (v == "k")
        return 0;
    if (v == "omega")
        return 0;
    return 0;
}

} // namespace simall::solver::bc
