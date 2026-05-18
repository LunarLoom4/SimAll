// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/AxisymmetricBc.cpp
// =============================================================================
#include "solver/bc/AxisymmetricBc.hpp"

#include <cctype>

namespace simall::solver::bc {

std::size_t AxisymmetricBc::apply(BcContext& ctx) {
    if (!ctx.mesh || !ctx.matrix || !ctx.rhs) return 0;
    const auto& F = ctx.mesh->faces();
    auto& A = *ctx.matrix;
    auto& b = *ctx.rhs;
    const std::string& v = ctx.variable;

    // Velocity radial-direction component is pinned to zero.
    // For axis=X, radial = (y, z); for axis=Y, radial = (x, z); etc.
    if (v.size() >= 3 && v[0] == 'U' && v[1] == '.') {
        const char comp = static_cast<char>(std::tolower(static_cast<unsigned char>(v[2])));
        bool isRadial = false;
        switch (p_.axis) {
            case AxisAlignment::X: isRadial = (comp == 'y' || comp == 'z'); break;
            case AxisAlignment::Y: isRadial = (comp == 'x' || comp == 'z'); break;
            case AxisAlignment::Z: isRadial = (comp == 'x' || comp == 'y'); break;
        }
        if (isRadial) {
            return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
                applyDirichlet(A, b, F.owner[f], 0.0);
            });
        }
        return 0;
    }
    // All scalars: zero gradient — no-op.
    return 0;
}

}  // namespace simall::solver::bc
