// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/SymmetryBc.cpp
// =============================================================================
#include "solver/bc/SymmetryBc.hpp"

#include <cmath>

namespace simall::solver::bc {

std::size_t SymmetryBc::apply(BcContext& ctx) {
    if (!ctx.mesh || !ctx.matrix || !ctx.rhs) return 0;
    const auto& F = ctx.mesh->faces();
    auto& A = *ctx.matrix;
    auto& b = *ctx.rhs;
    const std::string& v = ctx.variable;

    // Scalars: homogeneous Neumann — no contribution.
    if (v == "p" || v == "T" || v == "k" || v == "omega" || v == "epsilon")
        return 0;

    // Velocity component: project the cell velocity onto the face normal
    // and pin it to zero.  This is a simplified penalty form; the full
    // implementation is delegated to the projection step in solver code.
    if (v.size() >= 3 && v[0] == 'U' && v[1] == '.') {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
            const double mag = std::sqrt(Ax*Ax + Ay*Ay + Az*Az);
            if (mag <= 0.0) return;
            const double nx = Ax/mag, ny = Ay/mag, nz = Az/mag;
            // Penalty Dirichlet on the projected component is approximated
            // by partially raising the diagonal proportional to n_axis^2.
            double w = 0.0;
            switch (v[2]) {
                case 'x': case 'X': w = nx * nx; break;
                case 'y': case 'Y': w = ny * ny; break;
                case 'z': case 'Z': w = nz * nz; break;
                default: return;
            }
            const meshing::CellId c = F.owner[f];
            addToDiagonal(A, c, kPenalty * w);
            b[c] += 0.0;
        });
    }
    return 0;
}

}  // namespace simall::solver::bc
