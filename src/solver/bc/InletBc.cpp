// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/InletBc.cpp
// =============================================================================
#include "solver/bc/InletBc.hpp"

#include <cmath>

namespace simall::solver::bc {

namespace {

double velocityComponent(const std::string& v, const double U[3]) {
    if (v.size() >= 3 && v[0] == 'U' && v[1] == '.') {
        switch (v[2]) {
            case 'x': case 'X': return U[0];
            case 'y': case 'Y': return U[1];
            case 'z': case 'Z': return U[2];
        }
    }
    return 0.0;
}

double effectiveDensity(const BcContext& ctx, meshing::CellId c) {
    auto* rho = ctx.fields ? ctx.fields->find_scalar("rho") : nullptr;
    return (rho && c < rho->size()) ? (*rho)[c] : 1.225;
}

}  // anonymous

std::size_t InletBc::apply(BcContext& ctx) {
    if (!ctx.mesh || !ctx.matrix || !ctx.rhs) return 0;
    const auto& F = ctx.mesh->faces();
    auto& A = *ctx.matrix;
    auto& b = *ctx.rhs;
    const std::string& v = ctx.variable;

    if (v.size() >= 3 && v[0] == 'U' && v[1] == '.') {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            const meshing::CellId c = F.owner[f];
            double Uimp = 0.0;
            switch (p_.kind) {
                case InletKind::Velocity:
                    Uimp = velocityComponent(v, p_.velocity); break;
                case InletKind::MassFlow: {
                    const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
                    const double area = std::sqrt(Ax*Ax + Ay*Ay + Az*Az);
                    const double rho  = effectiveDensity(ctx, c);
                    if (area > 0 && rho > 0) {
                        const double Umag = p_.massFlowRate / (rho * area);
                        double U[3] = { Umag*p_.direction[0],
                                        Umag*p_.direction[1],
                                        Umag*p_.direction[2] };
                        Uimp = velocityComponent(v, U);
                    }
                    break;
                }
                case InletKind::PressureTotal:
                    // velocity comes from p_total - 0.5 rho U^2 = p_static.
                    // Solver-side coupling computes magnitude; here we leave
                    // velocity free.
                    return;
            }
            applyDirichlet(A, b, c, Uimp);
        });
    }

    if (v == "p") {
        if (p_.kind == InletKind::PressureTotal) {
            return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
                applyDirichlet(A, b, F.owner[f], p_.totalPressure);
            });
        }
        return 0;
    }
    if (v == "T") {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            applyDirichlet(A, b, F.owner[f], p_.staticTemperature);
        });
    }
    if (v == "k") {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            applyDirichlet(A, b, F.owner[f], p_.turbulence_k);
        });
    }
    if (v == "omega") {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            applyDirichlet(A, b, F.owner[f], p_.turbulence_omega);
        });
    }
    return 0;
}

}  // namespace simall::solver::bc
