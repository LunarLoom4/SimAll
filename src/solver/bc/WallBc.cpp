// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/WallBc.cpp
// =============================================================================
#include "solver/bc/WallBc.hpp"

#include <cmath>

namespace simall::solver::bc
{

namespace
{

// Pick the U_wall component matching variable "U.x"/"U.y"/"U.z".
double wallComponent(const std::string& var, const double U[3])
{
    if (var.size() >= 3 && var[0] == 'U' && var[1] == '.') {
        switch (var[2]) {
        case 'x':
        case 'X':
            return U[0];
        case 'y':
        case 'Y':
            return U[1];
        case 'z':
        case 'Z':
            return U[2];
        }
    }
    return 0.0;
}

} // namespace

std::size_t WallBc::apply(BcContext& ctx)
{
    if (!ctx.mesh || !ctx.matrix || !ctx.rhs)
        return 0;
    const auto& F = ctx.mesh->faces();
    auto& A = *ctx.matrix;
    auto& b = *ctx.rhs;

    const std::string& v = ctx.variable;

    // ---- Momentum components: U.x / U.y / U.z ----------------------------
    if (v.size() >= 3 && v[0] == 'U' && v[1] == '.') {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            const meshing::CellId c = F.owner[f];
            switch (p_.velocityMode) {
            case WallVelocityMode::NoSlip:
                applyDirichlet(A, b, c, 0.0);
                break;
            case WallVelocityMode::FreeSlip:
                /* tangential preserved; no contribution to normal eqn */ break;
            case WallVelocityMode::MovingWall:
                applyDirichlet(A, b, c, wallComponent(v, p_.U_wall));
                break;
            }
        });
    }

    // ---- Pressure: zero-gradient (Neumann homogeneous) -------------------
    if (v == "p" || v == "p_rgh") {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId) { /* no-op */ });
    }

    // ---- Temperature ----------------------------------------------------
    if (v == "T") {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            const meshing::CellId c = F.owner[f];
            const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
            const double area = std::sqrt(Ax * Ax + Ay * Ay + Az * Az);
            switch (p_.thermalMode) {
            case WallThermalMode::Adiabatic:
                /* zero gradient -> no contribution */ break;
            case WallThermalMode::Isothermal:
                applyDirichlet(A, b, c, p_.T_wall);
                break;
            case WallThermalMode::HeatFlux:
                applyNeumann(b, c, p_.heatFlux, area);
                break;
            case WallThermalMode::Convection:
                applyRobin(A, b, c, p_.h_conv, area, p_.T_infinity);
                break;
            }
        });
    }

    // ---- Turbulence variables (handled via wall functions) ---------------
    if (v == "k") {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            // Standard wall function: dk/dn = 0 at the wall.
            (void) f;
        });
    }
    if (v == "omega" || v == "epsilon") {
        return for_each_face(*ctx.mesh, [&](meshing::FaceId f) {
            // omega_wall = 60 nu / (beta1 y1^2)  — set Dirichlet to large
            // value; concrete y1, nu provided by TurbulenceModel.
            const meshing::CellId c = F.owner[f];
            applyDirichlet(A, b, c, 1.0e3);
        });
    }
    return 0;
}

} // namespace simall::solver::bc
