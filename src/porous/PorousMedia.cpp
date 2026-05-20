// =============================================================================
// SimAll Beta - Porous Media Subsystem
// File   : src/porous/PorousMedia.cpp
// =============================================================================
#include "porous/PorousMedia.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::porous
{

namespace
{

/// R = [axis1 | axis2 | axis3] columns; T = R diag(d) R^T
void build_aniso(const util::Vec3d& d,
                 const util::Vec3d& a1,
                 const util::Vec3d& a2,
                 const util::Vec3d& a3,
                 double T[9])
{
    const double R[9] = {a1.x, a2.x, a3.x, a1.y, a2.y, a3.y, a1.z, a2.z, a3.z};
    // M = R · diag(d)
    const double M[9] = {R[0] * d.x,
                         R[1] * d.y,
                         R[2] * d.z,
                         R[3] * d.x,
                         R[4] * d.y,
                         R[5] * d.z,
                         R[6] * d.x,
                         R[7] * d.y,
                         R[8] * d.z};
    // T = M · R^T
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0;
            for (int k = 0; k < 3; ++k)
                s += M[i * 3 + k] * R[j * 3 + k];
            T[i * 3 + j] = s;
        }
}

} // namespace

void PorousMedia::initialize(const meshing::Mesh& m, const std::vector<meshing::ZoneId>& cellZone)
{
    mesh_ = &m;
    const std::size_t nC = m.cells().size();
    D_.assign(nC * 9, 0.0);
    F_.assign(nC * 9, 0.0);
    active_.assign(nC, 0);
    for (const auto& z : zones_) {
        double Dz[9], Fz[9];
        build_aniso(z.dDiag, z.axis1, z.axis2, z.axis3, Dz);
        build_aniso(z.fDiag, z.axis1, z.axis2, z.axis3, Fz);
        for (std::size_t c = 0; c < nC; ++c) {
            if (cellZone.size() == nC && cellZone[c] == z.zone) {
                for (int k = 0; k < 9; ++k) {
                    D_[c * 9 + k] = Dz[k];
                    F_[c * 9 + k] = Fz[k];
                }
                active_[c] = 1;
            }
        }
    }
    std::size_t n = 0;
    for (auto a : active_)
        n += a;
    SIMALL_LOG_INFO("Porous", "initialized ", zones_.size(), " zones, ", n, " active cells");
}

void PorousMedia::apply(solver::FieldRegistry& F, double rho, double mu) const
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    auto& S = F.vector("S_porous", nC);
    const auto* U = F.find_vector("U");
    std::fill(S.x.begin(), S.x.end(), 0.0);
    std::fill(S.y.begin(), S.y.end(), 0.0);
    std::fill(S.z.begin(), S.z.end(), 0.0);
    if (!U)
        return;

    for (std::size_t c = 0; c < nC; ++c) {
        if (!active_[c])
            continue;
        const double ux = U->x[c], uy = U->y[c], uz = U->z[c];
        const double speed = std::sqrt(ux * ux + uy * uy + uz * uz);
        // M = (μ D + ½ ρ |U| F)
        const double* D = &D_[c * 9];
        const double* Fc = &F_[c * 9];
        const double a = mu;
        const double b = 0.5 * rho * speed;
        // s = -M · U
        S.x[c] = -((a * D[0] + b * Fc[0]) * ux + (a * D[1] + b * Fc[1]) * uy
                   + (a * D[2] + b * Fc[2]) * uz);
        S.y[c] = -((a * D[3] + b * Fc[3]) * ux + (a * D[4] + b * Fc[4]) * uy
                   + (a * D[5] + b * Fc[5]) * uz);
        S.z[c] = -((a * D[6] + b * Fc[6]) * ux + (a * D[7] + b * Fc[7]) * uy
                   + (a * D[8] + b * Fc[8]) * uz);
    }
}

} // namespace simall::porous
