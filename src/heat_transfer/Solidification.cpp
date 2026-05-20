// =============================================================================
// SimAll Beta - Heat-Transfer Subsystem
// File   : src/heat_transfer/Solidification.cpp
// =============================================================================
#include "heat_transfer/Solidification.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::heat_transfer
{

void Solidification::initialize(const meshing::Mesh& m,
                                solver::FieldRegistry& F,
                                SolidificationProps props)
{
    mesh_ = &m;
    p_ = props;
    const std::size_t nC = m.cells().size();
    F.scalar("f_l", nC);
    F.vector("S_DarcyMom", nC);
    F.scalar("S_LatentEn", nC);
    fl_prev_.assign(nC, 0.0);
    // Initialise f_l from T if available.
    if (auto* T = F.find_scalar("T")) {
        auto& fl = *F.find_scalar("f_l");
        for (std::size_t c = 0; c < nC; ++c) {
            const double t = (*T)[c];
            fl[c] = std::clamp(
                (t - p_.T_solidus) / std::max(1e-12, p_.T_liquidus - p_.T_solidus), 0.0, 1.0);
            fl_prev_[c] = fl[c];
        }
    }
}

void Solidification::update(double dt, solver::FieldRegistry& F)
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    const auto* T = F.find_scalar("T");
    const auto* U = F.find_vector("U");
    if (!T)
        return;
    auto& fl = *F.find_scalar("f_l");
    auto& Sd = *F.find_vector("S_DarcyMom");
    auto& SlH = *F.find_scalar("S_LatentEn");

    const double dT = std::max(1e-12, p_.T_liquidus - p_.T_solidus);
    const double invDt = (dt > 0) ? 1.0 / dt : 0.0;
    const double L_rho = p_.latentHeat * p_.rho;

    for (std::size_t c = 0; c < nC; ++c) {
        const double t = (*T)[c];
        const double f_new = std::clamp((t - p_.T_solidus) / dT, 0.0, 1.0);
        fl[c] = f_new;

        // Darcy momentum sink (Carman-Kozeny).
        const double f3 = f_new * f_new * f_new + p_.eps;
        const double Acoef = p_.A_mush * (1.0 - f_new) * (1.0 - f_new) / f3;
        if (U) {
            Sd.x[c] = -Acoef * U->x[c];
            Sd.y[c] = -Acoef * U->y[c];
            Sd.z[c] = -Acoef * U->z[c];
        } else {
            Sd.x[c] = Sd.y[c] = Sd.z[c] = 0.0;
        }

        // Latent heat source S_h = -ρ L ∂f_l/∂t (released on freezing).
        const double dfdt = (f_new - fl_prev_[c]) * invDt;
        SlH[c] = -L_rho * dfdt;
        fl_prev_[c] = f_new;
    }
}

} // namespace simall::heat_transfer
