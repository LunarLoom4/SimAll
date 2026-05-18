// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/CavitationZwart.cpp
// =============================================================================
#include "multiphase/CavitationZwart.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase {

void CavitationZwart::initialize(const meshing::Mesh& m, ZwartProps p) {
    mesh_ = &m;
    p_    = p;
    SIMALL_LOG_INFO("Cavitation",
        "Zwart-Gerber-Belamri initialized: p_sat=", p_.pSaturation,
        " r_nuc=", p_.r_nucleation, " R_B=", p_.R_bubble,
        " (", m.cells().size(), " cells)");
}

double CavitationZwart::apply(solver::FieldRegistry& F) {
    if (!mesh_) return 0.0;
    const std::size_t nC = mesh_->cells().size();
    const auto* pField = F.find_scalar("p");
    const auto* aField = F.find_scalar("alpha_v");
    auto& Sa = F.scalar("S_alpha_v", nC);
    auto& Sm = F.scalar("S_mass",    nC);
    std::fill(Sa.begin(), Sa.end(), 0.0);
    std::fill(Sm.begin(), Sm.end(), 0.0);
    if (!pField || !aField) return 0.0;

    const auto& Cc = mesh_->cells();
    const double rhoL = p_.rhoLiquid;
    const double rhoV = p_.rhoVapor;
    const double pV   = p_.pSaturation;
    const double R_B  = std::max(p_.R_bubble, 1e-12);
    const double r_n  = std::clamp(p_.r_nucleation, 0.0, 1.0);
    const double pref = 3.0 * rhoV / R_B;

    double total = 0.0;
    for (std::size_t c = 0; c < nC; ++c) {
        const double av = std::clamp((*aField)[c], 0.0, 1.0);
        const double dP = (*pField)[c] - pV;
        double mdot = 0.0;
        if (dP < 0.0) {
            // Vaporisation: m⁻ proportional to r_nuc (1-α_v)
            const double rate = std::sqrt((2.0/3.0) * (-dP) / rhoL);
            mdot = +p_.F_vap * pref * r_n * (1.0 - av) * rate;
        } else {
            // Condensation: m⁺ proportional to α_v
            const double rate = std::sqrt((2.0/3.0) * dP / rhoL);
            mdot = -p_.F_cond * pref * av * rate;
        }
        Sa[c] = mdot / std::max(rhoV, 1e-30);
        Sm[c] = mdot;
        total += mdot * Cc.volume[c];
    }
    return total;
}

}  // namespace simall::multiphase
