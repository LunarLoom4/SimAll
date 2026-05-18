// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/CavitationKunz.cpp
// =============================================================================
#include "multiphase/CavitationKunz.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase {

void CavitationKunz::initialize(const meshing::Mesh& m, KunzProps p) {
    mesh_ = &m;
    p_    = p;
    SIMALL_LOG_INFO("Cavitation",
        "Kunz initialized: p_sat=", p_.pSaturation,
        " C_dest=", p_.C_dest, " C_prod=", p_.C_prod,
        " (", m.cells().size(), " cells)");
}

double CavitationKunz::apply(solver::FieldRegistry& F) {
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
    const double tInf = p_.Lref / std::max(p_.Uref, 1e-12);
    const double dynP = 0.5 * rhoL * p_.Uref * p_.Uref;
    const double invDynPt = 1.0 / std::max(dynP * tInf, 1e-30);

    double total = 0.0;
    for (std::size_t c = 0; c < nC; ++c) {
        const double av = std::clamp((*aField)[c], 0.0, 1.0);
        const double al = 1.0 - av;
        const double dP = (*pField)[c] - pV;
        // Vaporisation (mass into vapor when p<p_v): ṁ⁺ on vapor side
        const double mDot_minus = -p_.C_dest * rhoV * al * std::min(0.0, dP) * invDynPt;
        // Condensation (vapor → liquid): ṁ⁻ on vapor side (negative)
        const double mDot_plus  =  p_.C_prod * rhoV * av * av * al / tInf;
        const double mdot = mDot_minus - mDot_plus;  // net into vapor
        Sa[c] = mdot / std::max(rhoV, 1e-30);
        Sm[c] = mdot;
        total += mdot * Cc.volume[c];
    }
    return total;
}

}  // namespace simall::multiphase
