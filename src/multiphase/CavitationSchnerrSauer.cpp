// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/CavitationSchnerrSauer.cpp
// =============================================================================
#include "multiphase/CavitationSchnerrSauer.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase
{

void CavitationSchnerrSauer::initialize(const meshing::Mesh& m, CavitationProps p)
{
    mesh_ = &m;
    p_ = p;
    const std::size_t nC = m.cells().size();
    SIMALL_LOG_INFO("Cavitation",
                    "Schnerr-Sauer initialized: p_sat=",
                    p_.pSaturation,
                    " ρ_l=",
                    p_.rhoLiquid,
                    " ρ_v=",
                    p_.rhoVapor,
                    " n0=",
                    p_.nucleiDensity,
                    " (",
                    nC,
                    " cells)");
}

double CavitationSchnerrSauer::apply(solver::FieldRegistry& F)
{
    if (!mesh_)
        return 0.0;
    const std::size_t nC = mesh_->cells().size();
    const auto* pField = F.find_scalar("p");
    const auto* aField = F.find_scalar("alpha_v");
    auto& Sa = F.scalar("S_alpha_v", nC);
    auto& Sm = F.scalar("S_mass", nC);
    std::fill(Sa.begin(), Sa.end(), 0.0);
    std::fill(Sm.begin(), Sm.end(), 0.0);
    if (!pField || !aField)
        return 0.0;

    const auto& Cc = mesh_->cells();
    const double rhoL = p_.rhoLiquid;
    const double rhoV = p_.rhoVapor;
    const double pV = p_.pSaturation;
    const double n0 = p_.nucleiDensity;
    const double Cevap = p_.evapCoeff;
    const double Ccond = p_.condCoeff;
    const double third = 1.0 / 3.0;
    const double prefV = (4.0 / 3.0) * M_PI * n0;

    double total = 0.0;
    for (std::size_t c = 0; c < nC; ++c) {
        double av = std::clamp((*aField)[c], 0.0, 1.0);
        // Avoid (1-α_v)=0 at fully vaporised cell — clamp away from extremes.
        const double avEff = std::clamp(av, 1e-12, 1.0 - 1e-12);
        // Bubble radius from local void fraction:
        // R_B = ( α_v / ( (1-α_v) · (4/3) π n₀ ) )^{1/3}
        const double R_B = std::cbrt(avEff / ((1.0 - avEff) * prefV));
        const double rhoM = avEff * rhoV + (1.0 - avEff) * rhoL;
        const double dP = (*pField)[c] - pV;
        const double mdot_coef = (3.0 * rhoV * rhoL / std::max(rhoM, 1e-30))
                                 * (avEff * (1.0 - avEff)) / std::max(R_B, 1e-30);
        double mdot = 0.0;
        if (dP < 0.0) {
            // Vaporisation (mass from liquid → vapor) : mdot > 0 on vapor side
            mdot = +Cevap * mdot_coef * std::sqrt(-(2.0 / 3.0) * dP / rhoL);
        } else {
            // Condensation (mass from vapor → liquid) : mdot < 0 on vapor side
            mdot = -Ccond * mdot_coef * std::sqrt((2.0 / 3.0) * dP / rhoL);
        }
        Sa[c] = mdot / std::max(rhoV, 1e-30); // [1/s]
        Sm[c] = mdot;                         // [kg/(m³·s)]
        total += mdot * Cc.volume[c];
    }
    return total;
}

} // namespace simall::multiphase
