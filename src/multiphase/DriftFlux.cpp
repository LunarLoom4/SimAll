// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/DriftFlux.cpp
// =============================================================================
#include "multiphase/DriftFlux.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase
{

namespace
{
double harmathy(double sigma, double g, double rhoL, double rhoG)
{
    const double drho = std::max(rhoL - rhoG, 1e-12);
    return 1.41 * std::pow(sigma * g * drho / (rhoL * rhoL), 0.25);
}
double ishii_chawla(double sigma, double g, double rhoL, double rhoG, double alphaG)
{
    const double drho = std::max(rhoL - rhoG, 1e-12);
    const double Ugj0 = 1.41 * std::pow(sigma * g * drho / (rhoL * rhoL), 0.25);
    // (1-α_g)^1.75 swarm correction
    return Ugj0 * std::pow(std::max(1.0 - alphaG, 0.0), 1.75);
}
} // namespace

void DriftFlux::initialize(const meshing::Mesh& m, DriftFluxProps props)
{
    mesh_ = &m;
    p_ = props;
    SIMALL_LOG_INFO("Multiphase",
                    "DriftFlux initialised: C0=",
                    p_.C0,
                    " ρ_l=",
                    p_.rhoLiquid,
                    " ρ_g=",
                    p_.rhoGas,
                    " σ=",
                    p_.surfaceTension,
                    " g=",
                    p_.gravity,
                    " (",
                    m.cells().size(),
                    " cells)");
}

double DriftFlux::apply(solver::FieldRegistry& F)
{
    if (!mesh_)
        return 0.0;
    const std::size_t nC = mesh_->cells().size();
    const auto* Umix = F.find_vector("U_mix");
    const auto* aGas = F.find_scalar("alpha_g");
    if (!Umix || !aGas)
        return 0.0;

    auto& Ugj = F.vector("U_drift_g", nC);
    auto& Ug = F.vector("U_g", nC);
    auto& Ul = F.vector("U_l", nC);

    double maxDrift = 0.0;
    for (std::size_t c = 0; c < nC; ++c) {
        const double aG = std::clamp((*aGas)[c], 0.0, 1.0);
        const double aL = 1.0 - aG;

        // Drift magnitude.
        double Ugj_mag = 0.0;
        switch (p_.driftLaw) {
        case DriftFluxProps::DriftLaw::Harmathy:
            Ugj_mag = harmathy(p_.surfaceTension, p_.gravity, p_.rhoLiquid, p_.rhoGas);
            break;
        case DriftFluxProps::DriftLaw::IshiiChawla:
            Ugj_mag = ishii_chawla(p_.surfaceTension, p_.gravity, p_.rhoLiquid, p_.rhoGas, aG);
            break;
        case DriftFluxProps::DriftLaw::UserConstant:
            Ugj_mag = p_.Ugj_const;
            break;
        }
        // Drift direction = -g/|g| (gas rises against gravity).
        // For now use +z = vertical convention.
        const double dirx = 0.0, diry = 0.0, dirz = 1.0;
        Ugj.x[c] = Ugj_mag * dirx;
        Ugj.y[c] = Ugj_mag * diry;
        Ugj.z[c] = Ugj_mag * dirz;

        // Zuber-Findlay reconstruction:
        //   U_g = C_0 j + U_gj,      j ≈ U_mix
        //   U_l = (j - α_g U_g) / α_l
        const double jx = Umix->x[c], jy = Umix->y[c], jz = Umix->z[c];
        Ug.x[c] = p_.C0 * jx + Ugj.x[c];
        Ug.y[c] = p_.C0 * jy + Ugj.y[c];
        Ug.z[c] = p_.C0 * jz + Ugj.z[c];
        const double alEff = std::max(aL, 1e-12);
        Ul.x[c] = (jx - aG * Ug.x[c]) / alEff;
        Ul.y[c] = (jy - aG * Ug.y[c]) / alEff;
        Ul.z[c] = (jz - aG * Ug.z[c]) / alEff;

        if (Ugj_mag > maxDrift)
            maxDrift = Ugj_mag;
    }
    return maxDrift;
}

} // namespace simall::multiphase
