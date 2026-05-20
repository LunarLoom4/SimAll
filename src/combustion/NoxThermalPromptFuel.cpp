// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/NoxThermalPromptFuel.cpp
// =============================================================================
#include "combustion/NoxThermalPromptFuel.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::combustion
{

namespace
{
constexpr double Ru = 8.314462618; // J/(mol·K)
constexpr double Avo = 6.02214076e23;
}

void NoxThermalPromptFuel::initialize(const meshing::Mesh& mesh, NoxProps props)
{
    mesh_ = &mesh;
    p_ = props;
    const std::size_t nC = mesh.cells().size();
    SIMALL_LOG_INFO("Combustion",
                    "NOx model initialised: thermal=",
                    p_.includeThermal,
                    " prompt=",
                    p_.includePrompt,
                    " fuel=",
                    p_.includeFuel,
                    " (",
                    nC,
                    " cells)");
}

double NoxThermalPromptFuel::apply(solver::FieldRegistry& F)
{
    if (!mesh_)
        return 0.0;
    const std::size_t nC = mesh_->cells().size();
    const auto* T = F.find_scalar("T");
    if (!T)
        return 0.0;

    const auto* Y_N2 = F.find_scalar("Y_N2");
    const auto* Y_O2 = F.find_scalar("Y_O2");
    const auto* Y_OH = F.find_scalar("Y_OH");
    const auto* Y_fuel = F.find_scalar("Y_fuel");
    const auto* Y_fuelN = F.find_scalar("Y_fuelN");
    const auto* rhoF = F.find_scalar("rho_mix");
    auto& S_NO = F.scalar("S_NO", nC);

    double total = 0.0;
    const auto& Cc = mesh_->cells();
    for (std::size_t i = 0; i < nC; ++i) {
        const double Tc = std::max((*T)[i], 300.0);
        const double rho = rhoF ? std::max((*rhoF)[i], 1e-6) : p_.rho_default;
        const double yN2 = Y_N2 ? (*Y_N2)[i] : 0.767; // air default
        const double yO2 = Y_O2 ? (*Y_O2)[i] : 0.233;
        const double yOH = Y_OH ? (*Y_OH)[i] : 0.0;
        const double yF = Y_fuel ? (*Y_fuel)[i] : 0.0;
        const double yFN = Y_fuelN ? (*Y_fuelN)[i] : 0.0;

        // Concentrations [mol/m³] = ρ Y / Mw.
        const double cN2 = rho * yN2 / p_.Mw_N2;
        const double cO2 = rho * yO2 / p_.Mw_O2;
        const double cFN = rho * yFN / p_.Mw_fuelN;

        double rate_mol = 0.0; // d[NO]/dt in mol/(m³·s)

        // ---------------- Thermal (extended Zeldovich) -----------------------
        if (p_.includeThermal && Tc > 1500.0) {
            // Quasi-steady N atom: [N] ≈ (k1f [N2][O]) / (k2f [O2] + k3f [OH])
            // Estimate [O] from O2 partial-equilibrium  [O] = K_O √[O2] (T)
            const double K_O = std::exp(-31090.0 / Tc) * 3.6e5; // Hanson-Salimian fit
            const double cO = K_O * std::sqrt(std::max(cO2, 0.0));
            const double k1f = 1.8e8 * std::exp(-318.0 / (Ru * Tc / 1000.0));
            const double k2f = 1.8e4 * Tc * std::exp(-4680.0 / Tc);
            const double k3f = 7.1e7 * std::exp(-450.0 / Tc);
            const double cOH = rho * yOH / 17.008e-3;
            const double denom = k2f * cO2 + k3f * cOH;
            const double cN = (denom > 1e-30) ? (k1f * cN2 * cO) / denom : 0.0;
            const double r_th = 2.0 * k1f * cN2 * cO; // factor 2: two NO per N consumed
            (void) cN;
            rate_mol += r_th;
        }

        // ---------------- Prompt (De Soete 1975) ------------------------------
        if (p_.includePrompt && Tc > 1000.0 && yF > 0.0) {
            const double k_pr = 6.4e6 * std::exp(-36510.0 / Tc);
            const double r_pr = k_pr * cN2 * cO2;
            rate_mol += r_pr;
        }

        // ---------------- Fuel-N (Bose-Wendt) ---------------------------------
        if (p_.includeFuel && cFN > 0.0) {
            const double k_f = 1.0e10 * std::exp(-33700.0 / Tc);
            const double r_f = k_f * cFN * cO2;
            rate_mol += r_f;
        }

        const double mdot = rate_mol * p_.Mw_NO; // kg/(m³·s)
        S_NO[i] += mdot;
        total += mdot * Cc.volume[i];
    }
    return total;
}

} // namespace simall::combustion
