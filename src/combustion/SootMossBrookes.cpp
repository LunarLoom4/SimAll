// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/SootMossBrookes.cpp
// =============================================================================
#include "combustion/SootMossBrookes.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::combustion
{

namespace
{
constexpr double Ru = 8.314462618; // J/(mol·K)
constexpr double Avo = 6.02214076e23;
constexpr double PI = 3.14159265358979323846;
constexpr double C_mw = 12.011e-3; // kg/mol carbon
} // namespace

void SootMossBrookes::initialize(const meshing::Mesh& mesh, MossBrookesProps props)
{
    mesh_ = &mesh;
    p_ = props;
    const std::size_t nC = mesh.cells().size();
    SIMALL_LOG_INFO("Combustion",
                    "SootMossBrookes initialised: C_alpha=",
                    p_.C_alpha,
                    " rho_s=",
                    p_.rho_soot,
                    " (",
                    nC,
                    " cells)");
}

double SootMossBrookes::apply(solver::FieldRegistry& F)
{
    if (!mesh_)
        return 0.0;
    const std::size_t nC = mesh_->cells().size();
    const auto* T = F.find_scalar("T");
    if (!T)
        return 0.0;
    const auto* rhoF = F.find_scalar("rho_mix");
    const auto* Yf = F.find_scalar("Y_fuel");
    const auto* Yoh = F.find_scalar("Y_OH");
    const auto* N = F.find_scalar("N_soot");
    const auto* Ys = F.find_scalar("Y_soot");
    auto& Sn = F.scalar("S_N_soot", nC);
    auto& Sm = F.scalar("S_Y_soot", nC);

    const auto& Cc = mesh_->cells();
    double total = 0.0;
    for (std::size_t i = 0; i < nC; ++i) {
        const double Tc = std::max((*T)[i], 300.0);
        const double rho = rhoF ? std::max((*rhoF)[i], 1e-6) : 1.0;
        const double yf = Yf ? std::max((*Yf)[i], 0.0) : 0.0;
        const double yoh = Yoh ? std::max((*Yoh)[i], 0.0) : 0.0;
        const double Ni = N ? std::max((*N)[i], 0.0) : 0.0;
        const double Yi = Ys ? std::max((*Ys)[i], 0.0) : 0.0;

        // Concentrations (mol/m³).
        const double cFuel = rho * yf / std::max(p_.Mw_fuel, 1e-30);

        // Soot particle diameter and surface area density.
        const double sootMass = rho * Yi; // kg/m³
        const double d_p = (Ni > 1e-30 && sootMass > 1e-30)
                               ? std::cbrt(6.0 * sootMass / (PI * p_.rho_soot * Ni))
                               : 0.0;
        const double A_s = PI * d_p * d_p * Ni; // m²/m³

        // 1) Nucleation (Lindstedt-style fuel-pyrolysis).
        const double R_nuc =
            p_.C_alpha * Avo * cFuel * std::exp(-p_.T_alpha / Tc); // particles/(m³·s)

        // 2) Coagulation sink for N.
        const double R_coag = p_.C_beta * std::sqrt(24.0 * Ru * Tc / (p_.rho_soot * Avo))
                              * std::sqrt(std::max(d_p, 1e-15)) * Ni * Ni;

        // 3) Surface growth.
        const double R_sg = p_.C_gamma * std::exp(-p_.T_gamma / Tc)
                            * std::sqrt(std::max(cFuel, 0.0)) * A_s; // kg/(m³·s)

        // 4) OH-oxidation (Neoh 1981 collision efficiency η).
        const double P_OH = rho * yoh / 17.008e-3 * Ru * Tc;          // partial pressure
        const double R_oxid = p_.C_oxid * A_s * P_OH / std::sqrt(Tc); // kg/(m³·s)

        // Source assembly.
        const double Sni = R_nuc - R_coag;                         // particles/(m³·s)
        const double Smi = (C_mw / Avo) * p_.C_atoms_per_p * R_nuc // nucleation mass
                           + R_sg - R_oxid;                        // kg/(m³·s)
        Sn[i] += Sni;
        Sm[i] += Smi;
        total += Smi * Cc.volume[i];
    }
    return total;
}

} // namespace simall::combustion
