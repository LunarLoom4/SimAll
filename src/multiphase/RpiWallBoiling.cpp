// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/RpiWallBoiling.cpp
// =============================================================================
#include "multiphase/RpiWallBoiling.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace simall::multiphase
{

namespace
{
constexpr double PI = 3.14159265358979323846;
}

void RpiWallBoiling::initialize(const meshing::Mesh& m,
                                const std::vector<solver::BoundarySpec>& bcs,
                                RpiProps p)
{
    mesh_ = &m;
    bcs_ = bcs;
    p_ = p;
    SIMALL_LOG_INFO("Multiphase",
                    "RPI wall boiling initialised: T_sat=",
                    p_.T_sat,
                    " h_fg=",
                    p_.h_fg,
                    " (",
                    m.cells().size(),
                    " cells, ",
                    bcs.size(),
                    " BCs)");
}

double RpiWallBoiling::apply(solver::FieldRegistry& F)
{
    if (!mesh_)
        return 0.0;
    const std::size_t nC = mesh_->cells().size();
    const auto* Tfield = F.find_scalar("T");
    if (!Tfield)
        return 0.0;

    auto& Sa = F.scalar("S_alpha_v", nC);
    auto& Sm = F.scalar("S_mass", nC);
    auto& Se = F.scalar("S_energy", nC);
    // Additive — do not zero existing sources here.

    // Build zone → wall temperature lookup from BC list.
    std::unordered_map<meshing::ZoneId, double> Twall;
    for (const auto& b : bcs_) {
        if (b.type == solver::BCType::Wall || b.type == solver::BCType::NoSlipWall) {
            Twall[b.zone] = b.scalarValue; // BoundarySpec.scalarValue = T_wall
        }
    }
    if (Twall.empty())
        return 0.0;

    const auto& Cc = mesh_->cells();
    const auto& Fc = mesh_->faces();
    const std::size_t nF = Fc.size();

    double total_q = 0.0;
    for (std::size_t f = 0; f < nF; ++f) {
        const auto n = Fc.neighbor[f];
        if (n != meshing::kBoundaryCell)
            continue; // interior face
        const auto zone = Fc.boundaryZone[f];
        auto it = Twall.find(zone);
        if (it == Twall.end())
            continue;
        const double T_w = it->second;
        const auto o = Fc.owner[f];
        if (o == meshing::kBoundaryCell)
            continue;

        const double T_l = (*Tfield)[o];
        const double dT_sup = T_w - p_.T_sat; // wall superheat
        const double dT_sub = p_.T_sat - T_l; // liquid subcooling
        if (dT_sup <= 0.0)
            continue; // no nucleate boiling

        const double areaMag = std::sqrt(Fc.areaX[f] * Fc.areaX[f] + Fc.areaY[f] * Fc.areaY[f]
                                         + Fc.areaZ[f] * Fc.areaZ[f]);
        const double Vc = std::max(Cc.volume[o], 1e-30);

        // 1) Bubble departure diameter d_w (Tolubinsky-Kostanchuk).
        const double d_w = p_.d_ref * std::exp(-std::max(dT_sub, 0.0) / p_.dT_ref);
        // 2) Nucleation site density N_w (Lemmert-Chawla).
        const double N_w = std::pow(185.0 * dT_sup, 1.805);
        // 3) Bubble departure frequency f_w (Cole).
        const double f_w = std::sqrt(4.0 * p_.g * (p_.rho_liquid - p_.rho_vapor)
                                     / (3.0 * std::max(d_w, 1e-12) * p_.rho_liquid));
        // 4) Nucleate-boiling area fraction A_b (Del Valle-Kenning), capped.
        const double A_b = std::min(PI * d_w * d_w * N_w * 0.25, 1.0);

        // Heat-flux partitioning.
        const double q_conv = (1.0 - A_b) * p_.h_conv * (T_w - T_l);
        const double q_quench = 2.0 * A_b
                                * std::sqrt(p_.lambda_l * p_.rho_liquid * p_.cp_liquid * f_w / PI)
                                * (T_w - T_l);
        const double q_evap = (PI / 6.0) * std::pow(d_w, 3) * N_w * f_w * p_.rho_vapor * p_.h_fg;
        const double q_total = q_conv + q_quench + q_evap;

        // Mass-transfer source per unit volume from the evaporative flux.
        const double mdot_per_vol = q_evap * areaMag / (std::max(p_.h_fg, 1e-30) * Vc);
        Sa[o] += mdot_per_vol / std::max(p_.rho_vapor, 1e-30);
        Sm[o] += mdot_per_vol;
        // Energy sink in the liquid cell (latent + sensible to vapour).
        Se[o] -= (q_evap * areaMag) / Vc;
        total_q += q_total * areaMag;
    }
    return total_q;
}

} // namespace simall::multiphase
