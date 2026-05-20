// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ParticleEvaporation.cpp
// =============================================================================
#include "particles/ParticleEvaporation.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::particles
{

namespace
{
constexpr double PI = 3.14159265358979323846;
constexpr double Ru = 8.314462618;
}

void ParticleEvaporation::initialize(const meshing::Mesh& mesh, EvaporationProps props)
{
    mesh_ = &mesh;
    p_ = props;
    SIMALL_LOG_INFO("Particles",
                    "Evaporation initialised (Tboil=",
                    p_.Tboil,
                    " L_vap=",
                    p_.L_vap,
                    " Abramzon=",
                    p_.abramzon,
                    ")");
}

double ParticleEvaporation::p_sat(double T) const
{
    // Clausius-Clapeyron with reference (Tboil, 1 atm).
    return 101325.0
           * std::exp(p_.L_vap * p_.Mw_vapor / Ru
                      * (1.0 / std::max(p_.Tboil, 1.0) - 1.0 / std::max(T, 1.0)));
}

double ParticleEvaporation::apply(double dt,
                                  LagrangianTracker& tracker,
                                  solver::FieldRegistry& F,
                                  double T_p)
{
    if (!mesh_)
        return 0.0;
    const std::size_t nC = mesh_->cells().size();
    const auto* Tg = F.find_scalar("T");
    const auto* Pg = F.find_scalar("p");
    const auto* Yv = F.find_scalar("Y_vapor");
    if (!Tg)
        return 0.0;
    auto& Sv = F.scalar("S_vapor", nC);
    auto& Se = F.scalar("S_energy", nC);

    auto& parts = tracker.mutable_particles();
    double m_evap_total = 0.0;
    const auto& Cc = mesh_->cells();
    for (auto& pt : parts) {
        if (!pt.active)
            continue;
        const auto c = pt.cell;
        if (c >= nC)
            continue;
        const double Tg_c = std::max((*Tg)[c], 280.0);
        const double Pg_c = Pg ? std::max((*Pg)[c], 1.0) : 101325.0;
        const double Yv_c = Yv ? std::clamp((*Yv)[c], 0.0, 0.99) : 0.0;

        // Droplet diameter from mass (assume single-component liquid).
        const double d_p = std::cbrt(6.0 * pt.mass / (PI * p_.rho_l));
        if (d_p < p_.d_min) {
            pt.active = false;
            continue;
        }

        // Vapour mass fraction at droplet surface (Raoult/Clausius-Clapeyron).
        const double ps = std::min(p_sat(T_p), 0.99 * Pg_c);
        const double Xs = ps / Pg_c; // mole fraction
        const double Ys = (Xs * p_.Mw_vapor) / (Xs * p_.Mw_vapor + (1 - Xs) * p_.Mw_gas);
        if (Ys <= Yv_c)
            continue; // no driving force
        const double BM = (Ys - Yv_c) / std::max(1.0 - Ys, 1e-6);

        // Reynolds, Schmidt, Sherwood (Ranz-Marshall + Abramzon-Sirignano).
        const double rho_g = p_.rho_g_ref * (293.15 / Tg_c); // approx ideal-gas
        const double Vrel =
            std::sqrt((pt.v.x) * (pt.v.x) + (pt.v.y) * (pt.v.y) + (pt.v.z) * (pt.v.z));
        const double Re = rho_g * Vrel * d_p / std::max(p_.mu_g, 1e-30);
        const double Sc = p_.mu_g / (rho_g * p_.Dv_ref);
        const double Sh0 =
            2.0 + 0.6 * std::sqrt(std::max(Re, 0.0)) * std::cbrt(std::max(Sc, 1e-30));
        double Sh = Sh0;
        if (p_.abramzon) {
            const double FM = std::pow(1.0 + BM, 0.7) * std::log(1.0 + BM) / std::max(BM, 1e-30);
            Sh = 2.0 + (Sh0 - 2.0) / std::max(FM, 1e-30);
        }

        // Evaporation mass flux from the droplet.
        const double mdot =
            PI * d_p * rho_g * p_.Dv_ref * Sh * std::log(1.0 + BM); // kg/s per parcel
        const double dm = std::min(mdot * dt, 0.95 * pt.mass);
        pt.mass -= dm;
        m_evap_total += dm;

        // Cell source contributions (parcel mass divided by cell volume).
        const double V = std::max(Cc.volume[c], 1e-30);
        Sv[c] += dm / (dt * V);
        Se[c] -= dm * p_.L_vap / (dt * V);

        // Disable parcel once below d_min after the mass update.
        const double d_new = std::cbrt(6.0 * pt.mass / (PI * p_.rho_l));
        if (d_new < p_.d_min)
            pt.active = false;
    }
    return m_evap_total;
}

} // namespace simall::particles
