// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ParticleHeating.cpp
// =============================================================================
#include "particles/ParticleHeating.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::particles {

namespace {
constexpr double PI    = 3.14159265358979323846;
constexpr double SIGMA = 5.670374419e-8;          // Stefan-Boltzmann
}

void ParticleHeating::initialize(const meshing::Mesh& mesh, HeatingProps props) {
    mesh_ = &mesh;
    p_    = props;
    SIMALL_LOG_INFO("Particles",
        "Heating initialised (cp_l=", p_.cp_l, " emiss=", p_.emissivity, ")");
}

void ParticleHeating::resize_to(std::size_t n) {
    if (T_p_.size() < n) T_p_.resize(n, p_.T_p_init);
}

double ParticleHeating::parcel_temperature(std::size_t i) const {
    return (i < T_p_.size()) ? T_p_[i] : p_.T_p_init;
}

void ParticleHeating::set_parcel_temperature(std::size_t i, double T) {
    if (i < T_p_.size()) T_p_[i] = std::clamp(T, p_.T_p_min, p_.T_p_max);
}

double ParticleHeating::apply(double dt, LagrangianTracker& tracker,
                              solver::FieldRegistry& F) {
    if (!mesh_) return 0.0;
    auto& parts = tracker.mutable_particles();
    resize_to(parts.size());

    const std::size_t nC = mesh_->cells().size();
    const auto* Tg = F.find_scalar("T");
    if (!Tg) return 0.0;
    const auto* Trad = p_.includeRadiation ? F.find_scalar("T_radiation") : nullptr;
    auto& Sq = F.scalar("S_particle_energy", nC);

    double Q_total = 0.0;
    const auto& Cc = mesh_->cells();
    for (std::size_t i = 0; i < parts.size(); ++i) {
        auto& pt = parts[i];
        if (!pt.active || pt.cell >= nC || pt.mass <= 0.0) continue;
        const double Tg_c = std::max((*Tg)[pt.cell], 1.0);
        // Diameter from mass (uses spec rho_l proxy: 1000 kg/m³ default).
        const double rho_l = 1000.0;
        const double d_p   = std::cbrt(6.0 * pt.mass / (PI * rho_l));
        const double A_p   = PI * d_p * d_p;
        const double rho_g = p_.rho_g_ref * (293.15 / Tg_c);
        const double Vrel  = std::sqrt(pt.v.x*pt.v.x + pt.v.y*pt.v.y + pt.v.z*pt.v.z);
        const double Re    = rho_g * Vrel * d_p / std::max(p_.mu_g, 1e-30);
        const double Pr    = p_.mu_g * p_.cp_g / std::max(p_.k_g, 1e-30);
        const double Nu    = 2.0 + 0.6 * std::sqrt(std::max(Re, 0.0)) * std::cbrt(Pr);
        const double h_c   = Nu * p_.k_g / std::max(d_p, 1e-30);
        const double Q_conv = h_c * A_p * (Tg_c - T_p_[i]);
        double Q_rad = 0.0;
        if (Trad) {
            const double Tr2 = (*Trad)[pt.cell] * (*Trad)[pt.cell];
            const double Tp2 = T_p_[i] * T_p_[i];
            Q_rad = p_.emissivity * SIGMA * A_p * (Tr2 * Tr2 - Tp2 * Tp2);
        }
        const double Q  = Q_conv + Q_rad;
        const double dT = dt * Q / std::max(pt.mass * p_.cp_l, 1e-30);
        T_p_[i] = std::clamp(T_p_[i] + dT, p_.T_p_min, p_.T_p_max);
        // Energy sink in the fluid cell (heat moves out of gas into droplet).
        Sq[pt.cell] -= Q / std::max(Cc.volume[pt.cell], 1e-30);
        Q_total += Q;
    }
    return Q_total;
}

}  // namespace simall::particles
