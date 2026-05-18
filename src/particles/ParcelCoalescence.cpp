// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ParcelCoalescence.cpp
// =============================================================================
#include "particles/ParcelCoalescence.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace simall::particles {

namespace {
constexpr double PI = 3.14159265358979323846;
}

void ParcelCoalescence::initialize(const meshing::Mesh& mesh,
                                   CoalescenceProps props) {
    mesh_ = &mesh;
    p_    = props;
    rng_.seed(p_.rngSeed);
    SIMALL_LOG_INFO("Particles",
        "Coalescence init: σ=", p_.sigma, " We_crit=", p_.We_crit);
}

std::size_t ParcelCoalescence::apply(double dt, LagrangianTracker& tracker) {
    if (!mesh_) return 0;
    auto& parts = tracker.mutable_particles();
    if (parts.size() < 2) return 0;

    // Bucket active parcels by host cell.
    std::unordered_map<meshing::CellId, std::vector<std::size_t>> bucket;
    bucket.reserve(parts.size());
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (parts[i].active && parts[i].cell >= 0)
            bucket[parts[i].cell].push_back(i);
    }
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::size_t events = 0;
    const auto& Cc = mesh_->cells();

    for (auto& [cid, ids] : bucket) {
        if (ids.size() < 2) continue;
        const double Vcell = std::max(Cc.volume[cid], 1e-30);
        // Random shuffle and sweep adjacent pairs.
        std::shuffle(ids.begin(), ids.end(), rng_);
        for (std::size_t k = 0; k + 1 < ids.size(); k += 2) {
            std::size_t i = ids[k], j = ids[k + 1];
            auto& pi = parts[i]; auto& pj = parts[j];
            if (!pi.active || !pj.active) continue;
            const double di = std::cbrt(6.0 * pi.mass / (PI * p_.rho_l));
            const double dj = std::cbrt(6.0 * pj.mass / (PI * p_.rho_l));
            const double ri = 0.5 * di, rj = 0.5 * dj;
            const double dvx = pi.v.x - pj.v.x;
            const double dvy = pi.v.y - pj.v.y;
            const double dvz = pi.v.z - pj.v.z;
            const double Vrel = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            if (Vrel < 1e-12) continue;

            // Parcel j number of droplets (assume each parcel = 1 droplet
            // multiplied by parcel count — fall back to mass ratio).
            const double Nj_phys = pj.mass / std::max(
                (PI / 6.0) * dj * dj * dj * p_.rho_l, 1e-30);
            const double nu  = PI * (ri + rj) * (ri + rj) * Vrel * Nj_phys / Vcell;
            const double P   = 1.0 - std::exp(-nu * dt);
            if (U(rng_) >= P) continue;

            const double r_small = std::min(ri, rj);
            const double We = p_.rho_l * Vrel * Vrel * r_small / std::max(p_.sigma, 1e-30);
            if (We < p_.We_crit) {
                // Coalesce: merge into the heavier parcel; deactivate lighter.
                std::size_t big = (pi.mass >= pj.mass) ? i : j;
                std::size_t sml = (big == i) ? j : i;
                const double m_new = parts[big].mass + parts[sml].mass;
                const double inv   = 1.0 / std::max(m_new, 1e-30);
                parts[big].v.x = (parts[big].mass * parts[big].v.x +
                                  parts[sml].mass * parts[sml].v.x) * inv;
                parts[big].v.y = (parts[big].mass * parts[big].v.y +
                                  parts[sml].mass * parts[sml].v.y) * inv;
                parts[big].v.z = (parts[big].mass * parts[big].v.z +
                                  parts[sml].mass * parts[sml].v.z) * inv;
                parts[big].mass = m_new;
                const double d_new = std::cbrt(6.0 * m_new / (PI * p_.rho_l));
                parts[big].inertia = 0.1 * m_new * d_new * d_new;
                parts[sml].active  = false;
                ++events;
            } else {
                // Grazing collision — elastic momentum swap (mass-conserving).
                const double m_tot = pi.mass + pj.mass;
                const double ux = (pi.mass * pi.v.x + pj.mass * pj.v.x) / m_tot;
                const double uy = (pi.mass * pi.v.y + pj.mass * pj.v.y) / m_tot;
                const double uz = (pi.mass * pi.v.z + pj.mass * pj.v.z) / m_tot;
                pi.v = {2*ux - pi.v.x, 2*uy - pi.v.y, 2*uz - pi.v.z};
                pj.v = {2*ux - pj.v.x, 2*uy - pj.v.y, 2*uz - pj.v.z};
            }
        }
    }
    return events;
}

}  // namespace simall::particles
