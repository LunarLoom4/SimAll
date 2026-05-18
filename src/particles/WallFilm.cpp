// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/WallFilm.cpp
// =============================================================================
#include "particles/WallFilm.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::particles {

namespace {
constexpr double PI = 3.14159265358979323846;
}

void WallFilm::initialize(const meshing::Mesh& mesh, WallFilmProps props) {
    mesh_ = &mesh;
    p_    = props;
    const std::size_t nF = mesh_->faces().size();
    h_.assign(nF, 0.0);
    m_.assign(nF, 0.0);
    u_.assign(nF, util::Vec3d{0,0,0});
    SIMALL_LOG_INFO("Particles",
        "Wall film init: faces=", nF, " zones=", p_.wallZones.size());
}

bool WallFilm::is_wall_face(std::int32_t zone) const {
    for (auto z : p_.wallZones) if (z == zone) return true;
    return false;
}

double WallFilm::film_thickness(std::size_t f) const {
    return (f < h_.size()) ? h_[f] : 0.0;
}
double WallFilm::film_mass(std::size_t f) const {
    return (f < m_.size()) ? m_[f] : 0.0;
}

double WallFilm::deposit(double dt, LagrangianTracker& tracker) {
    (void)dt;
    if (!mesh_) return 0.0;
    const auto& F = mesh_->faces();
    const auto& C = mesh_->cells();
    auto& parts   = tracker.mutable_particles();
    double m_dep_total = 0.0;

    for (auto& pt : parts) {
        if (!pt.active || pt.cell < 0) continue;
        const std::size_t off = C.faceOffsets[pt.cell];
        const std::size_t end = C.faceOffsets[pt.cell + 1];
        // Find the first wall face of the host cell whose outward normal
        // sees the parcel moving toward it.
        for (std::size_t k = off; k < end; ++k) {
            const std::size_t fi = C.faceIndices[k];
            if (!is_wall_face(F.boundaryZone[fi])) continue;
            const double aMag = std::sqrt(F.areaX[fi]*F.areaX[fi] +
                                          F.areaY[fi]*F.areaY[fi] +
                                          F.areaZ[fi]*F.areaZ[fi]);
            if (aMag < 1e-30) continue;
            const double nx = F.areaX[fi] / aMag;
            const double ny = F.areaY[fi] / aMag;
            const double nz = F.areaZ[fi] / aMag;
            const double vn = pt.v.x*nx + pt.v.y*ny + pt.v.z*nz;
            if (vn <= 0.0) continue;     // moving away from wall

            // Diameter from mass.
            const double d_p = std::cbrt(6.0 * pt.mass / (PI * p_.rho_l));
            const double We  = p_.rho_l * vn * vn * d_p / std::max(p_.sigma, 1e-30);

            double depositFrac = 0.0;
            if      (We < p_.Weks)  depositFrac = 1.0;            // stick
            else if (We < p_.Wespd) depositFrac = 1.0;            // spread
            else if (We < p_.Wereb) depositFrac = 0.0;            // rebound
            else                    depositFrac = 1.0 - p_.splash_ratio; // splash

            if (depositFrac <= 0.0) {
                // Pure rebound: reflect velocity, no mass added.
                pt.v.x -= 2.0 * vn * nx;
                pt.v.y -= 2.0 * vn * ny;
                pt.v.z -= 2.0 * vn * nz;
                break;
            }
            const double dm = depositFrac * pt.mass;
            m_[fi]    += dm;
            h_[fi]     = m_[fi] / (p_.rho_l * std::max(aMag, 1e-30));
            h_[fi]     = std::clamp(h_[fi], 0.0, p_.h_max);
            // Mass-weighted film tangential momentum (subtract normal comp).
            const double tvx = pt.v.x - vn * nx;
            const double tvy = pt.v.y - vn * ny;
            const double tvz = pt.v.z - vn * nz;
            const double w   = dm / std::max(m_[fi], 1e-30);
            u_[fi].x = (1 - w) * u_[fi].x + w * tvx;
            u_[fi].y = (1 - w) * u_[fi].y + w * tvy;
            u_[fi].z = (1 - w) * u_[fi].z + w * tvz;
            m_dep_total += dm;
            if (depositFrac >= 1.0) pt.active = false;
            else { pt.mass *= (1.0 - depositFrac); }
            break;
        }
    }
    return m_dep_total;
}

void WallFilm::advect(double dt) {
    if (!mesh_) return;
    // Simple shear-driven decay: ν h /δ² damping representing wall friction.
    for (std::size_t i = 0; i < h_.size(); ++i) {
        if (h_[i] < p_.h_min) { h_[i] = 0.0; m_[i] = 0.0; u_[i] = {0,0,0}; continue; }
        const double tauWall = p_.mu_l / std::max(h_[i], 1e-12);
        const double decay   = std::exp(-tauWall * dt / std::max(p_.rho_l * h_[i], 1e-12));
        u_[i].x *= decay; u_[i].y *= decay; u_[i].z *= decay;
    }
}

}  // namespace simall::particles
