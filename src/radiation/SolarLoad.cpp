// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/SolarLoad.cpp
// =============================================================================
#include "radiation/SolarLoad.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::radiation {

namespace {
constexpr double kPi = 3.14159265358979323846;

util::Vec3d normalise(util::Vec3d v) {
    const double m = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    return (m > 1e-30) ? util::Vec3d{v.x/m, v.y/m, v.z/m}
                       : util::Vec3d{0,0,-1};
}
}  // namespace

bool SolarLoad::initialize(const meshing::Mesh& mesh,
                           solver::FieldRegistry& F,
                           const SolarProps& props) {
    mesh_ = &mesh; F_ = &F; p_ = props;
    p_.sunDir = normalise(p_.sunDir);
    q_face_.assign(mesh.faces().size(), 0.0);
    SIMALL_LOG_INFO("Radiation",
        "SolarLoad init: I_dir=", p_.I_direct,
        " I_dif=", p_.I_diffuse, " useAshrae=", p_.useAshrae);
    return true;
}

const SolarFaceSpec* SolarLoad::spec_for_zone(meshing::ZoneId z) const {
    for (const auto& s : zones_) if (s.zone == z) return &s;
    return nullptr;
}

util::Vec3d SolarLoad::ashrae_sun_dir(double t) const {
    // Simplified ASHRAE solar position algorithm (Spencer 1971, ASHRAE 2009).
    const double dayFloat = t / 86400.0;
    const int    n_day    = 1 + (static_cast<int>(std::floor(dayFloat)) % 365);
    const double hour     = std::fmod(dayFloat, 1.0) * 24.0
                          + p_.timezone_hours;
    const double gamma    = 2.0 * kPi * (n_day - 1) / 365.0;
    const double decl     = 0.006918 - 0.399912*std::cos(gamma)
                          + 0.070257*std::sin(gamma)
                          - 0.006758*std::cos(2*gamma)
                          + 0.000907*std::sin(2*gamma);
    const double EoT      = 229.18 * (0.000075
                          + 0.001868*std::cos(gamma) - 0.032077*std::sin(gamma)
                          - 0.014615*std::cos(2*gamma) - 0.040849*std::sin(2*gamma));
    const double lstm     = 15.0 * p_.timezone_hours;
    const double tc       = 4.0 * (p_.longitude_deg - lstm) + EoT;
    const double solarT   = hour + tc / 60.0;
    const double Hrad     = (solarT - 12.0) * kPi / 12.0;
    const double Lrad     = p_.latitude_deg * kPi / 180.0;
    const double sinAlt   = std::sin(Lrad) * std::sin(decl)
                          + std::cos(Lrad) * std::cos(decl) * std::cos(Hrad);
    const double altitude = std::asin(std::clamp(sinAlt, -1.0, 1.0));
    const double cosAlt   = std::cos(altitude);
    const double sinAz    = -std::cos(decl) * std::sin(Hrad) / std::max(cosAlt, 1e-9);
    const double cosAz    = (std::sin(decl) - std::sin(Lrad) * sinAlt)
                          / std::max(std::cos(Lrad) * cosAlt, 1e-9);
    const double azimuth  = std::atan2(sinAz, cosAz);
    // ŝ FROM the sun toward the surface (negate observer direction).
    return normalise({ -cosAlt * std::sin(azimuth),
                       -cosAlt * std::cos(azimuth),
                       -sinAlt });
}

double SolarLoad::apply(double time_seconds) {
    if (!mesh_ || !F_) return 0.0;
    const auto& Ff = mesh_->faces();
    const util::Vec3d s = p_.useAshrae ? ashrae_sun_dir(time_seconds)
                                       : p_.sunDir;
    const double Idir = std::min(p_.I_direct, p_.I_total_max);
    const double Idif = std::max(p_.I_diffuse, 0.0);
    // Solar-zenith cosine for ground-reflected term (sun above horizon).
    const double cosZ = std::max(-s.z, 0.0);
    double total = 0.0;
    std::fill(q_face_.begin(), q_face_.end(), 0.0);
    for (std::size_t f = 0; f < Ff.size(); ++f) {
        if (Ff.neighbor[f] != meshing::kBoundaryCell) continue;
        const auto* spec = spec_for_zone(Ff.boundaryZone[f]);
        if (!spec) continue;
        const double ax = Ff.areaX[f], ay = Ff.areaY[f], az = Ff.areaZ[f];
        const double aMag = std::sqrt(ax*ax + ay*ay + az*az);
        if (aMag < 1e-30) continue;
        // Outward face normal.
        const double nx = ax / aMag, ny = ay / aMag, nz = az / aMag;
        const double cosI = std::max(-(s.x*nx + s.y*ny + s.z*nz), 0.0);
        const double Fsky = std::max(0.0, 0.5 * (1.0 + nz));
        const double Fgrd = std::max(0.0, 0.5 * (1.0 - nz));
        const double absorb = std::clamp(spec->absorptivity, 0.0, 1.0);
        const double q = absorb * (Idir * cosI
                                 + Idif * Fsky
                                 + p_.groundReflectance
                                   * (Idir * cosZ + Idif) * Fgrd);
        q_face_[f] = q;
        total += q * aMag;
    }
    return total;
}

double SolarLoad::face_flux(std::size_t f) const {
    return (f < q_face_.size()) ? q_face_[f] : 0.0;
}

}  // namespace simall::radiation
