// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ConeInjector.cpp
// =============================================================================
#include "particles/ConeInjector.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::particles
{

namespace
{
constexpr double PI = 3.14159265358979323846;
util::Vec3d normalise(util::Vec3d v)
{
    const double m = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return (m > 1e-30) ? util::Vec3d{v.x / m, v.y / m, v.z / m} : util::Vec3d{1, 0, 0};
}
void basis(const util::Vec3d& n, util::Vec3d& a, util::Vec3d& b)
{
    util::Vec3d ref = (std::abs(n.z) < 0.9) ? util::Vec3d{0, 0, 1} : util::Vec3d{1, 0, 0};
    a = normalise(
        {n.y * ref.z - n.z * ref.y, n.z * ref.x - n.x * ref.z, n.x * ref.y - n.y * ref.x});
    b = {n.y * a.z - n.z * a.y, n.z * a.x - n.x * a.z, n.x * a.y - n.y * a.x};
}
} // namespace

void ConeInjector::initialize(ConeInjectorProps props, DiameterSampler sampler)
{
    p_ = props;
    p_.axis = normalise(p_.axis);
    diaSampler_ = std::move(sampler);
    rng_.seed(p_.rngSeed);
    SIMALL_LOG_INFO("Particles",
                    "Cone injector kind=",
                    (p_.kind == ConeKind::Solid ? "Solid" : "Hollow"),
                    " outer=",
                    p_.theta_outer,
                    " inner=",
                    p_.theta_inner,
                    " swirl=",
                    p_.swirl_angle);
}

std::size_t ConeInjector::inject(double dt, std::size_t Nparcels, LagrangianTracker& tracker)
{
    (void) dt;
    if (Nparcels == 0)
        return 0;
    util::Vec3d a, b;
    basis(p_.axis, a, b);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::uniform_real_distribution<double> Uphi(0.0, 2.0 * PI);

    const double cos_outer = std::cos(p_.theta_outer);
    const double cos_inner = std::cos(p_.theta_inner);
    auto& parts = tracker.mutable_particles();
    for (std::size_t i = 0; i < Nparcels; ++i) {
        const double u = U(rng_);
        // Uniform sampling in cos θ on annulus [cos_outer, cos_inner].
        const double cmin = (p_.kind == ConeKind::Hollow) ? cos_outer : cos_outer;
        const double cmax = (p_.kind == ConeKind::Hollow) ? cos_inner : 1.0;
        const double ct = cmin + u * (cmax - cmin);
        const double st = std::sqrt(std::max(1.0 - ct * ct, 0.0));
        const double phi = Uphi(rng_) + p_.swirl_angle;
        const double cphi = std::cos(phi), sphi = std::sin(phi);

        util::Vec3d dir{ct * p_.axis.x + st * (cphi * a.x + sphi * b.x),
                        ct * p_.axis.y + st * (cphi * a.y + sphi * b.y),
                        ct * p_.axis.z + st * (cphi * a.z + sphi * b.z)};
        const double d = diaSampler_ ? diaSampler_() : p_.diameter_const;
        ParticleState s;
        s.x = p_.origin;
        s.v = {p_.speed * dir.x, p_.speed * dir.y, p_.speed * dir.z};
        s.mass = (PI / 6.0) * d * d * d * p_.rho_l;
        s.inertia = 0.1 * s.mass * d * d;
        s.cell = static_cast<meshing::CellId>(-1);
        s.active = true;
        parts.push_back(s);
    }
    return Nparcels;
}

} // namespace simall::particles
