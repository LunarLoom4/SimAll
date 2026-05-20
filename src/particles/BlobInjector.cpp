// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/BlobInjector.cpp
// =============================================================================
#include "particles/BlobInjector.hpp"

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

void BlobInjector::initialize(BlobInjectorProps props)
{
    p_ = props;
    p_.axis = normalise(p_.axis);
    rng_.seed(p_.rngSeed);
    SIMALL_LOG_INFO("Particles",
                    "Blob injector: d_nozzle=",
                    p_.nozzle_diameter,
                    " U_inj=",
                    injection_velocity(),
                    " duration=",
                    p_.duration);
}

double BlobInjector::injection_velocity() const noexcept
{
    return p_.Cd * std::sqrt(2.0 * std::max(p_.p_inj - p_.p_amb, 0.0) / std::max(p_.rho_l, 1e-30));
}

std::size_t BlobInjector::inject(double dt,
                                 std::size_t Nparcels,
                                 double elapsed,
                                 LagrangianTracker& tracker)
{
    if (elapsed < 0.0 || elapsed > p_.duration || Nparcels == 0)
        return 0;
    util::Vec3d a, b;
    basis(p_.axis, a, b);
    std::uniform_real_distribution<double> Ur(0.0, 1.0);
    std::uniform_real_distribution<double> Uth(0.0, 2.0 * PI);

    const double U_inj = injection_velocity();
    const double A_noz = PI * 0.25 * p_.nozzle_diameter * p_.nozzle_diameter;
    const double mdot = p_.Cd * p_.rho_l * U_inj * A_noz;
    const double m_target = mdot * dt;
    const double m_blob = (PI / 6.0) * std::pow(p_.nozzle_diameter, 3) * p_.rho_l;
    const double scale = (m_blob * Nparcels > 1e-30) ? m_target / (m_blob * Nparcels) : 1.0;

    auto& parts = tracker.mutable_particles();
    for (std::size_t i = 0; i < Nparcels; ++i) {
        ParticleState s;
        const double r = (0.5 * p_.nozzle_diameter) * std::sqrt(Ur(rng_));
        const double th = Uth(rng_);
        s.x = {p_.origin.x + r * (std::cos(th) * a.x + std::sin(th) * b.x),
               p_.origin.y + r * (std::cos(th) * a.y + std::sin(th) * b.y),
               p_.origin.z + r * (std::cos(th) * a.z + std::sin(th) * b.z)};
        s.v = {U_inj * p_.axis.x, U_inj * p_.axis.y, U_inj * p_.axis.z};
        s.mass = m_blob * scale;
        s.inertia = 0.1 * s.mass * p_.nozzle_diameter * p_.nozzle_diameter;
        s.cell = static_cast<meshing::CellId>(-1);
        s.active = true;
        parts.push_back(s);
    }
    return Nparcels;
}

} // namespace simall::particles
