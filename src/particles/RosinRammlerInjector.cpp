// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/RosinRammlerInjector.cpp
// =============================================================================
#include "particles/RosinRammlerInjector.hpp"

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
// Build an orthonormal basis (a, b, n) given the direction n.
void basis(const util::Vec3d& n, util::Vec3d& a, util::Vec3d& b)
{
    util::Vec3d ref = (std::abs(n.z) < 0.9) ? util::Vec3d{0, 0, 1} : util::Vec3d{1, 0, 0};
    a = normalise(
        {n.y * ref.z - n.z * ref.y, n.z * ref.x - n.x * ref.z, n.x * ref.y - n.y * ref.x});
    b = {n.y * a.z - n.z * a.y, n.z * a.x - n.x * a.z, n.x * a.y - n.y * a.x};
}
} // namespace

void RosinRammlerInjector::initialize(RosinRammlerProps props)
{
    p_ = props;
    p_.axis = normalise(p_.axis);
    rng_.seed(p_.rngSeed);
    SIMALL_LOG_INFO(
        "Particles", "Rosin-Rammler injector: X=", p_.X, " n=", p_.n, " mdot=", p_.mass_flowrate);
}

double RosinRammlerInjector::sample_diameter()
{
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int trial = 0; trial < 32; ++trial) {
        const double u = U(rng_);
        const double d = p_.X * std::pow(-std::log(std::max(1.0 - u, 1e-30)), 1.0 / p_.n);
        if (d >= p_.d_min && d <= p_.d_max)
            return d;
    }
    return std::clamp(p_.X, p_.d_min, p_.d_max);
}

std::size_t RosinRammlerInjector::inject(double dt,
                                         std::size_t Nparcels,
                                         LagrangianTracker& tracker)
{
    if (Nparcels == 0)
        return 0;
    util::Vec3d a, b;
    basis(p_.axis, a, b);
    std::uniform_real_distribution<double> Ur(0.0, 1.0);
    std::uniform_real_distribution<double> Uth(0.0, 2.0 * PI);

    // Sample diameters and compute parcel volumes / masses first.
    std::vector<double> dlist(Nparcels);
    double mass_per_parcel_sum = 0.0;
    for (auto& d : dlist) {
        d = sample_diameter();
        const double v_drop = (PI / 6.0) * d * d * d;
        mass_per_parcel_sum += v_drop * p_.rho_l;
    }
    // Required total mass this step (ṁ · dt) — distribute uniformly across
    // parcels by computing N_phys per parcel = (m_target / N) / m_drop.
    const double m_target = p_.mass_flowrate * dt;
    const double m_avg = mass_per_parcel_sum / static_cast<double>(Nparcels);
    const double scale = (m_avg > 1e-30) ? m_target / (m_avg * Nparcels) : 1.0;

    std::size_t injected = 0;
    auto& parts = tracker.mutable_particles();
    for (std::size_t i = 0; i < Nparcels; ++i) {
        ParticleState s;
        const double r = p_.nozzleRadius * std::sqrt(Ur(rng_));
        const double th = Uth(rng_);
        s.x = {p_.origin.x + r * (std::cos(th) * a.x + std::sin(th) * b.x),
               p_.origin.y + r * (std::cos(th) * a.y + std::sin(th) * b.y),
               p_.origin.z + r * (std::cos(th) * a.z + std::sin(th) * b.z)};
        s.v = {p_.speed * p_.axis.x, p_.speed * p_.axis.y, p_.speed * p_.axis.z};
        const double v_drop = (PI / 6.0) * dlist[i] * dlist[i] * dlist[i];
        s.mass = v_drop * p_.rho_l * scale;
        s.inertia = 0.1 * s.mass * dlist[i] * dlist[i];
        s.cell = static_cast<meshing::CellId>(-1);
        s.active = true;
        parts.push_back(s);
        ++injected;
    }
    return injected;
}

} // namespace simall::particles
