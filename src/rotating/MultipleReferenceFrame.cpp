// =============================================================================
// SimAll Beta - Rotating-Frame Subsystem
// File   : src/rotating/MultipleReferenceFrame.cpp
// =============================================================================
#include "rotating/MultipleReferenceFrame.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::rotating
{

void MultipleReferenceFrame::initialize(const meshing::Mesh& m,
                                        const std::vector<meshing::ZoneId>& cz)
{
    mesh_ = &m;
    const std::size_t nC = m.cells().size();
    wx_.assign(nC, 0.0);
    wy_.assign(nC, 0.0);
    wz_.assign(nC, 0.0);
    rx_.assign(nC, 0.0);
    ry_.assign(nC, 0.0);
    rz_.assign(nC, 0.0);
    active_.assign(nC, 0);
    const auto& C = m.cells();
    for (const auto& z : zones_) {
        const double an =
            std::sqrt(z.axis.x * z.axis.x + z.axis.y * z.axis.y + z.axis.z * z.axis.z);
        if (an < 1e-30)
            continue;
        const double ax = z.axis.x / an, ay = z.axis.y / an, az = z.axis.z / an;
        const double wx = ax * z.omega, wy = ay * z.omega, wz = az * z.omega;
        for (std::size_t c = 0; c < nC; ++c) {
            if (cz.size() == nC && cz[c] == z.zone) {
                wx_[c] = wx;
                wy_[c] = wy;
                wz_[c] = wz;
                rx_[c] = C.centroidX[c] - z.origin.x;
                ry_[c] = C.centroidY[c] - z.origin.y;
                rz_[c] = C.centroidZ[c] - z.origin.z;
                active_[c] = 1;
            }
        }
    }
    SIMALL_LOG_INFO(
        "MRF", "initialized ", zones_.size(), " zones, ", active_count(), " active cells");
}

void MultipleReferenceFrame::apply(solver::FieldRegistry& F, double rho) const
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    auto& S = F.vector("S_MRF", nC);
    std::fill(S.x.begin(), S.x.end(), 0.0);
    std::fill(S.y.begin(), S.y.end(), 0.0);
    std::fill(S.z.begin(), S.z.end(), 0.0);
    const auto* U = F.find_vector("U");
    if (!U)
        return;

    for (std::size_t c = 0; c < nC; ++c) {
        if (!active_[c])
            continue;
        const double ux = U->x[c], uy = U->y[c], uz = U->z[c];
        const double wx = wx_[c], wy = wy_[c], wz = wz_[c];
        // Coriolis: 2 ω × u
        const double cx = 2.0 * (wy * uz - wz * uy);
        const double cy = 2.0 * (wz * ux - wx * uz);
        const double cz = 2.0 * (wx * uy - wy * ux);
        // Centrifugal: ω × (ω × r)
        const double rx = rx_[c], ry = ry_[c], rz = rz_[c];
        const double ax = wy * rz - wz * ry;
        const double ay = wz * rx - wx * rz;
        const double az = wx * ry - wy * rx;
        const double dx = wy * az - wz * ay;
        const double dy = wz * ax - wx * az;
        const double dz = wx * ay - wy * ax;
        S.x[c] = -rho * (cx + dx);
        S.y[c] = -rho * (cy + dy);
        S.z[c] = -rho * (cz + dz);
    }
}

util::Vec3d MultipleReferenceFrame::frame_velocity(std::size_t c) const
{
    if (c >= active_.size() || !active_[c])
        return {0, 0, 0};
    return {wy_[c] * rz_[c] - wz_[c] * ry_[c],
            wz_[c] * rx_[c] - wx_[c] * rz_[c],
            wx_[c] * ry_[c] - wy_[c] * rx_[c]};
}

std::size_t MultipleReferenceFrame::active_count() const noexcept
{
    std::size_t n = 0;
    for (auto a : active_)
        n += a;
    return n;
}

} // namespace simall::rotating
