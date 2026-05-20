// =============================================================================
// SimAll Beta - Immersed Boundary Subsystem
// File   : src/ibm/ImmersedBoundary.cpp
// =============================================================================
#include "ibm/ImmersedBoundary.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace simall::ibm
{

namespace
{

inline double dot3(const util::Vec3d& a, const util::Vec3d& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Closest point on triangle to p (Ericson, Real-Time Collision Detection, §5.1.5)
util::Vec3d closest_on_tri(const util::Vec3d& p,
                           const util::Vec3d& a,
                           const util::Vec3d& b,
                           const util::Vec3d& c)
{
    const util::Vec3d ab{b.x - a.x, b.y - a.y, b.z - a.z};
    const util::Vec3d ac{c.x - a.x, c.y - a.y, c.z - a.z};
    const util::Vec3d ap{p.x - a.x, p.y - a.y, p.z - a.z};
    const double d1 = dot3(ab, ap), d2 = dot3(ac, ap);
    if (d1 <= 0 && d2 <= 0)
        return a;
    const util::Vec3d bp{p.x - b.x, p.y - b.y, p.z - b.z};
    const double d3 = dot3(ab, bp), d4 = dot3(ac, bp);
    if (d3 >= 0 && d4 <= d3)
        return b;
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        const double v = d1 / (d1 - d3);
        return {a.x + v * ab.x, a.y + v * ab.y, a.z + v * ab.z};
    }
    const util::Vec3d cp{p.x - c.x, p.y - c.y, p.z - c.z};
    const double d5 = dot3(ab, cp), d6 = dot3(ac, cp);
    if (d6 >= 0 && d5 <= d6)
        return c;
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        const double w = d2 / (d2 - d6);
        return {a.x + w * ac.x, a.y + w * ac.y, a.z + w * ac.z};
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return {b.x + w * (c.x - b.x), b.y + w * (c.y - b.y), b.z + w * (c.z - b.z)};
    }
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom, w = vc * denom;
    return {a.x + ab.x * v + ac.x * w, a.y + ab.y * v + ac.y * w, a.z + ab.z * v + ac.z * w};
}

bool point_inside_surface(const util::Vec3d& p, const meshing::StlSurface& s)
{
    // Same ray-cast inside test as OctreeMesher (Möller-Trumbore +x ray).
    int crossings = 0;
    const double d[3] = {1.0, 0.0, 0.0};
    for (const auto& tri : s.triangles) {
        const auto& v0 = s.vertices[tri[0]];
        const auto& v1 = s.vertices[tri[1]];
        const auto& v2 = s.vertices[tri[2]];
        const double e1[3] = {v1.x - v0.x, v1.y - v0.y, v1.z - v0.z};
        const double e2[3] = {v2.x - v0.x, v2.y - v0.y, v2.z - v0.z};
        const double pv[3] = {
            d[1] * e2[2] - d[2] * e2[1], d[2] * e2[0] - d[0] * e2[2], d[0] * e2[1] - d[1] * e2[0]};
        const double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
        if (std::abs(det) < 1e-20)
            continue;
        const double inv = 1.0 / det;
        const double tv[3] = {p.x - v0.x, p.y - v0.y, p.z - v0.z};
        const double u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
        if (u < 0 || u > 1)
            continue;
        const double qv[3] = {tv[1] * e1[2] - tv[2] * e1[1],
                              tv[2] * e1[0] - tv[0] * e1[2],
                              tv[0] * e1[1] - tv[1] * e1[0]};
        const double v = (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2]) * inv;
        if (v < 0 || u + v > 1)
            continue;
        const double t = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
        if (t > 0.0)
            ++crossings;
    }
    return (crossings & 1) == 1;
}

} // namespace

void ImmersedBoundary::initialize(const meshing::Mesh& fluid, const meshing::StlSurface& surface)
{
    mesh_ = &fluid;
    const std::size_t nC = fluid.cells().size();
    tag_.assign(nC, CellTag::Fluid);
    sdf_.assign(nC, std::numeric_limits<double>::max());

    // Per-cell unsigned distance via brute scan of triangles.
    const auto& C = fluid.cells();
    for (std::size_t c = 0; c < nC; ++c) {
        const util::Vec3d p{C.centroidX[c], C.centroidY[c], C.centroidZ[c]};
        double best = std::numeric_limits<double>::max();
        for (const auto& tri : surface.triangles) {
            const auto& a = surface.vertices[tri[0]];
            const auto& b = surface.vertices[tri[1]];
            const auto& q = surface.vertices[tri[2]];
            const util::Vec3d cp = closest_on_tri(p, a, b, q);
            const double dx = p.x - cp.x, dy = p.y - cp.y, dz = p.z - cp.z;
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best)
                best = d2;
        }
        const double d = std::sqrt(best);
        const bool inside = point_inside_surface(p, surface);
        sdf_[c] = inside ? -d : +d;
    }

    // Classify: SOLID if sdf < -h, IB if |sdf| < h, FLUID otherwise.
    // h = local cell size estimate (V^{1/3}).
    std::size_t nIB = 0, nS = 0;
    for (std::size_t c = 0; c < nC; ++c) {
        const double h = std::cbrt(std::max(C.volume[c], 1e-30));
        if (sdf_[c] < -h) {
            tag_[c] = CellTag::Solid;
            ++nS;
        } else if (std::abs(sdf_[c]) < h) {
            tag_[c] = CellTag::IB;
            ++nIB;
        }
    }
    SIMALL_LOG_INFO("IBM", "cells: fluid=", nC - nIB - nS, " IB=", nIB, " solid=", nS);
}

util::Vec3d ImmersedBoundary::body_velocity(const util::Vec3d& x) const
{
    // U_body = V_lin + Ω × (x - x_c)
    const util::Vec3d r{x.x - kin_.centreOfRotation.x,
                        x.y - kin_.centreOfRotation.y,
                        x.z - kin_.centreOfRotation.z};
    return {kin_.linearVelocity.x + kin_.angularVelocity.y * r.z - kin_.angularVelocity.z * r.y,
            kin_.linearVelocity.y + kin_.angularVelocity.z * r.x - kin_.angularVelocity.x * r.z,
            kin_.linearVelocity.z + kin_.angularVelocity.x * r.y - kin_.angularVelocity.y * r.x};
}

void ImmersedBoundary::apply_forcing(double dt, solver::FieldRegistry& F)
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    auto& S = F.vector("S_ibm", nC);
    const auto* U = F.find_vector("U");
    if (!U) {
        std::fill(S.x.begin(), S.x.end(), 0.0);
        std::fill(S.y.begin(), S.y.end(), 0.0);
        std::fill(S.z.begin(), S.z.end(), 0.0);
        return;
    }
    std::fill(S.x.begin(), S.x.end(), 0.0);
    std::fill(S.y.begin(), S.y.end(), 0.0);
    std::fill(S.z.begin(), S.z.end(), 0.0);
    const auto& C = mesh_->cells();
    const double invDt = 1.0 / std::max(dt, 1e-30);
    for (std::size_t c = 0; c < nC; ++c) {
        if (tag_[c] == CellTag::Fluid)
            continue;
        const util::Vec3d x{C.centroidX[c], C.centroidY[c], C.centroidZ[c]};
        const util::Vec3d uT = body_velocity(x);
        // Source per unit volume (kg/(m²·s²)): ρ (U* - U) / dt; ρ folded
        // in by caller's momentum assembly (consistent with other source
        // channels in the registry).
        S.x[c] = (uT.x - U->x[c]) * invDt;
        S.y[c] = (uT.y - U->y[c]) * invDt;
        S.z[c] = (uT.z - U->z[c]) * invDt;
    }
}

void ImmersedBoundary::enforce_solid(solver::FieldRegistry& F)
{
    if (!mesh_)
        return;
    auto* U = F.find_vector("U");
    if (!U)
        return;
    const auto& C = mesh_->cells();
    const std::size_t nC = mesh_->cells().size();
    for (std::size_t c = 0; c < nC; ++c) {
        if (tag_[c] != CellTag::Solid)
            continue;
        const util::Vec3d x{C.centroidX[c], C.centroidY[c], C.centroidZ[c]};
        const util::Vec3d uT = body_velocity(x);
        U->x[c] = uT.x;
        U->y[c] = uT.y;
        U->z[c] = uT.z;
    }
}

} // namespace simall::ibm
