// =============================================================================
// SimAll Beta - Morphing Subsystem
// File   : src/morphing/FfdBox.hpp
// Week   : 18
//
// Free-Form Deformation (Sederberg & Parry, 1986).  Embeds a region of
// the mesh in a trivariate Bernstein-Bezier lattice; moving control
// points of the lattice produces smooth, watertight, C¹-continuous
// deformations of the embedded geometry.  Industrial CFD usage:
// gradient-based shape optimisation, parametric what-if studies, and
// fluid-structure interaction body-fitted morphing.
//
// Mathematics:
//
//     Φ(s, t, u) = Σ_{i,j,k} B_i^l(s) B_j^m(t) B_k^n(u) · P_{ijk}
//
//   * (s, t, u) ∈ [0,1]³ are normalised lattice coordinates;
//   * P_{ijk} are (l+1)(m+1)(n+1) control points;
//   * deformed point = Φ(s, t, u) of the *current* control net.
//
// The FfdBox owns the original control-net layout (a regular grid in
// world coordinates), exposes a "displace" handle on each control point,
// and provides a `morph(pointsIn, pointsOut)` operation that re-evaluates
// every embedded vertex.
// =============================================================================
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace simall::morphing
{

struct Vec3
{
    double x = 0.0, y = 0.0, z = 0.0;
};

class FfdBox
{
public:
    /// Construct a box aligned with the axes spanning [origin, origin+size]
    /// with `nu × nv × nw` control points per dimension (degree = n-1).
    FfdBox(
        const Vec3& origin, const Vec3& size, std::uint32_t nu, std::uint32_t nv, std::uint32_t nw);

    [[nodiscard]] std::uint32_t nu() const noexcept { return nu_; }
    [[nodiscard]] std::uint32_t nv() const noexcept { return nv_; }
    [[nodiscard]] std::uint32_t nw() const noexcept { return nw_; }

    [[nodiscard]] const Vec3& control_point(std::uint32_t i,
                                            std::uint32_t j,
                                            std::uint32_t k) const;
    void set_control_point(std::uint32_t i, std::uint32_t j, std::uint32_t k, const Vec3& p);
    void displace_control_point(std::uint32_t i, std::uint32_t j, std::uint32_t k, const Vec3& d);
    void reset(); ///< restore the original (undeformed) control net

    /// Map a world-space point to lattice coordinates (s,t,u) ∈ [0,1]³.
    /// Returns false if the point lies outside the box.
    [[nodiscard]] bool world_to_lattice(const Vec3& world, Vec3& stu) const;

    /// Evaluate Φ(s,t,u) given the current control net.
    [[nodiscard]] Vec3 evaluate(double s, double t, double u) const;

    /// Apply the deformation to a batch of mesh vertices.  Vertices outside
    /// the box are copied unchanged.  In-place safe (pointsIn may equal pointsOut).
    void morph(const std::vector<Vec3>& pointsIn, std::vector<Vec3>& pointsOut) const;

private:
    Vec3 origin_;
    Vec3 size_;
    std::uint32_t nu_, nv_, nw_;
    std::vector<Vec3> P_;  // control points, row-major (i + nu*(j + nv*k))
    std::vector<Vec3> P0_; // initial control net (for reset())

    [[nodiscard]] std::size_t idx(std::uint32_t i, std::uint32_t j, std::uint32_t k) const noexcept
    {
        return std::size_t(i)
               + std::size_t(nu_) * (std::size_t(j) + std::size_t(nv_) * std::size_t(k));
    }
};

} // namespace simall::morphing
