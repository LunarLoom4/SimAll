// =============================================================================
// SimAll Beta - Utilities Subsystem
// File   : src/utilities/MathTypes.hpp
// Phase  : 24 / 6 (shared math primitives for CAD, mesh, solver)
//
// Lightweight POD vector/matrix types. We do NOT depend on Eigen here so that
// hot-loop kernels can include these freely without dragging template-heavy
// headers. Eigen is used for linear-algebra back-ends in solver/linalg.
// =============================================================================
#pragma once

#include <array>
#include <cmath>
#include <cstdint>

namespace simall::util
{

struct alignas(32) Vec3d
{
    double x{0.0}, y{0.0}, z{0.0};

    constexpr Vec3d() = default;
    constexpr Vec3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3d operator+(const Vec3d& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3d operator-(const Vec3d& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3d operator*(double s) const { return {x * s, y * s, z * s}; }
    constexpr double dot(const Vec3d& o) const { return x * o.x + y * o.y + z * o.z; }
    constexpr Vec3d cross(const Vec3d& o) const
    {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    double norm() const { return std::sqrt(x * x + y * y + z * z); }
    double norm2() const { return x * x + y * y + z * z; }
    Vec3d normalized() const
    {
        double n = norm();
        return n > 0 ? Vec3d{x / n, y / n, z / n} : Vec3d{};
    }
};

struct BoundingBox
{
    Vec3d min{1e300, 1e300, 1e300};
    Vec3d max{-1e300, -1e300, -1e300};

    void expand(const Vec3d& p)
    {
        if (p.x < min.x)
            min.x = p.x;
        if (p.y < min.y)
            min.y = p.y;
        if (p.z < min.z)
            min.z = p.z;
        if (p.x > max.x)
            max.x = p.x;
        if (p.y > max.y)
            max.y = p.y;
        if (p.z > max.z)
            max.z = p.z;
    }
    void expand(const BoundingBox& b)
    {
        expand(b.min);
        expand(b.max);
    }
    Vec3d center() const
    {
        return {(min.x + max.x) * 0.5, (min.y + max.y) * 0.5, (min.z + max.z) * 0.5};
    }
    Vec3d extent() const { return max - min; }
    bool valid() const { return min.x <= max.x && min.y <= max.y && min.z <= max.z; }
};

using PersistentId = std::uint64_t;
inline constexpr PersistentId kInvalidId = 0;

} // namespace simall::util
