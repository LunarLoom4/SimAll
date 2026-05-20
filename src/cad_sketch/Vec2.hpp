// =============================================================================
// SimAll Beta -- Parametric CAD / Sketcher
// File   : src/cad_sketch/Vec2.hpp
// Phase  : 23 Pass 23.1 (2D sketcher data model + constraint solver)
//
// Lightweight POD 2D-vector type used by the sketcher.  The sketcher is the
// only consumer of true 2D geometry in the codebase -- the rest of the CFD
// pipeline is 3D -- so we keep Vec2 local to cad_sketch instead of pushing
// it up into utilities/MathTypes.hpp.  Header is Eigen-free on purpose so
// it can be included by the public sketcher API without dragging in heavy
// linear-algebra templates.
// =============================================================================
#pragma once

#include <cmath>

namespace simall::cad::sketch {

struct Vec2 {
    double x{0.0};
    double y{0.0};

    constexpr Vec2() = default;
    constexpr Vec2(double xv, double yv) noexcept : x(xv), y(yv) {}

    constexpr Vec2 operator+(const Vec2& o) const noexcept { return {x + o.x, y + o.y}; }
    constexpr Vec2 operator-(const Vec2& o) const noexcept { return {x - o.x, y - o.y}; }
    constexpr Vec2 operator*(double s)      const noexcept { return {x * s, y * s}; }
    constexpr Vec2 operator/(double s)      const noexcept { return {x / s, y / s}; }
    constexpr Vec2 operator-()              const noexcept { return {-x, -y}; }

    constexpr double dot  (const Vec2& o) const noexcept { return x * o.x + y * o.y; }
    constexpr double cross(const Vec2& o) const noexcept { return x * o.y - y * o.x; }

    [[nodiscard]] double norm()  const noexcept { return std::sqrt(x * x + y * y); }
    [[nodiscard]] double norm2() const noexcept { return x * x + y * y; }

    [[nodiscard]] Vec2 normalized() const noexcept {
        const double n = norm();
        return n > 0.0 ? Vec2{x / n, y / n} : Vec2{};
    }
};

[[nodiscard]] inline double distance(const Vec2& a, const Vec2& b) noexcept {
    return (a - b).norm();
}

}  // namespace simall::cad::sketch
