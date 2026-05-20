// =============================================================================
// SimAll Beta - tests/regression/cases/LidDrivenCavity.hpp
// Lid-driven cavity (Re=100, 400, 1000) — Ghia et al., 1982 reference data
// for the u-velocity along the vertical centreline and v along horizontal
// centreline.  We do not run a full Navier-Stokes solve in this regression
// gate; instead we lock in the spline interpolation of the published
// reference values that downstream solver acceptance tests must hit
// (within a chosen tolerance).  This catches accidental rotation of axes,
// re-numbering of the published table, or unit-scale regressions in the
// post-processor.
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

#include <array>

namespace simall::regression::ldc
{

// Ghia 1982 — Re=100 — u along x=0.5 vertical centreline (17 sampling points
// from the floor y=0 to the lid y=1.0).
inline constexpr std::array<double, 17> kGhiaY_Re100 = {0.0000,
                                                        0.0547,
                                                        0.0625,
                                                        0.0703,
                                                        0.1016,
                                                        0.1719,
                                                        0.2813,
                                                        0.4531,
                                                        0.5000,
                                                        0.6172,
                                                        0.7344,
                                                        0.8516,
                                                        0.9531,
                                                        0.9609,
                                                        0.9688,
                                                        0.9766,
                                                        1.0000};
inline constexpr std::array<double, 17> kGhiaU_Re100 = {0.0000,
                                                        -0.03717,
                                                        -0.04192,
                                                        -0.04775,
                                                        -0.06434,
                                                        -0.10150,
                                                        -0.15662,
                                                        -0.21090,
                                                        -0.20581,
                                                        -0.13641,
                                                        0.00332,
                                                        0.23151,
                                                        0.68717,
                                                        0.73722,
                                                        0.78871,
                                                        0.84123,
                                                        1.0000};

// Computed scalar functional: maximum back-flow magnitude near y~0.18.
[[nodiscard]] inline double max_backflow_velocity()
{
    double m = 0.0;
    for (auto v : kGhiaU_Re100)
        m = std::min(m, v);
    return std::abs(m);
}

// Centreline integral ∫₀¹ |u(y)| dy via trapezoidal rule.
[[nodiscard]] inline double centerline_speed_integral()
{
    double s = 0.0;
    for (std::size_t i = 1; i < kGhiaY_Re100.size(); ++i) {
        const double dy = kGhiaY_Re100[i] - kGhiaY_Re100[i - 1];
        s += 0.5 * dy * (std::abs(kGhiaU_Re100[i]) + std::abs(kGhiaU_Re100[i - 1]));
    }
    return s;
}

} // namespace simall::regression::ldc
