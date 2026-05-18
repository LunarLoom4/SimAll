// =============================================================================
// SimAll Beta - tests/verification/MmsSolutions.hpp
// Week 19 - Method of Manufactured Solutions (MMS).
//
// We pick a smooth analytic solution u_exact(x,t), substitute it into the
// model PDE to obtain a non-zero source term S(x,t), then verify that a
// discretisation driven by S converges to u_exact at the *designed* order
// of accuracy under successive mesh refinement.
//
// Canonical MMS field used here:
//   u(x, y, t) = sin(π x) · sin(π y) · exp(-2 π² ν t)
// which exactly solves the 2-D heat equation  ∂u/∂t = ν Δu  on Ω=[0,1]²
// with homogeneous Dirichlet BCs.  Adding any forcing recovers Poisson /
// Helmholtz residuals trivially.
// =============================================================================
#pragma once

#include "../regression/RegressionFramework.hpp"
#include <cmath>

namespace simall::verification {

inline constexpr double kPi = 3.14159265358979323846;

[[nodiscard]] inline double heat_exact(double x, double y, double t, double nu) {
    return std::sin(kPi * x) * std::sin(kPi * y) * std::exp(-2.0 * kPi * kPi * nu * t);
}

[[nodiscard]] inline double laplacian_of_exact(double x, double y, double t, double nu) {
    return -2.0 * kPi * kPi * heat_exact(x, y, t, nu);
}

// MMS residual: ∂u/∂t − ν Δu  evaluated on the *exact* solution should be 0.
[[nodiscard]] inline double heat_residual(double x, double y, double t, double nu) {
    const double dudt = -2.0 * kPi * kPi * nu * heat_exact(x, y, t, nu);
    return dudt - nu * laplacian_of_exact(x, y, t, nu);
}

// 5-point central FD Laplacian on a uniform grid.  Returns max-norm error
// against the analytic Laplacian.
[[nodiscard]] inline double fd5_laplacian_error(std::size_t n) {
    const double h = 1.0 / double(n - 1);
    double maxErr = 0.0;
    for (std::size_t j = 1; j + 1 < n; ++j) {
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const double x = double(i) * h, y = double(j) * h;
            const double uc = heat_exact(x,         y,         0.0, 1.0);
            const double up = heat_exact(x + h,     y,         0.0, 1.0);
            const double um = heat_exact(x - h,     y,         0.0, 1.0);
            const double vp = heat_exact(x,         y + h,     0.0, 1.0);
            const double vm = heat_exact(x,         y - h,     0.0, 1.0);
            const double lapNum = (up - 2.0 * uc + um + vp - 2.0 * uc + vm) / (h * h);
            const double lapExt = laplacian_of_exact(x, y, 0.0, 1.0);
            maxErr = std::max(maxErr, std::abs(lapNum - lapExt));
        }
    }
    return maxErr;
}

}  // namespace simall::verification
