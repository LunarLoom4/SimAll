// =============================================================================
// SimAll Beta - tests/regression/cases/TaylorGreenVortex.hpp
// Taylor-Green vortex — analytic kinetic-energy decay for the 2-D base case:
//
//     E(t) = E_0 · exp(-2 ν k² t),   k = √2 (when L=2π domain)
//
// 3-D TGV reference is dissipation rate ε(t) tabulated by Brachet 1983; we
// reproduce the dimensionless peak ε_max ≈ 0.012 at t ≈ 9 for Re=1600.
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

namespace simall::regression::tgv
{

[[nodiscard]] inline double analytic_kinetic_energy_2d(double t,
                                                       double nu,
                                                       double k = std::sqrt(2.0),
                                                       double E0 = 1.0)
{
    return E0 * std::exp(-2.0 * nu * k * k * t);
}

// Brachet 1983 peak dissipation rate for the 3-D inviscid-limit TGV.
inline constexpr double kBrachetPeakDissipation = 0.0123;
inline constexpr double kBrachetPeakTime = 9.0;

} // namespace simall::regression::tgv
