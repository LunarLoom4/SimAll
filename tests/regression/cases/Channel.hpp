// =============================================================================
// SimAll Beta - tests/regression/cases/Channel.hpp
// Turbulent plane channel at Re_τ ≈ 180 (Kim, Moin, Moser 1987).  We use
// the classical Reichardt composite profile as the reference benchmark
// against which the RANS / LES post-processor velocity output is compared.
//
//     u+ = (1/κ) ln(1 + κ y+) + 7.8 · ( 1 - exp(-y+/11) − y+/11 · exp(-y+/3) )
//
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

namespace simall::regression::channel
{

[[nodiscard]] inline double reichardt_u_plus(double yPlus, double kappa = 0.41)
{
    if (yPlus <= 0.0)
        return 0.0;
    const double t1 = (1.0 / kappa) * std::log(1.0 + kappa * yPlus);
    const double t2 =
        7.8 * (1.0 - std::exp(-yPlus / 11.0) - (yPlus / 11.0) * std::exp(-yPlus / 3.0));
    return t1 + t2;
}

// Friction velocity prediction at the channel wall (h=1, u_b=1, Re_τ=180).
[[nodiscard]] inline double friction_velocity_estimate(double Re_tau, double Re_bulk)
{
    // Dean's correlation:  Cf ≈ 0.073 · Re_b^{-0.25}.  u_τ = u_b · sqrt(Cf/2).
    const double Cf = 0.073 * std::pow(Re_bulk, -0.25);
    (void) Re_tau;
    return std::sqrt(Cf / 2.0);
}

} // namespace simall::regression::channel
