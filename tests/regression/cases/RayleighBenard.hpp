// =============================================================================
// SimAll Beta - tests/regression/cases/RayleighBenard.hpp
// Rayleigh-Bénard convection: critical Rayleigh number for the onset of
// 2-D rolls between rigid no-slip plates with conducting top and bottom
// boundaries is Ra_c ≈ 1707.762 with critical wavenumber k_c ≈ 3.117.
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

namespace simall::regression::rb
{

inline constexpr double kRayleighCritical = 1707.762;
inline constexpr double kCriticalWavenumber = 3.117;

// Linear-theory growth rate (positive when Ra > Ra_c, negative below).
[[nodiscard]] inline double linear_growth_rate(double Ra)
{
    return (Ra - kRayleighCritical) / kRayleighCritical;
}

} // namespace simall::regression::rb
