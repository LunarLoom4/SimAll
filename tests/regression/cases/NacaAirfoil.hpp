// =============================================================================
// SimAll Beta - tests/regression/cases/NacaAirfoil.hpp
// Thin-airfoil theory benchmark: dCl/dα = 2π per radian for incompressible,
// inviscid, attached flow over a symmetric airfoil at small α.  We also
// provide the thickness distribution for NACA-4 series for the surface-mesh
// regression checker.
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

namespace simall::regression::naca
{

[[nodiscard]] inline double thin_airfoil_cl(double alphaRad)
{
    return 2.0 * 3.14159265358979323846 * alphaRad;
}

// NACA-4 half-thickness distribution: y_t(x) for max thickness t (fraction
// of chord), 0 ≤ x ≤ 1.
[[nodiscard]] inline double naca4_halfthickness(double x, double t)
{
    return 5.0 * t
           * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x
              - 0.1015 * x * x * x * x);
}

} // namespace simall::regression::naca
