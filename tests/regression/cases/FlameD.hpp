// =============================================================================
// SimAll Beta - tests/regression/cases/FlameD.hpp
// Sandia Flame D (Barlow & Frank, 1998).  Methane-air piloted jet, Re ≈ 22400.
// Reference centreline mixture-fraction decay  ξ(x/D) ≈ 5.4 / (x/D)  for
// far-field self-similar region (x/D > 30).  Centreline temperature peak
// ≈ 1900 K near x/D ≈ 45.  These two anchor numbers are used to gate
// combustion subsystem regressions.
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

namespace simall::regression::flameD {

[[nodiscard]] inline double centerline_mixture_fraction(double xOverD) {
    return (xOverD < 1.0) ? 1.0 : 5.4 / xOverD;
}

inline constexpr double kPeakTemperatureK = 1900.0;
inline constexpr double kPeakTemperatureXOverD = 45.0;

}  // namespace simall::regression::flameD
