// =============================================================================
// SimAll Beta - tests/regression/cases/Cylinder.hpp
// Flow past a circular cylinder.  Reference correlations:
//
//   * Strouhal-Reynolds (Roshko 1954):
//        St(Re) = 0.212 · (1 - 21.2 / Re)        for 50 ≤ Re ≤ 200
//
//   * Drag coefficient (Henderson 1995 fit):
//        Cd(Re) ≈ 1 + 10·Re^{-2/3}              for 1 ≤ Re ≤ 100
//
// Acceptance tolerance for production-grade DES regressions: ~3 % on St,
// ~5 % on Cd, when averaged over ≥30 shedding cycles.
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

namespace simall::regression::cylinder {

[[nodiscard]] inline double roshko_strouhal(double Re) {
    return 0.212 * (1.0 - 21.2 / std::max(Re, 1.0));
}

[[nodiscard]] inline double henderson_drag(double Re) {
    return 1.0 + 10.0 * std::pow(Re, -2.0 / 3.0);
}

}  // namespace simall::regression::cylinder
