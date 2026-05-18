// =============================================================================
// SimAll Beta - tests/regression/cases/Pipe.hpp
// Hagen-Poiseuille laminar pipe flow.  Used as the simplest closed-form
// check: parabolic velocity profile u(r) = 2·u_mean·(1 - r²/R²) and
// Darcy friction factor f = 64/Re.
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

namespace simall::regression::pipe {

[[nodiscard]] inline double darcy_friction_factor_laminar(double Re) {
    return 64.0 / std::max(Re, 1e-30);
}

// Colebrook turbulent friction factor (smooth pipe limit, Prandtl).
[[nodiscard]] inline double prandtl_friction_factor(double Re) {
    // Iterate 1/sqrt(f) = 2 log10(Re sqrt(f)) - 0.8
    double f = 0.02;
    for (int k = 0; k < 30; ++k) {
        const double lhs = 1.0 / std::sqrt(f);
        const double rhs = 2.0 * std::log10(Re * std::sqrt(f)) - 0.8;
        f = 1.0 / ((rhs + lhs) * 0.5) / ((rhs + lhs) * 0.5);
        if (std::abs(lhs - rhs) < 1e-8) break;
    }
    return f;
}

}  // namespace simall::regression::pipe
