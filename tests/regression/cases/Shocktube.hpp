// =============================================================================
// SimAll Beta - tests/regression/cases/Shocktube.hpp
// Sod shock tube — exact Riemann solution.  Reference test for any
// compressible scheme.  Standard initial condition:
//   left  : ρ=1.0, u=0.0, p=1.0
//   right : ρ=0.125, u=0.0, p=0.1
//   γ = 1.4
// At t=0.2 the exact post-shock pressure is ≈ 0.30313, post-shock density
// is ≈ 0.42632, contact velocity ≈ 0.92745, shock-front position ≈ 1.7522.
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

namespace simall::regression::shocktube {

struct SodReference {
    double pStar      = 0.30313;
    double uStar      = 0.92745;
    double rhoLeftStar  = 0.42632;
    double rhoRightStar = 0.26557;
    double shockSpeed = 1.7522;
};

[[nodiscard]] inline SodReference sod_reference() { return {}; }

// Helper: speed of sound for a perfect gas.
[[nodiscard]] inline double speed_of_sound(double gamma, double p, double rho) {
    return std::sqrt(gamma * p / std::max(rho, 1e-30));
}

}  // namespace simall::regression::shocktube
