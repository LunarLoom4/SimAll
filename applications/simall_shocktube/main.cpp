// =============================================================================
// SimAll Beta - Verification Application
// File   : applications/simall_shocktube/main.cpp
// Week   : 20 — Sod shock-tube reference printer.  Emits the exact
// Riemann-solution star-region values (p*, u*, ρ*) and the shock-front
// position at t = 0.2.  Acceptance: solver star-pressure within 1 %.
// =============================================================================
#include "../../tests/regression/cases/Shocktube.hpp"

#include <cstdio>

int main() {
    const auto r = simall::regression::shocktube::sod_reference();
    std::printf("SimAll Beta - simall_shocktube verification driver\n");
    std::printf("===================================================\n");
    std::printf("Sod problem (γ=1.4, x_diaphragm=0.5, t_final=0.2):\n");
    std::printf("  p*          = %.5f\n", r.pStar);
    std::printf("  u*          = %.5f\n", r.uStar);
    std::printf("  rho* (L)    = %.5f\n", r.rhoLeftStar);
    std::printf("  rho* (R)    = %.5f\n", r.rhoRightStar);
    std::printf("  shock speed = %.5f\n", r.shockSpeed);
    std::printf("\nAcceptance: |p_num − p*|/|p*| ≤ 0.01 on 200-cell uniform grid.\n");
    return 0;
}
