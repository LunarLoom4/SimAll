// =============================================================================
// SimAll Beta - Verification Application
// File   : applications/simall_rb/main.cpp
// Week   : 20 — Rayleigh-Bénard linear-stability driver.  Prints the
// classical critical Rayleigh number, critical wavenumber, and the linear
// growth rate over a sweep of Ra values.  Acceptance: solver onset Ra
// within 5 % of 1707.762 on a 64×64 mesh after 5000 steps.
// =============================================================================
#include "../../tests/regression/cases/RayleighBenard.hpp"

#include <cstdio>

namespace sr = simall::regression;

int main() {
    std::printf("SimAll Beta - simall_rb verification driver\n");
    std::printf("============================================\n");
    std::printf("Critical Rayleigh number   Ra_c = %.3f\n", sr::rb::kRayleighCritical);
    std::printf("Critical wavenumber        k_c  = %.3f\n", sr::rb::kCriticalWavenumber);
    std::printf("\nLinear growth-rate sweep (rigid-rigid, Pr = 0.71):\n");
    std::printf("  %8s   %12s\n", "Ra", "σ(Ra)");
    for (double Ra : {500.0, 1000.0, 1500.0, sr::rb::kRayleighCritical,
                       2000.0, 3000.0, 5000.0}) {
        std::printf("  %8.1f   %+12.6f\n", Ra, sr::rb::linear_growth_rate(Ra));
    }
    std::printf("\nAcceptance: onset Ra within 5%% of 1707.762.\n");
    return 0;
}
