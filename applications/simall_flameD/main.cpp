// =============================================================================
// SimAll Beta - Verification Application
// File   : applications/simall_flameD/main.cpp
// Week   : 20 — Sandia Flame D centreline reference printer.  Emits the
// far-field mixture-fraction scaling ξ(x/D) ≈ 5.4 / (x/D) and the
// reported peak temperature location.  Acceptance: solver ξ within 8 %
// at x/D ∈ {15, 30, 45, 60}.
// =============================================================================
#include "../../tests/regression/cases/FlameD.hpp"

#include <cstdio>

namespace sr = simall::regression;

int main() {
    std::printf("SimAll Beta - simall_flameD verification driver\n");
    std::printf("================================================\n");
    std::printf("Sandia Flame D (CH4/air piloted, Re ≈ 22400)\n\n");
    std::printf("  %6s   %14s\n", "x/D", "xi_centerline");
    for (double xd : {5.0, 15.0, 30.0, 45.0, 60.0, 75.0}) {
        std::printf("  %6.1f   %14.6f\n", xd,
                     sr::flameD::centerline_mixture_fraction(xd));
    }
    std::printf("\nPeak centreline temperature: %.0f K at x/D = %.1f\n",
                 sr::flameD::kPeakTemperatureK,
                 sr::flameD::kPeakTemperatureXOverD);
    std::printf("\nAcceptance: solver ξ within 8%% at x/D ∈ {15,30,45,60}.\n");
    return 0;
}
